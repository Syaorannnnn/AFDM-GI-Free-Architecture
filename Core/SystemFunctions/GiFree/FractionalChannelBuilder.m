classdef FractionalChannelBuilder < IChannelOperator
% FRACTIONALCHANNELBUILDER - 分数 Doppler 有效信道构造器
%
%   描述:
%   在 DAFT 等效模型下构建 GI-Free 多径有效信道矩阵，支持分数 Doppler，
%   核心使用 Dirichlet 核及其导数来构建路径响应和精细搜索基。
%
%   语法:
%   builder = FractionalChannelBuilder(cfg);
%   hPath = builder.buildPathMatrix(delay, doppler);
%   hEff  = builder.buildEffectiveChannel(pathParams);
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。
    properties (SetAccess = private)
        Config
    end

    methods

        % FRACTIONALCHANNELBUILDER 构造函数，绑定 GiFree 配置对象。
        function obj = FractionalChannelBuilder(cfg)
            arguments
                cfg (1,1) GiFreeConfig
            end
            obj.Config = cfg;
        end

        % BUILDPATHMATRIX 构造单路径有效信道矩阵（含分数 Doppler 展宽）。
        function hPath = buildPathMatrix(obj, delay, doppler)
            numSc    = obj.Config.NumSubcarriers;
            chirpC1  = obj.Config.ChirpParam1;
            chirpC2  = obj.Config.ChirpParam2;
            dirichletRadius = obj.Config.DirichletRadius;

            alphaInt = round(doppler);
            fracPart = doppler - alphaInt;
            locStep  = round(2 * numSc * chirpC1);
            locIndex = alphaInt + locStep * delay;

            kvRange  = -dirichletRadius:dirichletRadius;
            numKvPts = length(kvRange);

            dCoeffs = zeros(numKvPts, 1);
            for idx = 1:numKvPts
                dCoeffs(idx) = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvRange(idx), numSc);
            end

            pVec       = (0:numSc-1).';
            c2PSq      = chirpC2 * pVec .* pVec;
            phaseConst = numSc * chirpC1 * delay * delay;
            twoPiOverN = 2 * pi / numSc;

            maxNnz = numSc * numKvPts;
            rowAll = zeros(maxNnz, 1);
            colAll = zeros(maxNnz, 1);
            valAll = complex(zeros(maxNnz, 1));
            ptr    = 0;

            for idx = 1:numKvPts
                dc = dCoeffs(idx);
                if abs(dc) < 1e-15
                    continue;
                end

                qVec = mod(pVec + locIndex + kvRange(idx), numSc);
                phaseArg = twoPiOverN * ( ...
                    phaseConst - qVec * delay ...
                    + numSc * chirpC2 * (qVec .* qVec) ...
                    - numSc * c2PSq);
                phaseVals = exp(1j * phaseArg);

                blockRange = ptr+1 : ptr+numSc;
                rowAll(blockRange) = pVec + 1;
                colAll(blockRange) = qVec + 1;
                valAll(blockRange) = phaseVals * dc;
                ptr = ptr + numSc;
            end

            hPath = sparse(rowAll(1:ptr), colAll(1:ptr), valAll(1:ptr), numSc, numSc);
        end

        % BUILDEFFECTIVECHANNEL 将多径参数叠加得到总有效信道矩阵。
        function effectiveChannel = buildEffectiveChannel(obj, pathParams)
            numSc = obj.Config.NumSubcarriers;
            effectiveChannel  = sparse(numSc, numSc);
            for i = 1:size(pathParams, 1)
                effectiveChannel = effectiveChannel + pathParams(i,3) * ...
                    obj.buildPathMatrix(pathParams(i,1), pathParams(i,2));
            end
        end

        % BUILDPILOTRESPONSEVECTOR 构造指定导频索引处的单路径响应向量。
        function hCol = buildPilotResponseVector(obj, delay, doppler, pilotIdx)
            if nargin < 4
                pilotIdx = 0;
            end

            numSc    = obj.Config.NumSubcarriers;
            chirpC1  = obj.Config.ChirpParam1;
            chirpC2  = obj.Config.ChirpParam2;
            dirichletRadius = obj.Config.DirichletRadius;

            alphaInt = round(doppler);
            fracPart = doppler - alphaInt;
            locStep  = round(2 * numSc * chirpC1);
            locIndex = alphaInt + locStep * delay;

            hCol = zeros(numSc, 1);
            for kvShift = -dirichletRadius:dirichletRadius
                pIdx = mod(pilotIdx - locIndex - kvShift, numSc);
                phaseVal = exp(1j * 2 * pi * (chirpC1 * delay^2 ...
                    - pilotIdx * delay / numSc ...
                    + chirpC2 * pilotIdx^2 ...
                    - chirpC2 * pIdx^2));
                dCoeff = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvShift, numSc);
                hCol(pIdx + 1) = hCol(pIdx + 1) + phaseVal * dCoeff;
            end
        end

        % BUILDCOMPOSITEPILOTRESPONSE 构造包含导频幅度与序列的复合响应。
        function hComposite = buildCompositePilotResponse(obj, delay, doppler)
            localCfg = obj.Config;
            pilotResponse = obj.buildPilotResponseVector(delay, doppler, 0);
            hComposite = localCfg.PerPilotAmplitude * localCfg.PilotSequence * pilotResponse;
        end

        % APPLYPATHTOFRAME 将单路径算子作用到 DAFT 域发射帧。
        function y = applyPathToFrame(obj, delay, doppler, txFrame)
            y = obj.buildPathMatrix(delay, doppler) * txFrame;
        end

        % BUILDPATHRESPONSEDICTIONARY 构建多径响应字典矩阵。
        function dict = buildPathResponseDictionary(obj, pathParams, txFrame)
            numSc = obj.Config.NumSubcarriers;
            numPaths = size(pathParams, 1);
            dict = zeros(numSc, numPaths);
            for i = 1:numPaths
                dict(:, i) = obj.applyPathToFrame(pathParams(i,1), pathParams(i,2), txFrame);
            end
        end

        % BUILDDOPPLERBASIS 构造给定整数 Doppler 邻域的基函数集合。
        function [basisVecs, kvRange] = buildDopplerBasis(obj, delay, intDoppler, txFrame)
            numSc    = obj.Config.NumSubcarriers;
            dirichletRadius = obj.Config.DirichletRadius;
            chirpC1  = obj.Config.ChirpParam1;
            chirpC2  = obj.Config.ChirpParam2;

            kvRange = (-dirichletRadius:dirichletRadius).';
            numKv = length(kvRange);

            pVec = (0:numSc-1).';
            alphaInt = intDoppler;
            locStep = round(2 * numSc * chirpC1);
            locIndex = alphaInt + locStep * delay;
            phaseConst = numSc * chirpC1 * delay^2;
            c2PSq = chirpC2 * pVec .* pVec;
            twoPiOverN = 2 * pi / numSc;

            basisVecs = zeros(numSc, numKv);
            for idx = 1:numKv
                qVec = mod(pVec + locIndex + kvRange(idx), numSc);
                phaseArg = twoPiOverN * (phaseConst - qVec * delay ...
                    + numSc * chirpC2 * (qVec .* qVec) - numSc * c2PSq);
                basisVecs(:, idx) = exp(1j * phaseArg) .* txFrame(qVec + 1);
            end
        end

        % COMPUTEDIRICHLETWITHDERIV 计算 Dirichlet 系数及其对分数偏移的一阶导数。
        function [dVal, dDeriv] = computeDirichletWithDeriv(~, fracPart, shiftIndex, numSc)
            [dVal, dDeriv] = FractionalChannelBuilder.dirichletWithDeriv( ...
                fracPart, shiftIndex, numSc);
        end

    end

    methods (Static)

        % COMPUTEDIRICHLETCOEFF 计算分数 Doppler 下的 Dirichlet 系数。
        function dVal = computeDirichletCoeff(fracPart, shiftIndex, numSc)
            xVal = fracPart - shiftIndex;
            if abs(xVal) < 1e-12
                dVal = 1;
            elseif abs(fracPart) < 1e-12
                dVal = 0;
            else
                dVal = (1 / numSc) * (1 - exp(-1j * 2 * pi * fracPart)) / ...
                    (1 - exp(-1j * 2 * pi * xVal / numSc));
            end
        end

        % DIRICHLETWITHDERIV 计算 Dirichlet 系数及其导数（含数值稳定分支）。
        function [dVal, dDeriv] = dirichletWithDeriv(fracPart, shiftIndex, numSc)
            xVal = fracPart - shiftIndex;

            if abs(xVal) < 1e-10
                if abs(fracPart) < 1e-10
                    dVal = 1;
                    dDeriv = 0;
                else
                    dVal = 1;
                    h  = 1e-7;
                    dp = FractionalChannelBuilder.computeDirichletCoeff(fracPart+h, shiftIndex, numSc);
                    dm = FractionalChannelBuilder.computeDirichletCoeff(fracPart-h, shiftIndex, numSc);
                    dDeriv = (dp - dm) / (2*h);
                end
            elseif abs(fracPart) < 1e-10
                dVal = 0;
                expB = exp(-1j * 2 * pi * (-shiftIndex) / numSc);
                dDeriv = (1 / numSc) * (1j * 2 * pi) / (1 - expB);
            else
                expA = exp(-1j * 2 * pi * fracPart);
                expB = exp(-1j * 2 * pi * xVal / numSc);
                numeratorTerm = 1 - expA;
                denominatorTerm = 1 - expB;
                dVal = (1 / numSc) * numeratorTerm / denominatorTerm;
                numeratorDerivative = 1j * 2 * pi * expA;
                denominatorDerivative = (1j * 2 * pi / numSc) * expB;
                dDeriv = (1 / numSc) * (numeratorDerivative * denominatorTerm - numeratorTerm * denominatorDerivative) / (denominatorTerm * denominatorTerm);
            end
        end

    end

end


