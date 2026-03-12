classdef FractionalChannelBuilder < handle
    % FractionalChannelBuilder  分数多普勒信道矩阵构造工具
    %
    %   向量化实现: buildPathMatrix 按 kvShift 分块向量化,
    %   消除 N × numKvPts 的标量循环, 改为 numKvPts 次向量运算.

    properties (SetAccess = private)
        Config
    end

    methods

        function obj = FractionalChannelBuilder(cfg)
            arguments
                cfg (1,1) GiFreeConfig
            end
            obj.Config = cfg;
        end

        %% ---------- 单径信道矩阵 (向量化) ----------
        function hPath = buildPathMatrix(obj, delay, doppler)
            numSc    = obj.Config.NumSubcarriers;
            chirpC1  = obj.Config.ChirpParam1;
            chirpC2  = obj.Config.ChirpParam2;
            spreadKv = obj.Config.SpreadWidth;

            alphaInt = round(doppler);
            fracPart = doppler - alphaInt;
            locStep  = round(2 * numSc * chirpC1);
            locIndex = alphaInt + locStep * delay;

            kvRange  = -spreadKv:spreadKv;
            numKvPts = length(kvRange);

            % 预计算 Dirichlet 系数 (仅 numKvPts 个标量)
            dCoeffs = zeros(numKvPts, 1);
            for idx = 1:numKvPts
                dCoeffs(idx) = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvRange(idx), numSc);
            end

            pVec       = (0:numSc-1)';
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
                if abs(dc) < 1e-15, continue; end

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

        %% ---------- 多径有效信道 ----------
        function hEff = buildEffectiveChannel(obj, pathParams)
            numSc = obj.Config.NumSubcarriers;
            hEff  = sparse(numSc, numSc);
            for i = 1:size(pathParams, 1)
                hEff = hEff + pathParams(i,3) * ...
                    obj.buildPathMatrix(pathParams(i,1), pathParams(i,2));
            end
        end

        %% ---------- 导频响应向量 ----------
        function hCol = buildPilotResponseVector(obj, delay, doppler)
            numSc    = obj.Config.NumSubcarriers;
            chirpC1  = obj.Config.ChirpParam1;
            chirpC2  = obj.Config.ChirpParam2;
            spreadKv = obj.Config.SpreadWidth;

            alphaInt = round(doppler);
            fracPart = doppler - alphaInt;
            locStep  = round(2 * numSc * chirpC1);
            locIndex = alphaInt + locStep * delay;

            hCol = zeros(numSc, 1);
            for kvShift = -spreadKv:spreadKv
                pIdx     = mod(-locIndex - kvShift, numSc);
                phaseVal = exp(1j * 2 * pi * (chirpC1 * delay^2 - chirpC2 * pIdx^2));
                dCoeff   = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvShift, numSc);
                hCol(pIdx + 1) = phaseVal * dCoeff;
            end
        end

        %% ---------- 随机信道实现 ----------
        function [pathParams, hEff] = generateChannel(obj)
            cfg_      = obj.Config;
            numPaths  = cfg_.NumPaths;
            maxDelay  = cfg_.MaxDelaySpread;
            maxDopIdx = cfg_.MaxDopplerIndex;

            allDelays      = 0:maxDelay;
            selectedDelays = allDelays(randperm(length(allDelays), numPaths))';
            randomPhases   = -pi + 2 * pi * rand(numPaths, 1);
            dopplerShifts  = maxDopIdx * cos(randomPhases);

            if ~cfg_.UseFractionalDoppler
                dopplerShifts = round(dopplerShifts);
            end

            for p = 2:numPaths
                attempts = 0;
                while any(selectedDelays(1:p-1) == selectedDelays(p) & ...
                        round(dopplerShifts(1:p-1)) == round(dopplerShifts(p)))
                    randomPhases(p)  = -pi + 2 * pi * rand();
                    dopplerShifts(p) = maxDopIdx * cos(randomPhases(p));
                    if ~cfg_.UseFractionalDoppler
                        dopplerShifts(p) = round(dopplerShifts(p));
                    end
                    attempts = attempts + 1;
                    if attempts > 100, break; end
                end
            end

            complexGains = sqrt(1 / (2 * numPaths)) * ...
                (randn(numPaths, 1) + 1j * randn(numPaths, 1));
            pathParams = [selectedDelays, dopplerShifts, complexGains];
            hEff       = obj.buildEffectiveChannel(pathParams);
        end

    end

    methods (Static)

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

        function [dVal, dDeriv] = dirichletWithDeriv(fracPart, shiftIndex, numSc)
            % dirichletWithDeriv  Dirichlet 核系数及其对 fracPart 的解析导数
            xVal = fracPart - shiftIndex;

            if abs(xVal) < 1e-10
                if abs(fracPart) < 1e-10
                    dVal = 1; dDeriv = 0;
                else
                    dVal = 1;
                    h  = 1e-7;
                    dp = FractionalChannelBuilder.computeDirichletCoeff(fracPart+h, shiftIndex, numSc);
                    dm = FractionalChannelBuilder.computeDirichletCoeff(fracPart-h, shiftIndex, numSc);
                    dDeriv = (dp - dm) / (2*h);
                end
            elseif abs(fracPart) < 1e-10
                dVal   = 0;
                expB   = exp(-1j * 2 * pi * (-shiftIndex) / numSc);
                dDeriv = (1 / numSc) * (1j * 2 * pi) / (1 - expB);
            else
                expA = exp(-1j * 2 * pi * fracPart);
                expB = exp(-1j * 2 * pi * xVal / numSc);
                A = 1 - expA;  B = 1 - expB;
                dVal   = (1 / numSc) * A / B;
                Ap     = 1j * 2 * pi * expA;
                Bp     = (1j * 2 * pi / numSc) * expB;
                dDeriv = (1 / numSc) * (Ap * B - A * Bp) / (B * B);
            end
        end

    end

end
