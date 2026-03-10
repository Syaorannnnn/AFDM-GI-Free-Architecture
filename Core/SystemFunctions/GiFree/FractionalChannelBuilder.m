classdef FractionalChannelBuilder < handle
    % FractionalChannelBuilder  分数多普勒信道矩阵构造工具 (CAZAC 多导频扩展版)
    %
    %   新增/修改:
    %     buildPilotResponseVector : 接受可选 pilotIdx 参数 (默认 0, 向后兼容)
    %     buildCompositePilotResponse : K 个导频的相干合成响应向量

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

        %% ---------- 导频响应向量 (支持任意导频位置) ----------
        function hCol = buildPilotResponseVector(obj, delay, doppler, pilotIdx)
            % buildPilotResponseVector  计算单径 H_eff 的第 pilotIdx 列
            %
            %   pilotIdx : 0-based DAFT 域导频索引 (默认 0, 向后兼容)
            %
            %   数学推导:
            %     H_path(p, q) 中 q = mod(p + loc + kv, N), 其相位为
            %       exp(j*2π*(c1*l² - q*l/N + c2*q² - c2*p²))
            %     对第 pilotIdx 列 (q = pilotIdx), 对应行为
            %       p = mod(pilotIdx - loc - kv, N)
            %     相位变为
            %       exp(j*2π*(c1*l² - pilotIdx*l/N + c2*pilotIdx² - c2*p²))

            if nargin < 4, pilotIdx = 0; end

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
                pIdx     = mod(pilotIdx - locIndex - kvShift, numSc);
                phaseVal = exp(1j * 2 * pi * (chirpC1 * delay^2 ...
                    - pilotIdx * delay / numSc ...
                    + chirpC2 * pilotIdx^2 ...
                    - chirpC2 * pIdx^2));
                dCoeff   = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvShift, numSc);
                hCol(pIdx + 1) = phaseVal * dCoeff;
            end
        end

        %% ---------- CAZAC 多导频合成响应 ----------
        function hComposite = buildCompositePilotResponse(obj, delay, doppler)
            % buildCompositePilotResponse  K 个导频的相干合成响应向量
            %
            %   y_pilot = sum_k (A_p/sqrt(K)) * p_k * H_path(:, m_k)
            %
            %   K = 1 时等价于 PilotAmplitude * buildPilotResponseVector(delay, doppler, 0)

            cfg_     = obj.Config;
            K        = cfg_.NumPilots;
            positions = cfg_.PilotPositions;
            sequence = cfg_.PilotSequence;
            ampPP    = cfg_.PerPilotAmplitude;

            hComposite = zeros(cfg_.NumSubcarriers, 1);
            for k = 1:K
                hk = obj.buildPilotResponseVector(delay, doppler, positions(k));
                hComposite = hComposite + ampPP * sequence(k) * hk;
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
