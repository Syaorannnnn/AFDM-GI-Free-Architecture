classdef FractionalChannelBuilder
    % FractionalChannelBuilder  分数多普勒信道矩阵构造工具
    %   提供 Dirichlet 系数计算、单径矩阵构造、有效信道合成、
    %   导频响应向量生成，以及随机信道实现生成等静态方法。

    methods (Static)

        %% ---------- Dirichlet 插值系数 ----------
        function dVal = computeDirichletCoeff(fracPart, shiftIndex, numSc)
            % computeDirichletCoeff  计算分数多普勒的 Dirichlet 核系数
            %   fracPart   : 分数多普勒偏移量 (|fracPart| < 0.5)
            %   shiftIndex : 相对于中心格点的偏移索引
            %   numSc      : 子载波总数 N
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

        %% ---------- 单径信道矩阵 ----------
        function hPath = buildPathMatrix(delay, doppler, numSc, chirpC1, chirpC2, spreadKv)
            % buildPathMatrix  为单条路径 (delay, doppler) 构造 N x N 稀疏信道矩阵
            alphaInt   = round(doppler);
            fracPart   = doppler - alphaInt;
            locStep    = round(2 * numSc * chirpC1);
            locIndex   = alphaInt + locStep * delay;

            kvRange    = -spreadKv:spreadKv;
            numKvPts   = length(kvRange);

            % 预计算所有 Dirichlet 系数
            dCoeffs = zeros(numKvPts, 1);
            for idx = 1:numKvPts
                dCoeffs(idx) = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvRange(idx), numSc);
            end

            % 稀疏矩阵三元组预分配
            rowIdx = zeros(numSc * numKvPts, 1);
            colIdx = zeros(numSc * numKvPts, 1);
            vals   = complex(zeros(numSc * numKvPts, 1));
            ptr    = 0;

            for p = 0:(numSc - 1)
                for idx = 1:numKvPts
                    if abs(dCoeffs(idx)) < 1e-15, continue; end
                    kvShift = kvRange(idx);
                    qIdx = mod(p + locIndex + kvShift, numSc);
                    phaseVal = exp(1j * 2 * pi / numSc * ...
                        (numSc * chirpC1 * delay^2 - qIdx * delay + ...
                         numSc * chirpC2 * (qIdx^2 - p^2)));
                    ptr = ptr + 1;
                    rowIdx(ptr) = p + 1;
                    colIdx(ptr) = qIdx + 1;
                    vals(ptr)   = phaseVal * dCoeffs(idx);
                end
            end

            hPath = sparse(rowIdx(1:ptr), colIdx(1:ptr), vals(1:ptr), numSc, numSc);
        end

        %% ---------- 多径有效信道 ----------
        function hEff = buildEffectiveChannel(pathParams, numSc, chirpC1, chirpC2, spreadKv)
            % buildEffectiveChannel  叠加所有路径，构造总有效信道矩阵
            %   pathParams : P x 3 矩阵, 每行 [delay, doppler, complexGain]
            hEff = sparse(numSc, numSc);
            if isempty(pathParams), return; end
            for i = 1:size(pathParams, 1)
                hSingle = FractionalChannelBuilder.buildPathMatrix( ...
                    pathParams(i, 1), pathParams(i, 2), numSc, chirpC1, chirpC2, spreadKv);
                hEff = hEff + pathParams(i, 3) * hSingle;
            end
        end

        %% ---------- 导频响应向量 ----------
        function hCol = buildPilotResponseVector(delay, doppler, numSc, chirpC1, chirpC2, spreadKv)
            % buildPilotResponseVector  构造导频子载波 (q=0) 对应的信道列向量
            alphaInt   = round(doppler);
            fracPart   = doppler - alphaInt;
            locStep    = round(2 * numSc * chirpC1);
            locIndex   = alphaInt + locStep * delay;

            hCol = zeros(numSc, 1);
            for kvShift = -spreadKv:spreadKv
                pIdx = mod(-locIndex - kvShift, numSc);
                phaseVal = exp(1j * 2 * pi * (chirpC1 * delay^2 - chirpC2 * pIdx^2));
                dCoeff   = FractionalChannelBuilder.computeDirichletCoeff( ...
                    fracPart, kvShift, numSc);
                hCol(pIdx + 1) = phaseVal * dCoeff;
            end
        end

        %% ---------- 随机信道实现生成 ----------
        function [pathParams, hEff] = generateChannel(cfg)
            % generateChannel  按配置随机生成一次信道实现
            numPaths       = cfg.NumPaths;
            maxDelay       = cfg.MaxDelaySpread;
            maxDopplerIdx  = cfg.MaxDopplerIndex;
            numSc          = cfg.NumSubcarriers;
            chirpC1        = cfg.ChirpParam1;
            chirpC2        = cfg.ChirpParam2;
            spreadKv       = cfg.SpreadWidth;

            % 从 0:maxDelay 中不重复选取时延
            allDelays      = 0:maxDelay;
            selectedDelays = allDelays(randperm(length(allDelays), numPaths))';

            % 随机多普勒: 余弦模型 nu_i = k_max * cos(theta_i)
            randomPhases   = -pi + 2 * pi * rand(numPaths, 1);
            dopplerShifts  = maxDopplerIdx * cos(randomPhases);

            if ~cfg.UseFractionalDoppler
                dopplerShifts = round(dopplerShifts);
            end

            % 确保同一 (delay, round(doppler)) 格点不出现碰撞
            for pathIdx = 2:numPaths
                attempts = 0;
                while any(selectedDelays(1:pathIdx-1) == selectedDelays(pathIdx) & ...
                          round(dopplerShifts(1:pathIdx-1)) == round(dopplerShifts(pathIdx)))
                    randomPhases(pathIdx)  = -pi + 2 * pi * rand();
                    dopplerShifts(pathIdx) = maxDopplerIdx * cos(randomPhases(pathIdx));
                    if ~cfg.UseFractionalDoppler
                        dopplerShifts(pathIdx) = round(dopplerShifts(pathIdx));
                    end
                    attempts = attempts + 1;
                    if attempts > 100, break; end
                end
            end

            % 等功率 Rayleigh 增益
            complexGains = sqrt(1 / (2 * numPaths)) * ...
                (randn(numPaths, 1) + 1j * randn(numPaths, 1));
            pathParams   = [selectedDelays, dopplerShifts, complexGains];
            hEff         = FractionalChannelBuilder.buildEffectiveChannel( ...
                pathParams, numSc, chirpC1, chirpC2, spreadKv);
        end

    end
end
