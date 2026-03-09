classdef GiFreeEstimator
    % GiFreeEstimator  GI-Free AFDM 信道估计器
    %
    %   整合三大信道估计能力于一体:
    %     (1) 正则化 OMP: 从导频观测中检测路径位置并粗估增益
    %     (2) DD 多普勒精炼: 利用全帧信号精搜分数多普勒偏移
    %     (3) DD 增益重估: 利用全帧信号重估路径复增益
    %   并提供截断均值 (Trimmed Mean) 噪声功率估计工具。
    %
    %   调用链:
    %     Phase 1 使用 estimateByOmp (纯导频 OMP)
    %     Phase 2 使用 refineDopplerByDd + refineGainByDd (全帧 DD)
    %     噪声估计使用 estimateNoiseTrimmed

    methods (Static)

        %% ======================================================================
        %%  正则化 OMP 信道估计
        %% ======================================================================
        function [estPaths, hEffective] = estimateByOmp(rxSignal, cfg, ompRegParam)
            % estimateByOmp  干扰感知正则化 OMP 信道估计
            %   rxSignal    : 接收信号 (或经 SIC 清洁后的导频观测)
            %   cfg         : GiFreeConfig 配置对象
            %   ompRegParam : 增益估计正则化参数 (可选, 默认 = 0)

            if nargin < 3, ompRegParam = 0; end

            numSc          = cfg.NumSubcarriers;
            chirpC1        = cfg.ChirpParam1;
            chirpC2        = cfg.ChirpParam2;
            spreadKv       = cfg.SpreadWidth;
            numPathsUpper  = cfg.NumPathsUpper;
            pilotAmplitude = cfg.PilotAmplitude;

            estPaths        = zeros(numPathsUpper, 3);
            pilotDictionary = zeros(numSc, numPathsUpper);
            delayRecords    = zeros(numPathsUpper, 1);
            dopplerRecords  = zeros(numPathsUpper, 1);
            residual        = rxSignal;

            for pathIdx = 1:numPathsUpper

                % 步骤 1: 粗搜 — 整数格点峰值检测
                peakHit = GiFreeEstimator.findStrongestPeak(residual, cfg);
                delayVal   = peakHit(1);
                intDoppler = peakHit(2);

                % 步骤 2: 精搜 — fminbnd 连续优化分数多普勒
                if spreadKv == 0
                    fracDoppler = intDoppler;
                else
                    fracDoppler = GiFreeEstimator.refineFractionalDoppler( ...
                        residual, delayVal, intDoppler, cfg);
                end

                delayRecords(pathIdx)   = delayVal;
                dopplerRecords(pathIdx) = fracDoppler;

                % 步骤 3: 扩展导频字典
                pilotResponse = FractionalChannelBuilder.buildPilotResponseVector( ...
                    delayVal, fracDoppler, numSc, chirpC1, chirpC2, spreadKv);
                pilotDictionary(:, pathIdx) = pilotResponse * pilotAmplitude;

                % 步骤 4: 正则化联合 LS 增益估计
                activeDict   = pilotDictionary(:, 1:pathIdx);
                gramMatrix   = activeDict' * activeDict;

                if ompRegParam > 0
                    jointGains = (gramMatrix + ompRegParam * eye(pathIdx)) \ ...
                                 (activeDict' * rxSignal);
                else
                    jointGains = activeDict \ rxSignal;
                end

                % 步骤 5: 正交化残差
                residual = rxSignal - activeDict * jointGains;
            end

            % 组装输出
            for pathIdx = 1:numPathsUpper
                estPaths(pathIdx, :) = [delayRecords(pathIdx), ...
                                        dopplerRecords(pathIdx), ...
                                        jointGains(pathIdx)];
            end

            hEffective = FractionalChannelBuilder.buildEffectiveChannel( ...
                estPaths, numSc, chirpC1, chirpC2, spreadKv);
        end

        %% ======================================================================
        %%  DD 增益重估 (固定路径位置, 仅更新复增益)
        %% ======================================================================
        function [refinedPaths, hEffective] = refineGainByDd( ...
                estPaths, rxSignal, txFrameEst, cfg, ddRegParam)
            % refineGainByDd  利用全帧信号重估路径增益
            %   等效导频能量: ||x̂||^2 = Ap^2 + (N-1)*SNR >> Ap^2

            if nargin < 5, ddRegParam = 0; end

            numSc    = cfg.NumSubcarriers;
            chirpC1  = cfg.ChirpParam1;
            chirpC2  = cfg.ChirpParam2;
            spreadKv = cfg.SpreadWidth;
            numPaths = size(estPaths, 1);

            if numPaths == 0
                refinedPaths = estPaths;
                hEffective   = sparse(numSc, numSc);
                return;
            end

            % 构造全帧响应字典: fullDict(:,i) = H_i * x̂
            fullDict = zeros(numSc, numPaths);
            for i = 1:numPaths
                hSingle = FractionalChannelBuilder.buildPathMatrix( ...
                    estPaths(i, 1), estPaths(i, 2), numSc, chirpC1, chirpC2, spreadKv);
                fullDict(:, i) = hSingle * txFrameEst;
            end

            % 正则化 LS 增益重估
            gramFull = fullDict' * fullDict;
            if ddRegParam > 0
                newGains = (gramFull + ddRegParam * eye(numPaths)) \ (fullDict' * rxSignal);
            else
                newGains = fullDict \ rxSignal;
            end

            refinedPaths      = estPaths;
            refinedPaths(:,3) = newGains;
            hEffective = FractionalChannelBuilder.buildEffectiveChannel( ...
                refinedPaths, numSc, chirpC1, chirpC2, spreadKv);
        end

        %% ======================================================================
        %%  DD 多普勒精炼 + 增益重估
        %% ======================================================================
        function [refinedPaths, hEffective] = refineDopplerByDd( ...
                estPaths, rxSignal, txFrameEst, cfg, ddRegParam)
            % refineDopplerByDd  坐标下降精搜分数多普勒, 然后联合重估增益
            %   全帧匹配滤波的等效 SNR 远高于单导频, 多普勒精度大幅提升

            if nargin < 5, ddRegParam = 0; end

            numSc    = cfg.NumSubcarriers;
            chirpC1  = cfg.ChirpParam1;
            chirpC2  = cfg.ChirpParam2;
            spreadKv = cfg.SpreadWidth;
            numPaths = size(estPaths, 1);

            if numPaths == 0 || spreadKv == 0
                refinedPaths = estPaths;
                hEffective   = FractionalChannelBuilder.buildEffectiveChannel( ...
                    estPaths, numSc, chirpC1, chirpC2, spreadKv);
                return;
            end

            refinedPaths = estPaths;

            % 逐路径多普勒精炼 (坐标下降)
            for i = 1:numPaths
                % 计算排除第 i 条路径后的残差
                residualExcl = rxSignal;
                for j = 1:numPaths
                    if j == i, continue; end
                    hSingle = FractionalChannelBuilder.buildPathMatrix( ...
                        refinedPaths(j, 1), refinedPaths(j, 2), ...
                        numSc, chirpC1, chirpC2, spreadKv);
                    residualExcl = residualExcl - refinedPaths(j, 3) * hSingle * txFrameEst;
                end

                % 全帧匹配滤波搜索最优多普勒
                intDoppler = round(refinedPaths(i, 2));
                delayVal   = refinedPaths(i, 1);

                objFcn = @(nu) -GiFreeEstimator.computeFullFrameMetric( ...
                    delayVal, nu, residualExcl, txFrameEst, ...
                    numSc, chirpC1, chirpC2, spreadKv);

                optOpts = optimset('TolX', 1e-6, 'Display', 'off');
                bestNu  = fminbnd(objFcn, intDoppler - 0.499, intDoppler + 0.499, optOpts);
                refinedPaths(i, 2) = bestNu;
            end

            % 联合重估所有增益
            fullDict = zeros(numSc, numPaths);
            for i = 1:numPaths
                hSingle = FractionalChannelBuilder.buildPathMatrix( ...
                    refinedPaths(i, 1), refinedPaths(i, 2), ...
                    numSc, chirpC1, chirpC2, spreadKv);
                fullDict(:, i) = hSingle * txFrameEst;
            end

            gramFull = fullDict' * fullDict;
            if ddRegParam > 0
                newGains = (gramFull + ddRegParam * eye(numPaths)) \ (fullDict' * rxSignal);
            else
                newGains = fullDict \ rxSignal;
            end
            refinedPaths(:, 3) = newGains;

            hEffective = FractionalChannelBuilder.buildEffectiveChannel( ...
                refinedPaths, numSc, chirpC1, chirpC2, spreadKv);
        end

        %% ======================================================================
        %%  截断均值噪声功率估计
        %% ======================================================================
        function noisePowerEst = estimateNoiseTrimmed(rxSignal, hEff, txFrameEst, trimAlpha)
            % estimateNoiseTrimmed  使用截断均值估计等效噪声功率
            %   去掉最大的 trimAlpha 比例后取均值, 免疫错误符号的离群效应.
            %   比中位数高效 (利用更多正确样本), 仍然鲁棒.
            %
            %   trimAlpha : 截断比例, 如 0.1 表示去掉最大的 10%
            %   返回值    : 每元素平均等效噪声功率

            if nargin < 4, trimAlpha = 0.1; end

            residualPowers = abs(rxSignal - hEff * txFrameEst).^2;
            sortedPowers   = sort(residualPowers, 'ascend');
            numTotal       = length(sortedPowers);

            % 截断: 保留最小的 (1 - trimAlpha) 比例
            numKeep        = max(floor(numTotal * (1 - trimAlpha)), 1);
            noisePowerEst  = mean(sortedPowers(1:numKeep));
        end

        %% ======================================================================
        %%  内部辅助: 整数格点峰值检测
        %% ======================================================================
        function peakHit = findStrongestPeak(residualSignal, cfg)
            numSc          = cfg.NumSubcarriers;
            maxDelaySpread = cfg.MaxDelaySpread;
            maxDopplerIdx  = cfg.MaxDopplerIndex;
            locStep        = cfg.LocStep;

            maxMagnitude = -1;
            bestDelay    = 0;
            bestDoppler  = 0;

            for delayIdx = 0:maxDelaySpread
                for dopplerIdx = -maxDopplerIdx:maxDopplerIdx
                    locIndex  = dopplerIdx + locStep * delayIdx;
                    centerBin = mod(-locIndex, numSc);
                    gridVal   = abs(residualSignal(centerBin + 1));

                    if gridVal > maxMagnitude
                        maxMagnitude = gridVal;
                        bestDelay    = delayIdx;
                        bestDoppler  = dopplerIdx;
                    end
                end
            end

            peakHit = [bestDelay, bestDoppler];
        end

        %% ======================================================================
        %%  内部辅助: 分数多普勒连续精搜 (基于导频)
        %% ======================================================================
        function estDoppler = refineFractionalDoppler(residualSignal, delayVal, intDoppler, cfg)
            numSc    = cfg.NumSubcarriers;
            chirpC1  = cfg.ChirpParam1;
            chirpC2  = cfg.ChirpParam2;
            spreadKv = cfg.SpreadWidth;

            objFcn = @(fracShift) -GiFreeEstimator.computePilotMetric( ...
                fracShift, residualSignal, delayVal, intDoppler, ...
                numSc, chirpC1, chirpC2, spreadKv);

            optOpts     = optimset('TolX', 1e-4, 'Display', 'off');
            fracShiftEst = fminbnd(objFcn, -0.499, 0.499, optOpts);
            estDoppler   = intDoppler + fracShiftEst;
        end

        %% ======================================================================
        %%  内部辅助: 导频匹配滤波度量
        %% ======================================================================
        function metricVal = computePilotMetric( ...
                fracShift, signalVec, delayVal, intDoppler, numSc, chirpC1, chirpC2, spreadKv)
            pilotResp = FractionalChannelBuilder.buildPilotResponseVector( ...
                delayVal, intDoppler + fracShift, numSc, chirpC1, chirpC2, spreadKv);
            normSq    = real(pilotResp' * pilotResp);
            if normSq < 1e-15
                metricVal = 0;
            else
                metricVal = abs(pilotResp' * signalVec)^2 / normSq;
            end
        end

        %% ======================================================================
        %%  内部辅助: 全帧匹配滤波度量 (DD 多普勒搜索用)
        %% ======================================================================
        function metricVal = computeFullFrameMetric( ...
                delayVal, doppler, signalVec, txFrameEst, numSc, chirpC1, chirpC2, spreadKv)
            hSingle  = FractionalChannelBuilder.buildPathMatrix( ...
                delayVal, doppler, numSc, chirpC1, chirpC2, spreadKv);
            response = hSingle * txFrameEst;
            normSq   = real(response' * response);
            if normSq < 1e-15
                metricVal = 0;
            else
                metricVal = abs(response' * signalVec)^2 / normSq;
            end
        end

    end
end
