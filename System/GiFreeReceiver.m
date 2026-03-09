classdef GiFreeReceiver
    % GiFreeReceiver  最终版 GI-Free AFDM Turbo DD-CE 接收机
    %
    %   三阶段架构: Phase 1 粗调 → Phase 2 DD Turbo → 后判决精炼
    %
    %   =================================================================
    %   Phase 1: 干扰感知粗调 (建立可靠初始估计)
    %   =================================================================
    %     - 先验干扰功率初始正则化: lambda = (sigma_n^2 + sigma_ID2P^2) / SNR
    %     - 正则化 OMP 信道估计
    %     - MMSE 软符号用于 SIC 消除
    %     - 均值残差功率更新正则化 (安全, 适合高 BER 阶段)
    %
    %   =================================================================
    %   Phase 2: DD Turbo 精调 (利用数据提升信道估计)
    %   =================================================================
    %     - 输出 SINR 自适应软符号: 从 LMMSE 输出端估计实际检测噪声,
    %       使软符号收缩程度与检测可靠度精确匹配
    %     - DD 多普勒精炼 (前期) + DD 增益精炼 (后期)
    %     - 截断均值 (Trimmed Mean) 正则化:
    %         去掉最大的 alpha% 残差后取均值, 免疫错误符号离群效应
    %         alpha 随迭代递减: 初期 15% (保守) → 后期 5% (高效)
    %     - 正则化递减调度: 主动向 perfectCsiReg 逼近
    %
    %   =================================================================
    %   后判决精炼 (Post-Decision Refinement)
    %   =================================================================
    %     硬判决符号是最可靠的数据估计, 基于它做最后一轮:
    %       DD 增益精炼 → 紧凑截断均值 (alpha=3%) → LMMSE → 最终解调
    %     本质上是一步 "免费" 的额外迭代, 代价极低但收益显著。

    methods (Static)

        %% ========== 主接收入口 ==========
        function [detectedIndices, normalizedMse, estEffChannel] = receive( ...
                rxSignal, cfg, trueEffChannel, dataSnrLinear, noisePower)

            numSc      = cfg.NumSubcarriers;
            modOrder   = cfg.ModulationOrder;
            totalIter  = cfg.MaxSicIterations;
            pilotAmp   = cfg.PilotAmplitude;

            dataIndices = (2:numSc)';
            numData     = numSc - 1;
            pilotFrame  = zeros(numSc, 1);
            pilotFrame(1) = pilotAmp;

            % 预计算星座图
            qamConstellation    = qammod((0:modOrder-1)', modOrder, 'UnitAveragePower', true);
            scaledConstellation = qamConstellation * sqrt(dataSnrLinear);
            sparseEye           = speye(numData);

            % 理论最优正则化 (Perfect CSI 的值, 我们要逼近的目标)
            perfectCsiReg = noisePower / dataSnrLinear;

            % ================================================================
            %   先验干扰功率计算
            % ================================================================
            dataInterfPower     = (numData / numSc) * dataSnrLinear;
            effectiveNoisePower = noisePower + dataInterfPower;
            regParam            = effectiveNoisePower / dataSnrLinear;
            ompRegParam         = dataInterfPower / (pilotAmp^2);

            % ================================================================
            %   Phase 分配: 30% 粗调 + 60% DD Turbo + 10% 后判决精炼
            % ================================================================
            numPhase1       = max(round(totalIter * 0.30), 2);
            numPostDecision = 1;
            numPhase2       = max(totalIter - numPhase1 - numPostDecision, 3);

            % ================================================================
            %   初始 OMP + LMMSE
            % ================================================================
            [estPaths, hEff] = GiFreeEstimator.estimateByOmp( ...
                rxSignal, cfg, ompRegParam);

            cleanDataSig = rxSignal - hEff * pilotFrame;
            estFullSig   = GiFreeReceiver.lmmseDetect( ...
                cleanDataSig, hEff, regParam, numSc, sparseEye);

            % ================================================================
            %   Phase 1: 干扰感知粗调
            % ================================================================
            for iter = 1:numPhase1
                % MMSE 软符号 (用先验干扰功率作为 sigma^2)
                dataFrame = zeros(numSc, 1);
                rawDataEst  = estFullSig(dataIndices);
                softSymbols = GiFreeReceiver.computeSoftSymbols( ...
                    rawDataEst, scaledConstellation, effectiveNoisePower);
                dataFrame(dataIndices) = softSymbols;

                % ID2P 消除
                cleanPilotSig = rxSignal - hEff * dataFrame;

                % 均值残差 → 更新正则化 (Phase 1 用均值, 安全)
                pilotResidual     = cleanPilotSig - hEff * pilotFrame;
                residualPower     = real(pilotResidual' * pilotResidual) / numSc;
                residualInterf    = max(residualPower - noisePower, 0);
                effectiveNoisePower = noisePower + residualInterf;
                regParam    = max(effectiveNoisePower / dataSnrLinear, perfectCsiReg);
                ompRegParam = max(residualInterf / (pilotAmp^2), 0);

                % OMP 信道重估
                [estPaths, hEff] = GiFreeEstimator.estimateByOmp( ...
                    cleanPilotSig, cfg, ompRegParam);

                % LMMSE 数据检测
                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = GiFreeReceiver.lmmseDetect( ...
                    cleanDataSig, hEff, regParam, numSc, sparseEye);
            end

            % ================================================================
            %   Phase 2: DD Turbo 精调
            % ================================================================
            for iter = 1:numPhase2

                % ---- 输出 SINR 自适应软符号 ----
                rawDataEst    = estFullSig(dataIndices);
                outputNoiseVar = GiFreeReceiver.estimateOutputNoise( ...
                    rawDataEst, scaledConstellation, noisePower);

                softSymbols = GiFreeReceiver.computeSoftSymbols( ...
                    rawDataEst, scaledConstellation, outputNoiseVar);

                % 构造全帧估计
                txFrameEst    = zeros(numSc, 1);
                txFrameEst(1) = pilotAmp;
                txFrameEst(dataIndices) = softSymbols;

                % ---- ID2P 消除 → OMP 路径重检测 ----
                dataFrameDd = zeros(numSc, 1);
                dataFrameDd(dataIndices) = softSymbols;
                cleanPilotSig = rxSignal - hEff * dataFrameDd;

                pilotResidual  = cleanPilotSig - hEff * pilotFrame;
                residualPower  = real(pilotResidual' * pilotResidual) / numSc;
                residualInterf = max(residualPower - noisePower, 0);
                ompRegParam    = max(residualInterf / (pilotAmp^2), 0);

                [estPaths, ~] = GiFreeEstimator.estimateByOmp( ...
                    cleanPilotSig, cfg, ompRegParam);

                % ---- DD 精炼 ----
                txFramePower = real(txFrameEst' * txFrameEst) / numSc;
                ddRegParam   = max(noisePower / txFramePower, 1e-6);

                if iter <= 3
                    % 前期: 多普勒 + 增益精炼
                    [estPaths, hEff] = GiFreeEstimator.refineDopplerByDd( ...
                        estPaths, rxSignal, txFrameEst, cfg, ddRegParam);
                else
                    % 后期: 仅增益精炼 (多普勒已稳定)
                    [estPaths, hEff] = GiFreeEstimator.refineGainByDd( ...
                        estPaths, rxSignal, txFrameEst, cfg, ddRegParam);
                end

                % ---- 截断均值正则化 ----
                %   trimAlpha 随迭代递减: 初期保守, 后期高效
                trimAlpha        = max(0.15 - 0.02 * iter, 0.05);
                trimmedNoisePower = GiFreeEstimator.estimateNoiseTrimmed( ...
                    rxSignal, hEff, txFrameEst, trimAlpha);

                % 正则化: 截断均值估计, 带下限保护
                regParam = max(trimmedNoisePower / dataSnrLinear, perfectCsiReg);

                % 更新 effectiveNoisePower
                effectiveNoisePower = max(trimmedNoisePower, noisePower);

                % ---- LMMSE 数据检测 ----
                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = GiFreeReceiver.lmmseDetect( ...
                    cleanDataSig, hEff, regParam, numSc, sparseEye);
            end

            % ================================================================
            %   后判决精炼 (Post-Decision Refinement)
            % ================================================================
            %   硬判决符号是最可靠的数据估计:
            %     - 在 BER=1% 时, 253/255 个符号完全正确
            %     - 基于这些做 DD 增益重估, 信道估计最精确
            %     - 紧凑截断均值 (alpha=3%) 给出最接近理论值的正则化
            for iter = 1:numPostDecision
                % 硬判决
                normalizedEst = estFullSig(dataIndices) / sqrt(dataSnrLinear);
                hardDecIdx    = qamdemod(normalizedEst, modOrder, 'UnitAveragePower', true);
                hardDecSym    = qammod(hardDecIdx, modOrder, 'UnitAveragePower', true) * sqrt(dataSnrLinear);

                txFrameHard    = zeros(numSc, 1);
                txFrameHard(1) = pilotAmp;
                txFrameHard(dataIndices) = hardDecSym;

                % DD 增益精炼 (基于最可靠的硬判决)
                txFramePower = real(txFrameHard' * txFrameHard) / numSc;
                ddRegParam   = max(noisePower / txFramePower, 1e-6);

                [estPaths, hEff] = GiFreeEstimator.refineGainByDd( ...
                    estPaths, rxSignal, txFrameHard, cfg, ddRegParam);

                % 紧凑截断均值 (alpha=3%, 绝大部分符号正确时非常精确)
                trimmedNoisePower = GiFreeEstimator.estimateNoiseTrimmed( ...
                    rxSignal, hEff, txFrameHard, 0.03);

                regParam = max(trimmedNoisePower / dataSnrLinear, perfectCsiReg);

                % 最终 LMMSE
                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = GiFreeReceiver.lmmseDetect( ...
                    cleanDataSig, hEff, regParam, numSc, sparseEye);
            end

            % ================================================================
            %   最终解调
            % ================================================================
            normalizedDataSym = estFullSig(dataIndices) / sqrt(dataSnrLinear);
            detectedIndices   = qamdemod(normalizedDataSym, modOrder, 'UnitAveragePower', true);
            estEffChannel     = hEff;

            % 归一化 MSE
            normalizedMse = norm(full(hEff) - full(trueEffChannel), 'fro')^2 / ...
                            max(norm(full(trueEffChannel), 'fro')^2, 1e-20);
        end

        %% ========== LMMSE 数据检测 ==========
        function estSig = lmmseDetect(cleanSig, hEff, regParam, numSc, sparseEye)
            % lmmseDetect  线性最小均方误差均衡
            %   只对数据子载波 (2:N) 做均衡, 导频子载波 (1) 置零
            hData   = hEff(:, 2:numSc);
            estData = (hData' * hData + regParam * sparseEye) \ (hData' * cleanSig);
            estSig          = zeros(numSc, 1);
            estSig(2:numSc) = estData;
        end

        %% ========== MMSE 软符号估计 ==========
        function softSym = computeSoftSymbols(rawEst, constellation, noiseVar)
            % computeSoftSymbols  基于后验概率的 MMSE 软符号
            %   E[x|y] = sum_m x_m * P(x_m|y), P(x_m|y) ∝ exp(-|y-x_m|^2/sigma^2)
            %   干扰大时: 向原点收缩 (保守); 干扰小时: 趋近硬判决 (精确)
            numSymbols = length(rawEst);
            softSym    = zeros(numSymbols, 1);
            sigma2     = noiseVar;

            for k = 1:numSymbols
                distances  = abs(rawEst(k) - constellation).^2;
                logWeights = -distances / sigma2;
                logWeights = logWeights - max(logWeights);
                weights    = exp(logWeights);
                weights    = weights / sum(weights);
                softSym(k) = sum(weights .* constellation);
            end
        end

        %% ========== 输出 SINR 噪声估计 ==========
        function outputNoiseVar = estimateOutputNoise(rawEst, scaledConstellation, noisePower)
            % estimateOutputNoise  从 LMMSE 输出与最近星座点的距离估计检测噪声
            %   当大部分符号正确时, |x̂ - QAM_nearest(x̂)|^2 精确反映
            %   LMMSE 输出端的等效噪声方差 (包含热噪声 + 残余干扰)
            numSymbols    = length(rawEst);
            distToNearest = zeros(numSymbols, 1);

            for k = 1:numSymbols
                allDist = abs(rawEst(k) - scaledConstellation).^2;
                distToNearest(k) = min(allDist);
            end

            % 截断均值: 去掉最大 10% 的异常值 (错误判决导致的大距离)
            sortedDist    = sort(distToNearest, 'ascend');
            numKeep       = max(floor(numSymbols * 0.9), 1);
            outputNoiseVar = max(mean(sortedDist(1:numKeep)), noisePower);
        end

    end
end
