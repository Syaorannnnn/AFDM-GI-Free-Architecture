classdef GiFreeReceiver < handle
    % GiFreeReceiver  GI-Free AFDM Turbo DD-CE 接收机
    %
    %   三阶段架构: Phase 1 粗调 → Phase 2 DD Turbo → 后判决精炼
    %   Phase 1 使用直接 LMMSE, Phase 2 + 后判决使用 PCG 稀疏均衡.

    properties (SetAccess = private)
        Config
        Estimator
    end

    methods

        function obj = GiFreeReceiver(cfg, estimator)
            arguments
                cfg       (1,1) GiFreeConfig
                estimator (1,1) GiFreeEstimator
            end
            obj.Config    = cfg;
            obj.Estimator = estimator;
        end

        function [detectedIndices, normalizedMse, estEffChannel] = receive( ...
                obj, rxSignal, trueEffChannel, dataSnrLinear, noisePower)

            numSc      = obj.Config.NumSubcarriers;
            modOrder   = obj.Config.ModulationOrder;
            totalIter  = obj.Config.MaxSicIterations;
            pilotAmp   = obj.Config.PilotAmplitude;

            dataIndices = (2:numSc)';
            numData     = numSc - 1;
            pilotFrame  = zeros(numSc, 1);
            pilotFrame(1) = pilotAmp;

            qamConstellation    = qammod((0:modOrder-1)', modOrder, 'UnitAveragePower', true);
            scaledConstellation = qamConstellation * sqrt(dataSnrLinear);
            sparseEye           = speye(numData);

            perfectCsiReg       = noisePower / dataSnrLinear;
            dataInterfPower     = (numData / numSc) * dataSnrLinear;
            effectiveNoisePower = noisePower + dataInterfPower;
            regParam            = effectiveNoisePower / dataSnrLinear;
            ompRegParam         = dataInterfPower / (pilotAmp^2);

            % Phase 分配
            numPhase1       = max(round(totalIter * 0.30), 2);
            numPostDecision = 1;
            numPhase2       = max(totalIter - numPhase1 - numPostDecision, 3);

            % 初始 OMP + LMMSE
            [estPaths, hEff] = obj.Estimator.estimateByOmp(rxSignal, ompRegParam);
            cleanDataSig = rxSignal - hEff * pilotFrame;
            estFullSig   = GiFreeReceiver.lmmseDetect(cleanDataSig, hEff, regParam, numSc, sparseEye);

            % ============ Phase 1: 干扰感知粗调 (LMMSE) ============
            for iter = 1:numPhase1
                rawDataEst  = estFullSig(dataIndices);
                softSymbols = GiFreeReceiver.computeSoftSymbols( ...
                    rawDataEst, scaledConstellation, effectiveNoisePower);

                dataFrame = zeros(numSc, 1);
                dataFrame(dataIndices) = softSymbols;
                cleanPilotSig = rxSignal - hEff * dataFrame;

                pilotResidual       = cleanPilotSig - hEff * pilotFrame;
                residualPower       = real(pilotResidual' * pilotResidual) / numSc;
                residualInterf      = max(residualPower - noisePower, 0);
                effectiveNoisePower = noisePower + residualInterf;
                regParam    = max(effectiveNoisePower / dataSnrLinear, perfectCsiReg);
                ompRegParam = max(residualInterf / (pilotAmp^2), 0);

                [estPaths, hEff] = obj.Estimator.estimateByOmp(cleanPilotSig, ompRegParam);
                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = GiFreeReceiver.lmmseDetect(cleanDataSig, hEff, regParam, numSc, sparseEye);
            end

            % ============ Phase 2: DD Turbo 精调 (PCG) ============
            for iter = 1:numPhase2
                rawDataEst     = estFullSig(dataIndices);
                outputNoiseVar = GiFreeReceiver.estimateOutputNoise( ...
                    rawDataEst, scaledConstellation, noisePower);
                softSymbols    = GiFreeReceiver.computeSoftSymbols( ...
                    rawDataEst, scaledConstellation, outputNoiseVar);

                txFrameEst = zeros(numSc, 1);
                txFrameEst(1) = pilotAmp;
                txFrameEst(dataIndices) = softSymbols;

                % ID2P → OMP 路径重检测
                dataFrameDd = zeros(numSc, 1);
                dataFrameDd(dataIndices) = softSymbols;
                cleanPilotSig = rxSignal - hEff * dataFrameDd;

                pilotResidual  = cleanPilotSig - hEff * pilotFrame;
                residualPower  = real(pilotResidual' * pilotResidual) / numSc;
                residualInterf = max(residualPower - noisePower, 0);
                ompRegParam    = max(residualInterf / (pilotAmp^2), 0);
                [estPaths, ~]  = obj.Estimator.estimateByOmp(cleanPilotSig, ompRegParam);

                % DD 精炼
                txFramePower = real(txFrameEst' * txFrameEst) / numSc;
                ddRegParam   = max(noisePower / txFramePower, 1e-6);

                if iter <= 3
                    [estPaths, hEff] = obj.Estimator.refineDopplerByDd( ...
                        estPaths, rxSignal, txFrameEst, ddRegParam);
                else
                    [estPaths, hEff] = obj.Estimator.refineGainByDd( ...
                        estPaths, rxSignal, txFrameEst, ddRegParam);
                end

                trimAlpha         = max(0.15 - 0.02 * iter, 0.05);
                trimmedNoisePower = obj.Estimator.estimateNoiseTrimmed( ...
                    rxSignal, hEff, txFrameEst, trimAlpha);
                regParam = max(trimmedNoisePower / dataSnrLinear, perfectCsiReg);

                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = SparseLmmseSolver.solve(cleanDataSig, hEff, regParam, numSc, ...
                    'MaxIter', 20, 'Tolerance', 1e-6, 'WarmStart', estFullSig);
            end

            % ============ 后判决精炼 (PCG) ============
            for iter = 1:numPostDecision
                normalizedEst = estFullSig(dataIndices) / sqrt(dataSnrLinear);
                hardDecIdx    = qamdemod(normalizedEst, modOrder, 'UnitAveragePower', true);
                hardDecSym    = qammod(hardDecIdx, modOrder, 'UnitAveragePower', true) * sqrt(dataSnrLinear);

                txFrameHard    = zeros(numSc, 1);
                txFrameHard(1) = pilotAmp;
                txFrameHard(dataIndices) = hardDecSym;

                txFramePower = real(txFrameHard' * txFrameHard) / numSc;
                ddRegParam   = max(noisePower / txFramePower, 1e-6);
                [estPaths, hEff] = obj.Estimator.refineGainByDd( ...
                    estPaths, rxSignal, txFrameHard, ddRegParam);

                trimmedNoisePower = obj.Estimator.estimateNoiseTrimmed( ...
                    rxSignal, hEff, txFrameHard, 0.03);
                regParam = max(trimmedNoisePower / dataSnrLinear, perfectCsiReg);

                cleanDataSig = rxSignal - hEff * pilotFrame;
                estFullSig   = SparseLmmseSolver.solve(cleanDataSig, hEff, regParam, numSc, ...
                    'MaxIter', 15, 'Tolerance', 1e-6, 'WarmStart', estFullSig);
            end

            % 最终解调
            normalizedDataSym = estFullSig(dataIndices) / sqrt(dataSnrLinear);
            detectedIndices   = qamdemod(normalizedDataSym, modOrder, 'UnitAveragePower', true);
            estEffChannel     = hEff;
            normalizedMse     = norm(full(hEff) - full(trueEffChannel), 'fro')^2 / ...
                max(norm(full(trueEffChannel), 'fro')^2, 1e-20);
        end

    end

    methods (Static)

        function estSig = lmmseDetect(cleanSig, hEff, regParam, numSc, sparseEye)
            hData   = hEff(:, 2:numSc);
            estData = (hData' * hData + regParam * sparseEye) \ (hData' * cleanSig);
            estSig          = zeros(numSc, 1);
            estSig(2:numSc) = estData;
        end

        function softSym = computeSoftSymbols(rawEst, constellation, noiseVar)
            numSymbols = length(rawEst);
            softSym    = zeros(numSymbols, 1);
            for k = 1:numSymbols
                distances  = abs(rawEst(k) - constellation).^2;
                logWeights = -distances / noiseVar;
                logWeights = logWeights - max(logWeights);
                weights    = exp(logWeights);
                weights    = weights / sum(weights);
                softSym(k) = sum(weights .* constellation);
            end
        end

        function outputNoiseVar = estimateOutputNoise(rawEst, scaledConstellation, noisePower)
            numSymbols    = length(rawEst);
            distToNearest = zeros(numSymbols, 1);
            for k = 1:numSymbols
                distToNearest(k) = min(abs(rawEst(k) - scaledConstellation).^2);
            end
            sortedDist     = sort(distToNearest, 'ascend');
            numKeep        = max(floor(numSymbols * 0.9), 1);
            outputNoiseVar = max(mean(sortedDist(1:numKeep)), noisePower);
        end

    end

end
