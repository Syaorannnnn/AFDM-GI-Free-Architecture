classdef EpReceiver < handle
% EpReceiver: Embedded Pilot 接收端，含信道估计与数据检测。
    %
    %

    properties (Access = public)
        CsiMode
        EqualizerType
    end

    properties (Access = private)
        Config
        CurrentPilotPower
    end

    methods (Access = public)

        % EpReceiver: 函数实现见下方代码。
        function obj = EpReceiver(configObj)
            obj.Config = configObj;
        end

        % receive: 执行接收链路并输出检测结果。
        function rxData = receive(obj, rxSignal, noisePower, physicalChannelMatrix, pilotPower)
            obj.CurrentPilotPower = pilotPower;

            rxNoPrefix = rxSignal(obj.Config.PrefixLength + 1:end);

            % DAFT 鍙樻崲
            daftRx = AfdmTransforms.daft(rxNoPrefix, ...
                obj.Config.ChirpParam1, obj.Config.ChirpParam2);

            % 鑾峰彇绛夋晥淇￠亾鐭╅樀
            if obj.CsiMode == "Perfect"
                effectiveChannel = AfdmTransforms.generateEffectiveChannelMatrix( ...
                    physicalChannelMatrix, obj.Config);
            else
                [effectiveChannel, ~] = obj.runChannelEstimator(daftRx);
            end

            % ======================================================
            %
            % ======================================================
            pilotFrame = zeros(obj.Config.NumDataSubcarriers, 1);
            pilotFrame(obj.Config.PilotPos1) = sqrt(obj.CurrentPilotPower);
            pilotContribution = effectiveChannel * pilotFrame;
            cleanDataSignal = daftRx - pilotContribution;

            % MMSE 鍧囪　 (浠呬綔鐢ㄤ簬鏁版嵁瀛愯浇娉㈠垪)
            hData = effectiveChannel(:, obj.Config.DataPos1);
            estData = (hData' * hData + noisePower * eye(obj.Config.NumActiveCarriers)) ...
                \ (hData' * cleanDataSignal);

            % 硬判决
            rxData = qamdemod(estData, obj.Config.ModulationOrder, 'UnitAveragePower', true);
            rxData = rxData(:);
        end

    end

    methods (Access = private)

        % ================================================================
        % ================================================================
        % runChannelEstimator: 基于导频估计有效信道与路径参数。
        function [effectiveChannel, finalPaths] = runChannelEstimator(obj, rxSignalDaft)
            N = obj.Config.NumDataSubcarriers;
            c1 = obj.Config.ChirpParam1;
            c2 = obj.Config.ChirpParam2;

            mockTxFrame = zeros(N, 1);
            mockTxFrame(obj.Config.PilotPos1) = sqrt(obj.CurrentPilotPower);

            timePilot = AfdmTransforms.idaft(mockTxFrame, c1, c2);

            residual = rxSignalDaft;
            estPaths = zeros(obj.Config.NumPaths, 3);

            delayRange = 0:obj.Config.MaxDelaySamples;
            dopplerRange = -obj.Config.MaxDopplerIdx:obj.Config.MaxDopplerIdx;

            timeVec = (0:N - 1).';
            normDopplerVec = dopplerRange(:).' / N;
            dopplerPhaseMat = exp(-1j * 2 * pi * timeVec .* normDopplerVec);

            % 迭代处理：按当前索引更新状态。
            for pathIdx = 1:obj.Config.NumPaths
                bestCorr = -inf;
                bestDelay = 0;
                bestDopIdx = 1;

                for dIdx = 1:length(delayRange)
                    curDelay = delayRange(dIdx);
                    shifted = circshift(timePilot, curDelay);
                    candidates = shifted .* dopplerPhaseMat;  % 大小: N x numDoppler

                    daftCand = AfdmTransforms.daft(candidates, c1, c2);

                    for col = 1:size(daftCand, 2)
                        daftCand(:, col) = obj.applyPeakWindow(daftCand(:, col), 3);
                    end

                    energies = sum(abs(daftCand) .^ 2, 1);
                    corrs = abs(residual' * daftCand) ./ sqrt(energies + 1e-12);
                    [maxVal, maxIdx] = max(corrs);

                    if maxVal > bestCorr
                        bestCorr = maxVal;
                        bestDelay = curDelay;
                        bestDopIdx = maxIdx;
                    end
                end

                coarseDoppler = dopplerRange(bestDopIdx);

                costFn = @(frac) -obj.calcCorrelation( ...
                    bestDelay, coarseDoppler + frac, timePilot, residual, [c1, c2]);
                [fracOpt, ~] = fminbnd(costFn, -0.5, 0.5, optimset('TolX', 1e-3));
                finalDoppler = coarseDoppler + fracOpt;

                [cleanSig, pathEnergy] = obj.simulatePathEffect( ...
                    bestDelay, finalDoppler, timePilot, [c1, c2]);
                estGain = (cleanSig' * residual) / (pathEnergy + 1e-12);
                residual = residual - estGain * cleanSig;
                estPaths(pathIdx, :) = [bestDelay, finalDoppler, estGain];
            end

            % --- 鍏ㄥ眬 LS 鑱斿悎浼樺寲澧炵泭 ---
            basisMat = zeros(N, obj.Config.NumPaths);

            % 迭代处理：按当前索引更新状态。
            for pathIdx = 1:obj.Config.NumPaths
                basisMat(:, pathIdx) = obj.simulatePathEffect( ...
                    estPaths(pathIdx, 1), estPaths(pathIdx, 2), timePilot, [c1, c2]);
            end

            refinedGains = lsqminnorm(basisMat, rxSignalDaft);
            estPaths(:, 3) = refinedGains;
            finalPaths = estPaths;

            % --- 閲嶅缓绛夋晥淇￠亾 ---
            physChannel = obj.rebuildChannel(finalPaths);
            effectiveChannel = AfdmTransforms.generateEffectiveChannelMatrix(physChannel, obj.Config);
        end

        % ================================================================
        % ================================================================
        % simulatePathEffect: 仿真单路径对导频的响应。
        function [daftSig, energy] = simulatePathEffect(~, delay, doppler, timePilot, chirpParams)
            N = size(timePilot, 1);
            timeVec = (0:N - 1).';
            phaseRot = exp(-1j * 2 * pi * doppler / N .* timeVec);

            pathTimeDomain = circshift(timePilot, delay) .* phaseRot;
            daftSig = AfdmTransforms.daft(pathTimeDomain, chirpParams(1), chirpParams(2));


            if nargout > 1
                energy = sum(abs(daftSig) .^ 2);
            end
        end

        % ================================================================
        % ================================================================
        % applyPeakWindow: 仅保留峰值邻域抑制旁瓣干扰。
        function masked = applyPeakWindow(~, sig, radius)
            N = length(sig);
            [~, peakIdx] = max(abs(sig));
            winRaw = (peakIdx - radius:peakIdx + radius).';
            winIdx = mod(winRaw - 1, N) + 1;
            mask = zeros(N, 1);
            mask(winIdx) = 1;
            masked = sig .* mask;
        end

        % ================================================================
        % ================================================================
        % calcCorrelation: 计算候选路径与残差的相关性。
        function val = calcCorrelation(obj, delay, doppler, timePilot, rxSig, chirpParams)
            [daftSig, energy] = obj.simulatePathEffect(delay, doppler, timePilot, chirpParams);

            if energy > 1e-12
                val = abs(daftSig' * rxSig) / sqrt(energy);
            else
                val = 0;
            end
        end

        % ================================================================
        %
        %    H(r, c) = 危_i  h_i 路 exp(-j2蟺路谓_i路r/N_tot)
        %
        % ================================================================
        % rebuildChannel: 由估计路径重建时域物理信道矩阵。
        function physChannel = rebuildChannel(obj, pathParams)
            totalSc = obj.Config.TotalSubcarriers;
            numDataSc = obj.Config.NumDataSubcarriers;
            numPaths = size(pathParams, 1);

            % 多普勒缩放: 估计器用 N_data 归一化, 物理信道用 N_total 归一化
            dopplerScale = totalSc / numDataSc;

            allRows = [];
            allCols = [];
            allVals = [];

            colVec = (0:totalSc - 1).';
            for i = 1:numPaths
                delay = round(pathParams(i, 1));
                doppler = pathParams(i, 2) * dopplerScale;
                gain = pathParams(i, 3);

                rowVec = mod(colVec + delay, totalSc);

                % 多普勒相位: 作用在行索引上
                phases = exp(-1j * 2 * pi * doppler * rowVec / totalSc);
                vals = gain * phases;

                allRows = [allRows; rowVec + 1]; %#ok<AGROW>
                allCols = [allCols; colVec + 1]; %#ok<AGROW>
                allVals = [allVals; vals];        %#ok<AGROW>
            end

            physChannel = sparse(allRows, allCols, allVals, totalSc, totalSc);
        end

    end

end

