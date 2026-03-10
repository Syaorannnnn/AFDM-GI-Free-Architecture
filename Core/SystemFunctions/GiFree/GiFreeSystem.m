classdef GiFreeSystem < handle
    % GiFreeSystem  GI-Free AFDM 系统级编排器
    %
    %   Agile-AFDM 集成 (三模式):
    %     AgileQ > 0 时, runTrial 根据 AgileMode 自动选择 c2 优化策略:
    %       'standard'     — Q 次批量 IFFT, 标准 PAPR 度量
    %       'pilotAware'   — Q 次批量 IFFT, 导频感知交叉项度量
    %       'hierarchical' — ~3√Q 次 IFFT, 两阶段层级搜索
    %     AgileQ = 0 → 默认 c2, 完全向后兼容.

    properties (GetAccess = public, SetAccess = private)
        Config
        Transmitter
        ChannelBuilder
        Estimator
        Receiver
        C2Candidates
    end

    methods

        function obj = GiFreeSystem(cfg)
            arguments
                cfg (1,1) GiFreeConfig
            end
            obj.Config         = cfg;
            obj.Transmitter    = GiFreeTransmitter(cfg);
            obj.ChannelBuilder = FractionalChannelBuilder(cfg);
            obj.Estimator      = GiFreeEstimator(cfg, obj.ChannelBuilder);
            obj.Receiver       = GiFreeReceiver(cfg, obj.Estimator);

            if cfg.AgileQ > 0
                obj.C2Candidates = AgileC2Optimizer.generateCandidates( ...
                    cfg.NumSubcarriers, cfg.AgileQ);
            else
                obj.C2Candidates = [];
            end
        end

        function result = runTrial(obj, dataSnrLinear, noisePower)
            cfg   = obj.Config;
            numSc = cfg.NumSubcarriers;

            % ---- 1. 发射: DAFT 域组帧 (不依赖 c2) ----
            [txFrame, txDataIndices] = obj.Transmitter.transmit(dataSnrLinear);

            % ---- 2. Agile-c2: 按 AgileMode 选择最优 c2 ----
            agileIdx = 0;
            agileFFTs = 0;
            if cfg.AgileQ > 0
                switch lower(cfg.AgileMode)
                    case 'pilotaware'
                        [bestC2, agileIdx] = AgileC2Optimizer.selectPilotAware( ...
                            txFrame, obj.C2Candidates, cfg.PilotAmplitude);
                        agileFFTs = cfg.AgileQ;
                    case 'hierarchical'
                        [bestC2, agileIdx, hInfo] = AgileC2Optimizer.selectHierarchical( ...
                            txFrame, obj.C2Candidates);
                        agileFFTs = hInfo.totalFFTs;
                    otherwise  % 'standard'
                        [bestC2, agileIdx] = AgileC2Optimizer.selectForBlock( ...
                            txFrame, obj.C2Candidates);
                        agileFFTs = cfg.AgileQ;
                end
                cfg.C2Override = bestC2;
            else
                cfg.C2Override = [];
            end

            % ---- 3. PAPR 测量 ----
            tdSig   = AfdmTransforms.idaft(txFrame, cfg.ChirpParam1, cfg.ChirpParam2);
            paprLin = max(abs(tdSig).^2) / mean(abs(tdSig).^2);

            % ---- 4. 信道 + 噪声 ----
            [~, hEffTrue] = obj.ChannelBuilder.generateChannel();
            noiseVec = sqrt(noisePower / 2) * (randn(numSc,1) + 1j*randn(numSc,1));
            rxSignal = hEffTrue * txFrame + noiseVec;

            % ---- 5. 接收 ----
            [detIndicesSys, mseSys, ~] = obj.Receiver.receive( ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower);

            % ---- 6. Perfect CSI 对比 ----
            detIndicesRef = GiFreeSystem.perfectCsiDetect( ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower, cfg);

            % ---- 7. BER 统计 ----
            bitsPerSymbol = log2(cfg.ModulationOrder);
            [bitErrSys, ~, totalBits] = GiFreeSystem.countBitErrors( ...
                txDataIndices, detIndicesSys, bitsPerSymbol);
            [bitErrRef, ~, ~] = GiFreeSystem.countBitErrors( ...
                txDataIndices, detIndicesRef, bitsPerSymbol);

            % ---- 结果打包 ----
            result.bitErrorsSys = bitErrSys;
            result.bitErrorsRef = bitErrRef;
            result.totalBits    = totalBits;
            result.mseSystem    = mseSys;
            result.paprDb       = 10 * log10(paprLin);
            result.c2Used       = cfg.ChirpParam2;
            result.agileIdx     = agileIdx;
            result.agileFFTs    = agileFFTs;
        end

    end

    methods (Static)

        function detectedIndices = perfectCsiDetect( ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower, cfg)
            numSc    = cfg.NumSubcarriers;
            modOrder = cfg.ModulationOrder;

            pilotFrame = zeros(numSc, 1);
            pilotPos1  = cfg.PilotPositions + 1;
            pilotFrame(pilotPos1) = cfg.PerPilotAmplitude * cfg.PilotSequence;
            cleanDataSig = rxSignal - hEffTrue * pilotFrame;

            dataPos1 = cfg.DataPositions + 1;
            numData  = cfg.NumDataSymbols;
            regParam = noisePower / dataSnrLinear;
            hData    = hEffTrue(:, dataPos1);
            estData  = (hData' * hData + regParam * speye(numData)) \ ...
                       (hData' * cleanDataSig);

            detectedIndices = qamdemod( ...
                estData / sqrt(dataSnrLinear), modOrder, 'UnitAveragePower', true);
        end

        function [bitErrors, ber, totalBits] = countBitErrors(txIndices, rxIndices, bitsPerSymbol)
            txBits    = de2bi(txIndices, bitsPerSymbol, 'left-msb');
            rxBits    = de2bi(rxIndices, bitsPerSymbol, 'left-msb');
            bitErrors = sum(txBits(:) ~= rxBits(:));
            totalBits = numel(txBits);
            ber       = bitErrors / totalBits;
        end

    end

end
