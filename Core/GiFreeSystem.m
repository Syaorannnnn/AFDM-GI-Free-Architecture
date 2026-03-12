classdef GiFreeSystem < handle
    % GiFreeSystem  GI-Free AFDM 系统级编排器
    %
    %   顶层入口: 持有并组合 Transmitter, ChannelBuilder, Estimator, Receiver.
    %   所有子模块共享同一个 GiFreeConfig handle.
    %
    %   用法:
    %     cfg = GiFreeConfig();
    %     cfg.NumSubcarriers = 256;  ... % 外部初始化所有参数
    %     sys = GiFreeSystem(cfg);
    %     result = sys.runTrial(snrLin, noisePow);

    properties (SetAccess = private)
        Config
        Transmitter
        ChannelBuilder
        Estimator
        Receiver
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
        end

        function result = runTrial(obj, dataSnrLinear, noisePower)
            % runTrial  执行一次完整仿真 (内部随机信道)
            %
            %   输出:
            %     result.bitErrorsSys  : 系统 BER 比特错误数
            %     result.bitErrorsRef  : Perfect CSI 比特错误数
            %     result.totalBits     : 总比特数
            %     result.mseSystem     : 信道估计 NMSE

            numSc         = obj.Config.NumSubcarriers;
            modOrder      = obj.Config.ModulationOrder;
            bitsPerSymbol = log2(modOrder);

            % 发射
            [txFrame, txDataIndices] = obj.Transmitter.transmit(dataSnrLinear);

            % 信道 + 噪声
            [~, hEffTrue] = obj.ChannelBuilder.generateChannel();
            noiseVec = sqrt(noisePower / 2) * (randn(numSc,1) + 1j*randn(numSc,1));
            rxSignal = hEffTrue * txFrame + noiseVec;

            % 接收
            [detIndicesSys, mseSys, ~] = obj.Receiver.receive( ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower);

            % Perfect CSI 对比
            detIndicesRef = GiFreeSystem.perfectCsiDetect( ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower, obj.Config);

            % BER 统计 (MATLAB 内置 biterr)
            [bitErrSys, ~, totalBits] = GiFreeSystem.countBitErrors( ...
                txDataIndices, detIndicesSys, bitsPerSymbol);
            [bitErrRef, ~, ~] = GiFreeSystem.countBitErrors( ...
                txDataIndices, detIndicesRef, bitsPerSymbol);

            result.bitErrorsSys = bitErrSys;
            result.bitErrorsRef = bitErrRef;
            result.totalBits    = totalBits;
            result.mseSystem    = mseSys;
        end

    end

    methods (Static)

        function detectedIndices = perfectCsiDetect( ...
                rxSignal, hEffTrue, dataSnrLinear, noisePower, cfg)
            numSc    = cfg.NumSubcarriers;
            modOrder = cfg.ModulationOrder;

            pilotFrame    = zeros(numSc, 1);
            pilotFrame(1) = cfg.PilotAmplitude;
            cleanDataSig  = rxSignal - hEffTrue * pilotFrame;

            regParam = noisePower / dataSnrLinear;
            hData    = hEffTrue(:, 2:numSc);
            estData  = (hData' * hData + regParam * speye(numSc-1)) \ (hData' * cleanDataSig);

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
