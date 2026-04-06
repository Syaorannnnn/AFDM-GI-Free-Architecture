classdef GiFreeSystem < handle
% GIFREESYSTEM - GI-Free 端到端单次试验封装
%
%   描述:
%   该类串联 GI-Free 发射机、随机信道采样器、接收机与 Perfect-CSI 参考检测，
%   对外提供单次 Monte-Carlo 试验接口，并输出 BER/NMSE/PAPR 等关键指标。
%
%   语法:
%   systemObj = GiFreeSystem(cfg);
%   result = systemObj.runTrialBySnrDb(dataSnrDb, noisePowerLin);
%   result = systemObj.runTrial(dataSnrLin, noisePowerLin);
%
%   输入:
%   cfg - (GiFreeConfig) GI-Free 系统配置对象。
%
%   输出:
%   result - (struct) 单次试验结果，包含 bitErrorsSys/bitErrorsRef/totalBits/
%            mseSystem/paprDb 等字段。
%
%   另请参阅:
%   GiFreeConfig, GiFreeTransmitter, GiFreeReceiver, GiFreeEstimator
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。

    properties (GetAccess = public, SetAccess = private)
        Config
        Transmitter
        ChannelBuilder
        ChannelSampler
        Estimator
        Receiver
    end

    methods

        % GIFREESYSTEM 构造函数，初始化 GI-Free 收发链路组件。
        function obj = GiFreeSystem(cfg)
            arguments
                cfg (1,1) GiFreeConfig
            end
            cfg.validate();
            obj.Config = cfg;
            obj.Transmitter = GiFreeTransmitter(cfg);
            obj.ChannelBuilder = FractionalChannelBuilder(cfg);
            obj.ChannelSampler = GiFreeChannelSampler(cfg, obj.ChannelBuilder);
            obj.Estimator = GiFreeEstimator(cfg, obj.ChannelBuilder);
            obj.Receiver = GiFreeReceiver(cfg, obj.Estimator);
        end

        % RUNTRIALBYSNRDB 以 dB SNR 形式运行一次试验。
        function result = runTrialBySnrDb(obj, dataSnrDb, noisePowerLin)
            dataSnrLin = 10 ^ (dataSnrDb / 10);
            result = obj.runTrial(dataSnrLin, noisePowerLin);
        end

        % RUNTRIAL 执行一次完整 GI-Free 仿真试验（线性 SNR 入口）。
        function result = runTrial(obj, dataSnrLin, noisePowerLin)
            cfg = obj.Config;
            numSc = cfg.NumSubcarriers;

            % 动态导频模式：在发射前更新当前 SNR，使 PilotAmplitude 跟踪数据 SNR
            if cfg.UseDynamicPilot
                cfg.CurrentDataSnrLin = dataSnrLin;
            end

            [txFrame, txDataIndices] = obj.Transmitter.transmit(dataSnrLin);

            tdSig = AfdmTransforms.idaft(txFrame, cfg.ChirpParam1, cfg.ChirpParam2);
            paprLin = max(abs(tdSig).^2) / mean(abs(tdSig).^2);

            [~, trueEffectiveChannel] = obj.ChannelSampler.sampleChannel();
            noiseVec = sqrt(noisePowerLin / 2) * (randn(numSc,1) + 1j*randn(numSc,1));
            rxSignal = trueEffectiveChannel * txFrame + noiseVec;

            [detIndicesSys, mseSys, ~] = obj.Receiver.receive(rxSignal, trueEffectiveChannel, dataSnrLin, noisePowerLin);

            detIndicesRef = GiFreeSystem.perfectCsiDetect(rxSignal, trueEffectiveChannel, dataSnrLin, noisePowerLin, cfg);

            bitsPerSymbol = log2(cfg.ModulationOrder);
            [bitErrSys, ~, totalBits] = GiFreeSystem.countBitErrors(txDataIndices, detIndicesSys, bitsPerSymbol);
            [bitErrRef, ~, ~] = GiFreeSystem.countBitErrors(txDataIndices, detIndicesRef, bitsPerSymbol);

            result.bitErrorsSys = bitErrSys;
            result.bitErrorsRef = bitErrRef;
            result.totalBits = totalBits;
            result.mseSystem = mseSys;
            result.paprDb = 10 * log10(paprLin);
        end

        % RUNTRIALWITHCHANNEL 使用外部信道参数执行一次试验（用于公平对比）。
        %
        %   与 runTrial 的区别：接受外部 (delays, dopplers, gains) 向量构建信道，
        %   不调用内部 ChannelSampler，使 EP 与 GI-Free 可共享同一组信道参数。
        function result = runTrialWithChannel(obj, dataSnrLin, noisePowerLin, ...
                delays, dopplers, gains)
            cfg = obj.Config;
            numSc = cfg.NumSubcarriers;

            if cfg.UseDynamicPilot
                cfg.CurrentDataSnrLin = dataSnrLin;
            end

            [txFrame, txDataIndices] = obj.Transmitter.transmit(dataSnrLin);

            % 用外部信道参数构建有效信道矩阵
            pathParams = [delays(:), dopplers(:), gains(:)];
            trueEffectiveChannel = obj.ChannelBuilder.buildEffectiveChannel(pathParams);

            noiseVec = sqrt(noisePowerLin / 2) * (randn(numSc,1) + 1j*randn(numSc,1));
            rxSignal = trueEffectiveChannel * txFrame + noiseVec;

            [detIndicesSys, mseSys, ~] = obj.Receiver.receive( ...
                rxSignal, trueEffectiveChannel, dataSnrLin, noisePowerLin);

            detIndicesRef = GiFreeSystem.perfectCsiDetect( ...
                rxSignal, trueEffectiveChannel, dataSnrLin, noisePowerLin, cfg);

            bitsPerSymbol = log2(cfg.ModulationOrder);
            [bitErrSys, ~, totalBits] = GiFreeSystem.countBitErrors( ...
                txDataIndices, detIndicesSys, bitsPerSymbol);
            [bitErrRef, ~, ~] = GiFreeSystem.countBitErrors( ...
                txDataIndices, detIndicesRef, bitsPerSymbol);

            result.bitErrorsSys = bitErrSys;
            result.bitErrorsRef = bitErrRef;
            result.totalBits    = totalBits;
            result.mseSystem    = mseSys;
            result.paprDb       = 0;
        end

        % RUNTRIALORACLEID2P 以真实数据消除 ID2P，测试信道估计性能上界。
        %
        %   用已知发射数据完美消除 ID2P 干扰，然后仅做单次 OMP + LMMSE。
        %   该曲线表示 "若 CLTP 软判决完美准确，GI-Free 信道估计能达到的上界"。
        function result = runTrialOracleId2p(obj, dataSnrLin, noisePowerLin, ...
                delays, dopplers, gains)
            cfg = obj.Config;
            numSc = cfg.NumSubcarriers;

            if cfg.UseDynamicPilot
                cfg.CurrentDataSnrLin = dataSnrLin;
            end

            [txFrame, txDataIndices] = obj.Transmitter.transmit(dataSnrLin);

            pathParams = [delays(:), dopplers(:), gains(:)];
            trueEffectiveChannel = obj.ChannelBuilder.buildEffectiveChannel(pathParams);

            noiseVec = sqrt(noisePowerLin / 2) * (randn(numSc,1) + 1j*randn(numSc,1));
            rxSignal = trueEffectiveChannel * txFrame + noiseVec;

            % Oracle: 用真实数据帧完美消除 ID2P
            dataFrame = zeros(numSc, 1);
            dataFrame(cfg.DataPos1) = txFrame(cfg.DataPos1);
            cleanPilotSig = rxSignal - trueEffectiveChannel * dataFrame;

            % 无 ID2P 的单次 OMP 估计
            [~, estEffChannel] = obj.Estimator.estimateByOmp(cleanPilotSig, 0);

            % 用估计信道做 LMMSE 检测
            pilotFrame = zeros(numSc, 1);
            pilotFrame(cfg.PilotPos1) = cfg.PerPilotAmplitude * cfg.PilotSequence;
            cleanDataSig = rxSignal - estEffChannel * pilotFrame;
            regParam = noisePowerLin / dataSnrLin;
            hData = estEffChannel(:, cfg.DataPos1);
            numData = cfg.NumDataSymbols;
            estData = (hData' * hData + regParam * speye(numData)) \ ...
                      (hData' * cleanDataSig);
            detIndices = qamdemod(estData / sqrt(dataSnrLin), ...
                cfg.ModulationOrder, 'UnitAveragePower', true);

            bitsPerSymbol = log2(cfg.ModulationOrder);
            [bitErrSys, ~, totalBits] = GiFreeSystem.countBitErrors( ...
                txDataIndices, detIndices, bitsPerSymbol);

            result.bitErrorsSys = bitErrSys;
            result.bitErrorsRef = 0;
            result.totalBits    = totalBits;
            result.mseSystem    = norm(full(estEffChannel) - full(trueEffectiveChannel), 'fro')^2 / ...
                max(norm(full(trueEffectiveChannel), 'fro')^2, 1e-20);
            result.paprDb       = 0;
        end

    end

    methods (Static)

        % PERFECTCSIDETECT 使用真值信道进行 LMMSE 检测，作为性能参考上界。
        function detectedIndices = perfectCsiDetect(rxSignal, trueEffectiveChannel, dataSnrLin, noisePowerLin, cfg)
            numSc = cfg.NumSubcarriers;
            modOrder = cfg.ModulationOrder;

            pilotFrame = zeros(numSc, 1);
            pilotFrame(cfg.PilotPos1) = cfg.PerPilotAmplitude * cfg.PilotSequence;
            cleanDataSig = rxSignal - trueEffectiveChannel * pilotFrame;

            dataPos1 = cfg.DataPos1;
            numData = cfg.NumDataSymbols;
            regParam = noisePowerLin / dataSnrLin;
            hData = trueEffectiveChannel(:, dataPos1);
            estData = (hData' * hData + regParam * speye(numData)) \ ...
                      (hData' * cleanDataSig);

            detectedIndices = qamdemod(estData / sqrt(dataSnrLin), modOrder, 'UnitAveragePower', true);
        end

        % COUNTBITERRORS 统计误比特数并计算 BER。
        function [bitErrors, ber, totalBits] = countBitErrors(txIndices, rxIndices, bitsPerSymbol)
            txBits = de2bi(txIndices, bitsPerSymbol, 'left-msb');
            rxBits = de2bi(rxIndices, bitsPerSymbol, 'left-msb');
            bitErrors = sum(txBits(:) ~= rxBits(:));
            totalBits = numel(txBits);
            ber = bitErrors / totalBits;
        end

    end

end


