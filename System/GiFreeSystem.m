classdef GiFreeSystem
    % GiFreeSystem  GI-Free AFDM 系统级编排器
    %
    %   将发射、信道、接收三大模块解耦并串联:
    %     GiFreeTransmitter  → FractionalChannelBuilder → GiFreeReceiver
    %
    %   提供三个层级的仿真接口:
    %     runTrial      : 单次完整的 "发射→信道→接收→评估" 流程
    %     perfectCsiDetect : 已知真实信道的 LMMSE 检测 (性能上界)
    %     computeBer    : 符号索引 → 比特错误统计

    methods (Static)

        %% ========== 单次试验 ==========
        function trialResult = runTrial(cfg, dataSnrLinear, noisePower)
            % runTrial  执行一次完整的 GI-Free 系统仿真
            %
            %   输出结构体包含:
            %     berSystem    : 系统 BER
            %     berPerfCsi   : Perfect CSI BER (性能上界)
            %     mseSystem    : 信道估计 NMSE
            %     bitErrorsSys : 系统比特错误数
            %     bitErrorsRef : Perfect CSI 比特错误数
            %     totalBits    : 总比特数

            numSc         = cfg.NumSubcarriers;
            modOrder      = cfg.ModulationOrder;
            bitsPerSymbol = log2(modOrder);

            % ---- 发射 ----
            [txFrame, txDataIndices] = GiFreeTransmitter.transmit(cfg, dataSnrLinear);

            % ---- 信道 ----
            [~, hEffTrue] = FractionalChannelBuilder.generateChannel(cfg);
            noiseVec = sqrt(noisePower / 2) * ...
                (randn(numSc, 1) + 1j * randn(numSc, 1));
            rxSignal = hEffTrue * txFrame + noiseVec;

            % ---- 接收 (GiFreeReceiver) ----
            [detIndicesSys, mseSys, ~] = GiFreeReceiver.receive( ...
                rxSignal, cfg, hEffTrue, dataSnrLinear, noisePower);

            % ---- Perfect CSI 对比 ----
            detIndicesRef = GiFreeSystem.perfectCsiDetect( ...
                rxSignal, hEffTrue, cfg, dataSnrLinear, noisePower);

            % ---- 评估 ----
            [bitErrSys, totalBits] = GiFreeSystem.computeBer( ...
                txDataIndices, detIndicesSys, bitsPerSymbol);
            [bitErrRef, ~] = GiFreeSystem.computeBer( ...
                txDataIndices, detIndicesRef, bitsPerSymbol);

            % ---- 打包结果 ----
            trialResult.bitErrorsSys = bitErrSys;
            trialResult.bitErrorsRef = bitErrRef;
            trialResult.totalBits    = totalBits;
            trialResult.mseSystem    = mseSys;
        end

        %% ========== Perfect CSI 检测 ==========
        function detectedIndices = perfectCsiDetect( ...
                rxSignal, hEffTrue, cfg, dataSnrLinear, noisePower)
            % perfectCsiDetect  已知真实信道的 LMMSE 数据检测
            %   使用理论最优正则化 lambda = sigma_n^2 / SNR_data
            %   这是系统性能的理论上界 (仅受热噪声限制)

            numSc    = cfg.NumSubcarriers;
            modOrder = cfg.ModulationOrder;

            pilotFrame    = zeros(numSc, 1);
            pilotFrame(1) = cfg.PilotAmplitude;
            cleanDataSig  = rxSignal - hEffTrue * pilotFrame;

            regParam = noisePower / dataSnrLinear;
            hData    = hEffTrue(:, 2:numSc);
            estData  = (hData' * hData + regParam * speye(numSc - 1)) \ ...
                       (hData' * cleanDataSig);

            detectedIndices = qamdemod( ...
                estData / sqrt(dataSnrLinear), modOrder, 'UnitAveragePower', true);
        end

        %% ========== BER 计算 ==========
        function [bitErrors, totalBits] = computeBer(txIndices, rxIndices, bitsPerSymbol)
            % computeBer  比较发射与接收符号索引, 统计比特错误
            txBits    = de2bi(txIndices, bitsPerSymbol, 'left-msb');
            rxBits    = de2bi(rxIndices, bitsPerSymbol, 'left-msb');
            bitErrors = sum(txBits(:) ~= rxBits(:));
            totalBits = numel(txBits);
        end

    end
end
