classdef OfdmSystem < handle
% OfdmSystem: CP-OFDM 基线系统，含 LS/Perfect CSI 两种检测链路。
    % OfdmSystem  标准 CP-OFDM 梳状导频基线系统
    %
    % 用法:
    %   cfg = OfdmConfig();
    %   cfg.DftSize = 512;  cfg.CpLength = 4;
    %   cfg.PilotSpacing = 102;  cfg.ModulationOrder = 4;
    %   sys = OfdmSystem(cfg);
    %   result = sys.runTrial(dataSnrDb, pathDelays, pathDopplers, pathGains);
    %
    % 提供两种模式:
    %   runTrial          — 梳状导频 LS 估计 + DFT 插值 + 单抽头 MMSE 均衡
    %   runTrialPerfectCsi — 完美频域信道 + 全矩阵 MMSE 均衡

    properties (SetAccess = private)
        Config
    end

    methods

        % OfdmSystem: 函数实现见下方代码。
        function obj = OfdmSystem(cfg)
            arguments
                cfg (1, 1) OfdmConfig
            end
            obj.Config = cfg;
        end

        %% ========== 单次试验 (梳状导频 LS 估计) ==========
        % runTrial: 执行一次完整仿真试验。
        function result = runTrial(obj, dataSnrDb, pathDelays, pathDopplers, pathGains)
            bps = obj.Config.BitsPerSymbol;
            [rxFreq, ~, txDataIdx, pilotSymbols, noisePowerLin] = ...
                obj.prepareRxFreq(dataSnrDb, pathDelays, pathDopplers, pathGains);

            % LS + DFT 插值信道估计
            hEstDiag = obj.estimateChannelLsDft(rxFreq, pilotSymbols);

            % 单抽头 MMSE 均衡
            rxDataIdx = obj.equalizeOneTap(rxFreq, hEstDiag, noisePowerLin);

            % BER
            [result.bitErrors, result.totalBits] = OfdmSystem.computeBer(txDataIdx, rxDataIdx, bps);
        end

        %% ========== 单次试验 (完美 CSI, 全矩阵 MMSE) ==========
        % runTrialPerfectCsi: 在 Perfect CSI 假设下运行试验。
        function result = runTrialPerfectCsi(obj, dataSnrDb, pathDelays, pathDopplers, pathGains)
            N   = obj.Config.DftSize;
            bps = obj.Config.BitsPerSymbol;
            [rxFreq, physH, txDataIdx, pilotSymbols, noisePowerLin] = ...
                obj.prepareRxFreq(dataSnrDb, pathDelays, pathDopplers, pathGains);

            Ntotal = obj.Config.TotalFrameLength;
            cpLen  = obj.Config.CpLength;

            % 构造完美频域信道矩阵 H_freq = F · H_eff_time · F^H
            % CP 插入矩阵: [CP部分; I_N], 物理信道作用后去掉 CP 行
            insertMat = zeros(Ntotal, N);
            insertMat(cpLen + 1 : end, :) = eye(N);            % 主体: 单位阵
            insertMat(1 : cpLen, N - cpLen + 1 : N) = eye(cpLen); % CP: 复制末尾

            hEffTime  = physH(cpLen + 1 : end, :) * insertMat;   % N × N
            % H_freq = F · hEffTime · F^H  (O(N²) 实现)
            temp      = fft(hEffTime) / sqrt(N);                  % 左乘 F
            hFreqFull = (ifft(temp.') .* sqrt(N)).';              % 右乘 F^H

            % 全矩阵 MMSE 均衡
            rxDataIdx = obj.equalizeFullMmse(rxFreq, hFreqFull, noisePowerLin, pilotSymbols);

            % BER
            [result.bitErrors, result.totalBits] = OfdmSystem.computeBer(txDataIdx, rxDataIdx, bps);
        end

    end

    methods (Access = private)

        %% ========== 发射 ==========
        % transmit: 生成发射帧并返回发送符号。
        function [txTimeDomain, txDataIdx, pilotSymbols] = transmit(obj)
            N     = obj.Config.DftSize;
            cpLen = obj.Config.CpLength;
            modOrd = obj.Config.ModulationOrder;

            % 随机数据
            txDataIdx  = randi([0, modOrd - 1], obj.Config.NumDataCarriers, 1);
            dataSymbols = qammod(txDataIdx, modOrd, 'UnitAveragePower', true);

            % 导频: 全 1 BPSK (已知序列)
            pilotSymbols = ones(obj.Config.NumPilots, 1);

            % 组帧
            freqFrame = zeros(N, 1);
            freqFrame(obj.Config.PilotPos1) = pilotSymbols;
            freqFrame(obj.Config.DataPos1)  = dataSymbols;

            % IDFT → 时域 (与 AfdmTransforms.idft 归一化一致)
            timeFrame = ifft(freqFrame) * sqrt(N);

            % 标准 CP
            txTimeDomain = [timeFrame(end - cpLen + 1 : end); timeFrame];
        end

        %% ========== 公共接收前处理: 发射→信道→去CP→DFT ==========
        function [rxFreq, physH, txDataIdx, pilotSymbols, noisePowerLin] = ...
                prepareRxFreq(obj, dataSnrDb, pathDelays, pathDopplers, pathGains)
            N      = obj.Config.DftSize;
            Ntotal = obj.Config.TotalFrameLength;
            cpLen  = obj.Config.CpLength;

            dataSnrLin     = 10 ^ (dataSnrDb / 10);
            noisePowerLin = 1 / dataSnrLin;

            % 发射
            [txSig, txDataIdx, pilotSymbols] = obj.transmit();

            % 物理信道 (与 EpSystem 使用同一 LtvChannel)
            physH = LtvChannel(Ntotal, pathDelays, pathDopplers, pathGains);

            % 加噪
            noiseVec = sqrt(noisePowerLin / 2) * (randn(Ntotal, 1) + 1j * randn(Ntotal, 1));
            rxSig = physH * txSig + noiseVec;

            % 去 CP → DFT
            rxNoCP = rxSig(cpLen + 1 : end);
            rxFreq = fft(rxNoCP) / sqrt(N);
        end

        %% ========== LS + DFT 插值信道估计 ==========
        % estimateChannelLsDft: 函数实现见下方代码。
        function hEstDiag = estimateChannelLsDft(obj, rxFreq, pilotSymbols)
            N    = obj.Config.DftSize;
            Np   = obj.Config.NumPilots;
            Ltap = obj.Config.CpLength + 1;   % 最大可分辨抽头数

            % Step 1: 导频位置 LS 估计
            hPilotLs = rxFreq(obj.Config.PilotPos1) ./ pilotSymbols;

            % Step 2: Np 点 IFFT → 时域 CIR
            hCir = ifft(hPilotLs);

            % Step 3: 截断噪声抽头 (仅保留前 L+1 个有效抽头)
            if Ltap < Np
                hCir(Ltap + 1 : end) = 0;
            end

            % Step 4: 零填充到 N 点 + FFT → 全子载波频响插值
            hCirPadded = zeros(N, 1);
            hCirPadded(1 : Np) = hCir;
            hEstDiag = fft(hCirPadded) * (Np / N);
        end

        %% ========== 单抽头 MMSE 均衡 ==========
        % equalizeOneTap: 函数实现见下方代码。
        function rxDataIdx = equalizeOneTap(obj, rxFreq, hDiag, noisePowerLin)
            dataIdx = obj.Config.DataPos1;
            hData   = hDiag(dataIdx);
            yData   = rxFreq(dataIdx);

            % MMSE 权重: w[k] = conj(H[k]) / (|H[k]|² + σ²)
            wMmse = conj(hData) ./ (abs(hData).^2 + noisePowerLin);
            xHat  = wMmse .* yData;

            % 解调
            rxDataIdx = qamdemod(xHat, obj.Config.ModulationOrder, 'UnitAveragePower', true);
        end

        %% ========== 全矩阵 MMSE 均衡 (Perfect CSI) ==========
        % equalizeFullMmse: 函数实现见下方代码。
        function rxDataIdx = equalizeFullMmse(obj, rxFreq, hFreqFull, noisePowerLin, pilotSymbols)
            N        = obj.Config.DftSize;
            dataIdx  = obj.Config.DataPos1;
            pilotIdx = obj.Config.PilotPos1;

            % 去除导频贡献 (IP2D)
            pilotFrame = zeros(N, 1);
            pilotFrame(pilotIdx) = pilotSymbols;
            rxClean = rxFreq - hFreqFull * pilotFrame;

            % 数据列上的全矩阵 MMSE
            Hd   = hFreqFull(:, dataIdx);
            nD   = numel(dataIdx);
            xHat = (Hd' * Hd + noisePowerLin * eye(nD)) \ (Hd' * rxClean);

            % 解调
            rxDataIdx = qamdemod(xHat, obj.Config.ModulationOrder, 'UnitAveragePower', true);
        end

    end

    methods (Static)

        % computeBer: 统计误比特数与总比特数。
        function [bitErrors, totalBits] = computeBer(txIdx, rxIdx, bitsPerSymbol)
            txBits = de2bi(txIdx, bitsPerSymbol, 'left-msb');
            rxBits = de2bi(rxIdx, bitsPerSymbol, 'left-msb');
            bitErrors = sum(txBits(:) ~= rxBits(:));
            totalBits = numel(txBits);
        end

    end

end

