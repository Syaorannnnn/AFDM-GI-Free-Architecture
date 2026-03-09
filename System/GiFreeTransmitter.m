classdef GiFreeTransmitter
    % GiFreeTransmitter  GI-Free AFDM 发射机
    %   在 DAFT 域第 0 号子载波放置导频，其余子载波承载 QAM 数据符号。

    methods (Static)

        function [txFrame, txDataIndices] = transmit(cfg, dataSnrLinear)
            % transmit  生成一帧发射信号
            %   txFrame       : N x 1 复数发射向量
            %   txDataIndices : (N-1) x 1 QAM 符号索引 (用于后续 BER 计算)
            numSc    = cfg.NumSubcarriers;
            modOrder = cfg.ModulationOrder;

            % 随机数据符号
            txDataIndices = randi([0, modOrder - 1], numSc - 1, 1);
            qamSymbols    = qammod(txDataIndices, modOrder, 'UnitAveragePower', true);

            % 组帧: [导频; 数据]
            txFrame        = zeros(numSc, 1);
            txFrame(1)     = cfg.PilotAmplitude;
            txFrame(2:end) = qamSymbols * sqrt(dataSnrLinear);
        end

    end
end
