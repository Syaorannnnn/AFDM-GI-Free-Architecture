classdef GiFreeTransmitter < handle
    % GiFreeTransmitter  GI-Free AFDM 发射机 (CAZAC 多导频版)
    %
    %   K 个 DAFT 域子载波放置 ZC 导频序列, 其余 N-K 个承载 QAM 数据.
    %   K = 1 时退化为原始单导频行为.

    properties (SetAccess = private)
        Config
    end

    methods

        function obj = GiFreeTransmitter(cfg)
            arguments
                cfg (1,1) GiFreeConfig
            end
            obj.Config = cfg;
        end

        function [txFrame, txDataIndices] = transmit(obj, dataSnrLinear)
            numSc    = obj.Config.NumSubcarriers;
            modOrder = obj.Config.ModulationOrder;
            numData  = obj.Config.NumDataSymbols;

            % 数据符号 (N-K 个)
            txDataIndices = randi([0, modOrder-1], numData, 1);
            qamSymbols    = qammod(txDataIndices, modOrder, 'UnitAveragePower', true);

            % 组帧: 导频 + 数据
            txFrame = zeros(numSc, 1);
            pilotPos1 = obj.Config.PilotPositions + 1;       % 1-based
            dataPos1  = obj.Config.DataPositions + 1;         % 1-based
            txFrame(pilotPos1) = obj.Config.PerPilotAmplitude * obj.Config.PilotSequence;
            txFrame(dataPos1)  = qamSymbols * sqrt(dataSnrLinear);
        end

    end

end
