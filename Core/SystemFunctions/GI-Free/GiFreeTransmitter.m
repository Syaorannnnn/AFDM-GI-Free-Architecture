classdef GiFreeTransmitter < handle
    % GiFreeTransmitter  GI-Free AFDM 发射机
    %
    %   DAFT 域第 0 号子载波放置导频, 其余承载 QAM 数据.

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

            txDataIndices = randi([0, modOrder-1], numSc-1, 1);
            qamSymbols    = qammod(txDataIndices, modOrder, 'UnitAveragePower', true);

            txFrame        = zeros(numSc, 1);
            txFrame(1)     = obj.Config.PilotAmplitude;
            txFrame(2:end) = qamSymbols * sqrt(dataSnrLinear);
        end

    end

end
