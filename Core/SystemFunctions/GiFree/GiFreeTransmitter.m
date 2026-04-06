classdef GiFreeTransmitter < handle
% GIFREETRANSMITTER - GI-Free 发射端 DAFT 域组帧
%
%   描述:
%   在 DAFT 域生成一帧发送信号：固定导频位于 PilotPos1，其余 N-1 个子载波
%   承载 QAM 数据符号。数据符号按 dataSnrLin 做幅度缩放。
%
%   语法:
%   txObj = GiFreeTransmitter(cfg);
%   [txFrame, txDataIndices] = txObj.transmit(dataSnrLin);
%
%   输入:
%   cfg        - (GiFreeConfig) 系统配置对象。
%   dataSnrLin - (double) 数据符号线性 SNR。
%
%   输出:
%   txFrame       - (Nx1 complex) DAFT 域发射帧（含导频+数据）。
%   txDataIndices - ((N-1)x1 int) 发送符号索引。
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。

    properties (SetAccess = private)
        Config
    end

    methods

        % GIFREETRANSMITTER 构造函数，绑定配置对象。
        function obj = GiFreeTransmitter(cfg)
            arguments
                cfg (1,1) GiFreeConfig
            end
            obj.Config = cfg;
        end

        % TRANSMIT 生成 GI-Free 发射帧并返回真实发送的数据索引。
        function [txFrame, txDataIndices] = transmit(obj, dataSnrLin)
            numSubcarriers = obj.Config.NumSubcarriers;
            modulationOrder = obj.Config.ModulationOrder;
            numDataSymbols  = obj.Config.NumDataSymbols;

            txDataIndices = randi([0, modulationOrder-1], numDataSymbols, 1);
            qamSymbols    = qammod(txDataIndices, modulationOrder, 'UnitAveragePower', true);

            txFrame = zeros(numSubcarriers, 1);
            pilotPos1 = obj.Config.PilotPos1;
            dataPos1  = obj.Config.DataPos1;
            txFrame(pilotPos1) = obj.Config.PerPilotAmplitude * obj.Config.PilotSequence;
            txFrame(dataPos1)  = qamSymbols * sqrt(dataSnrLin);
        end

    end

end

