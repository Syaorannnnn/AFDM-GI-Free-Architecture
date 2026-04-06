classdef EpTransmitter < handle
% EpTransmitter: Embedded Pilot 发射端组帧与前缀插入。

    properties (Access = private)
        Config
    end

    methods (Access = public)

        % EpTransmitter: 函数实现见下方代码。
        function obj = EpTransmitter(configObj)
            obj.Config = configObj;
        end

        % transmit: 生成发射帧并返回发送符号。
        function [txSignal, originalDataIdx] = transmit(obj, pilotPowerLin)
            numDataSubcarriers = obj.Config.NumDataSubcarriers;

            originalDataIdx = randi([0, obj.Config.ModulationOrder - 1], ...
                obj.Config.NumActiveCarriers, 1);
            qamSymbols = qammod(originalDataIdx, obj.Config.ModulationOrder, ...
                'UnitAveragePower', true);

            daftFrame = zeros(numDataSubcarriers, 1);
            daftFrame(obj.Config.PilotPos1) = sqrt(pilotPowerLin);
            daftFrame(obj.Config.DataPos1) = qamSymbols;

            if strcmpi(obj.Config.WaveformType, "AFDM")
                timeFrame = AfdmTransforms.idaft(daftFrame, ...
                    obj.Config.ChirpParam1, obj.Config.ChirpParam2);
            else
                timeFrame = AfdmTransforms.idft(daftFrame);
            end

            txSignal = obj.addPrefix(timeFrame);
        end

    end

    methods (Access = private)

        % addPrefix: 函数实现见下方代码。
        function txSignal = addPrefix(obj, timeFrame)
            prefixLen = obj.Config.PrefixLength;
            numDataSubcarriers = obj.Config.NumDataSubcarriers;

            if strcmpi(obj.Config.WaveformType, "AFDM")
                phaseVec = exp(-1j * 2 * pi * obj.Config.ChirpParam1 * ...
                    (numDataSubcarriers ^ 2 + 2 * numDataSubcarriers * (-prefixLen:-1).'));
                prefix = timeFrame(end - prefixLen + 1:end) .* phaseVec;
            else
                prefix = timeFrame(end - prefixLen + 1:end);
            end

            txSignal = [prefix; timeFrame];
        end

    end

end

