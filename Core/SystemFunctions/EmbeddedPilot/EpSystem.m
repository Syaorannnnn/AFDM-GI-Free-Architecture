classdef EpSystem < handle
% EpSystem: Embedded Pilot 端到端单次试验封装。

    properties (SetAccess = private)
        Config
        Transmitter
        Receiver
    end

    methods

        % EpSystem: 函数实现见下方代码。
        function obj = EpSystem(cfg)
            arguments
                cfg (1, 1) EpConfig
            end

            obj.Config = cfg;
            obj.Transmitter = EpTransmitter(cfg);
            obj.Receiver = EpReceiver(cfg);
        end

        % runTrial: 执行一次完整仿真试验。
        function result = runTrial(obj, dataSnrDb, pathDelays, pathDopplers, pathGains)
            totalSc = obj.Config.TotalSubcarriers;
            bitsPerSymbol = obj.Config.BitsPerSymbol;

            dataSnrLin = 10 ^ (dataSnrDb / 10);
            pilotPowerLin = 10 ^ (obj.Config.PilotSnrDb / 10);
            noisePowerLin = 1 / dataSnrLin;

            [txSig, txData] = obj.Transmitter.transmit(pilotPowerLin);
            physH = LtvChannel(totalSc, pathDelays, pathDopplers, pathGains);

            noiseVec = sqrt(noisePowerLin / 2) * (randn(totalSc, 1) + 1j * randn(totalSc, 1));
            rxSig = physH * txSig + noiseVec;

            obj.Receiver.CsiMode = "Estimated";
            rxData = obj.Receiver.receive(rxSig, noisePowerLin, physH, pilotPowerLin);

            [result.bitErrors, result.totalBits] = EpSystem.computeBer(txData, rxData, bitsPerSymbol);
        end

        % runTrialPerfectCsi: 在 Perfect CSI 假设下运行试验。
        function result = runTrialPerfectCsi(obj, dataSnrDb, pathDelays, pathDopplers, pathGains)
            totalSc = obj.Config.TotalSubcarriers;
            bitsPerSymbol = obj.Config.BitsPerSymbol;

            dataSnrLin = 10 ^ (dataSnrDb / 10);
            pilotPowerLin = 10 ^ (obj.Config.PilotSnrDb / 10);
            noisePowerLin = 1 / dataSnrLin;

            [txSig, txData] = obj.Transmitter.transmit(pilotPowerLin);
            physH = LtvChannel(totalSc, pathDelays, pathDopplers, pathGains);

            noiseVec = sqrt(noisePowerLin / 2) * (randn(totalSc, 1) + 1j * randn(totalSc, 1));
            rxSig = physH * txSig + noiseVec;

            obj.Receiver.CsiMode = "Perfect";
            rxData = obj.Receiver.receive(rxSig, noisePowerLin, physH, pilotPowerLin);

            [result.bitErrors, result.totalBits] = EpSystem.computeBer(txData, rxData, bitsPerSymbol);
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

