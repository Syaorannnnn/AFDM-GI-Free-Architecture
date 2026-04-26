function [pilotFrame, pilotAmplitude] = buildGiFreePilotFrameForSnr(configObj, dataSnrLin)
% BUILDGIFREEPILOTFRAMEFORSNR 按当前数据 SNR 生成 GI-Free 导频帧与导频幅度。
%
%   语法:
%   [pilotFrame, pilotAmplitude] = buildGiFreePilotFrameForSnr(configObj, dataSnrLin)
%
%   输入:
%   configObj - (GiFreeConfig) GI-Free 配置对象。
%   dataSnrLin - (double) 当前数据 SNR 线性值。
%
%   输出:
%   pilotFrame - (Nx1 complex double) 当前 SNR 下的导频帧。
%   pilotAmplitude - (double) 当前 SNR 下的单导频幅度。
%
%   版本历史:
%   2026-04-21 - Aiden - 新增动态导频帧构造辅助函数。

    arguments
        configObj (1,1) GiFreeConfig
        dataSnrLin (1,1) double {mustBePositive}
    end

    if configObj.UseDynamicPilot
        configObj.CurrentDataSnrLin = dataSnrLin;
    end

    pilotAmplitude = configObj.PerPilotAmplitude;
    pilotFrame = zeros(configObj.NumSubcarriers, 1);
    pilotFrame(configObj.PilotPos1) = pilotAmplitude * configObj.PilotSequence;
end
