function [txFrame, txDataIndices] = transmitWithFixedSeed(systemObj, dataSnrLin, rngState)
% TRANSMITWITHFIXEDSEED 用固定随机态生成一帧，便于跨配置共享同一组数据。
%
%   语法:
%   [txFrame, txDataIndices] = transmitWithFixedSeed(systemObj, dataSnrLin, rngState)
%
%   输入:
%   systemObj - (GiFreeSystem) GI-Free 系统对象。
%   dataSnrLin - (double) 数据 SNR 线性值。
%   rngState - (struct) trial 开始时保存的随机态。
%
%   输出:
%   txFrame - (Nx1 double) 发射帧。
%   txDataIndices - (Kx1 double) 发送数据索引。
%
%   版本历史:
%   2026-04-19 - Aiden - 从 MainSimulation 提取为共享辅助函数。

    arguments
        systemObj (1,1) GiFreeSystem
        dataSnrLin (1,1) double {mustBePositive}
        rngState (1,1) struct
    end

    previousState = rng;
    cleanupObj = onCleanup(@() rng(previousState));
    rng(rngState);

    if systemObj.Config.UseDynamicPilot
        systemObj.Config.CurrentDataSnrLin = dataSnrLin;
    end
    [txFrame, txDataIndices] = systemObj.Transmitter.transmit(dataSnrLin);
    clear cleanupObj;
end
