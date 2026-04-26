function [fixedConfig, dynamicConfig] = createClipPilotConfigPair( ...
        numSubcarriers, modulationOrder, maxDelaySamples, maxDopplerIdx, ...
        dopplerGuard, numPaths, pilotSnrDb)
% CREATECLIPPILOTCONFIGPAIR 构造仅在导频策略上不同的一对 CLIP 配置。
%
%   语法:
%   [fixedConfig, dynamicConfig] = createClipPilotConfigPair( ...
%       numSubcarriers, modulationOrder, maxDelaySamples, maxDopplerIdx, ...
%       dopplerGuard, numPaths, pilotSnrDb)
%
%   输入:
%   numSubcarriers - (double) 子载波数量。
%   modulationOrder - (double) 调制阶数。
%   maxDelaySamples - (double) 最大时延索引。
%   maxDopplerIdx - (double) 最大整数 Doppler 索引。
%   dopplerGuard - (double) Doppler 保护宽度。
%   numPaths - (double) 路径数。
%   pilotSnrDb - (double) 固定导频 SNR，动态导频模式下同时作为基准偏置。
%
%   输出:
%   fixedConfig - (GiFreeConfig) 固定导频的最终 CLIP 配置。
%   dynamicConfig - (GiFreeConfig) 动态导频的最终 CLIP 配置。
%
%   版本历史:
%   2026-04-22 - 新增动态导频增益实验共享配置工厂。

    arguments
        numSubcarriers (1,1) double {mustBeInteger, mustBePositive}
        modulationOrder (1,1) double {mustBeInteger, mustBePositive}
        maxDelaySamples (1,1) double {mustBeInteger, mustBeNonnegative}
        maxDopplerIdx (1,1) double {mustBeInteger, mustBeNonnegative}
        dopplerGuard (1,1) double {mustBeInteger, mustBeNonnegative}
        numPaths (1,1) double {mustBeInteger, mustBePositive}
        pilotSnrDb (1,1) double {mustBeFinite}
    end

    baseConfig = GiFreeConfig();
    baseConfig.NumSubcarriers = numSubcarriers;
    baseConfig.ModulationOrder = modulationOrder;
    baseConfig.MaxDelaySamples = maxDelaySamples;
    baseConfig.MaxDopplerIdx = maxDopplerIdx;
    baseConfig.NumPaths = numPaths;
    baseConfig.DopplerGuard = dopplerGuard;
    baseConfig.DirichletRadius = 4;
    baseConfig.PilotSnrDb = pilotSnrDb;
    baseConfig.MaxSicIterations = 10;
    baseConfig.NumPathsUpper = 6;
    baseConfig.UseFractionalDoppler = true;
    baseConfig.EnableProgressiveCfar = true;
    baseConfig.EnablePathStabilityGate = true;
    baseConfig.UseWeightedPilotMetric = true;
    baseConfig.UseMatchedPilotMetric = true;

    fixedConfig = baseConfig.copy();
    fixedConfig.UseDynamicPilot = false;
    fixedConfig.DynamicPilotBaseDb = pilotSnrDb;

    dynamicConfig = baseConfig.copy();
    dynamicConfig.UseDynamicPilot = true;
    dynamicConfig.DynamicPilotBaseDb = pilotSnrDb;
end
