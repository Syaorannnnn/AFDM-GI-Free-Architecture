function configObj = createExp2GiFreeConfig( ...
        numSubcarriers, modulationOrder, maxDelaySamples, maxDopplerIdx, numPaths, pilotSnrDb)
% CREATEEXP2GIFREECONFIG 构造 Exp2 使用的 GI-Free 配置对象。
%
%   语法:
%   configObj = createExp2GiFreeConfig(numSubcarriers, modulationOrder, ...
%       maxDelaySamples, maxDopplerIdx, numPaths, pilotSnrDb)
%
%   输入:
%   numSubcarriers - (double) 子载波数量。
%   modulationOrder - (double) 调制阶数。
%   maxDelaySamples - (double) 最大时延索引。
%   maxDopplerIdx - (double) 最大整数 Doppler 索引。
%   numPaths - (double) 信道路径数。
%   pilotSnrDb - (double) 固定导频功率对应的 Pilot SNR dB。
%
%   输出:
%   configObj - (GiFreeConfig) Exp2 实验配置对象。
%
%   版本历史:
%   2026-04-21 - Aiden - 抽取 Exp2 GI-Free 配置构造函数。
%   2026-04-22 - 调整为固定导频版本，用于阶段一检测器归因实验。

    arguments
        numSubcarriers (1,1) double {mustBeInteger, mustBePositive}
        modulationOrder (1,1) double {mustBeInteger, mustBePositive}
        maxDelaySamples (1,1) double {mustBeInteger, mustBeNonnegative}
        maxDopplerIdx (1,1) double {mustBeInteger, mustBeNonnegative}
        numPaths (1,1) double {mustBeInteger, mustBePositive}
        pilotSnrDb (1,1) double {mustBeFinite}
    end

    configObj = GiFreeConfig();
    configObj.NumSubcarriers = numSubcarriers;
    configObj.ModulationOrder = modulationOrder;
    configObj.MaxDelaySamples = maxDelaySamples;
    configObj.MaxDopplerIdx = maxDopplerIdx;
    configObj.NumPaths = numPaths;
    configObj.DopplerGuard = 4;
    configObj.DirichletRadius = 0;
    configObj.PilotSnrDb = pilotSnrDb;
    configObj.MaxSicIterations = 10;
    configObj.NumPathsUpper = resolveExp2NumPathsUpper(numPaths);
    configObj.UseFractionalDoppler = false;
    configObj.UseDynamicPilot = false;
    configObj.DynamicPilotBaseDb = pilotSnrDb;
end
