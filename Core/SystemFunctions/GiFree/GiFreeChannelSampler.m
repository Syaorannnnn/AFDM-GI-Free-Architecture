classdef GiFreeChannelSampler < IChannelSampler
% GIFREECHANNELSAMPLER - GI-Free 随机信道参数采样器
%
%   描述:
%   负责生成随机路径参数 [delay, doppler, gain]，并调用 ChannelOperator
%   构造对应有效信道矩阵。支持整数/分数 Doppler 两种模式。
%
%   语法:
%   sampler = GiFreeChannelSampler(cfg, channelOperator);
%   pathParams = sampler.samplePathParams();
%   [pathParams, hEff] = sampler.sampleChannel();
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。

    properties (SetAccess = private)
        Config
        ChannelOperator
    end

    methods

        % GIFREECHANNELSAMPLER 构造函数，注入配置与信道算子。
        function obj = GiFreeChannelSampler(cfg, channelOperator)
            arguments
                cfg             (1,1) GiFreeConfig
                channelOperator (1,1) IChannelOperator
            end
            obj.Config = cfg;
            obj.ChannelOperator = channelOperator;
        end

        % SAMPLEPATHPARAMS 随机采样路径参数矩阵 [delay, doppler, gain]。
        function pathParams = samplePathParams(obj)
            localCfg   = obj.Config;
            numPaths   = localCfg.NumPaths;
            maxDelay   = localCfg.MaxDelaySamples;
            maxDopIdx  = localCfg.MaxDopplerIdx;

            allDelays = 0:maxDelay;
            selectedDelays = allDelays(randperm(length(allDelays), numPaths)).';
            randomPhases = -pi + 2 * pi * rand(numPaths, 1);
            dopplerShifts = maxDopIdx * cos(randomPhases);

            if ~localCfg.UseFractionalDoppler
                dopplerShifts = round(dopplerShifts);
            end

            % 为避免同一离散支撑重复，逐径检查 (delay, round(doppler)) 是否冲突。
            for pathIdx = 2:numPaths
                attempts = 0;
                while any(selectedDelays(1:pathIdx-1) == selectedDelays(pathIdx) & ...
                          round(dopplerShifts(1:pathIdx-1)) == round(dopplerShifts(pathIdx)))
                    randomPhases(pathIdx) = -pi + 2 * pi * rand();
                    dopplerShifts(pathIdx) = maxDopIdx * cos(randomPhases(pathIdx));
                    if ~localCfg.UseFractionalDoppler
                        dopplerShifts(pathIdx) = round(dopplerShifts(pathIdx));
                    end
                    attempts = attempts + 1;
                    if attempts > 100
                        break;
                    end
                end
            end

            complexGains = sqrt(1 / (2 * numPaths)) * ...
                (randn(numPaths, 1) + 1j * randn(numPaths, 1));

            pathParams = [selectedDelays, dopplerShifts, complexGains];
        end

        % SAMPLECHANNEL 采样路径并生成对应有效信道矩阵。
        function [pathParams, effectiveChannel] = sampleChannel(obj)
            pathParams = obj.samplePathParams();
            effectiveChannel = obj.ChannelOperator.buildEffectiveChannel(pathParams);
        end

    end

end


