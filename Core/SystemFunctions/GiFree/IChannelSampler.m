classdef (Abstract) IChannelSampler < handle
% ICHANNELSAMPLER - GI-Free 信道采样器接口
%
%   描述:
%   定义路径参数采样与整信道采样的统一入口。
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。

    methods (Abstract)
        % SAMPLEPATHPARAMS 采样多径参数 [delay, doppler, gain]。
        pathParams = samplePathParams(obj)
        % SAMPLECHANNEL 同时返回路径参数与有效信道矩阵。
        [pathParams, effectiveChannel] = sampleChannel(obj)
    end

end


