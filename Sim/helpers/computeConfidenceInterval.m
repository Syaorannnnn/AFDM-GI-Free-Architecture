function [ciLower, ciUpper] = computeConfidenceInterval(numErrors, totalBits, alpha)
% COMPUTECONFIDENCEINTERVAL 计算二项分布误码率的置信区间。
%
%   描述:
%   优先使用 Clopper-Pearson 精确区间；当错误数为 0 或 totalBits 时，
%   使用 Rule-of-3 近似，避免数值边界退化。
%
%   语法:
%   [ciLower, ciUpper] = computeConfidenceInterval(numErrors, totalBits)
%   [ciLower, ciUpper] = computeConfidenceInterval(numErrors, totalBits, alpha)
%
%   输入:
%   numErrors - (double) 错误事件个数，范围 [0, totalBits]。
%   totalBits - (double) 总观测比特数，必须为正整数。
%   alpha - (double) 显著性水平，默认 0.05，对应 95% 置信区间。
%
%   输出:
%   ciLower - (double) 置信区间下界。
%   ciUpper - (double) 置信区间上界。
%
%   版本历史:
%   2026-04-19 - Aiden - 新增 Clopper-Pearson 置信区间计算。

    if nargin < 3 || isempty(alpha)
        alpha = 0.05;
    end

    validateattributes(numErrors, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative'}, mfilename, 'numErrors');
    validateattributes(totalBits, {'numeric'}, ...
        {'scalar', 'integer', 'positive'}, mfilename, 'totalBits');
    validateattributes(alpha, {'numeric'}, ...
        {'scalar', '>', 0, '<', 1}, mfilename, 'alpha');

    if numErrors > totalBits
        error('computeConfidenceInterval:InvalidCount', ...
            'numErrors (%d) 不能大于 totalBits (%d)。', numErrors, totalBits);
    end

    if numErrors == 0
        ciLower = 0;
        ciUpper = min(3 / totalBits, 1);
        return;
    end

    if numErrors == totalBits
        ciLower = max(1 - 3 / totalBits, 0);
        ciUpper = 1;
        return;
    end

    ciLower = betaincinv(alpha / 2, numErrors, totalBits - numErrors + 1);
    ciUpper = betaincinv(1 - alpha / 2, numErrors + 1, totalBits - numErrors);
end
