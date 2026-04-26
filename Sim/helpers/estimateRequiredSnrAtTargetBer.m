function [requiredSnrDb, hasCrossing] = estimateRequiredSnrAtTargetBer( ...
        snrVecDb, berVec, targetBer)
% ESTIMATEREQUIREDSNRATTARGETBER 估计达到目标 BER 所需的 SNR。
%
%   描述:
%   在 BER 随 SNR 单调下降的前提下，使用 log10(BER) 线性插值提取
%   第一次达到目标 BER 的工作点。若扫描区间内没有交点，则返回 NaN。
%
%   语法:
%   [requiredSnrDb, hasCrossing] = estimateRequiredSnrAtTargetBer( ...
%       snrVecDb, berVec, targetBer)
%
%   输入:
%   snrVecDb - (Nx1 double) SNR 扫描向量。
%   berVec - (Nx1 double) 对应 BER 向量。
%   targetBer - (double) 目标 BER。
%
%   输出:
%   requiredSnrDb - (double) 达到目标 BER 所需的 SNR；若不存在则为 NaN。
%   hasCrossing - (logical) 扫描范围内是否存在交点。
%
%   版本历史:
%   2026-04-22 - 新增频谱效率权衡图工作点提取辅助函数。

    arguments
        snrVecDb (:,1) double
        berVec (:,1) double
        targetBer (1,1) double {mustBePositive, mustBeFinite}
    end

    if numel(snrVecDb) ~= numel(berVec)
        error('estimateRequiredSnrAtTargetBer:SizeMismatch', ...
            'snrVecDb 与 berVec 的长度必须一致。');
    end

    firstHitIdx = find(berVec <= targetBer, 1, 'first');
    if isempty(firstHitIdx)
        requiredSnrDb = NaN;
        hasCrossing = false;
        return;
    end

    hasCrossing = true;
    if firstHitIdx == 1
        requiredSnrDb = snrVecDb(1);
        return;
    end

    prevSnrDb = snrVecDb(firstHitIdx - 1);
    nextSnrDb = snrVecDb(firstHitIdx);
    prevLogBer = log10(max(berVec(firstHitIdx - 1), realmin));
    nextLogBer = log10(max(berVec(firstHitIdx), realmin));
    targetLogBer = log10(targetBer);

    if abs(nextLogBer - prevLogBer) < eps
        requiredSnrDb = nextSnrDb;
        return;
    end

    interpRatio = (targetLogBer - prevLogBer) / (nextLogBer - prevLogBer);
    interpRatio = min(max(interpRatio, 0), 1);
    requiredSnrDb = prevSnrDb + interpRatio * (nextSnrDb - prevSnrDb);
end
