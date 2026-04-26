function gapSummary = summarizeExp4SnrGap(resultStruct, targetBerVec)
% SUMMARIZEEXP4SNRGAP 计算 ProGuard-CLIP 相对 EP 的工作点 SNR 差值。
%
%   语法:
%   gapSummary = summarizeExp4SnrGap(resultStruct, targetBerVec)
%
%   输入:
%   resultStruct - (struct) Exp4 结果结构体，至少包含 snrVecDb、berEp、berAblation。
%   targetBerVec - (Nx1 double) 需要评估的目标 BER 向量。
%
%   输出:
%   gapSummary - (table) 每个目标 BER 下，ProGuard 与各 EP 曲线的工作点 SNR
%                以及对应差值（ProGuard - EP）。

    arguments
        resultStruct (1,1) struct
        targetBerVec (:,1) double {mustBePositive, mustBeFinite}
    end

    requiredFields = {'snrVecDb', 'berEp', 'berAblation'};
    for fieldIdx = 1:numel(requiredFields)
        if ~isfield(resultStruct, requiredFields{fieldIdx})
            error('summarizeExp4SnrGap:MissingField', ...
                'resultStruct 缺少字段 %s。', requiredFields{fieldIdx});
        end
    end

    % 修复: 支持 1 列或多列 berAblation，由 proGuardColumnIndex 指定列。
    % 若调用方未提供 proGuardColumnIndex，默认取最后一列（向后兼容）。
    numAblationCols = size(resultStruct.berAblation, 2);
    if isfield(resultStruct, 'proGuardColumnIndex')
        colIdx = resultStruct.proGuardColumnIndex;
    else
        colIdx = numAblationCols;
    end
    if colIdx < 1 || colIdx > numAblationCols
        error('summarizeExp4SnrGap:BadColumnIndex', ...
            'proGuardColumnIndex=%d 超出 berAblation 列数 %d。', ...
            colIdx, numAblationCols);
    end
    proGuardBer = resultStruct.berAblation(:, colIdx);
    gapCell = cell(numel(targetBerVec), 9);
    for targetIdx = 1:numel(targetBerVec)
        targetBer = targetBerVec(targetIdx);
        [proGuardSnrDb, proGuardHit] = estimateRequiredSnrAtTargetBer( ...
            resultStruct.snrVecDb(:), proGuardBer(:), targetBer);

        gapCell{targetIdx, 1} = targetBer;
        gapCell{targetIdx, 2} = proGuardSnrDb;
        gapCell{targetIdx, 3} = proGuardHit;

        for epIdx = 1:3
            [epSnrDb, epHit] = estimateRequiredSnrAtTargetBer( ...
                resultStruct.snrVecDb(:), resultStruct.berEp(:, epIdx), targetBer);
            baseColIdx = 3 + (epIdx - 1) * 2;
            gapCell{targetIdx, baseColIdx + 1} = epSnrDb;
            gapCell{targetIdx, baseColIdx + 2} = proGuardSnrDb - epSnrDb;
            if ~(proGuardHit && epHit)
                gapCell{targetIdx, baseColIdx + 2} = NaN;
            end
        end
    end

    gapSummary = cell2table(gapCell, 'VariableNames', { ...
        'TargetBer', ...
        'ProGuardRequiredSnrDb', 'ProGuardHasCrossing', ...
        'Ep1RequiredSnrDb', 'GapVsEp1Db', ...
        'Ep2RequiredSnrDb', 'GapVsEp2Db', ...
        'Ep4RequiredSnrDb', 'GapVsEp4Db'});
end
