function trialCountVec = buildExp2TrialCountVec(snrVecDb)
% BUILDEXP2TRIALCOUNTVEC 返回 Exp2 各 SNR 点的 trial 数设置。
%
%   描述:
%   低 SNR 区间保持较高样本量；高 SNR 区间在原有基础上小幅增加，以便更稳
%   地观察 P=5 的 BER 退化趋势。
%
%   输入:
%   snrVecDb - (Nx1 / 1xN double) SNR 采样点，单位 dB。
%
%   输出:
%   trialCountVec - (Nx1 double) 与 snrVecDb 一一对应的 trial 数。

    arguments
        snrVecDb (:,1) double
    end

    trialCountVec = 2000 * ones(size(snrVecDb));
    trialCountVec(snrVecDb <= 5) = 5000;
    trialCountVec(snrVecDb >= 15) = 4000;
end
