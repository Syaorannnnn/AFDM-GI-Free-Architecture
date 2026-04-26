function thresholdValue = paperCfarThreshold( ...
        rxSignal, numSubcarriers, pilotAmplitude, maxDelaySamples, maxDopplerIdx, targetFalseAlarmProb)
% PAPERCFARTHRESHOLD 计算相关搜索实验使用的 CFAR 度量门限。
%
%   描述:
%   返回与 OMP 检测度量同量纲的门限，即 residualPower × pilotAmplitude^2 × cfarCoeff。
%   若需要用于论文固定检测器的原始幅度比较，应再除以 pilotAmplitude^2 并开方。
%
%   语法:
%   thresholdValue = paperCfarThreshold(rxSignal, numSubcarriers, pilotAmplitude, ...
%       maxDelaySamples, maxDopplerIdx, targetFalseAlarmProb)
%
%   输入:
%   rxSignal - (Nx1 double) 当前接收信号或残差。
%   numSubcarriers - (double) 子载波数。
%   pilotAmplitude - (double) 导频幅度。
%   maxDelaySamples - (double) 最大时延索引。
%   maxDopplerIdx - (double) 最大整数 Doppler 索引。
%   targetFalseAlarmProb - (double) 目标虚警率。
%
%   输出:
%   thresholdValue - (double) 与 OMP 度量一致的门限。
%
%   版本历史:
%   2026-04-19 - Aiden - 新增相关搜索 CFAR 门限辅助函数。

    validateattributes(rxSignal, {'double'}, {'column'}, mfilename, 'rxSignal');
    validateattributes(numSubcarriers, {'numeric'}, {'scalar', 'integer', 'positive'}, ...
        mfilename, 'numSubcarriers');
    validateattributes(pilotAmplitude, {'numeric'}, {'scalar', 'positive'}, ...
        mfilename, 'pilotAmplitude');
    validateattributes(maxDelaySamples, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, ...
        mfilename, 'maxDelaySamples');
    validateattributes(maxDopplerIdx, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, ...
        mfilename, 'maxDopplerIdx');
    validateattributes(targetFalseAlarmProb, {'numeric'}, {'scalar', '>', 0, '<', 1}, ...
        mfilename, 'targetFalseAlarmProb');

    % 论文式 CFAR 系数：log(M/Pfa)，M 为 delay-doppler 搜索候选总数。
    numCandidates = (maxDelaySamples + 1) * (2 * maxDopplerIdx + 1);
    cfarCoeff = log(numCandidates / targetFalseAlarmProb);
    residualPower = real(rxSignal' * rxSignal) / numSubcarriers;
    thresholdValue = residualPower * pilotAmplitude^2 * cfarCoeff;
end
