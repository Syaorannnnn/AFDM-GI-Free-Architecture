function channelMatrixPhysical = buildTimeVaryingChannel( ...
        totalSubcarriers, pathDelays, pathDopplers, pathGains)
% BUILDTIMEVARYINGCHANNEL 构造离散时变物理信道矩阵。
%
% 输入:
%   totalSubcarriers - (double) 总子载波数量。
%   pathDelays       - (Px1 double) 各路径时延索引。
%   pathDopplers     - (Px1 double) 各路径 Doppler 索引。
%   pathGains        - (Px1 complex double) 各路径复增益。
%
% 输出:
%   channelMatrixPhysical - (NxN complex double) 物理信道矩阵。

    numPaths = length(pathDelays);

    % 构造循环移位矩阵基，用于模拟离散时延。
    baseVector = [zeros(1, totalSubcarriers - 1) 1];
    delayPermutationMatrix = toeplitz( ...
        [baseVector(1) fliplr(baseVector(2:end))], baseVector);
    channelMatrixPhysical = zeros(totalSubcarriers, totalSubcarriers);

    for pathIdx = 1:numPaths
        pathGain = pathGains(pathIdx);
        pathDelay = pathDelays(pathIdx);
        pathDoppler = pathDopplers(pathIdx);

        dopplerPhaseVec = exp(-1j * 2 * pi * pathDoppler * ...
            (0:totalSubcarriers - 1) / totalSubcarriers);
        dopplerMatrix = diag(dopplerPhaseVec);

        channelMatrixPhysical = channelMatrixPhysical + ...
            pathGain * dopplerMatrix * (delayPermutationMatrix ^ pathDelay);
    end
end
