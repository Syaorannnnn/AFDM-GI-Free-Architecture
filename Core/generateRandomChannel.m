function [delays, dopplers, gains] = generateRandomChannel( ...
        numPaths, maxDelay, maxDopplerIdx, isFractional)
% GENERATERANDOMCHANNEL 生成与主仿真口径一致的随机路径参数。
% 输入:
%   numPaths - 路径数。
%   maxDelay - 最大时延样本。
%   maxDopplerIdx - 最大 Doppler 索引。
%   isFractional - 是否生成分数 Doppler。
% 输出:
%   delays - 时延向量。
%   dopplers - Doppler 向量。
%   gains - 复增益向量。

    allDelays = (0:maxDelay).';
    delays = sort(allDelays(randperm(numel(allDelays), numPaths)));
    randomAngles = -pi + 2 * pi * rand(numPaths, 1);
    dopplers = maxDopplerIdx * cos(randomAngles);

    if ~isFractional
        dopplers = round(dopplers);
        for pathIdx = 2:numPaths
            attemptCount = 0;
            while any(delays(1:pathIdx-1) == delays(pathIdx) & ...
                    dopplers(1:pathIdx-1) == dopplers(pathIdx)) && attemptCount < 100
                randomAngles(pathIdx) = -pi + 2 * pi * rand();
                dopplers(pathIdx) = round(maxDopplerIdx * cos(randomAngles(pathIdx)));
                attemptCount = attemptCount + 1;
            end
        end
    end

    gains = sqrt(1 / (2 * numPaths)) * (randn(numPaths, 1) + 1j * randn(numPaths, 1));
    delays = delays(:);
    dopplers = dopplers(:);
    gains = gains(:);
end
