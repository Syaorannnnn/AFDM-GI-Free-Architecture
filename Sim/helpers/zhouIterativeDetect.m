function [estPaths, hEffective] = zhouIterativeDetect( ...
    rxSignal, configObj, channelBuilder, fixedThreshold, maxPaths)
%ZHOUITERATIVEDETECT Zhou 论文的迭代 SIC + 固定门限路径检测。
%
%   执行真正的迭代 SIC：每轮在残差上找最强峰，若峰值超过固定门限则记录，
%   构建原子加入字典，做 LS 增益估计，更新残差，直到峰值低于门限或达到路径上限。
%
%   输入:
%   rxSignal        - (Nx1 double) 接收信号（导频帧）
%   configObj       - (GiFreeConfig) 配置对象
%   channelBuilder  - (IChannelOperator) 信道构造算子
%   fixedThreshold  - (double) 固定门限值（与峰值度量同量纲）
%   maxPaths        - (double) 路径上限
%
%   输出:
%   estPaths    - (Px3) [delay, dopplerIdx, complexGain]
%   hEffective  - (NxN sparse/full) 有效信道矩阵
%
%   版本历史:
%   2026-04-21 - Aiden - 初始版本。

%% 参数验证与初始化
validateattributes(rxSignal, {'numeric'}, {'column', 'nonempty'});
validateattributes(fixedThreshold, {'numeric'}, {'scalar', 'positive'});
validateattributes(maxPaths, {'numeric'}, {'scalar', 'integer', 'positive'});

N = configObj.NumSubcarriers;
locStep = configObj.LocStep;
maxDelay = configObj.MaxDelaySamples;
maxDoppler = configObj.MaxDopplerIdx;
pilotPos0 = configObj.PilotPos0;
pilotSeq = configObj.PilotSequence;
pilotAmp = configObj.PerPilotAmplitude;
chirpC1 = configObj.ChirpParam1;
chirpC2 = configObj.ChirpParam2;

residual = rxSignal;
dictionary = zeros(N, maxPaths);
pathBuffer = zeros(maxPaths, 3);
numDetected = 0;

%% 迭代 SIC 主循环
for iterIdx = 1:maxPaths
    % Grid search: 在残差上找最强峰（度量与 OMP scoreCandidate 一致）
    bestMetric = -inf;
    for delayVal = 0:maxDelay
        for dopplerIdx = -maxDoppler:maxDoppler
            % 位置与相位计算（与 GiFreeEstimator.scoreCandidate 一致）
            locIndex = dopplerIdx + locStep * delayVal;
            responseColIdx = mod(pilotPos0 - locIndex, N);
            centerIdx = responseColIdx + 1;
            phaseVal = exp(1j * 2 * pi * (chirpC1 * delayVal^2 ...
                - pilotPos0 * delayVal / N + chirpC2 * pilotPos0^2 ...
                - chirpC2 * responseColIdx^2));
            % 度量：原子内积功率（单位：pilotAmp^2 * residual^2）
            metric = abs(conj(pilotAmp * pilotSeq * phaseVal) * residual(centerIdx))^2;
            if metric > bestMetric
                bestMetric = metric;
                bestDelay = delayVal;
                bestDoppler = dopplerIdx;
            end
        end
    end

    % 门限判决：若峰值低于固定门限，终止迭代
    if bestMetric < fixedThreshold
        break;
    end

    % 记录路径位置并构建原子
    numDetected = numDetected + 1;
    pathBuffer(numDetected, 1:2) = [bestDelay, bestDoppler];
    dictionary(:, numDetected) = channelBuilder.buildCompositePilotResponse(bestDelay, bestDoppler);

    % LS 估计增益并更新残差（对原信号做联合估计）
    activeDict = dictionary(:, 1:numDetected);
    gains = activeDict \ rxSignal;
    residual = rxSignal - activeDict * gains;
end

%% 输出
estPaths = pathBuffer(1:numDetected, :);
if numDetected > 0
    estPaths(:, 3) = gains;
    hEffective = channelBuilder.buildEffectiveChannel(estPaths);
else
    hEffective = sparse(N, N);
end

end