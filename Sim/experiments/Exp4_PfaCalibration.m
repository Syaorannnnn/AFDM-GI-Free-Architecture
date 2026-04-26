function [empiricalPfa, ciLower, ciUpper] = Exp4_PfaCalibration(configObj, numMonteCarlo)
% EXP4_PFACALIBRATION 验证定理 3.1 的虚警率校准曲线。
%
%   描述:
%   在纯噪声条件下，对目标 Pfa 扫描并统计经验虚警率，输出 95% 置信区间，
%   同时保存图像到 `Sim/Results/Figures/Fig4_PfaCalibration.fig`。
%
%   语法:
%   [empiricalPfa, ciLower, ciUpper] = Exp4_PfaCalibration(configObj)
%   [empiricalPfa, ciLower, ciUpper] = Exp4_PfaCalibration(configObj, numMonteCarlo)
%
%   输入:
%   configObj - (GiFreeConfig) Pfa 校准实验使用的配置对象。
%   numMonteCarlo - (double) 每个目标 Pfa 的蒙特卡洛次数，默认 1e5。
%
%   输出:
%   empiricalPfa - (6x1 double) 实测经验虚警率。
%   ciLower - (6x1 double) 95% 置信区间下界。
%   ciUpper - (6x1 double) 95% 置信区间上界。
%
%   版本历史:
%   2026-04-19 - Aiden - 新增 Pfa 校准实验。

    arguments
        configObj (1,1) GiFreeConfig
        numMonteCarlo (1,1) double {mustBePositive, mustBeInteger} = 1e5
    end

    targetPfaVec = [1e-3; 3e-3; 1e-2; 3e-2; 1e-1; 3e-1];
    numTargets = numel(targetPfaVec);
    numSubcarriers = configObj.NumSubcarriers;
    pilotAmplitude = configObj.PilotAmplitude;
    maxDelaySamples = configObj.MaxDelaySamples;
    maxDopplerIdx = configObj.MaxDopplerIdx;
    numCandidates = (maxDelaySamples + 1) * (2 * maxDopplerIdx + 1);

    empiricalPfa = zeros(numTargets, 1);
    ciLower = zeros(numTargets, 1);
    ciUpper = zeros(numTargets, 1);

    for targetIdx = 1:numTargets
        targetFalseAlarmProb = targetPfaVec(targetIdx);
        falseAlarmCount = 0;
        totalDecisionCount = numMonteCarlo * numCandidates;

        for trialIdx = 1:numMonteCarlo
            noiseResidual = sqrt(0.5) * ...
                (randn(numSubcarriers, 1) + 1j * randn(numSubcarriers, 1));
            thresholdValue = paperCfarThreshold( ...
                noiseResidual, numSubcarriers, pilotAmplitude, ...
                maxDelaySamples, maxDopplerIdx, targetFalseAlarmProb);

            for delayVal = 0:maxDelaySamples
                for dopplerIdx = -maxDopplerIdx:maxDopplerIdx
                    locIndex = dopplerIdx + configObj.LocStep * delayVal;
                    responsePos0 = mod(configObj.PilotPos0 - locIndex, numSubcarriers);
                    responsePos1 = responsePos0 + 1;
                    metricValue = abs(pilotAmplitude * noiseResidual(responsePos1))^2;
                    if metricValue >= thresholdValue
                        falseAlarmCount = falseAlarmCount + 1;
                    end
                end
            end
        end

        empiricalPfa(targetIdx) = falseAlarmCount / totalDecisionCount;
        [ciLower(targetIdx), ciUpper(targetIdx)] = computeConfidenceInterval( ...
            falseAlarmCount, totalDecisionCount);
    end

    plotPfaCalibrationFigure(targetPfaVec, empiricalPfa, ciLower, ciUpper);
end

function plotPfaCalibrationFigure(targetPfaVec, empiricalPfa, ciLower, ciUpper)
% PLOTPFACALIBRATIONFIGURE 绘制并保存 Pfa 校准曲线。
%
% 输入:
%   targetPfaVec - (Nx1 double) 目标虚警率。
%   empiricalPfa - (Nx1 double) 实测虚警率。
%   ciLower - (Nx1 double) 置信区间下界。
%   ciUpper - (Nx1 double) 置信区间上界。

    simDir = fileparts(fileparts(mfilename('fullpath')));
    figuresDir = fullfile(simDir, 'Results', 'Figures');
    if ~isfolder(figuresDir)
        mkdir(figuresDir);
    end

    minPositiveValue = min(targetPfaVec) / 10;
    safeCiLower = max(ciLower, minPositiveValue);
    safeCiUpper = max(ciUpper, minPositiveValue);

    figHandle = figure('Position', [80 80 560 420], 'Color', 'w');
    loglog(targetPfaVec, empiricalPfa, '-o', ...
        'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, ...
        'MarkerSize', 7, 'MarkerFaceColor', [0.85 0.33 0.10], ...
        'DisplayName', 'Measured');
    hold on;
    patch([targetPfaVec; flipud(targetPfaVec)], ...
        [safeCiLower; flipud(safeCiUpper)], ...
        [0.85 0.33 0.10], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'DisplayName', '95% CI');
    loglog(targetPfaVec, targetPfaVec, 'k--', ...
        'LineWidth', 1.4, 'DisplayName', 'y=x');
    hold off;

    applyIeeeStyle(gca);
    xlabel('Target P_{fa}');
    ylabel('Empirical P_{fa}');
    title('Fig.4: Pfa Calibration');
    legend('Location', 'northwest');
    savefig(figHandle, fullfile(figuresDir, 'Fig4_PfaCalibration.fig'));
end
