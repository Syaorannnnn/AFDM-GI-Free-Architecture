% MainGiFreeSim.m — GI-Free AFDM 最终系统仿真
%
% 对比两条核心曲线:
%   Perfect CSI    —— 已知真实信道的理论性能上界
%   GiFreeSystem   —— 最终版 Turbo DD-CE 接收机
%
% 最终版接收机的核心技术栈:
%   Phase 1  : 先验干扰感知正则化 + 正则化 OMP
%   Phase 2  : DD 多普勒/增益精炼 + 截断均值正则化 + 输出 SINR 软符号
%   PostDec  : 硬判决驱动的 DD 增益精炼 + 紧凑截断均值
%
% 图表配色: 深色莫兰迪 (Morandi) 色系, 纯白背景

clear; clc; close all; rng(42);

%% ======================== 1. 仿真参数 ========================

numSubcarriers = 256;
modulationOrder = 4;
maxDelaySpread = 2;
maxDopplerIndex = 2;
numPaths = 3;
pilotSnrDb = 45;
noisePower = 1;

dopplerGuard = 4;
spreadWidth = 4;

snrDbList = [0, 5, 10, 15, 20];
numTrials = 1e3;
numSnrPts = length(snrDbList);

totalSicIter = 10; % Phase1(3) + Phase2(6) + PostDecision(1)

% 构造配置
cfg = GiFreeConfig();
cfg.NumSubcarriers = numSubcarriers;
cfg.ModulationOrder = modulationOrder;
cfg.MaxDelaySpread = maxDelaySpread;
cfg.MaxDopplerIndex = maxDopplerIndex;
cfg.NumPaths = numPaths;
cfg.PilotSnrDb = pilotSnrDb;
cfg.MaxSicIterations = totalSicIter;
cfg.NumPathsUpper = numPaths + 1;
cfg.DopplerGuard = dopplerGuard;
cfg.SpreadWidth = spreadWidth;
cfg.UseFractionalDoppler = true;

% 结果存储: 列 1 = Perfect CSI, 列 2 = GiFreeSystem
berResults = zeros(numSnrPts, 2);
mseResults = zeros(numSnrPts, 2);

fprintf('============================================================\n');
fprintf('  GI-Free AFDM — Final System Simulation\n');
fprintf('  N = %d, %d-QAM, P = %d, PilotSNR = %ddB, SIC = %d iter\n', ...
    numSubcarriers, modulationOrder, numPaths, pilotSnrDb, totalSicIter);
fprintf('============================================================\n');

%% ======================== 2. 蒙特卡洛主循环 ========================

for snrIdx = 1:numSnrPts
    currentSnrDb = snrDbList(snrIdx);
    dataSnrLinear = 10 ^ (currentSnrDb / 10);
    bitsPerSymbol = log2(modulationOrder);

    cumBitErrSys = 0;
    cumBitErrRef = 0;
    cumTotalBits = 0;
    cumMse = 0;

    for trialIdx = 1:numTrials
        result = GiFreeSystem.runTrial(cfg, dataSnrLinear, noisePower);

        cumBitErrSys = cumBitErrSys + result.bitErrorsSys;
        cumBitErrRef = cumBitErrRef + result.bitErrorsRef;
        cumTotalBits = cumTotalBits + result.totalBits;
        cumMse = cumMse + result.mseSystem;
    end

    berResults(snrIdx, 1) = cumBitErrRef / max(cumTotalBits, 1);
    berResults(snrIdx, 2) = cumBitErrSys / max(cumTotalBits, 1);
    mseResults(snrIdx, 2) = cumMse / numTrials;

    berGap = berResults(snrIdx, 2) / max(berResults(snrIdx, 1), 1e-8);
    fprintf('SNR = %2d dB | BER: PerfCSI = %.2e  System = %.2e  (Gap = %.1f x) | NMSE = %.2e\n', ...
        currentSnrDb, berResults(snrIdx, 1), berResults(snrIdx, 2), ...
        berGap, mseResults(snrIdx, 2));
end

%% ======================== 3. 结果可视化 ========================

plotFinalResults(snrDbList, berResults, mseResults, cfg);
savefig(gcf, 'GiFreeFinalResults.fig');

%% ======================== 辅助函数: 绘图 ========================

function plotFinalResults(snrDbList, berResults, mseResults, cfg)

    colPerfCsi = [0.32, 0.30, 0.28]; % 深炭灰 (Charcoal)
    colSystem = [0.33, 0.47, 0.38]; % 深鼠尾草 (Deep Sage)
    colAccent = [0.55, 0.38, 0.40]; % 深灰玫瑰 (Dusty Mauve)
    colGridLine = [0.82, 0.80, 0.78]; % 浅暖灰 (Warm Silver)
    colFillArea = [0.85, 0.90, 0.86]; % 浅鼠尾草填充 (Light Sage Fill)
    colBarFace = [0.40, 0.55, 0.45]; % 中鼠尾草 (Medium Sage)
    colBarEdge = [0.28, 0.40, 0.32]; % 暗鼠尾草边框
    colDustyBlue = [0.40, 0.45, 0.58]; % 深灰蓝 (Dusty Steel)

    % 全局字体
    fontLabel = 11;
    fontTitle = 12;
    fontLeg = 9;

    numSc = cfg.NumSubcarriers;
    modOrder = cfg.ModulationOrder;
    numP = cfg.NumPaths;
    dopGuard = cfg.DopplerGuard;
    spreadKv = cfg.SpreadWidth;

    figure('Position', [80, 60, 1200, 900], 'Color', [1, 1, 1], ...
        'Name', 'GI-Free AFDM — Final System Performance');

    % ---- 子图 1: BER 核心对比 ----
    subplot(2, 2, 1);
    set(gca, 'Color', [1, 1, 1]);

    semilogy(snrDbList, berResults(:, 1), '-^', 'Color', colPerfCsi, ...
        'LineWidth', 2.2, 'MarkerSize', 9, 'MarkerFaceColor', colPerfCsi, ...
        'DisplayName', 'Perfect CSI'); hold on;
    semilogy(snrDbList, berResults(:, 2), '-s', 'Color', colSystem, ...
        'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colBarFace, ...
        'DisplayName', 'GiFreeSystem');

    grid on;
    set(gca, 'GridColor', colGridLine, 'GridAlpha', 0.6, 'MinorGridAlpha', 0.3);
    xlabel('SNR (dB)', 'FontSize', fontLabel);
    ylabel('BER', 'FontSize', fontLabel);
    title(sprintf('BER — N=%d, %d-QAM, P=%d', numSc, modOrder, numP), ...
        'FontSize', fontTitle);
    legend('Location', 'southwest', 'FontSize', fontLeg, 'Box', 'off');
    hold off;

    % ---- 子图 2: NMSE ----
    subplot(2, 2, 2);
    set(gca, 'Color', [1, 1, 1]);

    semilogy(snrDbList, mseResults(:, 2), '-s', 'Color', colSystem, ...
        'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colBarFace, ...
        'DisplayName', 'GiFreeSystem');

    grid on;
    set(gca, 'GridColor', colGridLine, 'GridAlpha', 0.6);
    xlabel('SNR (dB)', 'FontSize', fontLabel);
    ylabel('Normalized MSE', 'FontSize', fontLabel);
    title('Channel Estimation NMSE', 'FontSize', fontTitle);
    legend('Location', 'best', 'FontSize', fontLeg, 'Box', 'off');

    % ---- 子图 3: BER 差距分析 ----
    subplot(2, 2, 3);
    set(gca, 'Color', [1, 1, 1]);

    perfBer = max(berResults(:, 1), 1e-8);
    berGap = berResults(:, 2) ./ perfBer;

    % 填充区域: 差距可视化
    fill([snrDbList, fliplr(snrDbList)], ...
        [berGap', ones(1, length(snrDbList))], ...
        colFillArea, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off'); hold on;
    plot(snrDbList, berGap, '-s', 'Color', colSystem, ...
        'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colBarFace, ...
        'DisplayName', 'BER_{System} / BER_{PerfCSI}');
    yline(1, ':', 'Color', colGridLine, 'LineWidth', 1.8, ...
        'DisplayName', 'Perfect CSI = 1');

    grid on;
    set(gca, 'GridColor', colGridLine, 'GridAlpha', 0.6);
    xlabel('SNR (dB)', 'FontSize', fontLabel);
    ylabel('BER / BER_{PerfCSI}', 'FontSize', fontLabel);
    title('BER Gap to Perfect CSI', 'FontSize', fontTitle);
    legend('Location', 'best', 'FontSize', fontLeg, 'Box', 'off');
    hold off;

    % ---- 子图 4: 系统参数总结 (柱状图 + 注释) ----
    subplot(2, 2, 4);
    set(gca, 'Color', [1, 1, 1]);

    % 在选定 SNR 点的 BER 对比柱状图
    selectedSnrs = [5, 10, 15, 20, 25];
    selIdx = zeros(size(selectedSnrs));

    for i = 1:length(selectedSnrs)
        [~, selIdx(i)] = min(abs(snrDbList - selectedSnrs(i)));
    end

    berForBar = zeros(length(selectedSnrs), 2);

    for i = 1:length(selectedSnrs)
        berForBar(i, :) = berResults(selIdx(i), :);
    end

    barHandle = bar(1:length(selectedSnrs), berForBar, 'grouped');
    set(gca, 'YScale', 'log');
    barHandle(1).FaceColor = colPerfCsi;
    barHandle(1).EdgeColor = colPerfCsi;
    barHandle(2).FaceColor = colBarFace;
    barHandle(2).EdgeColor = colBarEdge;

    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%d dB', x), ...
        selectedSnrs, 'UniformOutput', false));

    % 标注差距倍数
    for i = 1:length(selectedSnrs)
        gapVal = berForBar(i, 2) / max(berForBar(i, 1), 1e-8);

        if gapVal < 100 && berForBar(i, 2) > 0
            barCenterX = i + 0.15;
            barTopY = berForBar(i, 2) * 1.5;
            text(barCenterX, barTopY, sprintf('%.1fx', gapVal), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'FontWeight', 'bold', 'Color', colAccent);
        end

    end

    grid on;
    set(gca, 'GridColor', colGridLine, 'GridAlpha', 0.6);
    xlabel('SNR', 'FontSize', fontLabel);
    ylabel('BER', 'FontSize', fontLabel);
    title('BER at Selected SNR Points', 'FontSize', fontTitle);
    legend('Perfect CSI', 'GiFreeSystem', ...
        'Location', 'northeast', 'FontSize', fontLeg, 'Box', 'off');

    % 全局标题
    sgtitle(sprintf( ...
        'GI-Free AFDM Final System — \\xi_\\nu=%d, k_\\nu=%d, Turbo DD-CE', ...
        dopGuard, spreadKv), 'FontSize', 14, 'FontWeight', 'bold', 'Color', colPerfCsi);
end
