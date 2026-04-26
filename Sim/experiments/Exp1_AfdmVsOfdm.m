function resultStruct = Exp1_AfdmVsOfdm( ...
        numSubcarriers, modulationOrder, bitsPerSymbol, maxDelaySamples, ...
        maxDopplerIdx, dopplerGuard, numPaths, pilotSnrDb, snrVecDb, figuresDir, trialCountOverride)
% EXP1_AFDMVSOFDM 执行 AFDM-EP 与 CP-OFDM 的 BER 基线对比。
%
%   描述:
%   该实验复用 Embedded Pilot AFDM 与 CP-OFDM 基线系统，在共享整数
%   Doppler 双色散信道下统计估计 CSI 与 Perfect CSI 两条链路的 BER。
%   该对照用于给出波形层面的基线，不承担保护带开销或 GI-Free 接收机归因。
%
%   输入:
%   numSubcarriers     - (double) 有效子载波数。
%   modulationOrder    - (double) QAM 调制阶数。
%   bitsPerSymbol      - (double) 每符号比特数。
%   maxDelaySamples    - (double) 最大时延样本。
%   maxDopplerIdx      - (double) 最大 Doppler 索引。
%   dopplerGuard       - (double) AFDM Doppler 保护参数。
%   numPaths           - (double) 信道路径数。
%   pilotSnrDb         - (double) AFDM-EP 导频 SNR。
%   snrVecDb           - (double) 数据 SNR 向量。
%   figuresDir         - (char/string) 图像输出目录。
%   trialCountOverride - (double) 可选，覆盖默认 trial 数。
%
%   输出:
%   resultStruct - (struct) BER、误比特数、总比特数、频谱效率与配置摘要。
%
%   版本历史:
%   2026-04-26 - Aiden - 恢复为 Part 1 主干实验函数。

    if nargin < 11 || isempty(trialCountOverride)
        trialCount = 1000;
    else
        trialCount = trialCountOverride;
    end

    snrVecDb = snrVecDb(:);
    afdmSystem = createAfdmEpSystem( ...
        numSubcarriers, modulationOrder, bitsPerSymbol, maxDelaySamples, ...
        maxDopplerIdx, dopplerGuard, numPaths, pilotSnrDb);
    ofdmSystem = createOfdmSystem( ...
        numSubcarriers, modulationOrder, maxDelaySamples, maxDopplerIdx, numPaths);

    resultStruct = initializeExp1Result(snrVecDb, trialCount, afdmSystem.Config, ofdmSystem.Config);

    fprintf('Exp1_AfdmVsOfdm: AFDM-EP vs CP-OFDM BER baseline\n');
    fprintf('  AFDM-EP SE: %.3f bps/Hz, OFDM SE: %.3f bps/Hz\n', ...
        resultStruct.spectralEfficiencyAfdm, resultStruct.spectralEfficiencyOfdm);
    fprintf('SNR   BER_AFDM_Est   BER_AFDM_PCSI   BER_OFDM_Est   BER_OFDM_PCSI\n');

    for snrIdx = 1:numel(snrVecDb)
        pointStats = runSingleAfdmOfdmPoint( ...
            afdmSystem, ofdmSystem, snrVecDb(snrIdx), ...
            trialCount, numPaths, maxDelaySamples, maxDopplerIdx);
        resultStruct = writeExp1Point(resultStruct, snrIdx, pointStats);
        fprintf(' %2d   %.2e      %.2e       %.2e      %.2e\n', ...
            snrVecDb(snrIdx), resultStruct.berAfdmEstimated(snrIdx), ...
            resultStruct.berAfdmPerfectCsi(snrIdx), resultStruct.berOfdmEstimated(snrIdx), ...
            resultStruct.berOfdmPerfectCsi(snrIdx));
    end

    plotAfdmVsOfdm(resultStruct, figuresDir);
end

function systemObj = createAfdmEpSystem( ...
        numSubcarriers, modulationOrder, bitsPerSymbol, maxDelaySamples, ...
        maxDopplerIdx, dopplerGuard, numPaths, pilotSnrDb)
% CREATEAFDMEPSYSTEM 构造 AFDM Embedded Pilot 基线系统。

    configObj = EpConfig();
    configObj.ModulationOrder = modulationOrder;
    configObj.BitsPerSymbol = bitsPerSymbol;
    configObj.MaxDopplerIdx = maxDopplerIdx;
    configObj.DopplerGuard = dopplerGuard;
    configObj.NumPaths = numPaths;
    configObj.PilotSnrDb = pilotSnrDb;
    configObj.MaxDelaySamples = maxDelaySamples;
    configObj.TotalSubcarriers = numSubcarriers + maxDelaySamples;
    systemObj = EpSystem(configObj);
end

function systemObj = createOfdmSystem( ...
        numSubcarriers, modulationOrder, maxDelaySamples, maxDopplerIdx, numPaths)
% CREATEOFDMSYSTEM 构造 CP-OFDM 梳状导频基线系统。

    configObj = OfdmConfig();
    configObj.ModulationOrder = modulationOrder;
    configObj.MaxDelaySamples = maxDelaySamples;
    configObj.NumPaths = numPaths;
    configObj.MaxDopplerIdx = maxDopplerIdx;
    configObj.CpLength = maxDelaySamples;
    configObj.PilotSpacing = floor(numSubcarriers / (maxDelaySamples + 1));
    configObj.DftSize = numSubcarriers;
    systemObj = OfdmSystem(configObj);
end

function resultStruct = initializeExp1Result(snrVecDb, trialCount, afdmConfig, ofdmConfig)
% INITIALIZEEXP1RESULT 初始化 Exp1 结果结构。

    numSnrPoints = numel(snrVecDb);
    resultStruct.snrVecDb = snrVecDb;
    resultStruct.numTrials = trialCount;
    resultStruct.berAfdmEstimated = zeros(numSnrPoints, 1);
    resultStruct.berAfdmPerfectCsi = zeros(numSnrPoints, 1);
    resultStruct.berOfdmEstimated = zeros(numSnrPoints, 1);
    resultStruct.berOfdmPerfectCsi = zeros(numSnrPoints, 1);
    resultStruct.bitErrorsAfdmEstimated = zeros(numSnrPoints, 1);
    resultStruct.bitErrorsAfdmPerfectCsi = zeros(numSnrPoints, 1);
    resultStruct.bitErrorsOfdmEstimated = zeros(numSnrPoints, 1);
    resultStruct.bitErrorsOfdmPerfectCsi = zeros(numSnrPoints, 1);
    resultStruct.totalBitsAfdm = zeros(numSnrPoints, 1);
    resultStruct.totalBitsOfdm = zeros(numSnrPoints, 1);
    resultStruct.spectralEfficiencyAfdm = ...
        afdmConfig.NumActiveCarriers * afdmConfig.BitsPerSymbol / afdmConfig.TotalSubcarriers;
    resultStruct.spectralEfficiencyOfdm = ...
        ofdmConfig.NumDataCarriers * ofdmConfig.BitsPerSymbol / ofdmConfig.TotalFrameLength;
end

function pointStats = runSingleAfdmOfdmPoint( ...
        afdmSystem, ofdmSystem, dataSnrDb, trialCount, ...
        numPaths, maxDelaySamples, maxDopplerIdx)
% RUNSINGLEAFDMOFDPOINT 统计一个 SNR 点的四条 BER 曲线。

    errorCounts = zeros(1, 4);
    totalBitsAfdm = 0;
    totalBitsOfdm = 0;

    for trialIdx = 1:trialCount
        [pathDelays, pathDopplers, pathGains] = generateRandomChannel( ...
            numPaths, maxDelaySamples, maxDopplerIdx, false);

        afdmEstimated = afdmSystem.runTrial(dataSnrDb, pathDelays, pathDopplers, pathGains);
        afdmPerfect = afdmSystem.runTrialPerfectCsi(dataSnrDb, pathDelays, pathDopplers, pathGains);
        ofdmEstimated = ofdmSystem.runTrial(dataSnrDb, pathDelays, pathDopplers, pathGains);
        ofdmPerfect = ofdmSystem.runTrialPerfectCsi(dataSnrDb, pathDelays, pathDopplers, pathGains);

        errorCounts = errorCounts + [afdmEstimated.bitErrors, afdmPerfect.bitErrors, ...
            ofdmEstimated.bitErrors, ofdmPerfect.bitErrors];
        totalBitsAfdm = totalBitsAfdm + afdmEstimated.totalBits;
        totalBitsOfdm = totalBitsOfdm + ofdmEstimated.totalBits;
        validatePerfectCsiBitCounts(afdmEstimated, afdmPerfect, ofdmEstimated, ofdmPerfect);
    end

    pointStats.ber = [ ...
        berFloor(errorCounts(1), totalBitsAfdm), berFloor(errorCounts(2), totalBitsAfdm), ...
        berFloor(errorCounts(3), totalBitsOfdm), berFloor(errorCounts(4), totalBitsOfdm)];
    pointStats.errorCounts = errorCounts;
    pointStats.totalBitsAfdm = totalBitsAfdm;
    pointStats.totalBitsOfdm = totalBitsOfdm;
end

function validatePerfectCsiBitCounts(afdmEstimated, afdmPerfect, ofdmEstimated, ofdmPerfect)
% VALIDATEPERFECTCSIBITCOUNTS 校验估计 CSI 与 Perfect CSI 的统计口径一致。

    if afdmEstimated.totalBits ~= afdmPerfect.totalBits
        error('Exp1_AfdmVsOfdm:AfdmBitCountMismatch', ...
            'AFDM estimated CSI and Perfect CSI trials returned different bit counts.');
    end
    if ofdmEstimated.totalBits ~= ofdmPerfect.totalBits
        error('Exp1_AfdmVsOfdm:OfdmBitCountMismatch', ...
            'OFDM estimated CSI and Perfect CSI trials returned different bit counts.');
    end
end

function resultStruct = writeExp1Point(resultStruct, snrIdx, pointStats)
% WRITEEXP1POINT 写入一个 SNR 点的统计结果。

    resultStruct.berAfdmEstimated(snrIdx) = pointStats.ber(1);
    resultStruct.berAfdmPerfectCsi(snrIdx) = pointStats.ber(2);
    resultStruct.berOfdmEstimated(snrIdx) = pointStats.ber(3);
    resultStruct.berOfdmPerfectCsi(snrIdx) = pointStats.ber(4);
    resultStruct.bitErrorsAfdmEstimated(snrIdx) = pointStats.errorCounts(1);
    resultStruct.bitErrorsAfdmPerfectCsi(snrIdx) = pointStats.errorCounts(2);
    resultStruct.bitErrorsOfdmEstimated(snrIdx) = pointStats.errorCounts(3);
    resultStruct.bitErrorsOfdmPerfectCsi(snrIdx) = pointStats.errorCounts(4);
    resultStruct.totalBitsAfdm(snrIdx) = pointStats.totalBitsAfdm;
    resultStruct.totalBitsOfdm(snrIdx) = pointStats.totalBitsOfdm;
end

function plotAfdmVsOfdm(resultStruct, figuresDir)
% PLOTAFDMVSOFDM 绘制 AFDM-EP 与 CP-OFDM 的 BER 曲线。

    plotStyle = PlotStyle();
    figureHandle = plotStyle.newFigure(640, 430);
    axHandle = axes('Parent', figureHandle);
    hold(axHandle, 'on');

    semilogy(axHandle, resultStruct.snrVecDb, resultStruct.berOfdmEstimated, '-o', ...
        'Color', PlotStyle.Blue, 'LineWidth', 1.8, 'MarkerSize', 7, ...
        'MarkerFaceColor', PlotStyle.Blue, 'DisplayName', 'OFDM Comb-LS');
    semilogy(axHandle, resultStruct.snrVecDb, resultStruct.berOfdmPerfectCsi, '--s', ...
        'Color', PlotStyle.Blue, 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', 'OFDM Perfect CSI');
    semilogy(axHandle, resultStruct.snrVecDb, resultStruct.berAfdmEstimated, '-^', ...
        'Color', PlotStyle.Red, 'LineWidth', 1.8, 'MarkerSize', 7, ...
        'MarkerFaceColor', PlotStyle.Red, 'DisplayName', 'AFDM-EP Estimated CSI');
    semilogy(axHandle, resultStruct.snrVecDb, resultStruct.berAfdmPerfectCsi, '--d', ...
        'Color', PlotStyle.Black, 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', 'AFDM-EP Perfect CSI');

    hold(axHandle, 'off');
    plotStyle.apply(axHandle, 10);
    set(axHandle, 'YScale', 'log');
    ylim(axHandle, [1e-5, 1]);
    xlabel(axHandle, 'Data SNR (dB)');
    ylabel(axHandle, 'BER');
    title(axHandle, 'Fig.1: AFDM-EP vs CP-OFDM');
    legend(axHandle, 'Location', 'southwest', 'FontSize', 8);
    savefig(figureHandle, fullfile(figuresDir, 'Fig1_AfdmVsOfdm.fig'));
end
