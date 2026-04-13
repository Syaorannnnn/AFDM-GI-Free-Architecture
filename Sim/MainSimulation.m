%% MAINSIMULATION - AFDM 全流程仿真主脚本（论文图复现）
%
% 描述:
% 本脚本统一生成 AFDM 系统核心对比图，覆盖以下实验部分：
% Part 1: AFDM vs OFDM（双色散信道）
% Part 2: CFAR 自适应门限 vs 固定门限
% Part 3: 分数 Doppler 场景下 EP / GI-Free 对比
%
% 依赖项:
% - MATLAB R2023b 或更高版本
% - Communications Toolbox
% - 仓库根目录已 addpath(genpath("./"))
%
% 如何运行:
% 1) 在 MATLAB 中将工作目录切换到仓库根目录。
% 2) 执行 MainSimulation。
%
% 输出:
% - 图文件: Sim/Results/Figures/*.fig
% - 日志:   Sim/Results/Logs/MainSimulation_*.txt
%
% 版本历史:
% 2026-04-01 - Aiden - 注释规范化。

clear; close all; clc;
addpath(genpath("./"));

fprintf('============================================================\n');
fprintf('     AFDM Full-Pipeline Simulation (Paper Figures)\n');
fprintf('============================================================\n\n');

simDir     = fileparts(mfilename('fullpath'));
figuresDir = fullfile(simDir, 'Results', 'Figures');
logsDir    = fullfile(simDir, 'Results', 'Logs');
if ~isfolder(figuresDir), mkdir(figuresDir); end
if ~isfolder(logsDir),    mkdir(logsDir);    end

% 日志镜像到 Logs，文件名附时间戳以保留历史记录。
simLogPath = fullfile(logsDir, sprintf('MainSimulation_%s.txt', ...
    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
diary(simLogPath);

%% ==================== 公共参数 ====================
numSc           = 512;      % DAFT 域子载波数
modOrder        = 4;        % QPSK
bps             = log2(modOrder);
maxDelay        = 4;        % 最大时延扩展样本数
maxDopplerIdx   = 2;        % 最大 Doppler 索引
defaultPaths    = 3;        % 默认路径数
pilotSnrEp      = 35;       % EP 导频 SNR (dB)
pilotSnrGf      = 45;       % GI-Free 导频 SNR (dB)

% 全局 SNR 扫描区间
snrVec = 0:5:20;
numSnr = numel(snrVec);

%% ==================== 配色/标记/样式 ====================
% 蓝色=基线/EP，红色=GI-Free，黑色虚线=理想参考，紫色=论文迭代
cBlue   = [0.00 0.45 0.74];    % 基线
cRed    = [0.85 0.33 0.10];    % GI-Free
cBlack  = [0.00 0.00 0.00];    % 理想参考
cPurple = [0.49 0.18 0.56];    % 论文迭代

%% ################################################################
%  Part 1: AFDM vs OFDM
%  ################################################################
fprintf('==================== Part 1: AFDM vs OFDM ====================\n');
fprintf('  AFDM: EP 单导频 + DAFT 域估计\n');
fprintf('  OFDM: 梳状导频 LS + DFT 插值 + 单抽头 MMSE (标准基线)\n\n');

dopplerGuard1 = 4;

% AFDM-EP 配置
cfgAfdm = EpConfig();
cfgAfdm.ModulationOrder  = modOrder;
cfgAfdm.BitsPerSymbol    = bps;
cfgAfdm.MaxDopplerIdx   = maxDopplerIdx;
cfgAfdm.DopplerGuard     = dopplerGuard1;
cfgAfdm.NumPaths         = defaultPaths;
cfgAfdm.PilotSnrDb         = pilotSnrEp;
cfgAfdm.MaxDelaySamples    = maxDelay;
cfgAfdm.TotalSubcarriers = numSc + maxDelay;   % 最后赋值触发 updateDerivedParams
sysAfdm = EpSystem(cfgAfdm);

% 标准 CP-OFDM 梳状导频基线
cfgOfdm = OfdmConfig();
cfgOfdm.ModulationOrder  = modOrder;  % BitsPerSymbol 由 Dependent 属性自动派生
cfgOfdm.MaxDelaySamples  = maxDelay;
cfgOfdm.NumPaths         = defaultPaths;
cfgOfdm.MaxDopplerIdx   = maxDopplerIdx;
cfgOfdm.CpLength         = maxDelay;
cfgOfdm.PilotSpacing     = floor(numSc / (maxDelay + 1));
cfgOfdm.DftSize          = numSc;
sysOfdm = OfdmSystem(cfgOfdm);

% 频谱效率对比
seAfdm = cfgAfdm.NumActiveCarriers * bps / cfgAfdm.TotalSubcarriers;
seOfdm = cfgOfdm.NumDataCarriers   * bps / cfgOfdm.TotalFrameLength;
fprintf('  AFDM-EP SE: %.3f bps/Hz  (数据载波 %d / 帧长 %d)\n', ...
    seAfdm, cfgAfdm.NumActiveCarriers, cfgAfdm.TotalSubcarriers);
fprintf('  OFDM   SE: %.3f bps/Hz  (数据载波 %d / 帧长 %d, 导频 %d 个, 间隔 %d)\n', ...
    seOfdm, cfgOfdm.NumDataCarriers, cfgOfdm.TotalFrameLength, ...
    cfgOfdm.NumPilots, cfgOfdm.PilotSpacing);
fprintf('\n');

nTrials1       = 1e3;
berAfdmEstVec  = zeros(numSnr, 1);
berAfdmPcsiVec = zeros(numSnr, 1);
berOfdmEstVec  = zeros(numSnr, 1);
berOfdmPcsiVec = zeros(numSnr, 1);

fprintf('SNR   BER_AFDM_Est  BER_AFDM_PCSI  BER_OFDM_Est  BER_OFDM_PCSI\n');
t1 = tic;

for si = 1:numSnr
    errAE = 0; bitsAE = 0;
    errAP = 0; bitsAP = 0;
    errOE = 0; bitsOE = 0;
    errOP = 0; bitsOP = 0;

    for tr = 1:nTrials1
        % 共享信道：分数 Doppler (Jakes)
        [dl, dp, gn] = generateRandomChannel(defaultPaths, maxDelay, maxDopplerIdx, true);

        rAE = sysAfdm.runTrial(snrVec(si), dl, dp, gn);
        rAP = sysAfdm.runTrialPerfectCsi(snrVec(si), dl, dp, gn);
        rOE = sysOfdm.runTrial(snrVec(si), dl, dp, gn);
        rOP = sysOfdm.runTrialPerfectCsi(snrVec(si), dl, dp, gn);

        errAE = errAE + rAE.bitErrors;  bitsAE = bitsAE + rAE.totalBits;
        errAP = errAP + rAP.bitErrors;  bitsAP = bitsAP + rAP.totalBits;
        errOE = errOE + rOE.bitErrors;  bitsOE = bitsOE + rOE.totalBits;
        errOP = errOP + rOP.bitErrors;  bitsOP = bitsOP + rOP.totalBits;
    end

    berAfdmEstVec(si)  = berFloor(errAE, bitsAE);
    berAfdmPcsiVec(si) = berFloor(errAP, bitsAP);
    berOfdmEstVec(si)  = berFloor(errOE, bitsOE);
    berOfdmPcsiVec(si) = berFloor(errOP, bitsOP);

    fprintf(' %2d   %.2e     %.2e      %.2e     %.2e\n', ...
        snrVec(si), berAfdmEstVec(si), berAfdmPcsiVec(si), ...
        berOfdmEstVec(si), berOfdmPcsiVec(si));
end
fprintf('Part 1 done (%.0f s)\n\n', toc(t1));

% ---- Fig.1 ----
fig1 = figure('Position', [80 80 560 420], 'Color', 'w');

% OFDM 梳状导频 LS（蓝色实线圆圈）
semilogy(snrVec, berOfdmEstVec, '-o', 'Color', cBlue, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cBlue, ...
    'DisplayName', 'OFDM Comb-LS (Est. CSI)');
hold on;

% OFDM Perfect CSI（蓝色虚线方块）
semilogy(snrVec, berOfdmPcsiVec, '--s', 'Color', cBlue, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'OFDM (Perfect CSI)');

% AFDM 浼拌 CSI (瀹炵嚎涓夎, 绾?
semilogy(snrVec, berAfdmEstVec, '-^', 'Color', cRed, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cRed, ...
    'DisplayName', 'AFDM-EP (Est. CSI)');

% AFDM Perfect CSI（黑色虚线菱形）
semilogy(snrVec, berAfdmPcsiVec, '--d', 'Color', cBlack, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'AFDM-EP (Perfect CSI)');

hold off;
applyIeeeStyle(gca);
ylim([1e-5 1]);
xlabel('Data SNR (dB)');
ylabel('BER');
title(sprintf('Fig.1: AFDM vs OFDM  N=%d P=%d Fractional Doppler', numSc, defaultPaths));
legend('Location', 'southwest');
savefig(fig1, fullfile(figuresDir, 'Fig1_AfdmVsOfdm.fig'));

%% ################################################################
%  Part 2: CFAR 自适应门限 vs 论文固定门限
%  ################################################################
fprintf('==================== Part 2: CFAR Threshold Benchmark ====================\n');

dopplerGuard3 = 4;

cfgGf3a = GiFreeConfig();
cfgGf3a.NumSubcarriers       = numSc;
cfgGf3a.ModulationOrder      = modOrder;
cfgGf3a.MaxDelaySamples       = maxDelay;
cfgGf3a.MaxDopplerIdx      = maxDopplerIdx;
cfgGf3a.NumPaths             = 5;
cfgGf3a.DopplerGuard         = dopplerGuard3;
cfgGf3a.DirichletRadius      = 0;
cfgGf3a.PilotSnrDb           = pilotSnrGf;
cfgGf3a.MaxSicIterations     = 10;
cfgGf3a.NumPathsUpper        = 7;
cfgGf3a.UseFractionalDoppler = false;
sysGf3a = GiFreeSystem(cfgGf3a);

% 缓存常用参数
pilotAmp  = cfgGf3a.PilotAmplitude;
c1Val     = cfgGf3a.ChirpParam1;
c2Val     = cfgGf3a.ChirpParam2;
locStep   = cfgGf3a.LocStep;
dataPos   = cfgGf3a.DataPos1;   % 1-based
nData3a   = cfgGf3a.NumDataSymbols;

pilotFrame3a    = zeros(numSc, 1);
pilotFrame3a(1) = cfgGf3a.PerPilotAmplitude * cfgGf3a.PilotSequence;

snrDb3     = 15;
snrLin3    = 10^(snrDb3/10);
paperThr   = 3 * sqrt(1 + snrLin3);      % 论文 eq.(17)
ompReg     = (nData3a/numSc) * snrLin3 / pilotAmp^2;
lmmseReg   = 1 / snrLin3;
nTrials3a  = 5e2;
pathCounts = 1:5;

% 统计量 [5 x 3]（论文首检、论文迭代、CFAR）
avgDet  = zeros(5,3);  avgPd   = zeros(5,3);
avgFa   = zeros(5,3);  avgNmse = zeros(5,3);
avgBer  = zeros(5,3);

t3a = tic;

for gi = 1:5
    nP = pathCounts(gi);
    cfgGf3a.NumPaths = nP;

    accDet = zeros(1,3); accPd = zeros(1,3); accFa = zeros(1,3);
    accNmse = zeros(1,3); accBerBits = zeros(1,3); accTotBits = 0;

    for tr = 1:nTrials3a
        % ---- 生成共享信号 ----
        [txFrame, txIdx] = sysGf3a.Transmitter.transmit(snrLin3);
        [trueParams, trueH] = sysGf3a.ChannelSampler.sampleChannel();
        nTrue = size(trueParams, 1);
        noise = sqrt(0.5) * (randn(numSc,1) + 1j*randn(numSc,1));
        rxSig = trueH * txFrame + noise;
        trueNorm = max(norm(full(trueH),'fro')^2, 1e-20);
        nBits = numel(txIdx) * bps;

        % (A) 论文首检: gamma = 3*sqrt(N0+Es)
        detA = paperDetect(rxSig, numSc, c1Val, c2Val, pilotAmp, locStep, maxDelay, maxDopplerIdx, paperThr);
        [cA, fA] = matchPaths(trueParams, detA);
        hA = buildHeff(sysGf3a.ChannelBuilder, detA);

        % (B) 论文迭代一次
        if size(detA,1) > 0
            cleanSig = rxSig - hA * pilotFrame3a;
            estSig = GiFreeReceiver.lmmseDetect(cleanSig, hA, lmmseReg, numSc, dataPos);
            dFrame = zeros(numSc,1);
            dFrame(dataPos) = estSig(dataPos);
            cleanPilot = rxSig - hA * dFrame;
            nDetA = size(detA,1);
            refinedThr = 3 * sqrt(1 + max((7-nDetA)/7, 0) * snrLin3);
            detB = paperDetect(cleanPilot, numSc, c1Val, c2Val, pilotAmp, locStep, maxDelay, maxDopplerIdx, refinedThr);
        else
            detB = detA;
        end
        [cB, fB] = matchPaths(trueParams, detB);
        hB = buildHeff(sysGf3a.ChannelBuilder, detB);

        % (C) OMP + CFAR
        [detC, hC] = sysGf3a.Estimator.estimateByOmp(rxSig, ompReg);
        [cC, fC] = matchPaths(trueParams, detC);

        % NMSE
        nA = norm(full(hA-trueH),'fro')^2/trueNorm;
        nB = norm(full(hB-trueH),'fro')^2/trueNorm;
        nC = norm(full(hC-trueH),'fro')^2/trueNorm;

        % BER（三种方法统一使用 LMMSE 检测，隔离门限贡献）
        bA = quickBer(rxSig, hA, pilotFrame3a, lmmseReg, numSc, dataPos, snrLin3, modOrder, txIdx);
        bB = quickBer(rxSig, hB, pilotFrame3a, lmmseReg, numSc, dataPos, snrLin3, modOrder, txIdx);
        bC = quickBer(rxSig, hC, pilotFrame3a, lmmseReg, numSc, dataPos, snrLin3, modOrder, txIdx);

        accDet  = accDet  + [size(detA,1), size(detB,1), size(detC,1)];
        accPd   = accPd   + [cA/nTrue, cB/nTrue, cC/nTrue];
        accFa   = accFa   + [fA, fB, fC];
        accNmse = accNmse + [nA, nB, nC];
        accBerBits = accBerBits + [bA, bB, bC]*nBits;
        accTotBits = accTotBits + nBits;
    end

    avgDet(gi,:)  = accDet/nTrials3a;
    avgPd(gi,:)   = accPd/nTrials3a;
    avgFa(gi,:)   = accFa/nTrials3a;
    avgNmse(gi,:) = accNmse/nTrials3a;
    avgBer(gi,:)  = accBerBits/accTotBits;

    fprintf('  P=%d: Pd=[%.0f%% %.0f%% %.0f%%]  NMSE=[%+.1f %+.1f %+.1f]dB  BER=[%.4f %.4f %.4f]\n', ...
        nP, avgPd(gi,:)*100, 10*log10(avgNmse(gi,:)+1e-20), avgBer(gi,:));
end
fprintf('Part 2 done (%.0f s)\n\n', toc(t3a));

% ---- Fig.2: 四子图 ----
legLbl  = {'Paper (first)', 'Paper (iter)', 'OMP+CFAR'};
pltClr  = {cBlue, cPurple, cRed};
pltMkr  = {'o', 's', '^'};
pltLine = {'-', '--', '-'};

fig3 = figure('Position', [60 60 1200 380], 'Color', 'w');

% (a) 检测数
subplot(1,4,1); hold on; grid on; box on;
for mi = 1:3
    plot(pathCounts, avgDet(:,mi), [pltLine{mi} pltMkr{mi}], ...
        'Color', pltClr{mi}, 'LineWidth', 1.6, 'MarkerSize', 6, ...
        'MarkerFaceColor', pltClr{mi});
end
plot(pathCounts, pathCounts, 'k:', 'LineWidth', 1);
applyIeeeStyle(gca, 9);
ylabel('Avg. detected paths'); xlabel('True path count P');
title('(a) Detection count');
legend([legLbl, {'Ideal'}], 'Location', 'northwest', 'FontSize', 7);

% (b) 虚警数
subplot(1,4,2); hold on; grid on; box on;
for mi = 1:3
    plot(pathCounts, avgFa(:,mi), [pltLine{mi} pltMkr{mi}], ...
        'Color', pltClr{mi}, 'LineWidth', 1.6, 'MarkerSize', 6, ...
        'MarkerFaceColor', pltClr{mi});
end
applyIeeeStyle(gca, 9);
ylabel('Avg. false alarms'); xlabel('P');
title('(b) False alarms');

% (c) NMSE
subplot(1,4,3); hold on; grid on; box on;
for mi = 1:3
    plot(pathCounts, 10*log10(avgNmse(:,mi)+1e-20), [pltLine{mi} pltMkr{mi}], ...
        'Color', pltClr{mi}, 'LineWidth', 1.6, 'MarkerSize', 6, ...
        'MarkerFaceColor', pltClr{mi});
end
applyIeeeStyle(gca, 9);
ylabel('NMSE (dB)'); xlabel('P');
title('(c) Channel estimation NMSE');

% (d) BER
subplot(1,4,4); hold on; grid on; box on;
for mi = 1:3
    plot(pathCounts, avgBer(:,mi), [pltLine{mi} pltMkr{mi}], ...
        'Color', pltClr{mi}, 'LineWidth', 1.6, 'MarkerSize', 6, ...
        'MarkerFaceColor', pltClr{mi});
end
applyIeeeStyle(gca, 9);
ylabel('BER'); xlabel('P');
title('(d) BER');

sgtitle(sprintf('Fig.2: CFAR Threshold Comparison  N=%d  DataSNR=%ddB  PilotSNR=%ddB', ...
    numSc, snrDb3, pilotSnrGf), 'FontSize', 12, 'FontName', 'Times New Roman');
savefig(fig3, fullfile(figuresDir, 'Fig2_CfarThreshold.fig'));

%% ################################################################
%  Part 3: 分数 Doppler：EP 保护开销 vs GI-Free 对比
%  ################################################################
fprintf('==================== Part 3: Fractional Doppler Focus ====================\n');

cfgInt = GiFreeConfig();
cfgInt.NumSubcarriers       = numSc;
cfgInt.ModulationOrder      = modOrder;
cfgInt.MaxDelaySamples       = maxDelay;
cfgInt.MaxDopplerIdx      = maxDopplerIdx;
cfgInt.NumPaths             = defaultPaths;
cfgInt.DopplerGuard         = dopplerGuard3;
cfgInt.DirichletRadius      = 0;        % 不建模 Dirichlet 半展宽
cfgInt.PilotSnrDb           = pilotSnrGf;
cfgInt.MaxSicIterations     = 10;
cfgInt.NumPathsUpper        = 6;
cfgInt.UseFractionalDoppler = true;      % 信道始终为分数 Doppler
cfgInt.UseDynamicPilot      = true;
cfgInt.DynamicPilotBaseDb   = 35;
cfgInt.EnableProgressiveCfar = false;
cfgInt.EnablePathStabilityGate = false;
sysInt = GiFreeSystem(cfgInt);

cfgFracBase = GiFreeConfig();
cfgFracBase.NumSubcarriers       = numSc;
cfgFracBase.ModulationOrder      = modOrder;
cfgFracBase.MaxDelaySamples      = maxDelay;
cfgFracBase.MaxDopplerIdx        = maxDopplerIdx;
cfgFracBase.NumPaths             = defaultPaths;
cfgFracBase.DopplerGuard         = dopplerGuard3;
cfgFracBase.DirichletRadius      = 4;       % 原始分数 Doppler 接收机
cfgFracBase.PilotSnrDb           = pilotSnrGf;
cfgFracBase.MaxSicIterations     = 10;
cfgFracBase.NumPathsUpper        = 6;
cfgFracBase.UseFractionalDoppler = true;
cfgFracBase.UseDynamicPilot      = false;
cfgFracBase.EnableProgressiveCfar = false;
cfgFracBase.EnablePathStabilityGate = false;
sysFracBase = GiFreeSystem(cfgFracBase);

cfgFrac = GiFreeConfig();
cfgFrac.NumSubcarriers       = numSc;
cfgFrac.ModulationOrder      = modOrder;
cfgFrac.MaxDelaySamples      = maxDelay;
cfgFrac.MaxDopplerIdx        = maxDopplerIdx;
cfgFrac.NumPaths             = defaultPaths;
cfgFrac.DopplerGuard         = dopplerGuard3;
cfgFrac.DirichletRadius      = 4;       % Dirichlet 核建模 + T1P 接收机优化
cfgFrac.PilotSnrDb           = pilotSnrGf;
cfgFrac.MaxSicIterations     = 10;
cfgFrac.NumPathsUpper        = 6;
cfgFrac.UseFractionalDoppler = true;
cfgFrac.UseDynamicPilot          = false;  % 固定导频模式（回退决策）
cfgFrac.DynamicPilotBaseDb       = 35;     % 备用参数（动态模式禁用时无效）
cfgFrac.EnableProgressiveCfar    = true;
cfgFrac.ProgressiveCfarInitScale = 3.0;
cfgFrac.ProgressiveCfarFinalScale = 1.0;
cfgFrac.EnablePathStabilityGate  = true;
cfgFrac.PathStabilityThreshold   = 2;
sysFrac = GiFreeSystem(cfgFrac);  % GI-Free 固定导频 + T1P 接收机（推荐配置）

% ---- EP 基准线：分数 Doppler 下不同保护带 ----
xiNuVec = [1, 2, 4];
nXiNu   = numel(xiNuVec);
sysEp3b = cell(nXiNu, 1);
seEp3b  = zeros(nXiNu, 1);

for xi = 1:nXiNu
    cfgEp = EpConfig();
    cfgEp.ModulationOrder  = modOrder;
    cfgEp.BitsPerSymbol    = bps;
    cfgEp.MaxDopplerIdx   = maxDopplerIdx;
    cfgEp.DopplerGuard     = xiNuVec(xi);
    cfgEp.NumPaths         = defaultPaths;
    cfgEp.PilotSnrDb         = pilotSnrEp;
    cfgEp.MaxDelaySamples    = maxDelay;
    cfgEp.TotalSubcarriers = numSc + maxDelay;
    sysEp3b{xi} = EpSystem(cfgEp);
    seEp3b(xi)  = cfgEp.NumActiveCarriers * bps / cfgEp.NumDataSubcarriers;
end

nTrials3b = 5e2;
berInt  = zeros(numSnr,1);  berFracBase  = zeros(numSnr,1);  berFrac  = zeros(numSnr,1);
berPcsi3b = zeros(numSnr,1);
nmseInt = zeros(numSnr,1);  nmseFracBase = zeros(numSnr,1);  nmseFrac = zeros(numSnr,1);
berEp3b = zeros(numSnr, nXiNu);

seGfInt  = cfgInt.NumDataSymbols * bps / numSc;
seGfFracBase = cfgFracBase.NumDataSymbols * bps / numSc;
seGfFrac = cfgFrac.NumDataSymbols * bps / numSc;
for xi = 1:nXiNu
    cfg = sysEp3b{xi}.Config;
    fprintf('  EP xi_nu=%d: %d data / %d sc -> SE = %.3f bits/s/Hz  (Q=%d)\n', ...
        xiNuVec(xi), cfg.NumActiveCarriers, cfg.NumDataSubcarriers, seEp3b(xi), ...
        cfg.ZeroPaddingLength * 2 + 1);
end
fprintf('  GI-Free Base: %d data / %d sc -> SE = %.3f bits/s/Hz\n', ...
    cfgFracBase.NumDataSymbols, numSc, seGfFracBase);
fprintf('  GI-Free T1P:  %d data / %d sc -> SE = %.3f bits/s/Hz\n', ...
    cfgFrac.NumDataSymbols, numSc, seGfFrac);
fprintf('\n');

fprintf(['SNR   BER_Int     BER_Base    BER_T1P     BER_PCSI    ', ...
         'NMSE_Int   NMSE_Base  NMSE_T1P   BER_EP1    BER_EP2    BER_EP4\n']);
t3b = tic;

for si = 1:numSnr
    snrLin = 10^(snrVec(si)/10);
    errInt = 0; errFracBase = 0; errFrac = 0; errPcsi = 0; totBits = 0;
    accNmseInt = 0; accNmseFracBase = 0; accNmseFrac = 0;
    errEp3b = zeros(1, nXiNu); bitsEp3b = zeros(1, nXiNu);

    for tr = 1:nTrials3b
        % 共享信道参数：GI-Free 与 EP 使用相同的 [delay, doppler, gain]
        trialRngState = rng;
        [txFrame, txIdx] = transmitWithFixedSeed(sysFrac, snrLin, trialRngState);
        [txFrameBase, txIdxBase] = transmitWithFixedSeed(sysFracBase, snrLin, trialRngState);
        if cfgInt.UseDynamicPilot
            cfgInt.CurrentDataSnrLin = snrLin;
        end
        [pathParams, trueH] = sysFrac.ChannelSampler.sampleChannel();
        dlShared = pathParams(:, 1);
        dpShared = pathParams(:, 2);
        gnShared = pathParams(:, 3);
        noise = sqrt(0.5) * (randn(numSc,1) + 1j*randn(numSc,1));
        rxSig = trueH * txFrame + noise;
        rxSigBase = trueH * txFrameBase + noise;
        nBits = numel(txIdx) * bps;

        % 整数接收机
        [detI, nmI, ~] = sysInt.Receiver.receive(rxSig, trueH, snrLin, 1);

        % 旧版分数接收机
        [detB, nmB, ~] = sysFracBase.Receiver.receive(rxSigBase, trueH, snrLin, 1);

        % 分数接收机
        [detF, nmF, ~] = sysFrac.Receiver.receive(rxSig, trueH, snrLin, 1);

        % Perfect CSI
        detP = GiFreeSystem.perfectCsiDetect(rxSig, trueH, snrLin, 1, cfgFrac);

        % EP 基准线（使用共享信道参数，实现严格公平对比）
        for xi = 1:nXiNu
            rEp = sysEp3b{xi}.runTrial(snrVec(si), dlShared, dpShared, gnShared);
            errEp3b(xi)  = errEp3b(xi)  + rEp.bitErrors;
            bitsEp3b(xi) = bitsEp3b(xi) + rEp.totalBits;
        end

        % BER 统计（使用 biterr）
        [neI, ~] = biterr(txIdx, detI, bps);
        [neB, ~] = biterr(txIdxBase, detB, bps);
        [neF, ~] = biterr(txIdx, detF, bps);
        [neP, ~] = biterr(txIdx, detP, bps);

        errInt  = errInt  + neI;
        errFracBase = errFracBase + neB;
        errFrac = errFrac + neF;
        errPcsi = errPcsi + neP;
        totBits = totBits + nBits;
        accNmseInt = accNmseInt + nmI;
        accNmseFracBase = accNmseFracBase + nmB;
        accNmseFrac = accNmseFrac + nmF;
    end

    berInt(si)    = berFloor(errInt,  totBits);
    berFracBase(si) = berFloor(errFracBase, totBits);
    berFrac(si)   = berFloor(errFrac, totBits);
    berPcsi3b(si) = berFloor(errPcsi, totBits);
    nmseInt(si)   = accNmseInt  / nTrials3b;
    nmseFracBase(si) = accNmseFracBase / nTrials3b;
    nmseFrac(si)  = accNmseFrac / nTrials3b;
    for xi = 1:nXiNu
        berEp3b(si, xi) = berFloor(errEp3b(xi), bitsEp3b(xi));
    end

    fprintf(' %2d   %.2e    %.2e    %.2e    %.2e    %+.1f      %+.1f      %+.1f      %.2e   %.2e   %.2e\n', ...
        snrVec(si), berInt(si), berFracBase(si), berFrac(si), berPcsi3b(si), ...
        10*log10(nmseInt(si)+1e-20), 10*log10(nmseFracBase(si)+1e-20), ...
        10*log10(nmseFrac(si)+1e-20), ...
        berEp3b(si,1), berEp3b(si,2), berEp3b(si,3));
end
fprintf('Part 3 done (%.0f s)\n\n', toc(t3b));

% ---- Fig.3: 三子图 ----
% EP 配色：浅蓝 -> 深蓝（随 xi_nu 递增）
cEpLight = [0.39 0.67 0.82];
cEpMid   = [0.00 0.45 0.74];
cEpDark  = [0.00 0.27 0.53];
cEp3     = {cEpLight, cEpMid, cEpDark};
mkEp3    = {'o', 's', 'd'};
lnEp3    = {'--', '-.', '-'};

fig4 = figure('Position', [60 60 1400 400], 'Color', 'w');

% (a) BER vs SNR
subplot(1,3,1); hold on; grid on; box on;
for xi = 1:nXiNu
    semilogy(snrVec, berEp3b(:,xi), [lnEp3{xi} mkEp3{xi}], ...
        'Color', cEp3{xi}, 'LineWidth', 1.6, 'MarkerSize', 6, ...
        'MarkerFaceColor', cEp3{xi}, ...
        'DisplayName', sprintf('EP \\xi_\\nu=%d (\\eta=%.2f)', xiNuVec(xi), seEp3b(xi)));
end
semilogy(snrVec, berInt, '-o', 'Color', cBlue, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cBlue, ...
    'DisplayName', sprintf('GI-Free IntRx (\\eta=%.2f)', seGfInt));
semilogy(snrVec, berFracBase, '-d', 'Color', cPurple, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cPurple, ...
    'DisplayName', sprintf('GI-Free FracRx Base (\\eta=%.2f)', seGfFracBase));
semilogy(snrVec, berFrac, '-^', 'Color', cRed, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cRed, ...
    'DisplayName', sprintf('GI-Free FracRx T1P (\\eta=%.2f)', seGfFrac));
semilogy(snrVec, berPcsi3b, '--s', 'Color', cBlack, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'Perfect CSI');
applyIeeeStyle(gca);
set(gca, 'YScale', 'log');
ylim([1e-4 1]);
xlabel('Data SNR (dB)'); ylabel('BER');
title('(a) BER vs SNR');
legend('Location', 'southwest', 'FontSize', 7);

% (b) NMSE vs SNR
subplot(1,3,2); hold on; grid on; box on;
plot(snrVec, 10*log10(nmseInt+1e-20), '-o', 'Color', cBlue, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cBlue, ...
    'DisplayName', 'GI-Free IntRx');
plot(snrVec, 10*log10(nmseFracBase+1e-20), '-d', 'Color', cPurple, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cPurple, ...
    'DisplayName', 'GI-Free FracRx Base');
plot(snrVec, 10*log10(nmseFrac+1e-20), '-^', 'Color', cRed, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cRed, ...
    'DisplayName', 'GI-Free FracRx T1P');
applyIeeeStyle(gca);
xlabel('Data SNR (dB)'); ylabel('NMSE (dB)');
title('(b) Channel Estimation NMSE');
legend('Location', 'northeast');

% (c) SE 柱状图
subplot(1,3,3); hold on; grid on; box on;
barLabels = cell(1, nXiNu + 1);
barValues = zeros(1, nXiNu + 1);
barColors = [cell2mat(cEp3'); cRed];
for xi = 1:nXiNu
    barLabels{xi} = sprintf('EP \\xi_\\nu=%d', xiNuVec(xi));
    barValues(xi) = seEp3b(xi);
end
barLabels{end} = 'GI-Free T1P';
barValues(end) = seGfFrac;
bh = bar(barValues, 0.6, 'FaceColor', 'flat');
for bi = 1:numel(barValues)
    bh.CData(bi,:) = barColors(bi,:);
end
set(gca, 'XTickLabel', barLabels, 'XTick', 1:numel(barValues));
applyIeeeStyle(gca, 9);
ylabel('SE (bits/s/Hz)');
title('(c) Spectral Efficiency');

sgtitle(sprintf('Fig.3: Fractional Doppler  N=%d P=%d k_{max}=%d  EP Overhead vs GI-Free', ...
    numSc, defaultPaths, maxDopplerIdx), ...
    'FontSize', 12, 'FontName', 'Times New Roman');
savefig(fig4, fullfile(figuresDir, 'Fig3_FractionalDoppler.fig'));

%% ==================== 汇总 ====================
fprintf('============================================================\n');
fprintf('  All simulations complete. 3 figures generated.\n');
fprintf('  Fig.1: AFDM vs OFDM  (Fractional Doppler)\n');
fprintf('  Fig.2: CFAR Threshold (Variable path count)\n');
fprintf('  Fig.3: Fractional Doppler (EP xi_nu=1/2/4 vs GI-Free Base/T1P)\n');
fprintf('  Figures -> Sim/Results/Figures/\n');
fprintf('  Log    -> Sim/Results/Logs/');
fprintf('============================================================\n');
diary off;

%% ==================== 辅助函数 ====================

function [txFrame, txDataIndices] = transmitWithFixedSeed(systemObj, dataSnrLin, rngState)
% TRANSMITWITHFIXEDSEED 用固定随机态生成一帧，便于跨配置共享同一组数据。
% 输入:
%   systemObj - GI-Free 系统对象。
%   dataSnrLin - 数据 SNR 线性值。
%   rngState - trial 开始时保存的随机态。
% 输出:
%   txFrame - 发射帧。
%   txDataIndices - 发送数据索引。

    previousState = rng;
    cleanupObj = onCleanup(@() rng(previousState));
    rng(rngState);

    if systemObj.Config.UseDynamicPilot
        systemObj.Config.CurrentDataSnrLin = dataSnrLin;
    end
    [txFrame, txDataIndices] = systemObj.Transmitter.transmit(dataSnrLin);
    clear cleanupObj;
end

function detected = paperDetect(rxSig, numSc, c1, c2, amp, locStep, maxDelay, maxDop, thr)
% PAPERDETECT 论文固定门限路径检测实现 [Zhou et al. 2024, eq.(17)-(18)]
    maxPossible = (maxDelay + 1) * (2 * maxDop + 1);
    buf   = zeros(maxPossible, 3);
    count = 0;
    for l = 0:maxDelay
        for k = -maxDop:maxDop
            pos = mod(-(k + locStep*l), numSc);
            if abs(rxSig(pos+1)) >= thr
                phase        = exp(1j*2*pi*(c1*l^2 - c2*pos^2));
                gain         = rxSig(pos+1) / (phase * amp);
                count        = count + 1;
                buf(count, :) = [l, k, gain];
            end
        end
    end
    detected = buf(1:count, :);
end

function hEff = buildHeff(builder, pathParams)
% BUILDHEFF 从路径参数构造有效信道矩阵
    numSc = builder.Config.NumSubcarriers;
    if isempty(pathParams)
        hEff = sparse(numSc, numSc);
    else
        hEff = builder.buildEffectiveChannel(pathParams);
    end
end

function [nCorrect, nFalse] = matchPaths(trueP, estP)
% MATCHPATHS 路径匹配统计（按 delay + round(doppler)）
    nT = size(trueP, 1);  nE = size(estP, 1);
    matched = false(nT, 1);  falseAlarm = true(nE, 1);
    for ei = 1:nE
        for ti = 1:nT
            if ~matched(ti) && estP(ei,1)==trueP(ti,1) && round(estP(ei,2))==round(trueP(ti,2))
                matched(ti) = true;
                falseAlarm(ei) = false;
                break;
            end
        end
    end
    nCorrect = sum(matched);
    nFalse   = sum(falseAlarm);
end

function val = quickBer(rxSig, hEst, pilotFrame, regParam, numSc, dataPos, snrLin, modOrd, txIdx)
% QUICKBER 计算一次 LMMSE 检测 BER（用于快速诊断）
    cleanSig = rxSig - hEst * pilotFrame;
    estSig   = GiFreeReceiver.lmmseDetect(cleanSig, hEst, regParam, numSc, dataPos);
    rxIdx    = qamdemod(estSig(dataPos)/sqrt(snrLin), modOrd, 'UnitAveragePower', true);
    [nErr, ~] = biterr(txIdx, rxIdx, log2(modOrd));
    val = nErr / (numel(txIdx) * log2(modOrd));
end


