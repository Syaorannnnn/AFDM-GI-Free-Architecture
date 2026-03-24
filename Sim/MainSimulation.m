%% MainSimulation.m
%
%  叙事主线: AFDM 能否同时实现高频谱效率和高可靠性?
%
%  Part 1 — AFDM vs OFDM (分数多普勒, EP 架构, Perfect CSI 对比)
%  Part 2 — EP vs GI-Free (整数多普勒, 频谱效率对比)
%  Part 3a — CFAR 自适应门限 vs 论文固定门限 (路径数不确定)
%  Part 3b — 分数多普勒: 整数接收机 vs 分数接收机 (error floor)
%  Part 4 — Agile-c2 PAPR CCDF (完整收发链路)
%
%  配色方案 (对齐 GI-Free 论文 2404.01088 Fig.5 风格):
%    蓝色实线圆圈  — 基线/论文方法
%    红色实线三角  — 本文方法
%    黑色虚线方块  — 理想参考 (Perfect CSI / 理想值)
%    紫色虚线方块  — 论文迭代 (仅 Part 3a)

clear all; close all; clc;
addpath(genpath("./"));

fprintf('============================================================\n');
fprintf('     AFDM Full-Pipeline Simulation (Paper Figures)\n');
fprintf('============================================================\n\n');

% --- 确保输出目录存在 ---
if ~isfolder('./Sim/Results')
    mkdir('./Sim/Results');
end

%% ==================== 公共参数 ====================
numSc           = 512;      % DAFT 域子载波数
modOrder        = 4;        % QPSK
bps             = log2(modOrder);
maxDelay        = 4;        % 最大时延扩展 (支持最多 5 条路径)
maxDopplerIdx   = 2;        % 最大多普勒指数
defaultPaths    = 3;        % 默认路径数
pilotSnrEp      = 35;       % EP 导频 SNR (dB)
pilotSnrGf      = 45;       % GI-Free 导频 SNR (dB)

% 全局 SNR 扫描
snrVec = 0:5:20;
numSnr = numel(snrVec);

%% ==================== 配色/标记/样式====================
% 蓝圈=EP/论文, 红三角=GI-Free/本文, 黑虚方=理想
cBlue   = [0.00 0.45 0.74];    % 基线
cRed    = [0.85 0.33 0.10];    % 本文方法
cBlack  = [0.00 0.00 0.00];    % 理想参考
cPurple = [0.49 0.18 0.56];    % 论文迭代

%% ################################################################
%  Part 1: AFDM vs OFDM
%  ################################################################
fprintf('==================== Part 1: AFDM vs OFDM ====================\n');
fprintf('  Fractional Doppler, EP architecture\n');
fprintf('  Only difference: modulation basis (chirp vs sinusoid)\n\n');

dopplerGuard1 = 4;

cfgAfdm = EpConfig();
cfgAfdm.TotalSubcarriers = numSc + maxDelay;
cfgAfdm.MaxPathDelays    = maxDelay;
cfgAfdm.MaxNormDoppler   = maxDopplerIdx;
cfgAfdm.DopplerGuard     = dopplerGuard1;
cfgAfdm.NumPaths         = defaultPaths;
cfgAfdm.PilotSnr         = pilotSnrEp;
cfgAfdm.WaveformType     = "AFDM";

cfgOfdm = EpConfig();
cfgOfdm.TotalSubcarriers = numSc + maxDelay;
cfgOfdm.MaxPathDelays    = maxDelay;
cfgOfdm.MaxNormDoppler   = maxDopplerIdx;
cfgOfdm.DopplerGuard     = dopplerGuard1;
cfgOfdm.NumPaths         = defaultPaths;
cfgOfdm.PilotSnr         = pilotSnrEp;
cfgOfdm.WaveformType     = "OFDM";

sysAfdm = EpSystem(cfgAfdm);
sysOfdm = EpSystem(cfgOfdm);

nTrials1      = 500;
berAfdmVec    = zeros(numSnr, 1);
berOfdmVec    = zeros(numSnr, 1);
berPcsi1Vec   = zeros(numSnr, 1);

fprintf('SNR   BER_AFDM    BER_OFDM    BER_PCSI\n');
t1 = tic;

for si = 1:numSnr
    errA = 0; bitsA = 0;
    errO = 0; bitsO = 0;
    errP = 0; bitsP = 0;

    for tr = 1:nTrials1
        % 共享信道: 分数多普勒 (Jakes)
        [dl, dp, gn] = generateRandomChannel(defaultPaths, maxDelay, maxDopplerIdx, true);

        rA = sysAfdm.runTrial(snrVec(si), dl, dp, gn);
        rO = sysOfdm.runTrialPerfectCsi(snrVec(si), dl, dp, gn);
        rP = sysAfdm.runTrialPerfectCsi(snrVec(si), dl, dp, gn);

        errA  = errA  + rA.bitErrors;   bitsA = bitsA + rA.totalBits;
        errO  = errO  + rO.bitErrors;   bitsO = bitsO + rO.totalBits;
        errP  = errP  + rP.bitErrors;   bitsP = bitsP + rP.totalBits;
    end

    berAfdmVec(si)  = berFloor(errA, bitsA);
    berOfdmVec(si)  = berFloor(errO, bitsO);
    berPcsi1Vec(si) = berFloor(errP, bitsP);

    fprintf(' %2d   %.2e    %.2e    %.2e\n', ...
        snrVec(si), berAfdmVec(si), berOfdmVec(si), berPcsi1Vec(si));
end
fprintf('Part 1 done (%.0f s)\n\n', toc(t1));

% ---- Fig.1 ----
fig1 = figure('Position', [80 80 560 420], 'Color', 'w');
semilogy(snrVec, berOfdmVec, '-s', 'Color', cBlue, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cBlue, ...
    'DisplayName', 'OFDM-EP (Perfect CSI)');
hold on;
semilogy(snrVec, berAfdmVec, '-^', 'Color', cRed, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cRed, ...
    'DisplayName', 'AFDM-EP (Estimated CSI)');
semilogy(snrVec, berPcsi1Vec, '--o', 'Color', cBlack, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'AFDM-EP (Perfect CSI)');
hold off;
applyIeeeStyle(gca);
ylim([1e-5 1]);
xlabel('Data SNR (dB)');
ylabel('BER');
title(sprintf('Fig.1: AFDM vs OFDM  N=%d P=%d Fractional Doppler', numSc, defaultPaths));
legend('Location', 'southwest');
savefig(fig1, './Sim/Results/Fig1_AfdmVsOfdm.fig');

%% ################################################################
%  Part 2: EP vs GI-Free
%  ################################################################
fprintf('==================== Part 2: EP vs GI-Free ====================\n');
fprintf('  Integer Doppler (fair to EP), same c1\n');

dopplerGuard2 = 0;   % 整数多普勒不需要额外保护

cfgEp2 = EpConfig();
cfgEp2.TotalSubcarriers = numSc + maxDelay;
cfgEp2.MaxPathDelays    = maxDelay;
cfgEp2.MaxNormDoppler   = maxDopplerIdx;
cfgEp2.DopplerGuard     = dopplerGuard2;
cfgEp2.NumPaths         = defaultPaths;
cfgEp2.PilotSnr         = pilotSnrEp;
sysEp2 = EpSystem(cfgEp2);

cfgGf2 = GiFreeConfig();
cfgGf2.NumSubcarriers       = numSc;
cfgGf2.ModulationOrder      = modOrder;
cfgGf2.MaxDelaySpread       = maxDelay;
cfgGf2.MaxDopplerIndex      = maxDopplerIdx;
cfgGf2.NumPaths             = defaultPaths;
cfgGf2.DopplerGuard         = dopplerGuard2;
cfgGf2.SpreadWidth          = 0;
cfgGf2.PilotSnrDb           = pilotSnrGf;
cfgGf2.MaxSicIterations     = 10;
cfgGf2.NumPathsUpper        = 6;
cfgGf2.UseFractionalDoppler = false;
sysGf2 = GiFreeSystem(cfgGf2);

seEp = cfgEp2.NumActiveCarriers * bps / cfgEp2.NumDataSubcarriers;
seGf = cfgGf2.NumDataSymbols    * bps / numSc;
fprintf('  EP:      %d data / %d sc -> SE = %.3f bits/s/Hz\n', ...
    cfgEp2.NumActiveCarriers, cfgEp2.NumDataSubcarriers, seEp);
fprintf('  GI-Free: %d data / %d sc -> SE = %.3f bits/s/Hz\n', ...
    cfgGf2.NumDataSymbols, numSc, seGf);
fprintf('  SE gain: %.1f%%\n\n', (seGf/seEp - 1)*100);

% Part 2: 独立 SNR 范围 (避免高 SNR 下 PDR 崩塌)
snrVec2  = 0:2:20;
nSnr2    = numel(snrVec2);
nTrials2 = 1e3;

berEp2   = zeros(nSnr2, 1);
berGf2   = zeros(nSnr2, 1);
berPcsi2 = zeros(nSnr2, 1);

fprintf('SNR   BER_EP      BER_GF      BER_PCSI\n');
t2 = tic;

for si = 1:nSnr2
    snrLin = 10^(snrVec2(si)/10);
    errEp = 0; bitsEp = 0;
    errGf = 0; errPcsi = 0; bitsGf = 0;

    for tr = 1:nTrials2
        % EP: 整数多普勒
        [dl, dp, gn] = generateRandomChannel(defaultPaths, maxDelay, maxDopplerIdx, false);
        rE = sysEp2.runTrial(snrVec2(si), dl, dp, gn);
        errEp = errEp + rE.bitErrors;  bitsEp = bitsEp + rE.totalBits;

        % GI-Free: 独立信道实现, 同统计参数
        rG = sysGf2.runTrial(snrLin, 1);
        errGf   = errGf   + rG.bitErrorsSys;
        errPcsi = errPcsi + rG.bitErrorsRef;
        bitsGf  = bitsGf  + rG.totalBits;
    end

    berEp2(si)   = berFloor(errEp,   bitsEp);
    berGf2(si)   = berFloor(errGf,   bitsGf);
    berPcsi2(si) = berFloor(errPcsi, bitsGf);

    fprintf(' %2d   %.2e    %.2e    %.2e\n', snrVec2(si), berEp2(si), berGf2(si), berPcsi2(si));
end
fprintf('Part 2 done (%.0f s)\n\n', toc(t2));

% ---- Fig.2 ----
fig2 = figure('Position', [100 100 560 420], 'Color', 'w');
semilogy(snrVec2, berEp2, '-o', 'Color', cBlue, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cBlue, ...
    'DisplayName', sprintf('EP (\\eta=%.2f, %d data)', seEp, cfgEp2.NumActiveCarriers));
hold on;
semilogy(snrVec2, berGf2, '-^', 'Color', cRed, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cRed, ...
    'DisplayName', sprintf('GI-Free (\\eta=%.2f, %d data)', seGf, cfgGf2.NumDataSymbols));
semilogy(snrVec2, berPcsi2, '--s', 'Color', cBlack, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'Perfect CSI');
hold off;
applyIeeeStyle(gca);
ylim([1e-5 1]);
xlabel('Data SNR (dB)');
ylabel('BER');
title(sprintf('Fig.2: EP vs GI-Free  N=%d P=%d Integer Doppler', numSc, defaultPaths));
legend('Location', 'southwest');
savefig(fig2, './Sim/Results/Fig2_EpVsGiFree.fig');

%% ################################################################
%  Part 3a: CFAR 门限 vs 论文固定门限
%  ################################################################
fprintf('==================== Part 3a: CFAR Threshold Benchmark ====================\n');

dopplerGuard3 = 4;

cfgGf3a = GiFreeConfig();
cfgGf3a.NumSubcarriers       = numSc;
cfgGf3a.ModulationOrder      = modOrder;
cfgGf3a.MaxDelaySpread       = maxDelay;
cfgGf3a.MaxDopplerIndex      = maxDopplerIdx;
cfgGf3a.NumPaths             = 5;
cfgGf3a.DopplerGuard         = dopplerGuard3;
cfgGf3a.SpreadWidth          = 0;
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
dataPos   = cfgGf3a.DataPositions + 1;   % 1-based
nData3a   = cfgGf3a.NumDataSymbols;

pilotFrame3a    = zeros(numSc, 1);
pilotFrame3a(1) = cfgGf3a.PerPilotAmplitude * cfgGf3a.PilotSequence;

snrDb3     = 15;
snrLin3    = 10^(snrDb3/10);
paperThr   = 3 * sqrt(1 + snrLin3);      % 论文 eq.(17)
ompReg     = (nData3a/numSc) * snrLin3 / pilotAmp^2;
lmmseReg   = 1 / snrLin3;
nTrials3a  = 300;
pathCounts = 1:5;

% 统计量: [5 x 3] (论文首检, 论文迭代, CFAR)
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
        [trueParams, trueH] = sysGf3a.ChannelBuilder.generateChannel();
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

        % BER (三种统一 LMMSE 检测, 隔离门限贡献)
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
fprintf('Part 3a done (%.0f s)\n\n', toc(t3a));

% ---- Fig.3: 四子图 ----
legLbl  = {'Paper (first)', 'Paper (iter)', 'OMP+CFAR'};
pltClr  = {cBlue, cPurple, cRed};
pltMkr  = {'o', 's', '^'};
pltLine = {'-', '--', '-'};

fig3 = figure('Position', [60 60 1200 380], 'Color', 'w');

% (a) 检出数
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

sgtitle(sprintf('Fig.3: CFAR Threshold Comparison  N=%d  DataSNR=%ddB  PilotSNR=%ddB', ...
    numSc, snrDb3, pilotSnrGf), 'FontSize', 12, 'FontName', 'Times New Roman');
savefig(fig3, './Sim/Results/Fig3_CfarThreshold.fig');

%% ################################################################
%  Part 3b: 分数多普勒 — 整数接收机 vs 分数接收机
%  ################################################################
fprintf('==================== Part 3b: Fractional Doppler ====================\n');

cfgInt = GiFreeConfig();
cfgInt.NumSubcarriers       = numSc;
cfgInt.ModulationOrder      = modOrder;
cfgInt.MaxDelaySpread       = maxDelay;
cfgInt.MaxDopplerIndex      = maxDopplerIdx;
cfgInt.NumPaths             = defaultPaths;
cfgInt.DopplerGuard         = dopplerGuard3;
cfgInt.SpreadWidth          = 0;        % 不建模 Dirichlet 扩展
cfgInt.PilotSnrDb           = pilotSnrGf;
cfgInt.MaxSicIterations     = 10;
cfgInt.NumPathsUpper        = 6;
cfgInt.UseFractionalDoppler = true;      % 信道始终为分数多普勒
sysInt = GiFreeSystem(cfgInt);

cfgFrac = GiFreeConfig();
cfgFrac.NumSubcarriers       = numSc;
cfgFrac.ModulationOrder      = modOrder;
cfgFrac.MaxDelaySpread       = maxDelay;
cfgFrac.MaxDopplerIndex      = maxDopplerIdx;
cfgFrac.NumPaths             = defaultPaths;
cfgFrac.DopplerGuard         = dopplerGuard3;
cfgFrac.SpreadWidth          = 4;       % Dirichlet 核建模 + 割线法精炼
cfgFrac.PilotSnrDb           = pilotSnrGf;
cfgFrac.MaxSicIterations     = 10;
cfgFrac.NumPathsUpper        = 6;
cfgFrac.UseFractionalDoppler = true;
sysFrac = GiFreeSystem(cfgFrac);

nTrials3b = 80;
berInt  = zeros(numSnr,1);  berFrac  = zeros(numSnr,1);
berPcsi3b = zeros(numSnr,1);
nmseInt = zeros(numSnr,1);  nmseFrac = zeros(numSnr,1);

fprintf('SNR   BER_Int     BER_Frac    BER_PCSI    NMSE_Int   NMSE_Frac\n');
t3b = tic;

for si = 1:numSnr
    snrLin = 10^(snrVec(si)/10);
    errInt = 0; errFrac = 0; errPcsi = 0; totBits = 0;
    accNmseInt = 0; accNmseFrac = 0;

    for tr = 1:nTrials3b
        % 共享信号: 分数多普勒信道
        [txFrame, txIdx] = sysFrac.Transmitter.transmit(snrLin);
        [~, trueH] = sysFrac.ChannelBuilder.generateChannel();
        noise = sqrt(0.5) * (randn(numSc,1) + 1j*randn(numSc,1));
        rxSig = trueH * txFrame + noise;
        nBits = numel(txIdx) * bps;

        % 整数接收机
        [detI, nmI, ~] = sysInt.Receiver.receive(rxSig, trueH, snrLin, 1);

        % 分数接收机
        [detF, nmF, ~] = sysFrac.Receiver.receive(rxSig, trueH, snrLin, 1);

        % Perfect CSI
        detP = GiFreeSystem.perfectCsiDetect(rxSig, trueH, snrLin, 1, cfgFrac);

        % BER 统计 (使用 biterr)
        [neI, ~] = biterr(txIdx, detI, bps);
        [neF, ~] = biterr(txIdx, detF, bps);
        [neP, ~] = biterr(txIdx, detP, bps);

        errInt  = errInt  + neI;  errFrac = errFrac + neF;  errPcsi = errPcsi + neP;
        totBits = totBits + nBits;
        accNmseInt = accNmseInt + nmI;  accNmseFrac = accNmseFrac + nmF;
    end

    berInt(si)    = berFloor(errInt,  totBits);
    berFrac(si)   = berFloor(errFrac, totBits);
    berPcsi3b(si) = berFloor(errPcsi, totBits);
    nmseInt(si)   = accNmseInt  / nTrials3b;
    nmseFrac(si)  = accNmseFrac / nTrials3b;

    fprintf(' %2d   %.2e    %.2e    %.2e    %+.1f      %+.1f\n', ...
        snrVec(si), berInt(si), berFrac(si), berPcsi3b(si), ...
        10*log10(nmseInt(si)+1e-20), 10*log10(nmseFrac(si)+1e-20));
end
fprintf('Part 3b done (%.0f s)\n\n', toc(t3b));

% ---- Fig.4: 双子图 ----
fig4 = figure('Position', [80 80 960 380], 'Color', 'w');

subplot(1,2,1); hold on; grid on; box on;
semilogy(snrVec, berInt, '-o', 'Color', cBlue, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cBlue, ...
    'DisplayName', 'Integer Rx (SpreadWidth=0)');
semilogy(snrVec, berFrac, '-^', 'Color', cRed, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cRed, ...
    'DisplayName', 'Fractional Rx (SpreadWidth=4)');
semilogy(snrVec, berPcsi3b, '--s', 'Color', cBlack, ...
    'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'Perfect CSI');
applyIeeeStyle(gca);
ylim([1e-4 1]);
xlabel('Data SNR (dB)'); ylabel('BER');
title('(a) BER vs SNR');
legend('Location', 'southwest', 'FontSize', 9);

subplot(1,2,2); hold on; grid on; box on;
plot(snrVec, 10*log10(nmseInt+1e-20), '-o', 'Color', cBlue, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cBlue, ...
    'DisplayName', 'Integer Rx');
plot(snrVec, 10*log10(nmseFrac+1e-20), '-^', 'Color', cRed, ...
    'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cRed, ...
    'DisplayName', 'Fractional Rx');
applyIeeeStyle(gca);
xlabel('Data SNR (dB)'); ylabel('NMSE (dB)');
title('(b) Channel Estimation NMSE');
legend('Location', 'northeast');

sgtitle(sprintf('Fig.4: Fractional Doppler Channel  N=%d P=%d k_{max}=%d', ...
    numSc, defaultPaths, maxDopplerIdx), ...
    'FontSize', 12, 'FontName', 'Times New Roman');
savefig(fig4, './Sim/Results/Fig4_FractionalDoppler.fig');

%% ################################################################
%  Part 4: Agile-c2 PAPR
%  ################################################################
fprintf('==================== Part 4: Agile-c2 PAPR ====================\n');

cfgPapr = GiFreeConfig();
cfgPapr.NumSubcarriers       = numSc;
cfgPapr.ModulationOrder      = modOrder;
cfgPapr.MaxDelaySpread       = maxDelay;
cfgPapr.MaxDopplerIndex      = maxDopplerIdx;
cfgPapr.NumPaths             = defaultPaths;
cfgPapr.DopplerGuard         = dopplerGuard3;
cfgPapr.SpreadWidth          = 0;
cfgPapr.PilotSnrDb           = pilotSnrGf;
cfgPapr.MaxSicIterations     = 10;
cfgPapr.NumPathsUpper        = 6;
cfgPapr.UseFractionalDoppler = false;
cfgPapr.AgileQ               = 0;
txPapr = GiFreeTransmitter(cfgPapr);

agileQ = 32;
candidates = AgileC2Optimizer.generateCandidates(numSc, agileQ);

nFrames4  = 2000;
snrLin4   = 10^(15/10);
paprBase  = zeros(nFrames4, 1);
paprAgile = zeros(nFrames4, 1);

t4 = tic;
idxVec = (0:numSc-1).';
chirpPhaseMat = (2*pi*idxVec.^2) * candidates.';

for fi = 1:nFrames4
    [daftFrame, ~] = txPapr.transmit(snrLin4);

    % 基线 PAPR
    c2base = cfgPapr.C2BaseValue;
    tdBase = ifft(daftFrame .* exp(1j * (2*pi*idxVec.^2) * c2base)) * sqrt(numSc);
    pwBase = abs(tdBase).^2;
    paprBase(fi) = 10*log10(max(pwBase)/mean(pwBase));

    % 敏捷 PAPR: 批量 IFFT
    preChirped = daftFrame .* exp(1j * chirpPhaseMat);
    tdAll = ifft(preChirped) * sqrt(numSc);
    pwAll = abs(tdAll).^2;
    paprAllLin = max(pwAll,[],1) ./ mean(pwAll,1);
    [~, bestIdx] = min(paprAllLin);
    paprAgile(fi) = 10*log10(paprAllLin(bestIdx));
end

fprintf('Part 4 done (%.1f s)\n', toc(t4));

p90base  = prctile(paprBase, 90);
p90agile = prctile(paprAgile, 90);
fprintf('  PAPR@90%%ile: Baseline=%.2f dB  Agile(Q=%d)=%.2f dB  Reduction=%.2f dB\n\n', ...
    p90base, agileQ, p90agile, p90base - p90agile);

% ---- Fig.5 ----
fig5 = figure('Position', [120 120 560 420], 'Color', 'w');
ccdfProb = (nFrames4:-1:1) / nFrames4;
semilogy(sort(paprBase), ccdfProb, '-', 'Color', cBlue, 'LineWidth', 2, ...
    'DisplayName', sprintf('Baseline c_2 (P90=%.1f dB)', p90base));
hold on;
semilogy(sort(paprAgile), ccdfProb, '-', 'Color', cRed, 'LineWidth', 2, ...
    'DisplayName', sprintf('Agile Q=%d (P90=%.1f dB)', agileQ, p90agile));
yline(0.1, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
hold off;
applyIeeeStyle(gca);
ylim([1e-2 1]);
xlabel('PAPR (dB)'); ylabel('CCDF = Pr(PAPR > x)');
title(sprintf('Fig.5: Agile-c_2 PAPR CCDF  N=%d  Q=%d', numSc, agileQ));
legend('Location', 'southwest');
savefig(fig5, './Sim/Results/Fig5_AgilePapr.fig');

%% ==================== 汇总 ====================
fprintf('============================================================\n');
fprintf('  All simulations complete. 5 figures generated.\n');
fprintf('  Fig.1: AFDM vs OFDM  (Fractional Doppler)\n');
fprintf('  Fig.2: EP vs GI-Free (Integer Doppler + SE)\n');
fprintf('  Fig.3: CFAR Threshold (Variable path count)\n');
fprintf('  Fig.4: Fractional Doppler (Int Rx vs Frac Rx)\n');
fprintf('  Fig.5: Agile-c2 PAPR CCDF\n');
fprintf('============================================================\n');

%% ==================== 辅助函数 ====================

function val = berFloor(numErrors, totalBits)
% berFloor  BER 计算 (含零错误时的 Clopper-Pearson 95% 置信上限)
%
%   Rule of 3: P(X=0|p) = (1-p)^N <= 0.05  =>  p <= 3/N
%   含义: 95% 置信度下真实 BER 不超过此值, 防止 semilogy 断裂
    if numErrors > 0
        val = numErrors / totalBits;
    else
        val = 3 / totalBits;
    end
end

function [delays, dopplers, gains] = generateRandomChannel(numPaths, maxDelay, maxDoppler, isFractional)
% generateRandomChannel  随机信道参数生成 (Jakes 多普勒频谱)
    allDelays = (0:maxDelay).';
    delays    = sort(allDelays(randperm(numel(allDelays), numPaths)));
    thetas    = -pi + 2*pi*rand(numPaths, 1);
    dopplers  = maxDoppler * cos(thetas);

    if ~isFractional
        dopplers = round(dopplers);
        for k = 2:numPaths
            attempts = 0;
            while any(delays(1:k-1) == delays(k) & ...
                      dopplers(1:k-1) == dopplers(k)) && attempts < 100
                thetas(k)   = -pi + 2*pi*rand();
                dopplers(k) = round(maxDoppler * cos(thetas(k)));
                attempts    = attempts + 1;
            end
        end
    end

    gains = sqrt(1/(2*numPaths)) * (randn(numPaths,1) + 1j*randn(numPaths,1));
    delays   = delays(:);
    dopplers = dopplers(:);
    gains    = gains(:);
end

function detected = paperDetect(rxSig, numSc, c1, c2, amp, locStep, maxDelay, maxDop, thr)
% paperDetect  论文固定门限路径检测 [Zhou et al. 2024, eq.(17)-(18)]
    detected = zeros(0, 3);
    for l = 0:maxDelay
        for k = -maxDop:maxDop
            pos = mod(-(k + locStep*l), numSc);
            if abs(rxSig(pos+1)) >= thr
                phase = exp(1j*2*pi*(c1*l^2 - c2*pos^2));
                gain  = rxSig(pos+1) / (phase * amp);
                detected(end+1, :) = [l, k, gain]; %#ok<AGROW>
            end
        end
    end
end

function hEff = buildHeff(builder, pathParams)
% buildHeff  从检测路径参数构造等效信道矩阵
    numSc = builder.Config.NumSubcarriers;
    if isempty(pathParams)
        hEff = sparse(numSc, numSc);
    else
        hEff = builder.buildEffectiveChannel(pathParams);
    end
end

function [nCorrect, nFalse] = matchPaths(trueP, estP)
% matchPaths  路径匹配 (按 delay + round(doppler))
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
% quickBer  LMMSE 检测 BER (隔离门限贡献, 使用 biterr)
    cleanSig = rxSig - hEst * pilotFrame;
    estSig   = GiFreeReceiver.lmmseDetect(cleanSig, hEst, regParam, numSc, dataPos);
    rxIdx    = qamdemod(estSig(dataPos)/sqrt(snrLin), modOrd, 'UnitAveragePower', true);
    [nErr, ~] = biterr(txIdx, rxIdx, log2(modOrd));
    val = nErr / (numel(txIdx) * log2(modOrd));
end

function applyIeeeStyle(ax, fontSize)
% applyIeeeStyle  IEEE 论文风格: Times New Roman, 合适字号, 紧凑布局
%   参考 GI-Free 论文 (2404.01088v1) Fig.5 配色和排版
    if nargin < 2, fontSize = 11; end
    set(ax, ...
        'FontName',   'Times New Roman', ...
        'FontSize',   fontSize, ...
        'LineWidth',  0.8, ...
        'TickDir',    'in', ...
        'TickLength', [0.015 0.015], ...
        'Box',        'on');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
end