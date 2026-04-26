function tests = testSimulationRefactor
% TESTSIMULATIONREFACTOR - 第四章仿真重构的关键行为测试
%
%   描述:
%   覆盖 GI-Free 配置深拷贝、置信区间、Exp1 主入口、Exp2/Exp3/Exp4 helper、CFAR 门限、
%   固定门限 OMP 以及 Pfa 校准实验的轻量级冒烟验证。
%
%   语法:
%   results = runtests("Tests/testSimulationRefactor.m");
%
%   输出:
%   tests - (matlab.unittest.Test) 函数式测试集合。
%
%   版本历史:
%   2026-04-19 - Aiden - 新增第四章仿真重构测试。

    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% SETUPONCE 为测试加入仓库路径。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(repoRoot));
    testCase.TestData.RepoRoot = repoRoot;
end

function testComputeConfidenceIntervalZeroErrors(testCase)
% TESTCOMPUTECONFIDENCEINTERVALZEROERRORS 验证零错误事件使用 Rule-of-3。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    [ciLower, ciUpper] = computeConfidenceInterval(0, 5000);

    verifyEqual(testCase, ciLower, 0);
    verifyEqual(testCase, ciUpper, 3 / 5000, 'AbsTol', 1e-15);
end

function testGiFreeConfigCopyDoesNotShareHandle(testCase)
% TESTGIFREECONFIGCOPYDOESNOTSHAREHANDLE 验证配置复制不会共享 handle 状态。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    sourceConfig = createMinimalConfig();
    sourceConfig.EnableProgressiveCfar = true;

    copiedConfig = sourceConfig.copy();
    copiedConfig.EnableProgressiveCfar = false;
    copiedConfig.NumSubcarriers = 128;

    verifyTrue(testCase, sourceConfig.EnableProgressiveCfar);
    verifyFalse(testCase, copiedConfig.EnableProgressiveCfar);
    verifyEqual(testCase, sourceConfig.NumSubcarriers, 64);
    verifyEqual(testCase, copiedConfig.NumSubcarriers, 128);
end

function testCreateExp2GiFreeConfigUsesFixedPilot(testCase)
% TESTCREATEEXP2GIFREECONFIGUSESFIXEDPILOT 验证 Exp2 配置使用固定导频。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    pilotSnrDb = 35;
    configObj = createExp2GiFreeConfig(512, 4, 4, 2, 5, pilotSnrDb);
    configObj.CurrentDataSnrLin = 100;

    expectedPilotPower = 10 ^ (pilotSnrDb / 10);
    verifyFalse(testCase, configObj.UseDynamicPilot);
    verifyEqual(testCase, configObj.PilotAmplitude ^ 2, expectedPilotPower, 'RelTol', 1e-12);
end

function testCreateExp2GiFreeConfigLoosensNumPathsUpperForP5(testCase)
% TESTCREATEEXP2GIFREECONFIGLOOSENSNUMPATHSUPPERFORP5 验证 Exp2 为高路径数放宽路径上限。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    configObj = createExp2GiFreeConfig(512, 4, 4, 2, 5, 35);

    verifyEqual(testCase, configObj.NumPathsUpper, 10);
end

function testBuildExp2TrialCountVecRaisesHighSnrSamples(testCase)
% TESTBUILDEXP2TRIALCOUNTVECRAISESHIGHSNRSAMPLES 验证 Exp2 高 SNR 使用更高 trial 数。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    snrVecDb = (0:5:20).';
    trialCountVec = buildExp2TrialCountVec(snrVecDb);

    verifyEqual(testCase, trialCountVec, [5000; 5000; 2000; 4000; 4000]);
end

function testResolveExp2PathCountVecSupportsP5Only(testCase)
% TESTRESOLVEEXP2PATHCOUNTVECSUPPORTSP5ONLY 验证 Exp2 可只运行 P=5。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    verifyEqual(testCase, resolveExp2PathCountVec([]), [3, 5]);
    verifyEqual(testCase, resolveExp2PathCountVec(5), 5);
end

function testBuildGiFreePilotFrameForSnrRefreshesDynamicPilot(testCase)
% TESTBUILDGIFREEPILOTFRAMEFORSNRREFRESHESDYNAMICPILOT 验证导频帧按当前 SNR 构造。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    configObj = createMinimalConfig();
    configObj.UseDynamicPilot = true;
    configObj.DynamicPilotBaseDb = 20;
    configObj.CurrentDataSnrLin = 1;
    dataSnrLin = 100;

    [pilotFrame, pilotAmplitude] = buildGiFreePilotFrameForSnr(configObj, dataSnrLin);

    expectedPilotPower = dataSnrLin * 10 ^ (configObj.DynamicPilotBaseDb / 10);
    verifyEqual(testCase, configObj.CurrentDataSnrLin, dataSnrLin);
    verifyEqual(testCase, pilotAmplitude ^ 2, expectedPilotPower, 'RelTol', 1e-12);
    verifyEqual(testCase, abs(pilotFrame(configObj.PilotPos1)) ^ 2, expectedPilotPower, 'RelTol', 1e-12);
end

function testCreateClipPilotConfigPairOnlyDiffersInPilotPolicy(testCase)
% TESTCREATECLIPPILOTCONFIGPAIRONLYDIFFERSINPILOTPOLICY 验证 Exp3 配置对仅区分导频策略。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    [fixedConfig, dynamicConfig] = createClipPilotConfigPair( ...
        128, 4, 3, 2, 4, 3, 35);

    verifyEqual(testCase, fixedConfig.NumSubcarriers, dynamicConfig.NumSubcarriers);
    verifyEqual(testCase, fixedConfig.MaxDelaySamples, dynamicConfig.MaxDelaySamples);
    verifyEqual(testCase, fixedConfig.MaxDopplerIdx, dynamicConfig.MaxDopplerIdx);
    verifyEqual(testCase, fixedConfig.NumPaths, dynamicConfig.NumPaths);
    verifyEqual(testCase, fixedConfig.DirichletRadius, dynamicConfig.DirichletRadius);
    verifyTrue(testCase, fixedConfig.UseFractionalDoppler);
    verifyTrue(testCase, dynamicConfig.UseFractionalDoppler);
    verifyTrue(testCase, fixedConfig.UseMatchedPilotMetric);
    verifyTrue(testCase, dynamicConfig.UseMatchedPilotMetric);
    verifyFalse(testCase, fixedConfig.UseDynamicPilot);
    verifyTrue(testCase, dynamicConfig.UseDynamicPilot);

    fixedConfig.CurrentDataSnrLin = 1;
    dynamicConfig.CurrentDataSnrLin = 100;
    verifyEqual(testCase, fixedConfig.PilotAmplitude ^ 2, 10 ^ (35 / 10), 'RelTol', 1e-12);
    verifyEqual(testCase, dynamicConfig.PilotAmplitude ^ 2, 100 * 10 ^ (35 / 10), 'RelTol', 1e-12);
end

function testEstimateRequiredSnrAtTargetBerInterpolatesLogScale(testCase)
% TESTESTIMATEREQUIREDSNRATTARGETBERINTERPOLATESLOGSCALE 验证目标 BER 工作点插值。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    snrVecDb = [0; 5; 10];
    berVec = [1e-1; 1e-2; 1e-3];

    [requiredSnrDb, hasCrossing] = estimateRequiredSnrAtTargetBer( ...
        snrVecDb, berVec, 10 ^ (-1.5));

    verifyTrue(testCase, hasCrossing);
    verifyEqual(testCase, requiredSnrDb, 2.5, 'AbsTol', 1e-12);
end

function testSummarizeExp4SnrGapReturnsExpectedGapTable(testCase)
% TESTSUMMARIZEEXP4SNRGAPRETURNSEXPECTEDGAPTABLE 验证 Exp4 工作点差值摘要。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    resultStruct = struct();
    resultStruct.snrVecDb = [0; 5; 10; 15];
    resultStruct.berEp = [ ...
        2e-1, 2e-1, 2e-1; ...
        6e-2, 5e-2, 4e-2; ...
        8e-3, 6e-3, 4e-3; ...
        8e-4, 6e-4, 4e-4];
    resultStruct.berAblation = zeros(4, 4);
    resultStruct.berAblation(:, 4) = [2e-1; 7e-2; 1.2e-2; 1e-3];

    gapSummary = summarizeExp4SnrGap(resultStruct, [1e-2; 1e-3]);

    verifySize(testCase, gapSummary, [2 9]);
    verifyTrue(testCase, gapSummary.ProGuardHasCrossing(1));
    verifyGreaterThan(testCase, gapSummary.GapVsEp1Db(1), 0);
    verifyGreaterThan(testCase, gapSummary.GapVsEp4Db(1), 0);
    verifyFalse(testCase, isnan(gapSummary.GapVsEp4Db(2)));
end

function testMainSimulationPassesSnrVecToExp4FractionalDoppler(testCase)
% TESTMAINSIMULATIONPASSESSNRVECTOEXP4FRACTIONALDOPPLER 验证 Part 4 调用传入 SNR 向量。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    mainSimulationPath = fullfile(testCase.TestData.RepoRoot, 'Sim', 'MainSimulation.m');
    mainSimulationText = fileread(mainSimulationPath);

    expectedCallPattern = ['params\.DefaultPaths,\s*params\.PilotSnrEpDb,\s*', ...
        'params\.PilotSnrGiFreeDb,\s*params\.SnrVecDb,\s*figuresDir'];
    verifyNotEmpty(testCase, regexp(mainSimulationText, expectedCallPattern, 'once'));
end

function testMainSimulationRunsExp1AfdmVsOfdm(testCase)
% TESTMAINSIMULATIONRUNSEXP1AFDMVSOFDM 验证主入口保留 AFDM vs OFDM 实验。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    mainSimulationPath = fullfile(testCase.TestData.RepoRoot, 'Sim', 'MainSimulation.m');
    exp1Path = fullfile(testCase.TestData.RepoRoot, 'Sim', 'experiments', 'Exp1_AfdmVsOfdm.m');
    mainSimulationText = fileread(mainSimulationPath);

    verifyTrue(testCase, isfile(exp1Path));
    verifyNotEmpty(testCase, regexp(mainSimulationText, 'runPart1AfdmVsOfdm', 'once'));
    verifyNotEmpty(testCase, regexp(mainSimulationText, 'Exp1_AfdmVsOfdm', 'once'));
end

function testPublicCfarThresholdMatchesFormula(testCase)
% TESTPUBLICCFARTHRESHOLDMATCHESFORMULA 验证公开 CFAR 门限接口的公式一致性。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    configObj = createMinimalConfig();
    systemObj = GiFreeSystem(configObj);
    weightedResidual = ones(configObj.NumSubcarriers, 1);
    thresholdScale = 1.2;

    thresholdValue = systemObj.Estimator.computeCfarThreshold( ...
        weightedResidual, configObj.NumSubcarriers, thresholdScale);

    numCandidates = (configObj.MaxDelaySamples + 1) * ...
        (2 * configObj.MaxDopplerIdx + 1);
    residualPower = real(weightedResidual' * weightedResidual) / configObj.NumSubcarriers;
    expectedThreshold = residualPower * configObj.PilotAmplitude^2 * ...
        log(numCandidates / systemObj.Estimator.cfarFalseAlarmProb) * thresholdScale;

    verifyEqual(testCase, thresholdValue, expectedThreshold, 'RelTol', 1e-12);
end

function testConfigCfarFalseAlarmOverridesEstimatorDefault(testCase)
% TESTCONFIGCFARFALSEALARMOVERRIDESESTIMATORDEFAULT 验证配置级 Pfa 可覆盖估计器默认值。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    configObj = createMinimalConfig();
    configObj.CfarFalseAlarmProb = 1e-2;
    systemObj = GiFreeSystem(configObj);
    weightedResidual = ones(configObj.NumSubcarriers, 1);

    thresholdValue = systemObj.Estimator.computeCfarThreshold( ...
        weightedResidual, configObj.NumSubcarriers, 1.0);

    numCandidates = (configObj.MaxDelaySamples + 1) * ...
        (2 * configObj.MaxDopplerIdx + 1);
    expectedThreshold = configObj.PilotAmplitude^2 * ...
        log(numCandidates / configObj.CfarFalseAlarmProb);

    verifyEqual(testCase, thresholdValue, expectedThreshold, 'RelTol', 1e-12);
end

function testMatchedPilotMetricCfarThresholdOmitsPilotPower(testCase)
% TESTMATCHEDPILOTMETRICCFARTHRESHOLDOMITSPILOTPOWER 验证匹配统计量门限量纲。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    configObj = createMinimalConfig();
    configObj.UseMatchedPilotMetric = true;
    systemObj = GiFreeSystem(configObj);
    weightedResidual = ones(configObj.NumSubcarriers, 1);

    thresholdValue = systemObj.Estimator.computeCfarThreshold( ...
        weightedResidual, configObj.NumSubcarriers, 1.0);

    numCandidates = (configObj.MaxDelaySamples + 1) * ...
        (2 * configObj.MaxDopplerIdx + 1);
    expectedThreshold = log(numCandidates / configObj.CfarFalseAlarmProb);

    verifyEqual(testCase, thresholdValue, expectedThreshold, 'RelTol', 1e-12);
end

function testEstimateByOmpWithFixedThresholdSuppressesPeak(testCase)
% TESTESTIMATEBYOMPWITHFIXEDTHRESHOLDSUPPRESSESPEAK 验证固定高门限抑制 OMP 检测。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    configObj = createMinimalConfig();
    systemObj = GiFreeSystem(configObj);
    rxSignal = 2 * systemObj.ChannelBuilder.buildCompositePilotResponse(0, 0);

    [estimatedPaths, effectiveChannel, ompDiag] = ...
        systemObj.Estimator.estimateByOmpWithFixedThreshold( ...
            rxSignal, 0, ones(configObj.NumSubcarriers, 1), realmax);

    verifyEmpty(testCase, estimatedPaths);
    verifySize(testCase, effectiveChannel, [configObj.NumSubcarriers, configObj.NumSubcarriers]);
    verifyEqual(testCase, ompDiag.modeUsed, 'fixed_threshold');
end

function testPaperCfarThresholdMatchesMetricScale(testCase)
% TESTPAPERCFARTHRESHOLDMATCHESMETRICSCALE 验证相关搜索 CFAR 门限的量纲。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    rxSignal = ones(64, 1);
    pilotAmplitude = 10;
    targetFalseAlarmProb = 0.03;

    thresholdValue = paperCfarThreshold(rxSignal, 64, pilotAmplitude, 2, 1, targetFalseAlarmProb);

    expectedThreshold = pilotAmplitude^2 * log(9 / targetFalseAlarmProb);
    verifyEqual(testCase, thresholdValue, expectedThreshold, 'RelTol', 1e-12);
end

function testExp4PfaCalibrationSmoke(testCase)
% TESTEXP4PFACALIBRATIONSMMOKE 验证 Pfa 校准实验可返回有限概率向量。
%
% 输入:
%   testCase - (matlab.unittest.FunctionTestCase) 测试上下文。

    rng(42);
    configObj = createMinimalConfig();
    [empiricalPfa, ciLower, ciUpper] = Exp4_PfaCalibration(configObj, 20);

    verifySize(testCase, empiricalPfa, [6 1]);
    verifySize(testCase, ciLower, [6 1]);
    verifySize(testCase, ciUpper, [6 1]);
    verifyGreaterThanOrEqual(testCase, empiricalPfa, zeros(6, 1));
    verifyLessThanOrEqual(testCase, empiricalPfa, ones(6, 1));
end

function configObj = createMinimalConfig()
% CREATEMINIMALCONFIG 构造满足正交约束的轻量 GI-Free 配置。
%
% 输出:
%   configObj - (GiFreeConfig) 用于单元测试的最小配置。

    configObj = GiFreeConfig();
    configObj.NumSubcarriers = 64;
    configObj.ModulationOrder = 4;
    configObj.MaxDelaySamples = 2;
    configObj.MaxDopplerIdx = 1;
    configObj.NumPaths = 1;
    configObj.DopplerGuard = 2;
    configObj.DirichletRadius = 0;
    configObj.PilotSnrDb = 20;
    configObj.MaxSicIterations = 2;
    configObj.NumPathsUpper = 2;
    configObj.UseFractionalDoppler = false;
    configObj.EnableProgressiveCfar = true;
    configObj.EnablePathStabilityGate = true;
end
