function MainSimulation()
% MAINSIMULATION 运行第四章主干仿真实验并生成论文图。
%
% 描述:
%   主入口只负责路径、日志、公共参数和 Part 1 至 Part 4 实验分派。
%
% 输出:
%   图文件: Sim/Results/Figures/*.fig
%   日志:   Sim/Results/Logs/MainSimulation_*.txt

    close all; clc;
    simDir = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(simDir);
    addpath(genpath(repoRoot));

    [figuresDir, logsDir] = initializeSimulationOutputs(simDir);
    logPath = fullfile(logsDir, sprintf('MainSimulation_%s.txt', ...
        char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
    diary(logPath);
    diaryCleanup = onCleanup(@() diary('off'));

    fprintf('============================================================\n');
    fprintf('     AFDM Chapter-4 Simulation Pipeline\n');
    fprintf('============================================================\n\n');

    params = buildMainSimulationParams();
    runPart1AfdmVsOfdm(params, figuresDir);
    part2Result = runPart2DetectorAttribution(params, figuresDir);
    printPart2Tables(part2Result);
    runPart3DynamicPilot(params, figuresDir);
    runPart4FractionalDoppler(params, figuresDir);

    fprintf('============================================================\n');
    fprintf('  Chapter-4 simulations complete.\n');
    fprintf('  Fig.1: AFDM-EP vs CP-OFDM BER baseline\n');
    fprintf('  Fig.2: Phase-1 detector attribution under fixed pilot\n');
    fprintf('  Fig.3: Dynamic pilot gain for final CLIP receiver\n');
    fprintf('  Fig.4: EP vs GI-Free fractional-Doppler comparison\n');
    fprintf('  Figures -> %s\n', figuresDir);
    fprintf('  Log     -> %s\n', logPath);
    fprintf('============================================================\n');
    clear diaryCleanup;
end

function [figuresDir, logsDir] = initializeSimulationOutputs(simDir)
% INITIALIZESIMULATIONOUTPUTS 创建图像与日志输出目录。

    figuresDir = fullfile(simDir, 'Results', 'Figures');
    logsDir = fullfile(simDir, 'Results', 'Logs');
    if ~isfolder(figuresDir)
        mkdir(figuresDir);
    end
    if ~isfolder(logsDir)
        mkdir(logsDir);
    end
end

function params = buildMainSimulationParams()
% BUILDMAINSIMULATIONPARAMS 返回第四章主干实验公共参数。

    params.NumSubcarriers = 512;
    params.ModulationOrder = 4;
    params.BitsPerSymbol = log2(params.ModulationOrder);
    params.MaxDelaySamples = 4;
    params.MaxDopplerIdx = 2;
    params.DopplerGuard = 4;
    params.DefaultPaths = 3;
    params.Exp3PathCount = 5;
    params.PilotSnrEpDb = 35;
    params.PilotSnrGiFreeDb = 35;
    params.SnrVecDb = (0:5:20).';
    params.Exp1TrialCount = 1000;
end

function runPart1AfdmVsOfdm(params, figuresDir)
% RUNPART1AFDMVSOFDM 执行 AFDM-EP 与 CP-OFDM 基线对比实验。

    fprintf('==================== Part 1: AFDM vs OFDM ====================\n');
    partTimer = tic;
    Exp1_AfdmVsOfdm( ...
        params.NumSubcarriers, params.ModulationOrder, params.BitsPerSymbol, ...
        params.MaxDelaySamples, params.MaxDopplerIdx, params.DopplerGuard, ...
        params.DefaultPaths, params.PilotSnrEpDb, params.SnrVecDb, figuresDir, ...
        params.Exp1TrialCount);
    fprintf('Part 1 done (%.0f s)\n\n', toc(partTimer));
end

function resultStruct = runPart2DetectorAttribution(params, figuresDir)
% RUNPART2DETECTORATTRIBUTION 执行固定导频阶段一检测器归因实验。

    fprintf('==================== Part 2: Phase-1 Detector Attribution ====================\n');
    configObj = createExp2GiFreeConfig( ...
        params.NumSubcarriers, params.ModulationOrder, params.MaxDelaySamples, ...
        params.MaxDopplerIdx, 5, params.PilotSnrGiFreeDb);
    systemObj = GiFreeSystem(configObj);

    partTimer = tic;
    resultStruct = Exp2_CfarThreshold( ...
        systemObj, params.ModulationOrder, params.BitsPerSymbol, figuresDir);
    fprintf('Part 2 done (%.0f s)\n\n', toc(partTimer));
end

function runPart3DynamicPilot(params, figuresDir)
% RUNPART3DYNAMICPILOT 执行动、静态导频对比实验。

    fprintf('==================== Part 3: Dynamic Pilot Gain ====================\n');
    partTimer = tic;
    Exp3_DynamicPilotGain( ...
        params.NumSubcarriers, params.ModulationOrder, params.BitsPerSymbol, ...
        params.MaxDelaySamples, params.MaxDopplerIdx, params.DopplerGuard, ...
        params.Exp3PathCount, params.PilotSnrGiFreeDb, params.SnrVecDb, figuresDir);
    fprintf('Part 3 done (%.0f s)\n\n', toc(partTimer));
end

function runPart4FractionalDoppler(params, figuresDir)
% RUNPART4FRACTIONALDOPPLER 执行 EP 与 GI-Free 的分数 Doppler 对比实验。

    fprintf('==================== Part 4: EP vs GI-Free in Fractional Doppler ====================\n');
    partTimer = tic;
    Exp4_FractionalDoppler( ...
        params.NumSubcarriers, params.ModulationOrder, params.BitsPerSymbol, ...
        params.MaxDelaySamples, params.MaxDopplerIdx, params.DopplerGuard, ...
        params.DefaultPaths, params.PilotSnrEpDb, params.PilotSnrGiFreeDb, ...
        params.SnrVecDb, figuresDir);
    fprintf('Part 4 done (%.0f s)\n\n', toc(partTimer));
end

function printPart2Tables(resultStruct)
% PRINTPART2TABLES 打印 Part 2 的 P=3 与 P=5 消融表格。

    snrPrintVec = [10, 15, 20];
    fprintf('\n=== Table 4.X: P=3 Ablation ===\n');
    printAblationTable(resultStruct, 3, snrPrintVec);
    fprintf('\n=== Table 4.Y: P=5 Ablation ===\n');
    printAblationTable(resultStruct, 5, snrPrintVec);
end

function printAblationTable(resultStruct, pathCount, snrPrintVec)
% PRINTABLATIONTABLE 打印指定路径数下的 BER、Pd 与 NMSE 摘要。

    pathIdx = find(resultStruct.pathCountVec == pathCount, 1);
    if isempty(pathIdx)
        fprintf('  No P=%d data.\n', pathCount);
        return;
    end

    fprintf('Method                      | Integer@10/15/20 dB | Fractional@10/15/20 dB\n');
    fprintf('                            | BER  Pd   NMSE      | BER  Pd   NMSE\n');
    for methodIdx = 1:numel(resultStruct.methodLabels)
        fprintf('%-27s | ', resultStruct.methodLabels{methodIdx});
        for scenarioIdx = 1:numel(resultStruct.scenarioLabels)
            printScenarioCells(resultStruct, pathIdx, scenarioIdx, methodIdx, snrPrintVec);
            if scenarioIdx == 1
                fprintf('| ');
            end
        end
        fprintf('\n');
    end
end

function printScenarioCells(resultStruct, pathIdx, scenarioIdx, methodIdx, snrPrintVec)
% PRINTSCENARIOCELLS 打印单个场景的三个 SNR 摘要单元。

    for snrDb = snrPrintVec
        snrIdx = find(resultStruct.snrVecDb == snrDb, 1);
        if isempty(snrIdx)
            fprintf('N/A  N/A  N/A       ');
            continue;
        end
        berValue = resultStruct.avgBer(snrIdx, pathIdx, scenarioIdx, methodIdx);
        pdValue = resultStruct.avgPd(snrIdx, pathIdx, scenarioIdx, methodIdx);
        nmseValue = resultStruct.avgNmse(snrIdx, pathIdx, scenarioIdx, methodIdx);
        fprintf('%.2e %.2f %.2e ', berValue, pdValue, nmseValue);
    end
end
