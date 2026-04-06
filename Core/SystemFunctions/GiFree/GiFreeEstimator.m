classdef GiFreeEstimator < handle
% GIFREEESTIMATOR - GI-Free 稀疏路径估计与参数精炼模块
%
%   描述:
%   提供 OMP+CFAR 初始化、数据导向 Doppler/增益精炼、残差路径注入候选生成、
%   以及噪声功率截尾估计等能力。输出路径矩阵统一为 [delay, doppler, gain]。
%
%   语法:
%   estimator = GiFreeEstimator(cfg, channelBuilder);
%   [estPaths, hEff] = estimator.estimateByOmp(rxSignal, ompRegParam, obsWeights);
%   [estPaths, hEff, ompDiag] = estimator.estimateByOmp(rxSignal, ompRegParam, ...
%       obsWeights, searchModeOverride);
%
%   输出:
%   estPaths - (Px3) [delay, doppler, complexGain]。
%   hEff     - (NxN sparse/full) 估计有效信道矩阵。
%   ompDiag  - (struct) OMP 诊断信息，包含搜索模式、候选统计、残差与最终支撑签名。
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。

    properties (SetAccess = private)
        Config
        ChannelBuilder
    end

    properties (Constant)
        % 虚警率3%，在高压力场景下降低漏检真实路径的风险
        cfarFalseAlarmProb = 0.03
    end

    methods

        % GIFREEESTIMATOR 构造函数，绑定配置与信道构造算子。
        function obj = GiFreeEstimator(cfg, chBuilder)
            arguments
                cfg       (1,1) GiFreeConfig
                chBuilder (1,1) IChannelOperator
            end
            obj.Config = cfg;
            obj.ChannelBuilder = chBuilder;
        end

        % ESTIMATEBYOMP 通过 OMP + CFAR 执行稀疏路径初始化估计。
        function [estPaths, hEffective, ompDiag] = estimateByOmp( ...
                obj, rxSignal, ompRegParam, observationWeights, searchModeOverride)
            if nargin < 3
                ompRegParam = 0;
            end
            if nargin < 4 || isempty(observationWeights)
                observationWeights = ones(size(rxSignal));
            end
            if nargin < 5
                searchModeOverride = '';
            end

            searchMode = obj.resolveSearchMode(searchModeOverride);
            if strcmpi(searchMode, 'beam')
                [estPaths, hEffective, ompDiag] = obj.estimateByBeamOmp( ...
                    rxSignal, ompRegParam, observationWeights);
                return;
            end

            [estPaths, hEffective, ompDiag] = obj.estimateByGreedyOmp( ...
                rxSignal, ompRegParam, observationWeights, zeros(0, 3), 'greedy');
        end

        % ESTIMATEBYGREEDYOMP 执行标准 OMP，并可从已有支撑集继续扩展。
        function [estPaths, hEffective, ompDiag] = estimateByGreedyOmp( ...
                obj, rxSignal, ompRegParam, observationWeights, seedPaths, modeLabel)
            if nargin < 5
                seedPaths = zeros(0, 3);
            end
            if nargin < 6 || isempty(modeLabel)
                modeLabel = 'greedy';
            end

            numSubcarriers = obj.Config.NumSubcarriers;
            numPathsUpper  = obj.Config.NumPathsUpper;
            numCandidates = (obj.Config.MaxDelaySamples + 1) * ...
                            (2 * obj.Config.MaxDopplerIdx + 1);
            cfarCoeff = log(numCandidates / obj.cfarFalseAlarmProb);
            pilotAmpSq = obj.Config.PilotAmplitude^2;
            fracRefineCallCount = 0;
            ompDiag = obj.createEmptyOmpDiag(modeLabel);

            estPaths = seedPaths;
            numDetected = size(seedPaths, 1);
            pilotDictionary = zeros(numSubcarriers, numPathsUpper);
            if numDetected > 0
                for pathIdx = 1:numDetected
                    pilotDictionary(:, pathIdx) = obj.ChannelBuilder.buildCompositePilotResponse( ...
                        seedPaths(pathIdx, 1), seedPaths(pathIdx, 2));
                end
                activeDict = pilotDictionary(:, 1:numDetected);
                gramMatrix = activeDict' * activeDict;
                if ompRegParam > 0
                    jointGains = (gramMatrix + ompRegParam * eye(numDetected)) \ (activeDict' * rxSignal);
                else
                    jointGains = activeDict \ rxSignal;
                end
                estPaths(:, 3) = jointGains;
                residual = rxSignal - activeDict * jointGains;
            else
                residual = rxSignal;
            end
            ompDiag.weightedResidualNormPerDepth = obj.appendResidualNorm( ...
                ompDiag.weightedResidualNormPerDepth, residual, observationWeights, numSubcarriers);

            for pathIdx = (numDetected + 1):numPathsUpper
                [peakHit, peakMetric] = obj.findStrongestPeakExcluding(residual, estPaths, observationWeights);

                weightedResidual = sqrt(observationWeights) .* residual;
                residualPower = real(weightedResidual' * weightedResidual) / numSubcarriers;
                threshold = residualPower * pilotAmpSq * cfarCoeff;
                if isempty(peakHit) || peakMetric < threshold
                    break;
                end

                delayVal = peakHit(1);
                dopplerIdx = peakHit(2);
                if obj.Config.SpreadWidth == 0
                    fracDoppler = dopplerIdx;
                else
                    fracDoppler = obj.refineFractionalDoppler(residual, delayVal, dopplerIdx);
                    fracRefineCallCount = fracRefineCallCount + 1;
                end

                compositeResp = obj.ChannelBuilder.buildCompositePilotResponse(delayVal, fracDoppler);
                pilotDictionary(:, pathIdx) = compositeResp;
                numDetected = pathIdx;

                activeDict = pilotDictionary(:, 1:pathIdx);
                gramMatrix = activeDict' * activeDict;
                if ompRegParam > 0
                    jointGains = (gramMatrix + ompRegParam * eye(pathIdx)) \ (activeDict' * rxSignal);
                else
                    jointGains = activeDict \ rxSignal;
                end
                residual = rxSignal - activeDict * jointGains;
                ompDiag.weightedResidualNormPerDepth = obj.appendResidualNorm( ...
                    ompDiag.weightedResidualNormPerDepth, residual, observationWeights, numSubcarriers);

                estPaths(pathIdx, :) = [delayVal, fracDoppler, jointGains(pathIdx)];
            end

            estPaths = estPaths(1:numDetected, :);
            if numDetected > 0
                estPaths(:, 3) = jointGains;
                hEffective = obj.ChannelBuilder.buildEffectiveChannel(estPaths);
            else
                hEffective = sparse(numSubcarriers, numSubcarriers);
            end
            ompDiag.fracRefineCallCount = fracRefineCallCount;
            ompDiag.finalSupportSignature = obj.buildSupportSignature(estPaths);
        end

        % ESTIMATEBYBEAMOMP 在前几层使用 beam search，之后回退到 greedy OMP。
        function [estPaths, hEffective, ompDiag] = estimateByBeamOmp(obj, rxSignal, ompRegParam, observationWeights)
            beamStrategyMode = obj.resolveBeamStrategyMode();
            [beamWidth, beamDepth] = obj.resolveBeamBaseDims(beamStrategyMode);
            probeBeamWidth = obj.resolveBeamProbeWidth(beamWidth, beamStrategyMode);
            expandFactor = max(round(obj.Config.BeamExpandFactor), 1);
            probeExpandWidth = probeBeamWidth * expandFactor;
            numSubcarriers = obj.Config.NumSubcarriers;
            numCandidates = (obj.Config.MaxDelaySamples + 1) * ...
                            (2 * obj.Config.MaxDopplerIdx + 1);
            cfarCoeff = log(numCandidates / obj.cfarFalseAlarmProb);
            pilotAmpSq = obj.Config.PilotAmplitude^2;
            minImproveRatio = max(obj.Config.BeamMinImproveRatio, 0);
            adaptiveImproveRatio = max(obj.Config.BeamAdaptiveImproveRatio, 0);
            adaptiveSpreadThreshold = max(round(obj.Config.BeamAdaptiveSpreadThreshold), 0);
            uncertainMetricRatio = max(obj.Config.BeamUncertainMetricRatio, 1);
            enableSignaturePrune = logical(obj.Config.EnableBeamSignaturePrune);
            fracRefineCallCount = 0;

            beamStates = struct('supportPaths', zeros(0,3), ...
                                'dictionaryMatrix', zeros(numSubcarriers,0), ...
                                'weightedResidualNorm', real((sqrt(observationWeights).*rxSignal)' * (sqrt(observationWeights).*rxSignal)) / numSubcarriers, ...
                                'residual', rxSignal);
            ompDiag = obj.createEmptyOmpDiag('beam');
            bestResidualPrev = beamStates.weightedResidualNorm;
            prevDepthImproveRatio = inf;

            for depthIdx = 1:beamDepth
                expandedStates = repmat(beamStates(1), 0, 1);
                coarseProposalTable = zeros(0, 5);
                uncertainStateIdx = zeros(0, 1);
                for stateIdx = 1:numel(beamStates)
                    state = beamStates(stateIdx);
                    weightedResidual = sqrt(observationWeights) .* state.residual;
                    residualPower = real(weightedResidual' * weightedResidual) / numSubcarriers;
                    threshold = residualPower * pilotAmpSq * cfarCoeff;

                    [candidateList, ~] = obj.selectTopCandidates( ...
                        state.residual, state.supportPaths, observationWeights, probeExpandWidth, threshold);
                    if isempty(candidateList)
                        continue;
                    end

                    if obj.isStateBeamUncertain(candidateList, uncertainMetricRatio)
                        uncertainStateIdx(end+1, 1) = stateIdx; %#ok<AGROW>
                    end
                    localRank = (1:size(candidateList, 1)).';
                    stateCol = stateIdx * ones(size(candidateList, 1), 1);
                    coarseProposalTable = [coarseProposalTable; ...
                        [stateCol, candidateList(:, 1), candidateList(:, 2), candidateList(:, 3), localRank]]; %#ok<AGROW>
                end
                currentBeamWidth = obj.resolveBeamWidthForDepth( ...
                    beamWidth, probeBeamWidth, beamStrategyMode, prevDepthImproveRatio, ...
                    ~isempty(uncertainStateIdx), adaptiveImproveRatio, adaptiveSpreadThreshold);
                currentExpandWidth = currentBeamWidth * expandFactor;
                selectedProposalTable = obj.selectBeamExpansionProposals( ...
                    coarseProposalTable, currentBeamWidth, currentExpandWidth, uncertainStateIdx);
                ompDiag.candidateCountPerDepth(end+1, 1) = size(selectedProposalTable, 1);

                for proposalIdx = 1:size(selectedProposalTable, 1)
                    state = beamStates(selectedProposalTable(proposalIdx, 1));
                    delayVal = selectedProposalTable(proposalIdx, 2);
                    dopplerIdx = selectedProposalTable(proposalIdx, 3);
                    if obj.Config.SpreadWidth == 0
                        fracDoppler = dopplerIdx;
                    else
                        fracDoppler = obj.refineFractionalDoppler(state.residual, delayVal, dopplerIdx);
                        fracRefineCallCount = fracRefineCallCount + 1;
                    end

                    compositeResp = obj.ChannelBuilder.buildCompositePilotResponse(delayVal, fracDoppler);
                    expandedDict = [state.dictionaryMatrix, compositeResp];
                    pathCount = size(expandedDict, 2);
                    gramMatrix = expandedDict' * expandedDict;
                    if ompRegParam > 0
                        jointGains = (gramMatrix + ompRegParam * eye(pathCount)) \ (expandedDict' * rxSignal);
                    else
                        jointGains = expandedDict \ rxSignal;
                    end
                    residual = rxSignal - expandedDict * jointGains;
                    supportPaths = [state.supportPaths; [delayVal, fracDoppler, jointGains(end)]];
                    supportPaths(:, 3) = jointGains;

                    beamState.supportPaths = supportPaths;
                    beamState.dictionaryMatrix = expandedDict;
                    beamState.residual = residual;
                    weightedResidual = sqrt(observationWeights) .* residual;
                    beamState.weightedResidualNorm = real(weightedResidual' * weightedResidual) / numSubcarriers;
                    expandedStates(end+1) = beamState; %#ok<AGROW>
                end

                if isempty(expandedStates)
                    break;
                end
                if enableSignaturePrune
                    expandedStates = obj.pruneBeamStatesBySignature(expandedStates);
                end

                residualNormVec = [expandedStates.weightedResidualNorm].';
                [~, sortIdx] = sort(residualNormVec, 'ascend');
                keepCount = min(currentBeamWidth, numel(sortIdx));
                beamStates = expandedStates(sortIdx(1:keepCount));
                currentBestResidual = beamStates(1).weightedResidualNorm;
                ompDiag.beamStateCountPerDepth(end+1, 1) = numel(beamStates);
                ompDiag.weightedResidualNormPerDepth(end+1, 1) = currentBestResidual;
                ompDiag.beamDepthUsed = depthIdx;

                improveRatio = max(bestResidualPrev - currentBestResidual, 0) / ...
                    max(bestResidualPrev, 1e-12);
                bestResidualPrev = currentBestResidual;
                prevDepthImproveRatio = improveRatio;
                if improveRatio < minImproveRatio
                    break;
                end
            end

            if isempty(beamStates)
                [estPaths, hEffective, greedyDiag] = obj.estimateByGreedyOmp( ...
                    rxSignal, ompRegParam, observationWeights, zeros(0,3), 'beam');
                ompDiag.fracRefineCallCount = fracRefineCallCount + greedyDiag.fracRefineCallCount;
                ompDiag.weightedResidualNormPerDepth = [ ...
                    ompDiag.weightedResidualNormPerDepth; greedyDiag.weightedResidualNormPerDepth];
                ompDiag.finalSupportSignature = greedyDiag.finalSupportSignature;
                return;
            end

            residualNormVec = [beamStates.weightedResidualNorm].';
            [~, bestStateIdx] = min(residualNormVec);
            bestState = beamStates(bestStateIdx);
            [estPaths, hEffective, greedyDiag] = obj.estimateByGreedyOmp( ...
                rxSignal, ompRegParam, observationWeights, bestState.supportPaths, 'beam');
            ompDiag.fracRefineCallCount = fracRefineCallCount + greedyDiag.fracRefineCallCount;
            ompDiag.weightedResidualNormPerDepth = [ ...
                ompDiag.weightedResidualNormPerDepth; greedyDiag.weightedResidualNormPerDepth];
            ompDiag.finalSupportSignature = greedyDiag.finalSupportSignature;
        end

        % SELECTTOPCANDIDATES 返回当前残差下 top-B 个未占用候选。
        function [candidateList, totalPassingCount] = selectTopCandidates( ...
                obj, residualSignal, supportPaths, observationWeights, beamWidth, threshold)
            delayVec = 0:obj.Config.MaxDelaySamples;
            dopplerIdxVec = -obj.Config.MaxDopplerIdx:obj.Config.MaxDopplerIdx;
            candidateTable = zeros(0, 3);

            for delayVal = delayVec
                for dopplerIdx = dopplerIdxVec
                    if ~isempty(supportPaths)
                        supportSig = [round(supportPaths(:,1)), round(supportPaths(:,2))];
                        if any(supportSig(:,1) == delayVal & supportSig(:,2) == dopplerIdx)
                            continue;
                        end
                    end
                    metric = obj.scoreCandidate( ...
                        residualSignal, observationWeights, delayVal, dopplerIdx, obj.Config.NumSubcarriers, ...
                        obj.Config.LocStep, obj.Config.ChirpParam1, obj.Config.ChirpParam2, ...
                        obj.Config.PilotPos0, obj.Config.PilotSequence, obj.Config.PerPilotAmplitude);
                    if metric >= threshold
                        candidateTable(end+1, :) = [delayVal, dopplerIdx, metric]; %#ok<AGROW>
                    end
                end
            end

            if isempty(candidateTable)
                candidateList = zeros(0, 3);
                totalPassingCount = 0;
                return;
            end

            [~, sortIdx] = sort(candidateTable(:, 3), 'descend');
            keepCount = min(beamWidth, length(sortIdx));
            candidateList = candidateTable(sortIdx(1:keepCount), :);
            totalPassingCount = size(candidateTable, 1);
        end

        % REFINEDOPPLERBYDD 通过数据导向策略精炼分数 Doppler，并重估增益。
        function [refinedPaths, hEffective] = refineDopplerByDd(obj, ...
                estPaths, rxSignal, txFrameEst, ddRegParam)

            if nargin < 5
                ddRegParam = 0;
            end

            numSc = obj.Config.NumSubcarriers;
            spreadKv = obj.Config.SpreadWidth;
            numPaths = size(estPaths, 1);

            if numPaths == 0 || spreadKv == 0
                refinedPaths = estPaths;
                hEffective = obj.ChannelBuilder.buildEffectiveChannel(estPaths);
                return;
            end

            refinedPaths = estPaths;
            pathResponses = obj.ChannelBuilder.buildPathResponseDictionary(refinedPaths, txFrameEst);
            totalContrib = pathResponses * refinedPaths(:,3);

            for i = 1:numPaths
                delayVal = refinedPaths(i, 1);
                intDoppler = round(refinedPaths(i, 2));

                contribI = refinedPaths(i, 3) * pathResponses(:, i);
                residualExcl = rxSignal - (totalContrib - contribI);

                [basisVecs, kvRange] = obj.ChannelBuilder.buildDopplerBasis( ...
                    delayVal, intDoppler, txFrameEst);
                frac0 = refinedPaths(i, 2) - intDoppler;
                fracOpt = obj.secantDopplerSearch( ...
                    residualExcl, basisVecs, kvRange, numSc, frac0, -0.499, 0.499);
                refinedPaths(i, 2) = intDoppler + fracOpt;

                pathResponses(:, i) = obj.ChannelBuilder.applyPathToFrame( ...
                    delayVal, refinedPaths(i,2), txFrameEst);
                totalContrib = totalContrib - contribI + refinedPaths(i, 3) * pathResponses(:, i);
            end

            gramFull = pathResponses' * pathResponses;
            if ddRegParam > 0
                newGains = (gramFull + ddRegParam * eye(numPaths)) \ ...
                    (pathResponses' * rxSignal);
            else
                newGains = pathResponses \ rxSignal;
            end
            refinedPaths(:, 3) = newGains;
            hEffective = obj.ChannelBuilder.buildEffectiveChannel(refinedPaths);
        end

        % REFINEGAINBYDD 在固定 delay/doppler 支撑下重估路径增益。
        function [refinedPaths, hEffective] = refineGainByDd(obj, ...
                estPaths, rxSignal, txFrameEst, ddRegParam)
            if nargin < 5
                ddRegParam = 0;
            end

            numSc = obj.Config.NumSubcarriers;
            numPaths = size(estPaths, 1);

            if numPaths == 0
                refinedPaths = estPaths;
                hEffective = sparse(numSc, numSc);
                return;
            end

            fullDict = obj.ChannelBuilder.buildPathResponseDictionary(estPaths, txFrameEst);

            gramFull = fullDict' * fullDict;
            if ddRegParam > 0
                newGains = (gramFull + ddRegParam * eye(numPaths)) \ ...
                    (fullDict' * rxSignal);
            else
                newGains = fullDict \ rxSignal;
            end

            refinedPaths = estPaths;
            refinedPaths(:, 3) = newGains;
            hEffective = obj.ChannelBuilder.buildEffectiveChannel(refinedPaths);
        end

        % ESTIMATENOISETRIMMED 基于截尾残差统计估计噪声功率。
        function noisePowerEst = estimateNoiseTrimmed(~, rxSignal, effectiveChannel, txFrameEst, trimAlpha)
            if nargin < 5
                trimAlpha = 0.1;
            end
            residualPowers = abs(rxSignal - effectiveChannel * txFrameEst).^2;
            sortedPowers = sort(residualPowers, 'ascend');
            numKeep = max(floor(length(sortedPowers) * (1 - trimAlpha)), 1);
            noisePowerEst = mean(sortedPowers(1:numKeep));
        end

        % FINDSTRONGESTPEAK 在 delay-doppler 候选网格上搜索最强峰值。
        function [peakHit, peakMetric] = findStrongestPeak(obj, residualSignal, observationWeights)
            numSubcarriers = obj.Config.NumSubcarriers;
            locStep = obj.Config.LocStep;
            chirpC1 = obj.Config.ChirpParam1;
            chirpC2 = obj.Config.ChirpParam2;
            pilotPos0 = obj.Config.PilotPos0;
            pilotSeq = obj.Config.PilotSequence;
            pilotAmp = obj.Config.PerPilotAmplitude;

            delayVec = 0:obj.Config.MaxDelaySamples;
            dopplerIdxVec = -obj.Config.MaxDopplerIdx:obj.Config.MaxDopplerIdx;
            [delayGrid, dopplerGrid] = meshgrid(delayVec, dopplerIdxVec);
            delayVecFlat = delayGrid(:);
            dopplerVecFlat = dopplerGrid(:);

            metrics = zeros(length(delayVecFlat), 1);
            for candidateIdx = 1:length(delayVecFlat)
                delayVal = delayVecFlat(candidateIdx);
                dopplerIdx = dopplerVecFlat(candidateIdx);
                metrics(candidateIdx) = obj.scoreCandidate( ...
                    residualSignal, observationWeights, delayVal, dopplerIdx, ...
                    numSubcarriers, locStep, chirpC1, chirpC2, pilotPos0, pilotSeq, pilotAmp);
            end
            [peakMetric, bestIdx] = max(metrics);
            peakHit = [delayVecFlat(bestIdx), dopplerVecFlat(bestIdx)];
        end

        % PROPOSERESIDUALPATH 在当前支撑外提出一个残差路径候选。
        function [newPath, peakMetric, threshold] = proposeResidualPath( ...
                obj, residualSignal, excludedPaths, ompRegParam, thresholdScale, observationWeights)
            if nargin < 4
                ompRegParam = 0;
            end
            if nargin < 5
                thresholdScale = 0.85;
            end
            if nargin < 6 || isempty(observationWeights)
                observationWeights = ones(size(residualSignal));
            end

            numSubcarriers = obj.Config.NumSubcarriers;
            numCandidates = (obj.Config.MaxDelaySamples + 1) * ...
                (2 * obj.Config.MaxDopplerIdx + 1);
            cfarCoeff = log(numCandidates / obj.cfarFalseAlarmProb);
            pilotAmpSq = obj.Config.PilotAmplitude^2;

            weightedResidual = sqrt(observationWeights) .* residualSignal;
            residualPower = real(weightedResidual' * weightedResidual) / numSubcarriers;
            threshold = residualPower * pilotAmpSq * cfarCoeff * thresholdScale;

            [peakHit, peakMetric] = obj.findStrongestPeakExcluding( ...
                residualSignal, excludedPaths, observationWeights);
            if isempty(peakHit) || peakMetric < threshold
                newPath = [];
                return;
            end

            delayVal = peakHit(1);
            intDoppler = peakHit(2);
            if obj.Config.SpreadWidth == 0
                fracDoppler = intDoppler;
            else
                fracDoppler = obj.refineFractionalDoppler(residualSignal, delayVal, intDoppler);
            end

            compositeResp = obj.ChannelBuilder.buildCompositePilotResponse(delayVal, fracDoppler);
            respNormSq = real(compositeResp' * compositeResp);
            if respNormSq < 1e-15
                newPath = [];
                return;
            end

            if ompRegParam > 0
                gain = (respNormSq + ompRegParam) \ (compositeResp' * residualSignal);
            else
                gain = (compositeResp' * residualSignal) / respNormSq;
            end

            newPath = [delayVal, fracDoppler, gain];
        end

        % FINDSTRONGESTPEAKEXCLUDING 在排除已有支撑后搜索最强峰值。
        function [peakHit, peakMetric] = findStrongestPeakExcluding( ...
                obj, residualSignal, excludedPaths, observationWeights)
            numSubcarriers = obj.Config.NumSubcarriers;
            locStep = obj.Config.LocStep;
            chirpC1 = obj.Config.ChirpParam1;
            chirpC2 = obj.Config.ChirpParam2;
            pilotPos0 = obj.Config.PilotPos0;
            pilotSeq = obj.Config.PilotSequence;
            pilotAmp = obj.Config.PerPilotAmplitude;

            delayVec = 0:obj.Config.MaxDelaySamples;
            dopplerIdxVec = -obj.Config.MaxDopplerIdx:obj.Config.MaxDopplerIdx;

            if isempty(excludedPaths)
                excludedSig = zeros(0, 2);
            else
                excludedSig = [round(excludedPaths(:, 1)), round(excludedPaths(:, 2))];
            end

            bestMetric = -inf;
            bestHit = [];
            for delayVal = delayVec
                for dopplerIdx = dopplerIdxVec
                    if ~isempty(excludedSig) && any(excludedSig(:,1) == delayVal & excludedSig(:,2) == dopplerIdx)
                        continue;
                    end
                    metric = obj.scoreCandidate( ...
                        residualSignal, observationWeights, delayVal, dopplerIdx, ...
                        numSubcarriers, locStep, chirpC1, chirpC2, pilotPos0, pilotSeq, pilotAmp);
                    if metric > bestMetric
                        bestMetric = metric;
                        bestHit = [delayVal, dopplerIdx];
                    end
                end
            end

            peakHit = bestHit;
            peakMetric = bestMetric;
        end

        % SCORECANDIDATE 计算单个 delay-doppler 候选的加权匹配得分。
        function metric = scoreCandidate(obj, residualSignal, obsWeights, ...
                delayVal, dopplerIdx, numSubcarriers, locStep, chirpC1, chirpC2, ...
                pilotPos0, pilotSeq, pilotAmp)
            locIndex = dopplerIdx + locStep * delayVal;
            clusterWidth = max(round(obj.Config.PilotClusterWidth), 0);

            zCenter = 0;
            clusterEnergy = 0;
            responseColIdx = mod(pilotPos0 - locIndex, numSubcarriers);
            phaseVal = exp(1j * 2 * pi * (chirpC1 * delayVal^2 ...
                - pilotPos0 * delayVal / numSubcarriers + chirpC2 * pilotPos0^2 - chirpC2 * responseColIdx^2));
            centerIdx = responseColIdx + 1;
            zCenter = obsWeights(centerIdx) .* conj(pilotAmp * pilotSeq * phaseVal) .* residualSignal(centerIdx);

            if clusterWidth > 0
                clusterIdx = mod(responseColIdx + (-clusterWidth:clusterWidth), numSubcarriers) + 1;
                clusterEnergy = sum(obsWeights(clusterIdx) .* abs(residualSignal(clusterIdx)).^2);
            end

            metric = abs(zCenter)^2;
            if obj.Config.UseWeightedPilotMetric && clusterWidth > 0
                metric = metric + 0.35 * clusterEnergy / (2 * clusterWidth + 1);
            end
        end

        % REFINEFRACTIONALDOPPLER 在整数 Doppler 附近做一维分数偏移搜索。
        function estDoppler = refineFractionalDoppler(obj, residualSignal, delayVal, intDoppler)
            objFcn = @(fracShift) -obj.computePilotMetric(fracShift, residualSignal, delayVal, intDoppler);
            optOpts = optimset('TolX', 1e-4, 'Display', 'off');
            fracShiftEst = fminbnd(objFcn, -0.499, 0.499, optOpts);
            estDoppler = intDoppler + fracShiftEst;
        end

        % COMPUTEPILOTMETRIC 计算给定分数偏移下的归一化匹配度量。
        function metricVal = computePilotMetric(obj, fracShift, signalVec, delayVal, intDoppler)
            compositeResp = obj.ChannelBuilder.buildCompositePilotResponse(delayVal, intDoppler + fracShift);
            normSq = real(compositeResp' * compositeResp);
            if normSq < 1e-15
                metricVal = 0;
            else
                metricVal = abs(compositeResp' * signalVec)^2 / normSq;
            end
        end

    end

    methods (Access = private)

        % RESOLVESEARCHMODE 解析单次 OMP 调用应使用的搜索模式。
        function searchMode = resolveSearchMode(obj, searchModeOverride)
            if nargin < 2 || isempty(searchModeOverride)
                searchMode = obj.Config.EarlyPathSearchMode;
            else
                searchMode = searchModeOverride;
            end

            if ~(strcmpi(searchMode, 'greedy') || strcmpi(searchMode, 'beam'))
                error('GiFreeEstimator:InvalidEarlyPathSearchMode', ...
                    'EarlyPathSearchMode 必须为 ''greedy'' 或 ''beam''。');
            end
            searchMode = lower(searchMode);
        end

        % RESOLVEBEAMSTRATEGYMODE 解析 Beam 搜索策略模式。
        function strategyMode = resolveBeamStrategyMode(obj)
            strategyMode = lower(obj.Config.BeamStrategyMode);
            validModes = {'fastadaptive', 'performance', 'fixed'};
            if ~any(strcmp(strategyMode, validModes))
                error('GiFreeEstimator:InvalidBeamStrategyMode', ...
                    'BeamStrategyMode 必须为 ''fastAdaptive''、''performance'' 或 ''fixed''。');
            end
        end

        % RESOLVEBEAMPROBEWIDTH 返回当前策略下用于粗筛的最大 Beam 宽度。
        function probeBeamWidth = resolveBeamProbeWidth(obj, baseBeamWidth, strategyMode)
            if strcmp(strategyMode, 'fastadaptive')
                probeBeamWidth = min(baseBeamWidth + 1, obj.Config.NumPathsUpper);
            else
                probeBeamWidth = baseBeamWidth;
            end
        end

        % RESOLVEBEAMBASEDIMS 根据策略模式解析基准 Beam 宽度与深度。
        function [beamWidth, beamDepth] = resolveBeamBaseDims(obj, strategyMode)
            if strcmp(strategyMode, 'performance')
                beamWidth = max(round(obj.Config.PerformanceBeamWidth), 1);
                beamDepth = max(round(obj.Config.PerformanceBeamDepth), 1);
            else
                beamWidth = max(round(obj.Config.BeamWidth), 1);
                beamDepth = max(round(obj.Config.BeamDepth), 1);
            end
            beamDepth = min(beamDepth, obj.Config.NumPathsUpper);
        end

        % RESOLVEBEAMWIDTHFORDEPTH 根据场景难度与上一层收益决定当前层 Beam 宽度。
        function currentBeamWidth = resolveBeamWidthForDepth(obj, ...
                baseBeamWidth, probeBeamWidth, strategyMode, prevDepthImproveRatio, ...
                hasUncertainState, adaptiveImproveRatio, adaptiveSpreadThreshold)
            currentBeamWidth = baseBeamWidth;
            if ~strcmp(strategyMode, 'fastadaptive')
                return;
            end

            if obj.Config.SpreadWidth >= adaptiveSpreadThreshold || ...
                    prevDepthImproveRatio < adaptiveImproveRatio || hasUncertainState
                currentBeamWidth = probeBeamWidth;
            end
        end

        % CREATEEMPTYOMPDIAG 构造字段完整的 OMP 诊断结构体。
        function ompDiag = createEmptyOmpDiag(~, modeLabel)
            ompDiag = struct( ...
                'modeUsed', lower(modeLabel), ...
                'beamDepthUsed', 0, ...
                'beamStateCountPerDepth', zeros(0, 1), ...
                'candidateCountPerDepth', zeros(0, 1), ...
                'fracRefineCallCount', 0, ...
                'weightedResidualNormPerDepth', zeros(0, 1), ...
                'finalSupportSignature', '');
        end

        % APPENDRESIDUALNORM 记录当前加权残差范数。
        function residualVec = appendResidualNorm(~, residualVec, residualSignal, observationWeights, numSubcarriers)
            weightedResidual = sqrt(observationWeights) .* residualSignal;
            residualNorm = real(weightedResidual' * weightedResidual) / numSubcarriers;
            residualVec(end+1, 1) = residualNorm; %#ok<AGROW>
        end

        % BUILDSUPPORTSIGNATURE 生成排序后的 delay-doppler 支撑签名。
        function signature = buildSupportSignature(~, pathSet)
            if isempty(pathSet)
                signature = '';
                return;
            end

            roundedSupport = [round(pathSet(:, 1)), round(pathSet(:, 2))];
            roundedSupport = sortrows(roundedSupport, [1 2]);
            tokenList = arrayfun(@(rowIdx) sprintf('%d:%d', ...
                roundedSupport(rowIdx, 1), roundedSupport(rowIdx, 2)), ...
                1:size(roundedSupport, 1), 'UniformOutput', false);
            signature = strjoin(tokenList, '|');
        end

        % PRUNEBEAMSTATESBYSIGNATURE 对重复支撑签名仅保留残差最优分支。
        function prunedStates = pruneBeamStatesBySignature(obj, beamStates)
            if isempty(beamStates)
                prunedStates = beamStates;
                return;
            end

            prunedStates = repmat(beamStates(1), 0, 1);
            signatureList = cell(0, 1);
            for stateIdx = 1:numel(beamStates)
                stateSig = obj.buildSupportSignature(beamStates(stateIdx).supportPaths);
                sigIdx = find(strcmp(signatureList, stateSig), 1);
                if isempty(sigIdx)
                    signatureList{end+1, 1} = stateSig; %#ok<AGROW>
                    prunedStates(end+1, 1) = beamStates(stateIdx); %#ok<AGROW>
                    continue;
                end

                if beamStates(stateIdx).weightedResidualNorm < prunedStates(sigIdx).weightedResidualNorm
                    prunedStates(sigIdx) = beamStates(stateIdx);
                end
            end
        end

        % ISSTATEBEAMUNCERTAIN 判断某个 state 的前两名粗候选是否过于接近。
        function isUncertain = isStateBeamUncertain(~, candidateList, uncertainMetricRatio)
            if size(candidateList, 1) < 2
                isUncertain = false;
                return;
            end

            topMetric = candidateList(1, 3);
            secondMetric = candidateList(2, 3);
            metricRatio = topMetric / max(secondMetric, 1e-12);
            isUncertain = metricRatio < uncertainMetricRatio;
        end

        % SELECTBEAMEXPANSIONPROPOSALS 从全局粗筛候选中选取少量待精炼分支。
        function selectedProposalTable = selectBeamExpansionProposals( ...
                ~, coarseProposalTable, beamWidth, refineBudget, uncertainStateIdx)
            if isempty(coarseProposalTable)
                selectedProposalTable = zeros(0, 5);
                return;
            end

            refineBudget = max(round(refineBudget), beamWidth);
            [~, sortIdx] = sort(coarseProposalTable(:, 4), 'descend');
            sortedProposalTable = coarseProposalTable(sortIdx, :);
            uncertainStateIdx = unique(uncertainStateIdx(:));

            primaryMask = sortedProposalTable(:, 5) == 1;
            primaryProposalTable = sortedProposalTable(primaryMask, :);
            keepPrimaryCount = min(size(primaryProposalTable, 1), refineBudget);
            selectedProposalTable = primaryProposalTable(1:keepPrimaryCount, :);

            if keepPrimaryCount >= refineBudget
                return;
            end

            remainingMask = ~primaryMask;
            secondChoiceMask = remainingMask & sortedProposalTable(:, 5) == 2 & ...
                ismember(sortedProposalTable(:, 1), uncertainStateIdx);
            secondChoiceProposalTable = sortedProposalTable(secondChoiceMask, :);
            keepSecondCount = min(refineBudget - size(selectedProposalTable, 1), ...
                size(secondChoiceProposalTable, 1));
            if keepSecondCount > 0
                selectedProposalTable = [selectedProposalTable; ...
                    secondChoiceProposalTable(1:keepSecondCount, :)]; %#ok<AGROW>
            end

            if size(selectedProposalTable, 1) >= refineBudget
                return;
            end

            remainingProposalTable = sortedProposalTable(remainingMask & ~secondChoiceMask, :);
            keepRemainCount = min(refineBudget - size(selectedProposalTable, 1), ...
                size(remainingProposalTable, 1));
            if keepRemainCount > 0
                selectedProposalTable = [selectedProposalTable; remainingProposalTable(1:keepRemainCount, :)]; %#ok<AGROW>
            end
        end

        % SECANTDOPPLERSEARCH 以割线法在分数 Doppler 区间内优化度量。
        function fracOpt = secantDopplerSearch(obj, ...
                residual, basisVecs, kvRange, numSc, frac0, lb, ub)
            maxIter = 6;
            tolGrad = 1e-8;

            [f0, g0] = obj.evalMetricAndGrad(frac0, residual, basisVecs, kvRange, numSc);
            fracBest = frac0;
            fBest = f0;

            delta = 0.02;
            frac1 = frac0 + delta;
            if frac1 > ub
                frac1 = frac0 - delta;
            end
            [f1, g1] = obj.evalMetricAndGrad(frac1, residual, basisVecs, kvRange, numSc);
            if f1 > fBest
                fracBest = frac1;
                fBest = f1;
            end

            % NOTE: 以梯度零点为目标迭代，兼顾边界裁剪和最优值追踪。
            for iter = 1:maxIter
                dg = g1 - g0;
                if abs(dg) < 1e-20
                    break;
                end
                fracNew = max(lb+1e-6, min(ub-1e-6, frac1 - g1*(frac1-frac0)/dg));
                [fNew, gNew] = obj.evalMetricAndGrad(fracNew, residual, basisVecs, kvRange, numSc);
                if fNew > fBest
                    fracBest = fracNew;
                    fBest = fNew;
                end
                if abs(gNew) < tolGrad * max(abs(fNew), 1e-10)
                    break;
                end
                frac0 = frac1;
                g0 = g1;
                frac1 = fracNew;
                g1 = gNew;
            end

            fracOpt = fracBest;
        end

        % EVALMETRICANDGRAD 计算分数 Doppler 度量值及其梯度。
        function [fval, fgrad] = evalMetricAndGrad(obj, ...
                fracPart, residual, basisVecs, kvRange, numSc)
            numKv = length(kvRange);
            dCoeffs = zeros(numKv, 1);
            dDerivs = zeros(numKv, 1);
            for idx = 1:numKv
                [dCoeffs(idx), dDerivs(idx)] = obj.ChannelBuilder.computeDirichletWithDeriv( ...
                    fracPart, kvRange(idx), numSc);
            end
            s  = basisVecs * dCoeffs;
            sp = basisVecs * dDerivs;
            a  = residual' * s;
            ap = residual' * sp;
            b  = real(s' * s);
            bp = 2 * real(s' * sp);
            absA2 = real(a * conj(a));
            bSafe = max(b, 1e-15);
            fval = absA2 / bSafe;
            fgrad = (2*real(conj(a)*ap)*b - absA2*bp) / (bSafe*bSafe);
        end

    end

end


