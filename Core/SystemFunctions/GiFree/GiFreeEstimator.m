classdef GiFreeEstimator < handle
    % GiFreeEstimator  GI-Free AFDM 信道估计器 (CAZAC 多导频 + 自适应门限版)
    %
    %   核心变更:
    %     (1) OMP 峰值检测: K 个导频相干合并, 含相位校正的匹配滤波
    %     (2) 自适应 CFAR 停止门限: 基于残差噪声底估计, 数学推导如下
    %     (3) 导频域分数多普勒搜索: 使用合成导频响应
    %
    %   CFAR 门限推导:
    %     相干合并度量 M = |sum_k conj(a_k·φ_k) · r[bin_k]|²
    %     H0 下 (无路径): r[bin_k] ~ CN(0, σ²_res), 各 bin 独立
    %       z = sum_k conj(a_k·φ_k)·r[bin_k] ~ CN(0, σ²_res · ||a||²)
    %       其中 ||a||² = sum_k |a_k|² = A_p²
    %       M = |z|² ~ Exp(σ²_res · A_p²)
    %     虚警率控制 (Bonferroni 校正):
    %       P_fa_per_test = P_fa_total / numCandidates
    %       γ = σ²_res · A_p² · ln(numCandidates / P_fa_total)
    %     其中 σ²_res = ||residual||² / N 从残差自适应估计

    properties (SetAccess = private)
        Config
        ChannelBuilder
    end

    properties (Constant)
        CfarPfa = 0.01     % 总虚警率
    end

    methods

        function obj = GiFreeEstimator(cfg, chBuilder)
            arguments
                cfg       (1,1) GiFreeConfig
                chBuilder (1,1) FractionalChannelBuilder
            end
            obj.Config         = cfg;
            obj.ChannelBuilder = chBuilder;
        end

        %% ==================== OMP 信道估计 (CAZAC + CFAR) ====================
        function [estPaths, hEffective] = estimateByOmp(obj, rxSignal, ompRegParam)
            if nargin < 3, ompRegParam = 0; end

            numSc          = obj.Config.NumSubcarriers;
            numPathsUpper  = obj.Config.NumPathsUpper;

            % CFAR 门限系数: ln(numCandidates / P_fa)
            numCand   = (obj.Config.MaxDelaySpread + 1) * ...
                        (2 * obj.Config.MaxDopplerIndex + 1);
            cfarCoeff = log(numCand / obj.CfarPfa);
            Ap2       = obj.Config.PilotAmplitude^2;

            estPaths        = zeros(numPathsUpper, 3);
            pilotDictionary = zeros(numSc, numPathsUpper);
            residual        = rxSignal;
            numDetected     = 0;

            for pathIdx = 1:numPathsUpper
                % 多导频相干合并峰值检测
                [peakHit, peakMetric] = obj.findStrongestPeak(residual);

                % 自适应 CFAR 停止判据
                residualPower = real(residual' * residual) / numSc;
                threshold     = residualPower * Ap2 * cfarCoeff;
                if peakMetric < threshold
                    break;
                end

                delayVal   = peakHit(1);
                intDoppler = peakHit(2);

                if obj.Config.SpreadWidth == 0
                    fracDoppler = intDoppler;
                else
                    fracDoppler = obj.refineFractionalDoppler( ...
                        residual, delayVal, intDoppler);
                end

                % 使用合成导频响应构造字典原子
                compositeResp = obj.ChannelBuilder.buildCompositePilotResponse( ...
                    delayVal, fracDoppler);
                pilotDictionary(:, pathIdx) = compositeResp;
                numDetected = pathIdx;

                % 联合 LS 增益估计
                activeDict = pilotDictionary(:, 1:pathIdx);
                gramMatrix = activeDict' * activeDict;
                if ompRegParam > 0
                    jointGains = (gramMatrix + ompRegParam * eye(pathIdx)) \ ...
                        (activeDict' * rxSignal);
                else
                    jointGains = activeDict \ rxSignal;
                end
                residual = rxSignal - activeDict * jointGains;

                estPaths(pathIdx, :) = [delayVal, fracDoppler, jointGains(pathIdx)];
            end

            % 修剪至实际检测到的路径数
            estPaths = estPaths(1:numDetected, :);

            if numDetected > 0
                estPaths(:, 3) = jointGains;
                hEffective = obj.ChannelBuilder.buildEffectiveChannel(estPaths);
            else
                hEffective = sparse(numSc, numSc);
            end
        end

        %% ==================== DD 多普勒精炼 ====================
        function [refinedPaths, hEffective] = refineDopplerByDd(obj, ...
                estPaths, rxSignal, txFrameEst, ddRegParam)

            if nargin < 5, ddRegParam = 0; end

            numSc    = obj.Config.NumSubcarriers;
            spreadKv = obj.Config.SpreadWidth;
            chirpC1  = obj.Config.ChirpParam1;
            chirpC2  = obj.Config.ChirpParam2;
            numPaths = size(estPaths, 1);

            if numPaths == 0 || spreadKv == 0
                refinedPaths = estPaths;
                hEffective   = obj.ChannelBuilder.buildEffectiveChannel(estPaths);
                return;
            end

            refinedPaths = estPaths;
            kvRange      = (-spreadKv:spreadKv)';
            numKv        = length(kvRange);
            pVec         = (0:numSc-1)';
            twoPiOverN   = 2 * pi / numSc;

            pathResponses = zeros(numSc, numPaths);
            for i = 1:numPaths
                hSingle = obj.ChannelBuilder.buildPathMatrix( ...
                    refinedPaths(i,1), refinedPaths(i,2));
                pathResponses(:,i) = hSingle * txFrameEst;
            end
            totalContrib = pathResponses * refinedPaths(:,3);

            for i = 1:numPaths
                delayVal   = refinedPaths(i, 1);
                intDoppler = round(refinedPaths(i, 2));

                contribI     = refinedPaths(i, 3) * pathResponses(:, i);
                residualExcl = rxSignal - (totalContrib - contribI);

                alphaInt   = intDoppler;
                locStep    = round(2 * numSc * chirpC1);
                locIndex   = alphaInt + locStep * delayVal;
                phaseConst = numSc * chirpC1 * delayVal^2;
                c2PSq      = chirpC2 * pVec .* pVec;

                basisVecs = zeros(numSc, numKv);
                for idx = 1:numKv
                    qVec = mod(pVec + locIndex + kvRange(idx), numSc);
                    phaseArg = twoPiOverN * (phaseConst - qVec * delayVal ...
                        + numSc * chirpC2 * (qVec .* qVec) - numSc * c2PSq);
                    basisVecs(:, idx) = exp(1j * phaseArg) .* txFrameEst(qVec + 1);
                end

                frac0   = refinedPaths(i, 2) - intDoppler;
                fracOpt = GiFreeEstimator.secantDopplerSearch( ...
                    residualExcl, basisVecs, kvRange, numSc, frac0, -0.499, 0.499);
                refinedPaths(i, 2) = intDoppler + fracOpt;

                hNew = obj.ChannelBuilder.buildPathMatrix(delayVal, refinedPaths(i,2));
                pathResponses(:, i) = hNew * txFrameEst;
                totalContrib = totalContrib - contribI ...
                    + refinedPaths(i, 3) * pathResponses(:, i);
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

        %% ==================== DD 增益重估 ====================
        function [refinedPaths, hEffective] = refineGainByDd(obj, ...
                estPaths, rxSignal, txFrameEst, ddRegParam)
            if nargin < 5, ddRegParam = 0; end

            numSc    = obj.Config.NumSubcarriers;
            numPaths = size(estPaths, 1);

            if numPaths == 0
                refinedPaths = estPaths;
                hEffective   = sparse(numSc, numSc);
                return;
            end

            fullDict = zeros(numSc, numPaths);
            for i = 1:numPaths
                hSingle = obj.ChannelBuilder.buildPathMatrix( ...
                    estPaths(i,1), estPaths(i,2));
                fullDict(:, i) = hSingle * txFrameEst;
            end

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

        %% ==================== 截断均值噪声估计 ====================
        function noisePowerEst = estimateNoiseTrimmed(~, rxSignal, hEff, txFrameEst, trimAlpha)
            if nargin < 5, trimAlpha = 0.1; end
            residualPowers = abs(rxSignal - hEff * txFrameEst).^2;
            sortedPowers   = sort(residualPowers, 'ascend');
            numKeep        = max(floor(length(sortedPowers) * (1 - trimAlpha)), 1);
            noisePowerEst  = mean(sortedPowers(1:numKeep));
        end

    end

    methods (Access = private)

        function [peakHit, peakMetric] = findStrongestPeak(obj, residualSignal)
            % 多导频相干合并峰值检测 (含相位校正的匹配滤波)
            %
            %   对每个候选 (l, α), 计算:
            %     z = sum_k conj(a_k · φ_k(l,α)) · residual[bin_k(l,α)]
            %   其中 φ_k = exp(j·2π·(c1·l² - m_k·l/N + c2·m_k² - c2·p_k²))
            %   为第 k 个导频在 (l,α) 下的信道响应相位.

            numSc   = obj.Config.NumSubcarriers;
            locStep = obj.Config.LocStep;
            chirpC1 = obj.Config.ChirpParam1;
            chirpC2 = obj.Config.ChirpParam2;
            K       = obj.Config.NumPilots;
            pilotPos = obj.Config.PilotPositions;
            pilotSeq = obj.Config.PilotSequence;
            ampPP   = obj.Config.PerPilotAmplitude;

            delays   = 0:obj.Config.MaxDelaySpread;
            dopplers = -obj.Config.MaxDopplerIndex:obj.Config.MaxDopplerIndex;
            [dGrid, nuGrid] = meshgrid(delays, dopplers);
            dVec  = dGrid(:);
            nuVec = nuGrid(:);

            metrics = zeros(length(dVec), 1);
            for g = 1:length(dVec)
                l = dVec(g);
                alpha = nuVec(g);
                locIndex = alpha + locStep * l;
                z = 0;
                for k = 1:K
                    mk = pilotPos(k);
                    pk = mod(mk - locIndex, numSc);
                    phaseK = exp(1j * 2 * pi * (chirpC1 * l^2 ...
                        - mk * l / numSc + chirpC2 * mk^2 - chirpC2 * pk^2));
                    z = z + conj(ampPP * pilotSeq(k) * phaseK) ...
                        * residualSignal(pk + 1);
                end
                metrics(g) = abs(z)^2;
            end
            [peakMetric, bestIdx] = max(metrics);
            peakHit = [dVec(bestIdx), nuVec(bestIdx)];
        end

        function estDoppler = refineFractionalDoppler(obj, residualSignal, delayVal, intDoppler)
            % 使用合成导频响应进行分数多普勒搜索
            objFcn = @(fracShift) -obj.computePilotMetric( ...
                fracShift, residualSignal, delayVal, intDoppler);
            optOpts      = optimset('TolX', 1e-4, 'Display', 'off');
            fracShiftEst = fminbnd(objFcn, -0.499, 0.499, optOpts);
            estDoppler   = intDoppler + fracShiftEst;
        end

        function metricVal = computePilotMetric(obj, fracShift, signalVec, delayVal, intDoppler)
            % 合成导频响应度量 (自动适配 K 个导频)
            compositeResp = obj.ChannelBuilder.buildCompositePilotResponse( ...
                delayVal, intDoppler + fracShift);
            normSq = real(compositeResp' * compositeResp);
            if normSq < 1e-15
                metricVal = 0;
            else
                metricVal = abs(compositeResp' * signalVec)^2 / normSq;
            end
        end

    end

    %% ==================== 静态: 割线法 + 解析梯度 ====================
    methods (Static)

        function fracOpt = secantDopplerSearch( ...
                residual, basisVecs, kvRange, numSc, frac0, lb, ub)
            maxIter = 6;
            tolGrad = 1e-8;

            [f0, g0] = GiFreeEstimator.evalMetricAndGrad( ...
                frac0, residual, basisVecs, kvRange, numSc);
            fracBest = frac0;  fBest = f0;

            delta = 0.02;
            frac1 = frac0 + delta;
            if frac1 > ub, frac1 = frac0 - delta; end
            [f1, g1] = GiFreeEstimator.evalMetricAndGrad( ...
                frac1, residual, basisVecs, kvRange, numSc);
            if f1 > fBest, fracBest = frac1; fBest = f1; end

            for iter = 1:maxIter
                dg = g1 - g0;
                if abs(dg) < 1e-20, break; end
                fracNew = max(lb+1e-6, min(ub-1e-6, frac1 - g1*(frac1-frac0)/dg));
                [fNew, gNew] = GiFreeEstimator.evalMetricAndGrad( ...
                    fracNew, residual, basisVecs, kvRange, numSc);
                if fNew > fBest, fracBest = fracNew; fBest = fNew; end
                if abs(gNew) < tolGrad * max(abs(fNew), 1e-10), break; end
                frac0 = frac1; g0 = g1;
                frac1 = fracNew; g1 = gNew;
            end

            fracOpt = fracBest;
        end

        function [fval, fgrad] = evalMetricAndGrad( ...
                fracPart, residual, basisVecs, kvRange, numSc)
            numKv   = length(kvRange);
            dCoeffs = zeros(numKv, 1);
            dDerivs = zeros(numKv, 1);
            for idx = 1:numKv
                [dCoeffs(idx), dDerivs(idx)] = ...
                    FractionalChannelBuilder.dirichletWithDeriv( ...
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
            fval  = absA2 / bSafe;
            fgrad = (2*real(conj(a)*ap)*b - absA2*bp) / (bSafe*bSafe);
        end

    end

end
