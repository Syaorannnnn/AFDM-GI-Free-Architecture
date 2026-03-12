classdef GiFreeEstimator < handle
    % GiFreeEstimator  GI-Free AFDM 信道估计器 (统一优化版)
    %
    %   整合三大信道估计能力:
    %     (1) 正则化 OMP: 导频域路径检测 + 分数多普勒粗搜
    %     (2) DD 多普勒精炼: 基预计算 + 解析梯度割线法 (5x 加速)
    %     (3) DD 增益重估: 全帧信号联合 LS
    %
    %   优化技术:
    %     - 向量化整数格点峰值检测 (消除双层 for)
    %     - 基响应预计算 (消除 fminbnd 内层 buildPathMatrix)
    %     - Dirichlet 核解析导数 + 割线法 (3-5 步替代 fminbnd 8-12 次)
    %     - 增量残差更新 (O(N) 替代 P-1 次稀疏 matvec)

    properties (SetAccess = private)
        Config
        ChannelBuilder
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

        %% ==================== OMP 信道估计 ====================
        function [estPaths, hEffective] = estimateByOmp(obj, rxSignal, ompRegParam)
            if nargin < 3, ompRegParam = 0; end

            numSc          = obj.Config.NumSubcarriers;
            numPathsUpper  = obj.Config.NumPathsUpper;
            pilotAmplitude = obj.Config.PilotAmplitude;

            estPaths        = zeros(numPathsUpper, 3);
            pilotDictionary = zeros(numSc, numPathsUpper);
            residual        = rxSignal;

            for pathIdx = 1:numPathsUpper
                peakHit    = obj.findStrongestPeak(residual);
                delayVal   = peakHit(1);
                intDoppler = peakHit(2);

                if obj.Config.SpreadWidth == 0
                    fracDoppler = intDoppler;
                else
                    fracDoppler = obj.refineFractionalDoppler( ...
                        residual, delayVal, intDoppler);
                end

                pilotResponse = obj.ChannelBuilder.buildPilotResponseVector( ...
                    delayVal, fracDoppler);
                pilotDictionary(:, pathIdx) = pilotResponse * pilotAmplitude;

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

            % 最终增益由联合 LS 确定
            estPaths(:, 3) = jointGains;
            hEffective = obj.ChannelBuilder.buildEffectiveChannel(estPaths);
        end

        %% ==================== DD 多普勒精炼 ====================
        function [refinedPaths, hEffective] = refineDopplerByDd(obj, ...
                estPaths, rxSignal, txFrameEst, ddRegParam)
            % refineDopplerByDd  基预计算 + 割线法版多普勒 DD 精炼
            %
            %   核心优化:
            %     H(d,ν)*x̂ = Σ_{kv} d(ν,kv) · basis_kv
            %     basis_kv 在搜索区间内固定, 仅 Dirichlet 系数随 ν 变化.
            %     割线法利用解析梯度在 3-5 步内收敛.

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

            % 预计算所有路径的响应向量 (增量残差用)
            pathResponses = zeros(numSc, numPaths);
            for i = 1:numPaths
                hSingle = obj.ChannelBuilder.buildPathMatrix( ...
                    refinedPaths(i,1), refinedPaths(i,2));
                pathResponses(:,i) = hSingle * txFrameEst;
            end
            totalContrib = pathResponses * refinedPaths(:,3);

            % 逐路径多普勒精炼
            for i = 1:numPaths
                delayVal   = refinedPaths(i, 1);
                intDoppler = round(refinedPaths(i, 2));

                % 增量残差 (O(N))
                contribI     = refinedPaths(i, 3) * pathResponses(:, i);
                residualExcl = rxSignal - (totalContrib - contribI);

                % 基响应预计算: 对每个 kvShift 计算固定的 Φ_kv * txFrameEst
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

                % 割线法搜索最优分数多普勒
                frac0   = refinedPaths(i, 2) - intDoppler;
                fracOpt = GiFreeEstimator.secantDopplerSearch( ...
                    residualExcl, basisVecs, kvRange, numSc, frac0, -0.499, 0.499);
                refinedPaths(i, 2) = intDoppler + fracOpt;

                % 更新缓存
                hNew = obj.ChannelBuilder.buildPathMatrix(delayVal, refinedPaths(i,2));
                pathResponses(:, i) = hNew * txFrameEst;
                totalContrib = totalContrib - contribI ...
                    + refinedPaths(i, 3) * pathResponses(:, i);
            end

            % 联合重估增益
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

        function peakHit = findStrongestPeak(obj, residualSignal)
            % 向量化格点搜索 (消除双层 for)
            numSc    = obj.Config.NumSubcarriers;
            locStep  = obj.Config.LocStep;
            delays   = 0:obj.Config.MaxDelaySpread;
            dopplers = -obj.Config.MaxDopplerIndex:obj.Config.MaxDopplerIndex;
            [dGrid, nuGrid] = meshgrid(delays, dopplers);
            dVec  = dGrid(:);
            nuVec = nuGrid(:);
            centerBins = mod(-(nuVec + locStep * dVec), numSc) + 1;
            [~, bestIdx] = max(abs(residualSignal(centerBins)));
            peakHit = [dVec(bestIdx), nuVec(bestIdx)];
        end

        function estDoppler = refineFractionalDoppler(obj, residualSignal, delayVal, intDoppler)
            objFcn = @(fracShift) -obj.computePilotMetric( ...
                fracShift, residualSignal, delayVal, intDoppler);
            optOpts      = optimset('TolX', 1e-4, 'Display', 'off');
            fracShiftEst = fminbnd(objFcn, -0.499, 0.499, optOpts);
            estDoppler   = intDoppler + fracShiftEst;
        end

        function metricVal = computePilotMetric(obj, fracShift, signalVec, delayVal, intDoppler)
            pilotResp = obj.ChannelBuilder.buildPilotResponseVector( ...
                delayVal, intDoppler + fracShift);
            normSq = real(pilotResp' * pilotResp);
            if normSq < 1e-15
                metricVal = 0;
            else
                metricVal = abs(pilotResp' * signalVec)^2 / normSq;
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
