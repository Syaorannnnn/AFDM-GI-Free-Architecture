classdef AgileC2Optimizer
    % AgileC2Optimizer  Agile-AFDM c2 敏捷优化器 (三模式向量化版)
    %
    %   三种选择模式:
    %
    %   (1) selectForBlock — 标准 PAPR 度量
    %       min_{c2} max_n |s[n]|² / mean_n |s[n]|²
    %       直接优化完整时域信号的 PAPR.
    %
    %   (2) selectPilotAware — 导频感知交叉项度量 (本文贡献)
    %       利用导频位于 m=0 的结构性质:
    %         IDAFT 前 chirp 为 exp(+j2πc2·m²), 当 m=0 时恒等于 1.
    %         → 经 IFFT 归一化后, 导频时域贡献 y_p = A_p/√N (实数常量, 不依赖 c2)
    %         → 完整信号 y[n] = A_p/√N + y_d[n; c2]
    %         → |s[n]|² = (A_p/√N)² + Δ[n], 其中 Δ[n] = |y_d|² + 2(A_p/√N)·Re(y_d)
    %       优化目标: min_{c2} max_n Δ[n] (最小化导频功率底噪之上的峰值)
    %
    %       理论优势:
    %         标准 PAPR = (A_p²/N + max Δ) / (A_p²/N + σ²_d)
    %           当 A_p²/N >> σ²_d 时, PAPR → 1 + max Δ/A_p²·N, 灵敏度 → 0
    %         导频感知度量 Δ 中交叉项 2(A_p/√N)·Re(y_d) 随 A_p 线性增长
    %           → 在高导频功率区间, 交叉项主导 Δ, 导频感知的候选区分度更高
    %         交叉点: 当 A_p/√N ≈ σ_d (数据 RMS) 时, 两种模式开始分化
    %           以 N=512, QPSK, dataSnr=15dB 为例: σ_d ≈ 5.6, 故交叉点 ≈ 42 dB
    %
    %   (3) selectHierarchical — 两阶段层级搜索 (本文贡献)
    %       Stage 1: 从 Q 个候选中等距抽取 Q₁ ≈ √Q 个做粗搜
    %       Stage 2: 在最优粗候选的邻域中取 ~2√Q 个做精搜
    %       总 FFT 次数: ≈ 3√Q  (vs 穷举 Q)
    %       → Q=64 时 24 次 vs 64 次 (2.7×加速)
    %       → Q=256 时 48 次 vs 256 次 (5.3×加速)
    %
    %   向量化:
    %     所有模式均使用 N×Q 矩阵批量 IFFT, 不含循环.

    methods (Static)

        function candidates = generateCandidates(N, Q)
            % generateCandidates  生成 Q 个等间距 c2 候选值
            %   c2_q = 1/(2N²π) + q/(Q·N),  q = 0,...,Q-1
            %   覆盖 [c2_base, c2_base + 1/N) 一个 PAPR 完整周期
            if nargin < 2, Q = 32; end
            c2Base     = 1 / (2 * N^2 * pi);
            candidates = c2Base + (0:Q-1).' / (Q * N);
        end

        %% ============ 模式 1: 标准 PAPR ============
        function [bestC2, bestIdx, allPapr] = selectForBlock(daftFrame, candidates)
            % selectForBlock  标准 PAPR 度量 (向量化批量 IFFT)
            %
            %   度量: PAPR = max|s|² / mean|s|²
            %   后 chirp c1 不影响 |s[n]|², 无需计算.

            N = length(daftFrame);

            % 前 chirp: conj(chirpC2) = exp(+j·2π·c2_q·m²)
            chirpPhases = (2 * pi * (0:N-1).'.^2) * candidates.';    % N×Q
            preChirped  = daftFrame .* exp(1j * chirpPhases);          % N×Q

            % 批量 IFFT
            td = ifft(preChirped) * sqrt(N);                           % N×Q

            % PAPR (逐列)
            pw = abs(td).^2;
            allPapr = 10 * log10(max(pw,[],1) ./ mean(pw,1)).';

            [~, bestIdx] = min(allPapr);
            bestC2 = candidates(bestIdx);
        end

        %% ============ 模式 2: 导频感知交叉项度量 ============
        function [bestC2, bestIdx, allDelta] = selectPilotAware(daftFrame, candidates, pilotAmp)
            % selectPilotAware  导频感知度量 (向量化批量 IFFT)
            %
            %   数学推导:
            %     导频位于 m=0, 前 chirp 因子 exp(+j2πc2·0²)=1.
            %     经 IFFT (含 √N 归一化) 后, 导频时域贡献为:
            %       y_p[n] = A_p/√N  (实数常量, 不依赖 c2, 不依赖 n)
            %     完整信号 y[n] = A_p/√N + y_d[n; c2], 于是:
            %       |y[n]|² = (A_p/√N)² + Δ[n]
            %       Δ[n] = |y_d|² + 2·(A_p/√N)·Re(y_d)
            %
            %   优化目标: min_{c2} max_n Δ[n]
            %     → 最小化导频功率底噪之上的峰值超额功率
            %
            %   灵敏度分析:
            %     标准 PAPR ≈ 1 + max Δ/(A_p²/N + σ²_d), 灵敏度 → 0 当 A_p → ∞
            %     Δ 中交叉项 2(A_p/√N)·Re(y_d) 随 A_p 线性增长
            %     → 高导频功率时, 导频感知度量对 c2 的区分度优于标准 PAPR
            %
            %   输入:
            %     daftFrame  : N×1 DAFT 域帧 (含导频)
            %     candidates : Q×1 候选 c2
            %     pilotAmp   : 导频幅度 A_p (实数, = cfg.PilotAmplitude)

            N = length(daftFrame);

            % 提取数据帧 (导频位于 1-based 索引 1, 置零)
            dataFrame    = daftFrame;
            dataFrame(1) = 0;

            % 前 chirp + 批量 IFFT → y_d[n; c2] (仅数据部分)
            chirpPhases = (2 * pi * (0:N-1).'.^2) * candidates.';
            yd = ifft(dataFrame .* exp(1j * chirpPhases)) * sqrt(N);   % N×Q

            % 导频时域幅度: A_p/√N (IFFT 归一化)
            pilotTd = pilotAmp / sqrt(N);

            % 导频感知交叉项度量: Δ[n] = |y_d|² + 2·(A_p/√N)·Re(y_d)
            delta = abs(yd).^2 + 2 * pilotTd * real(yd);              % N×Q

            % 峰值超额功率 (逐列)
            allDelta = max(delta, [], 1).';                            % Q×1

            [~, bestIdx] = min(allDelta);
            bestC2 = candidates(bestIdx);
        end

        %% ============ 模式 3: 两阶段层级搜索 ============
        function [bestC2, bestIdx, info] = selectHierarchical(daftFrame, candidates)
            % selectHierarchical  两阶段层级搜索
            %
            %   原理:
            %     PAPR 关于 c2 在一个周期内平滑变化, 相邻候选强相关.
            %     因此用粗搜定位最优区域 + 精搜定位最优点, 可以跳过大量中间候选.
            %
            %   复杂度:
            %     Stage 1: √Q 次 FFT (粗搜)
            %     Stage 2: ~2√Q 次 FFT (精搜)
            %     总计 ≈ 3√Q 次 FFT  (vs 穷举 Q 次)

            Q  = length(candidates);
            Q1 = ceil(sqrt(Q));
            coarseStep = max(floor(Q / Q1), 1);

            % ---- Stage 1: 粗搜 (等距子集) ----
            coarseIdx = 1:coarseStep:Q;
            [~, bestCoarseLocal] = AgileC2Optimizer.selectForBlock( ...
                daftFrame, candidates(coarseIdx));
            bestCoarseGlobal = coarseIdx(bestCoarseLocal);

            % ---- Stage 2: 精搜 (最优粗候选 ± coarseStep 邻域) ----
            fineStart = max(1, bestCoarseGlobal - coarseStep);
            fineEnd   = min(Q, bestCoarseGlobal + coarseStep);
            fineIdx   = fineStart:fineEnd;
            [bestC2, bestFineLocal] = AgileC2Optimizer.selectForBlock( ...
                daftFrame, candidates(fineIdx));
            bestIdx = fineIdx(bestFineLocal);

            info.coarseFFTs = length(coarseIdx);
            info.fineFFTs   = length(fineIdx);
            info.totalFFTs  = info.coarseFFTs + info.fineFFTs;
        end

    end

end