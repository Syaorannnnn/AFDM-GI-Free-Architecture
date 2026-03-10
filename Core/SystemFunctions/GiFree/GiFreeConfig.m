classdef GiFreeConfig < handle
    % GiFreeConfig  GI-Free AFDM 系统参数配置
    %
    %   Agile-AFDM 优化:
    %     AgileQ > 0 时启用 c2 敏捷优化, 模式由 AgileMode 控制:
    %       'standard'     — 标准 PAPR 度量 (min max|s|²/mean|s|²)
    %       'pilotAware'   — 导频感知交叉项度量 (min max Δ[n])
    %       'hierarchical' — 两阶段层级搜索 (3√Q 次 FFT 替代 Q 次)
    %     AgileQ = 0 时使用默认 c2, 完全向后兼容.

    properties
        NumSubcarriers
        ModulationOrder
        MaxDelaySpread
        MaxDopplerIndex
        NumPaths
        DopplerGuard
        SpreadWidth
        PilotSnrDb
        MaxSicIterations
        NumPathsUpper
        UseFractionalDoppler
        NumPilots   = 1
        ZcRootIndex = 1
        AgileQ      = 0            % 候选数, 0=禁用, 推荐 32
        AgileMode   = 'standard'   % 'standard' | 'pilotAware' | 'hierarchical'
        C2Override  = []            % 内部使用, 由 GiFreeSystem 自动写入
    end

    properties (Dependent)
        ChirpParam1
        ChirpParam2
        C2BaseValue
        PilotAmplitude
        LocStep
        PilotZoneSpan
        PilotPositions
        DataPositions
        NumDataSymbols
        PilotSequence
        PerPilotAmplitude
    end

    methods

        function val = get.ChirpParam1(obj)
            val = (2 * (obj.MaxDopplerIndex + obj.DopplerGuard) + 1) ...
                / (2 * obj.NumSubcarriers);
        end

        function val = get.ChirpParam2(obj)
            if ~isempty(obj.C2Override)
                val = obj.C2Override;
            else
                val = 1 / (2 * obj.NumSubcarriers^2 * pi);
            end
        end

        function val = get.C2BaseValue(obj)
            val = 1 / (2 * obj.NumSubcarriers^2 * pi);
        end

        function val = get.PilotAmplitude(obj)
            val = sqrt(10^(obj.PilotSnrDb / 10));
        end

        function val = get.LocStep(obj)
            val = 2 * (obj.MaxDopplerIndex + obj.DopplerGuard) + 1;
        end

        function val = get.PilotZoneSpan(obj)
            val = obj.LocStep * (obj.MaxDelaySpread + 1);
        end

        function val = get.PilotPositions(obj)
            val = (0:obj.NumPilots - 1).' * obj.PilotZoneSpan;
        end

        function val = get.DataPositions(obj)
            allIdx  = (0:obj.NumSubcarriers - 1).';
            isPilot = ismember(allIdx, obj.PilotPositions);
            val     = allIdx(~isPilot);
        end

        function val = get.NumDataSymbols(obj)
            val = obj.NumSubcarriers - obj.NumPilots;
        end

        function val = get.PilotSequence(obj)
            K = obj.NumPilots;
            if K <= 1
                val = 1;
            else
                n = (0:K-1).';
                u = obj.ZcRootIndex;
                if mod(K, 2) == 1
                    val = exp(-1j * pi * u * n .* (n + 1) / K);
                else
                    val = exp(-1j * pi * u * n .* n / K);
                end
            end
        end

        function val = get.PerPilotAmplitude(obj)
            val = obj.PilotAmplitude / sqrt(obj.NumPilots);
        end

    end

end
