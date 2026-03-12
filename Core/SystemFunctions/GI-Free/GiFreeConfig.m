classdef GiFreeConfig < handle
    % GiFreeConfig  GI-Free AFDM 系统参数配置
    %
    %   handle 类: 所有持有引用的模块共享同一份配置.
    %   所有必要参数不做内部初始化, 一律由外部赋值.
    %   Dependent 属性根据基本参数自动派生.
    %
    %   用法:
    %     cfg = GiFreeConfig();
    %     cfg.NumSubcarriers       = 256;
    %     cfg.ModulationOrder      = 4;
    %     cfg.MaxDelaySpread       = 2;
    %     cfg.MaxDopplerIndex      = 2;
    %     cfg.NumPaths             = 3;
    %     cfg.DopplerGuard         = 4;
    %     cfg.SpreadWidth          = 4;
    %     cfg.PilotSnrDb           = 45;
    %     cfg.MaxSicIterations     = 10;
    %     cfg.NumPathsUpper        = 4;
    %     cfg.UseFractionalDoppler = true;

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
    end

    properties (Dependent)
        ChirpParam1
        ChirpParam2
        PilotAmplitude
        LocStep
    end

    methods

        function val = get.ChirpParam1(obj)
            val = (2 * (obj.MaxDopplerIndex + obj.DopplerGuard) + 1) ...
                / (2 * obj.NumSubcarriers);
        end

        function val = get.ChirpParam2(obj)
            val = 1 / (2 * obj.NumSubcarriers^2 * pi);
        end

        function val = get.PilotAmplitude(obj)
            val = sqrt(10^(obj.PilotSnrDb / 10));
        end

        function val = get.LocStep(obj)
            val = 2 * (obj.MaxDopplerIndex + obj.DopplerGuard) + 1;
        end

    end

end
