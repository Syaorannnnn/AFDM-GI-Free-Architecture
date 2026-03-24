classdef OfdmConfig < handle
    % OfdmConfig: 标准 CP-OFDM 梳状导频系统参数配置
    %
    % 设计理念:
    %   以 DftSize (DFT长度) 和 CpLength (循环前缀) 为顶层输入,
    %   TotalFrameLength = DftSize + CpLength 与 EpConfig.TotalSubcarriers 一致,
    %   确保与 AFDM-EP 的空口帧长公平对比.
    %
    % 导频方案: 等间距梳状导频 (comb-type), 间距 PilotSpacing

    properties (Access = public)
        DftSize                 % N: DFT 长度 (匹配 EpConfig.NumDataSubcarriers)
        CpLength                % 循环前缀长度 (= MaxPathDelays)
        ModulationOrder         % QAM 阶数 (4, 16, 64, ...)
        PilotSpacing            % D_p: 梳状导频间距
        MaxPathDelays           % 最大信道时延扩展 (样本数)
        NumPaths                % 信道路径数
        MaxNormDoppler          % 最大归一化多普勒指数
    end

    properties (SetAccess = private)
        TotalFrameLength        % N + CpLength (空口帧长)
        NumPilots               % 导频数量
        PilotIndices            % 导频子载波索引 (1-based)
        DataIndices             % 数据子载波索引 (1-based)
        NumDataCarriers         % 数据子载波数量
    end

    properties (Dependent)
        BitsPerSymbol           % log2(ModulationOrder), 派生属性
    end

    methods

        function set.DftSize(obj, val)
            obj.DftSize = val; obj.updateDerivedParams();
        end

        function set.CpLength(obj, val)
            obj.CpLength = val; obj.updateDerivedParams();
        end

        function set.ModulationOrder(obj, val)
            obj.ModulationOrder = val;
            obj.updateDerivedParams();
        end

        function set.PilotSpacing(obj, val)
            obj.PilotSpacing = val; obj.updateDerivedParams();
        end

        function val = get.BitsPerSymbol(obj)
            val = log2(obj.ModulationOrder);
        end

        function updateDerivedParams(obj)
            if isempty(obj.DftSize) || isempty(obj.CpLength) || ...
               isempty(obj.PilotSpacing)
                return;
            end

            N  = obj.DftSize;
            Dp = obj.PilotSpacing;

            obj.TotalFrameLength = N + obj.CpLength;

            % 梳状导频: 子载波 1, 1+Dp, 1+2*Dp, ...
            obj.PilotIndices  = (1:Dp:N).';
            obj.NumPilots     = numel(obj.PilotIndices);

            allIdx = (1:N).';
            obj.DataIndices     = setdiff(allIdx, obj.PilotIndices);
            obj.NumDataCarriers = numel(obj.DataIndices);

            % 导频间距可分辨验证: D_p <= N / (L_cp + 1)
            maxAllowed = floor(N / (obj.CpLength + 1));
            if Dp > maxAllowed
                warning('OfdmConfig:PilotSpacing', ...
                    'D_p=%d 超过延迟可分辨极限 N/(L+1)=%d, 信道估计会产生混叠', ...
                    Dp, maxAllowed);
            end
        end

    end

end
