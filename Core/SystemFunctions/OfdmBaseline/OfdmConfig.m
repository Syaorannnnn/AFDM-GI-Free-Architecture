classdef OfdmConfig < handle
% OfdmConfig: OFDM 基线配置对象与索引派生。
    %
    %                PilotIndices/DataIndices

    properties (Access = public)
        DftSize
        CpLength
        ModulationOrder
        PilotSpacing
        MaxDelaySamples
        NumPaths
        MaxDopplerIdx
    end

    properties (SetAccess = private)
        TotalFrameLength
        NumPilots
        PilotPos1Internal
        DataPos1Internal
        NumDataCarriers
    end

    properties (Dependent)
        BitsPerSymbol

        PilotPos1
        DataPos1
        PilotPos0
        DataPos0

        MaxPathDelays
        MaxNormDoppler
        PilotIndices
        DataIndices
    end

    methods

        function set.DftSize(obj, val)
            obj.DftSize = val;
            obj.updateDerivedParams();
        end

        function set.CpLength(obj, val)
            obj.CpLength = val;
            obj.updateDerivedParams();
        end

        function set.ModulationOrder(obj, val)
            obj.ModulationOrder = val;
            obj.updateDerivedParams();
        end

        function set.PilotSpacing(obj, val)
            obj.PilotSpacing = val;
            obj.updateDerivedParams();
        end

        function set.MaxPathDelays(obj, val)
            OfdmConfig.warnDeprecatedAlias('MaxPathDelays', 'MaxDelaySamples');
            obj.MaxDelaySamples = val;
        end

        function val = get.MaxPathDelays(obj)
            val = obj.MaxDelaySamples;
        end

        function set.MaxNormDoppler(obj, val)
            OfdmConfig.warnDeprecatedAlias('MaxNormDoppler', 'MaxDopplerIdx');
            obj.MaxDopplerIdx = val;
        end

        function val = get.MaxNormDoppler(obj)
            val = obj.MaxDopplerIdx;
        end

        function val = get.BitsPerSymbol(obj)
            val = log2(obj.ModulationOrder);
        end

        function val = get.PilotPos1(obj)
            val = obj.PilotPos1Internal;
        end

        function val = get.DataPos1(obj)
            val = obj.DataPos1Internal;
        end

        function val = get.PilotPos0(obj)
            val = obj.PilotPos1Internal - 1;
        end

        function val = get.DataPos0(obj)
            val = obj.DataPos1Internal - 1;
        end

        function val = get.PilotIndices(obj)
            OfdmConfig.warnDeprecatedAlias('PilotIndices', 'PilotPos1');
            val = obj.PilotPos1;
        end

        function val = get.DataIndices(obj)
            OfdmConfig.warnDeprecatedAlias('DataIndices', 'DataPos1');
            val = obj.DataPos1;
        end

        % updateDerivedParams: 根据基础参数更新所有派生量。
        function updateDerivedParams(obj)
            if isempty(obj.DftSize) || isempty(obj.CpLength) || isempty(obj.PilotSpacing)
                return;
            end

            dftSize = obj.DftSize;
            pilotSpacing = obj.PilotSpacing;

            obj.TotalFrameLength = dftSize + obj.CpLength;

            obj.PilotPos1Internal = (1:pilotSpacing:dftSize).';
            obj.NumPilots = numel(obj.PilotPos1Internal);

            allPos1 = (1:dftSize).';
            obj.DataPos1Internal = setdiff(allPos1, obj.PilotPos1Internal);
            obj.NumDataCarriers = numel(obj.DataPos1Internal);

            maxAllowed = floor(dftSize / (obj.CpLength + 1));
            if pilotSpacing > maxAllowed
                warning('OfdmConfig:PilotSpacing', ...
                    'PilotSpacing=%d is too large for current CP length. Suggested upper bound is %d.', ...
                    pilotSpacing, maxAllowed);
            end
        end

    end

    methods (Static, Access = private)

        % warnDeprecatedAlias: 兼容旧字段名并给出一次性弃用告警。
        function warnDeprecatedAlias(oldName, newName)
            persistent warnedMap;
            if isempty(warnedMap)
                warnedMap = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            end

            warnKey = [oldName '->' newName];
            if ~isKey(warnedMap, warnKey)
                warning('OfdmConfig:DeprecatedAlias', ...
                    'Property "%s" is deprecated. Use "%s" instead.', oldName, newName);
                warnedMap(warnKey) = true;
            end
        end

    end

end

