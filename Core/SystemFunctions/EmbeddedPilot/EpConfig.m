classdef EpConfig < handle
% EpConfig: Embedded Pilot 链路的配置对象及派生参数。
    %
    %

    properties (Access = public)
        BitsPerSymbol
        ModulationOrder

        TotalSubcarriers
        MaxDelaySamples
        MaxDopplerIdx
        DopplerGuard
        NumPaths

        PilotSnrDb

        ManualZeroPaddingLength
    end

    properties (Dependent)
        MaxPathDelays
        MaxNormDoppler
        PilotSnr

        PilotPos1
        PilotPos0
        DataPos1
        DataPos0
    end

    properties (SetAccess = private)
        NumDataSubcarriers
        PrefixLength
        ChirpParam1
        ChirpParam2
        ZeroPaddingLength
        TheoreticalZeroPaddingLength
        PilotSequenceLength
        PilotIndex
        PilotSequence
        NumActiveCarriers
        ActiveIndices
    end

    methods

        function set.TotalSubcarriers(obj, val)
            obj.TotalSubcarriers = val;
            obj.updateDerivedParams();
        end

        function set.MaxDelaySamples(obj, val)
            obj.MaxDelaySamples = val;
            obj.updateDerivedParams();
        end

        function set.MaxDopplerIdx(obj, val)
            obj.MaxDopplerIdx = val;
            obj.updateDerivedParams();
        end

        function set.ManualZeroPaddingLength(obj, val)
            obj.ManualZeroPaddingLength = val;
            obj.updateDerivedParams();
        end

        function set.MaxPathDelays(obj, val)
            EpConfig.warnDeprecatedAlias('MaxPathDelays', 'MaxDelaySamples');
            obj.MaxDelaySamples = val;
        end

        function val = get.MaxPathDelays(obj)
            val = obj.MaxDelaySamples;
        end

        function set.MaxNormDoppler(obj, val)
            EpConfig.warnDeprecatedAlias('MaxNormDoppler', 'MaxDopplerIdx');
            obj.MaxDopplerIdx = val;
        end

        function val = get.MaxNormDoppler(obj)
            val = obj.MaxDopplerIdx;
        end

        function set.PilotSnr(obj, val)
            EpConfig.warnDeprecatedAlias('PilotSnr', 'PilotSnrDb');
            obj.PilotSnrDb = val;
        end

        function val = get.PilotSnr(obj)
            val = obj.PilotSnrDb;
        end

        function val = get.PilotPos1(obj)
            val = obj.PilotIndex;
        end

        function val = get.PilotPos0(obj)
            val = obj.PilotIndex - 1;
        end

        function val = get.DataPos1(obj)
            val = obj.ActiveIndices(:);
        end

        function val = get.DataPos0(obj)
            val = obj.ActiveIndices(:) - 1;
        end

        % updateDerivedParams: 根据基础参数更新所有派生量。
        function updateDerivedParams(obj)
            if isempty(obj.TotalSubcarriers) || isempty(obj.MaxDelaySamples) || ...
               isempty(obj.MaxDopplerIdx)    || isempty(obj.DopplerGuard)
                return;
            end

            % N_data = N_total - PrefixLength
            obj.PrefixLength = obj.MaxDelaySamples;
            numDataSubcarriers = obj.TotalSubcarriers - obj.PrefixLength;
            obj.NumDataSubcarriers = numDataSubcarriers;

            if numDataSubcarriers <= 0
                error('EpConfig:BadN', ...
                    'TotalSubcarriers (%d) must be larger than PrefixLength (%d).', ...
                    obj.TotalSubcarriers, obj.PrefixLength);
            end

            obj.ChirpParam1 = ...
                (2 * (obj.MaxDopplerIdx + obj.DopplerGuard) + 1) / (2 * numDataSubcarriers);
            obj.ChirpParam2 = 1 / (numDataSubcarriers ^ 2 * 2 * pi);

            if (2 * (obj.MaxDopplerIdx + obj.DopplerGuard) * (obj.MaxDelaySamples + 1)) ...
                    + obj.MaxDelaySamples > numDataSubcarriers
                error('EpConfig:Orth', 'Subcarriers do not satisfy orthogonality condition.');
            end

            zpTheory = ceil((obj.MaxDelaySamples + 1) * ...
                (2 * (obj.MaxDopplerIdx + obj.DopplerGuard) + 1) - 1);
            obj.TheoreticalZeroPaddingLength = zpTheory;

            if obj.ManualZeroPaddingLength > 0
                obj.ZeroPaddingLength = obj.ManualZeroPaddingLength;
            else
                obj.ZeroPaddingLength = zpTheory;
            end

            obj.PilotSequenceLength = 1;
            obj.PilotIndex = obj.ZeroPaddingLength + 1;
            obj.PilotSequence = 1;

            pilotSpan = 2 * obj.ZeroPaddingLength + obj.PilotSequenceLength;
            obj.NumActiveCarriers = numDataSubcarriers - pilotSpan;
            obj.ActiveIndices = (pilotSpan + 1):numDataSubcarriers;

            if obj.NumActiveCarriers <= 0
                error('EpConfig:NoData', ...
                    'Pilot/ZP overhead (%d) exceeds available subcarriers (%d).', ...
                    pilotSpan, numDataSubcarriers);
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
                warning('EpConfig:DeprecatedAlias', ...
                    'Property "%s" is deprecated. Use "%s" instead.', oldName, newName);
                warnedMap(warnKey) = true;
            end
        end

    end

end

