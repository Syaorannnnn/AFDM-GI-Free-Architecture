classdef SparseLmmseSolver
% SPARSELMMSESOLVER - 稀疏 LMMSE 的 PCG 求解器
%
%   描述:
%   求解 (H_D^H H_D + lambda I)x = H_D^H y 的数据子载波估计问题，
%   通过对角预条件 PCG 降低大规模稀疏系统求解开销。
%
%   语法:
%   [estSignal, info] = SparseLmmseSolver.solve(cleanSignal, effectiveChannel, ...
%       regularization, numSubcarriers, 'MaxIter', 30, 'Tolerance', 1e-6);
%
%   输出:
%   estSignal - (Nx1 complex) 仅数据位置填充值，其余位置为零。
%   info      - (struct) 迭代次数与相对残差。
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。

    methods (Static)

        % SOLVE 使用预条件共轭梯度 (PCG) 求解稀疏 LMMSE 系统。
        function [estSignal, info] = solve(cleanSignal, effectiveChannel, regularization, numSubcarriers, options)
            arguments
                cleanSignal  (:,1)
                effectiveChannel
                regularization  (1,1) double
                numSubcarriers  (1,1) double
                options.MaxIter     (1,1) double = 30
                options.Tolerance   (1,1) double = 1e-6
                options.WarmStart   (:,1) double = []
                options.DataPos1    (:,1) double = []
                options.DataIndices (:,1) double = []
            end

            maxIter = options.MaxIter;
            residualTolerance = options.Tolerance;

            if ~isempty(options.DataPos1)
                dataPos1 = options.DataPos1;
            elseif ~isempty(options.DataIndices)
                SparseLmmseSolver.warnDeprecatedAlias('DataIndices', 'DataPos1');
                dataPos1 = options.DataIndices;
            else
                dataPos1 = (2:numSubcarriers)';
            end

            effectiveChannelData = effectiveChannel(:, dataPos1);
            numData = length(dataPos1);

            rhsVector = effectiveChannelData' * cleanSignal;
            rhsNorm = norm(rhsVector);

            if rhsNorm < 1e-15
                estSignal = zeros(numSubcarriers, 1);
                info.iterations = 0;
                info.relResidual = 0;
                return;
            end

            diagonalPreconditioner = zeros(numData, 1);
            for k = 1:numData
                [~, ~, nonzeroValues] = find(effectiveChannelData(:, k));
                diagonalPreconditioner(k) = real(nonzeroValues' * nonzeroValues) + regularization;
            end
            diagonalPreconditioner = max(diagonalPreconditioner, 1e-15);

            if ~isempty(options.WarmStart) && length(options.WarmStart) == numSubcarriers
                solutionData = options.WarmStart(dataPos1);
            elseif ~isempty(options.WarmStart) && length(options.WarmStart) == numData
                solutionData = options.WarmStart;
            else
                solutionData = zeros(numData, 1);
            end

            applyNormalMatrix = @(vector) effectiveChannelData' * (effectiveChannelData * vector) + regularization * vector;

            residualVector = rhsVector - applyNormalMatrix(solutionData);
            preconditionedResidual = residualVector ./ diagonalPreconditioner;
            searchDirection = preconditionedResidual;
            residualInnerProductOld = real(residualVector' * preconditionedResidual);
            relativeResidual = norm(residualVector) / rhsNorm;

            finalIter = maxIter;
            % NOTE: 主循环跟踪相对残差，满足阈值后提前停止以降低计算开销。
            for iter = 1:maxIter
                normalTimesDirection = applyNormalMatrix(searchDirection);
                directionInnerProduct = real(searchDirection' * normalTimesDirection);
                if abs(directionInnerProduct) < 1e-15
                    finalIter = iter;
                    break;
                end

                stepSize = residualInnerProductOld / directionInnerProduct;
                solutionData = solutionData + stepSize * searchDirection;
                residualVector = residualVector - stepSize * normalTimesDirection;
                relativeResidual = norm(residualVector) / rhsNorm;
                if relativeResidual < residualTolerance
                    finalIter = iter;
                    break;
                end

                preconditionedResidual = residualVector ./ diagonalPreconditioner;
                residualInnerProductNew = real(residualVector' * preconditionedResidual);
                betaCoeff = residualInnerProductNew / max(residualInnerProductOld, 1e-15);
                searchDirection = preconditionedResidual + betaCoeff * searchDirection;
                residualInnerProductOld = residualInnerProductNew;
            end

            estSignal = zeros(numSubcarriers, 1);
            estSignal(dataPos1) = solutionData;

            if nargout >= 2
                info.iterations = finalIter;
                info.relResidual = relativeResidual;
            end
        end

    end

    methods (Static, Access = private)

        % WARNDEPRECATEDALIAS 对旧字段名访问给出一次性弃用告警。
        function warnDeprecatedAlias(oldName, newName)
            persistent warnedMap;
            if isempty(warnedMap)
                warnedMap = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            end

            warnKey = [oldName '->' newName];
            if ~isKey(warnedMap, warnKey)
                warning('SparseLmmseSolver:DeprecatedAlias', ...
                    'Property "%s" is deprecated. Use "%s" instead.', oldName, newName);
                warnedMap(warnKey) = true;
            end
        end

    end

end


