classdef SparseLmmseSolver
    % SparseLmmseSolver  基于 PCG 的稀疏 LMMSE 均衡器 (通用数据索引版)
    %
    %   新增 DataIndices 选项: 指定数据列索引 (1-based).
    %   未指定时默认 2:numSc (向后兼容单导频模式).

    methods (Static)

        function [estSig, info] = solve(cleanSig, hEff, regParam, numSc, options)
            arguments
                cleanSig  (:,1)
                hEff
                regParam  (1,1) double
                numSc     (1,1) double
                options.MaxIter     (1,1) double = 30
                options.Tolerance   (1,1) double = 1e-6
                options.WarmStart   (:,1) double = []
                options.DataIndices (:,1) double = []
            end

            maxIter = options.MaxIter;
            tol     = options.Tolerance;

            % 数据列索引: 默认 2:numSc (单导频兼容)
            if isempty(options.DataIndices)
                dataIndices = (2:numSc)';
            else
                dataIndices = options.DataIndices;
            end
            hData   = hEff(:, dataIndices);
            numData = length(dataIndices);

            % 右端向量
            b     = hData' * cleanSig;
            bnorm = norm(b);

            if bnorm < 1e-15
                estSig = zeros(numSc, 1);
                info.iterations = 0; info.relResidual = 0;
                return;
            end

            % 对角预条件
            diagPrecond = zeros(numData, 1);
            for k = 1:numData
                [~, ~, vals] = find(hData(:, k));
                diagPrecond(k) = real(vals' * vals) + regParam;
            end
            diagPrecond = max(diagPrecond, 1e-15);

            % 初始化 (支持 warm start)
            if ~isempty(options.WarmStart) && length(options.WarmStart) == numSc
                x = options.WarmStart(dataIndices);
            elseif ~isempty(options.WarmStart) && length(options.WarmStart) == numData
                x = options.WarmStart;
            else
                x = zeros(numData, 1);
            end

            applyA = @(v) hData' * (hData * v) + regParam * v;

            % PCG 迭代
            r      = b - applyA(x);
            z      = r ./ diagPrecond;
            p      = z;
            rzOld  = real(r' * z);
            relRes = norm(r) / bnorm;

            finalIter = maxIter;
            for iter = 1:maxIter
                Ap  = applyA(p);
                pAp = real(p' * Ap);
                if abs(pAp) < 1e-15, finalIter = iter; break; end

                alpha = rzOld / pAp;
                x     = x + alpha * p;
                r     = r - alpha * Ap;
                relRes = norm(r) / bnorm;
                if relRes < tol, finalIter = iter; break; end

                z     = r ./ diagPrecond;
                rzNew = real(r' * z);
                p     = z + (rzNew / max(rzOld, 1e-15)) * p;
                rzOld = rzNew;
            end

            % 组装全长输出
            estSig = zeros(numSc, 1);
            estSig(dataIndices) = x;

            if nargout >= 2
                info.iterations  = finalIter;
                info.relResidual = relRes;
            end
        end

    end

end
