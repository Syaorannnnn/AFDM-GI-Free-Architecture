classdef SparseLmmseSolver
    % SparseLmmseSolver  基于预条件共轭梯度 (PCG) 的稀疏 LMMSE 均衡器
    %
    %   求解: (H_d^H H_d + lambda I) x = H_d^H y
    %   绝不显式构造 Gram 矩阵, 仅通过两次稀疏 matvec 实现 A*x.

    methods (Static)

        function [estSig, info] = solve(cleanSig, hEff, regParam, numSc, options)
            arguments
                cleanSig  (:,1)
                hEff
                regParam  (1,1) double
                numSc     (1,1) double
                options.MaxIter   (1,1) double = 30
                options.Tolerance (1,1) double = 1e-6
                options.WarmStart (:,1) double = []
            end

            maxIter = options.MaxIter;
            tol     = options.Tolerance;
            hData   = hEff(:, 2:numSc);
            numData = numSc - 1;

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

            % 初始化
            if ~isempty(options.WarmStart) && length(options.WarmStart) == numSc
                x = options.WarmStart(2:numSc);
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

            estSig          = zeros(numSc, 1);
            estSig(2:numSc) = x;

            if nargout >= 2
                info.iterations  = finalIter;
                info.relResidual = relRes;
            end
        end

    end

end
