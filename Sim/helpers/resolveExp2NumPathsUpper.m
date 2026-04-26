function numPathsUpper = resolveExp2NumPathsUpper(numPaths)
% RESOLVEEXP2NUMPATHSUPPER 返回 Exp2 的路径搜索上限。
%
%   描述:
%   为避免高路径数场景下搜索空间过紧，Exp2 对路径上限采用与真实路径数
%   相关的放宽规则。当前规则保持 P=3 时与历史配置一致，同时放宽 P=5。
%
%   输入:
%   numPaths - (double) 真实路径数。
%
%   输出:
%   numPathsUpper - (double) OMP 与 Zhou 检测共用的路径搜索上限。

    arguments
        numPaths (1,1) double {mustBeInteger, mustBePositive}
    end

    numPathsUpper = max(8, 2 * numPaths);
end
