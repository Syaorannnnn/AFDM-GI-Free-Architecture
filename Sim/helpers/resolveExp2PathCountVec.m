function pathCountVec = resolveExp2PathCountVec(pathCountVecOverride)
% RESOLVEEXP2PATHCOUNTVEC 返回 Exp2 需要扫描的路径数集合。
%
%   描述:
%   默认运行 P=3 与 P=5 两组路径数。若提供 override，则只运行指定路径数，
%   用于定向复现实验。
%
%   输入:
%   pathCountVecOverride - (double) 可选，路径数标量或向量。
%
%   输出:
%   pathCountVec - (1xN double) 去重后的路径数向量。

    if nargin < 1 || isempty(pathCountVecOverride)
        pathCountVec = [3, 5];
        return;
    end

    validateattributes(pathCountVecOverride, {'numeric'}, ...
        {'vector', 'integer', 'positive'}, mfilename, 'pathCountVecOverride');
    pathCountVec = unique(pathCountVecOverride(:).', 'stable');
end
