function applyIeeeStyle(axHandle, fontSize)
% APPLYIEEESTYLE 应用仓库统一绘图风格。
% 输入:
%   axHandle - 目标坐标轴句柄。
%   fontSize - 可选字号，默认 11。
% 输出:
%   无。

    if nargin < 2
        fontSize = 11;
    end

    plotStyle = PlotStyle();
    plotStyle.apply(axHandle, fontSize);
end
