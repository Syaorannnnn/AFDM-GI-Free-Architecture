function optionValue = validateTextOption(optionValue, validOptions, optionName, useUpper)
% VALIDATETEXTOPTION 统一文本选项并校验取值。
% 输入:
%   optionValue - 待校验选项。
%   validOptions - 允许取值列表。
%   optionName - 选项名，用于报错信息。
%   useUpper - 是否转为大写后比较，默认 false。
% 输出:
%   optionValue - 统一大小写后的 string 标量。

    if nargin < 4
        useUpper = false;
    end

    if useUpper
        optionValue = upper(string(optionValue));
    else
        optionValue = lower(string(optionValue));
    end

    if ~isscalar(optionValue) || ~any(optionValue == validOptions)
        error('Simulation:InvalidOption', ...
            '%s must be one of: %s.', optionName, strjoin(cellstr(validOptions), ', '));
    end
end
