function fieldValue = getFieldWithDefault(inputStruct, fieldName, defaultValue)
% GETFIELDWITHDEFAULT 从结构体读取字段，不存在时返回默认值。
% 输入:
%   inputStruct - 输入结构体。
%   fieldName - 待读取字段名。
%   defaultValue - 字段不存在时的默认值。
% 输出:
%   fieldValue - 读取结果或默认值。

    if isfield(inputStruct, fieldName)
        fieldValue = inputStruct.(fieldName);
    else
        fieldValue = defaultValue;
    end
end
