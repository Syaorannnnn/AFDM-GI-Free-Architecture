function structArray = structAppend(structArray, structItem)
% STRUCTAPPEND 追加结构体元素，并检查字段集合一致性。
% 输入:
%   structArray - 原结构体数组。
%   structItem - 待追加结构体。
% 输出:
%   structArray - 追加后的结构体数组。

    if isempty(structArray)
        structArray = structItem;
        return;
    end

    if ~isequal(sort(fieldnames(structArray)), sort(fieldnames(structItem)))
        error('Simulation:StructFieldMismatch', ...
            'Struct append failed because field sets are inconsistent.');
    end

    structArray(end + 1) = structItem;
end
