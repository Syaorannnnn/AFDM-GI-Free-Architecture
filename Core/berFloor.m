function berValue = berFloor(numErrors, totalBits)
% BERFLOOR 计算 BER，零错误时采用 Rule-of-3 下界。
% 输入:
%   numErrors - 误比特数。
%   totalBits - 总比特数。
% 输出:
%   berValue - BER 或 Rule-of-3 下界。

    if totalBits <= 0
        error('Simulation:InvalidBitCount', ...
            'totalBits must be positive.');
    end

    if numErrors > 0
        berValue = numErrors / totalBits;
    else
        berValue = 3 / totalBits;
    end
end
