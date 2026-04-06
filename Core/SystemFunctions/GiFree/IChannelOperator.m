classdef (Abstract) IChannelOperator < handle
% ICHANNELOPERATOR - GI-Free 有效信道算子接口
%
%   描述:
%   定义路径级/整体有效信道构建与 Doppler 基函数生成的统一接口，
%   便于估计器与信道构造模块解耦。
%
%   版本历史:
%   2026-04-01 - Aiden - 注释规范化。

    methods (Abstract)
        % BUILDPATHMATRIX 构造单一路径对应的有效信道矩阵。
        hPath = buildPathMatrix(obj, delay, doppler)
        % BUILDEFFECTIVECHANNEL 将多径参数叠加为总有效信道矩阵。
        effectiveChannel = buildEffectiveChannel(obj, pathParams)
        % BUILDPILOTRESPONSEVECTOR 构造某导频索引对应的路径响应列向量。
        hCol = buildPilotResponseVector(obj, delay, doppler, pilotIdx)
        % BUILDCOMPOSITEPILOTRESPONSE 构造导频复合响应（含导频幅值与序列）。
        hComposite = buildCompositePilotResponse(obj, delay, doppler)
        % APPLYPATHTOFRAME 将单一路径作用到给定发射帧。
        y = applyPathToFrame(obj, delay, doppler, txFrame)
        % BUILDPATHRESPONSEDICTIONARY 构造多径响应字典矩阵。
        dict = buildPathResponseDictionary(obj, pathParams, txFrame)
        % BUILDDOPPLERBASIS 构造整数 Doppler 邻域基函数集合。
        [basisVecs, kvRange] = buildDopplerBasis(obj, delay, intDoppler, txFrame)
        % COMPUTEDIRICHLETWITHDERIV 计算 Dirichlet 系数与一阶导数。
        [dVal, dDeriv] = computeDirichletWithDeriv(obj, fracPart, shiftIndex, numSc)
    end

end


