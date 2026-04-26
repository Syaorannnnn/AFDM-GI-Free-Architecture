function [rxPilotSignal, truePaths, cfg] = synthesizeThreePathSignal(cfg)
% SYNTHESIZETHREEPATHSIGNAL 构造已知 3 条路径的合成导频信号，供 OMP 一致性测试使用。
%
%   描述:
%   基于给定配置生成确定性 3 路径信道下的导频接收信号，路径参数固定（不随机），
%   使测试可重现并能验证真实路径的候选排名。
%   信号构造步骤：导频帧 → FractionalChannelBuilder 施加信道 → 无噪声（高 SNR 等效）。
%
%   语法:
%   [rxPilotSignal, truePaths, cfg] = synthesizeThreePathSignal(cfg)
%   [rxPilotSignal, truePaths, cfg] = synthesizeThreePathSignal()
%
%   输入:
%   cfg - (可选, GiFreeConfig) 测试用配置对象；缺省时内部构造默认测试配置。
%
%   输出:
%   rxPilotSignal - (Nx1 complex) 导频域接收信号（无噪声）。
%   truePaths     - (1×3 struct，含 delay/dopplerIdx) 三条真实路径描述符。
%   cfg           - (GiFreeConfig) 实际使用的配置对象（入参原样返回或新建）。
%
%   路径参数（固定）:
%   路径 1: delay=0, dopplerIdx=1,  增益=0.8+0j
%   路径 2: delay=2, dopplerIdx=-1, 增益=0.5+0.3j
%   路径 3: delay=1, dopplerIdx=2,  增益=0.3-0.2j
%
%   版本历史:
%   2026-04-24 - 为 OMP 一致性测试新增。

    if nargin < 1 || isempty(cfg)
        cfg = buildDefaultTestGiFreeConfig();
    end

    % ---- 固定路径参数（确定性，可重现）----
    truePaths = struct( ...
        'delay',      {0,  2,  1}, ...
        'dopplerIdx', {1, -1,  2});

    pathGains = [0.8 + 0j; 0.5 + 0.3j; 0.3 - 0.2j];

    % ---- 构造信道构造算子 ----
    chBuilder = FractionalChannelBuilder(cfg);

    % ---- 构造导频帧（仅导频位置非零）----
    numSc      = cfg.NumSubcarriers;
    pilotFrame = zeros(numSc, 1);
    pilotAmp   = cfg.PerPilotAmplitude;
    pilotSeq   = cfg.PilotSequence;
    pilotPos1  = cfg.PilotPos1;            % 1-based 索引
    pilotFrame(pilotPos1) = pilotAmp * pilotSeq;

    % ---- 叠加三条路径响应（无噪声）----
    rxPilotSignal = zeros(numSc, 1);
    for pathIdx = 1:numel(truePaths)
        dv = truePaths(pathIdx).delay;
        av = truePaths(pathIdx).dopplerIdx;
        pathMat = chBuilder.buildPathMatrix(dv, av);
        rxPilotSignal = rxPilotSignal + pathGains(pathIdx) * (pathMat * pilotFrame);
    end
end

% ============================================================
%  内部 Helper: 构造默认测试配置
% ============================================================

function cfg = buildDefaultTestGiFreeConfig()
% BUILDDEFAULTTESTGIFREECONFIG 返回满足 Theorem-1 的最小测试用 GiFreeConfig。
%
%   参数选择原则：
%   - NumSubcarriers=128 以降低运行时间
%   - MaxDelaySamples=2, MaxDopplerIdx=2 保证 truePaths 中最大 delay/doppler 均合法
%   - DopplerGuard=1 以使 Theorem-1 验证通过
%   - DirichletRadius=4 与主实验一致

    cfg = GiFreeConfig();
    cfg.NumSubcarriers   = 128;
    cfg.ModulationOrder  = 4;
    cfg.MaxDelaySamples  = 2;
    cfg.MaxDopplerIdx    = 2;
    cfg.DopplerGuard     = 1;
    cfg.NumPaths         = 3;
    cfg.NumPathsUpper    = 6;
    cfg.DirichletRadius  = 4;
    cfg.PilotSnrDb       = 35;
    cfg.MaxSicIterations = 3;
    cfg.UseFractionalDoppler = false;
    cfg.validate();
end
