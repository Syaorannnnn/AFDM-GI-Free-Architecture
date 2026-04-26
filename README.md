# AFDM 系统性能研究

本仓库是基于 MATLAB 的 AFDM 系统仿真与论文实验工程，当前主线已完成第四章实验重构与两个 worktree 成果整合。代码覆盖 OFDM 基线、嵌入式导频（EP）基线、GI-Free 架构、分数 Doppler 信道建模，以及最终 CLIP 接收机的主干仿真流程。

> 状态：阶段性完成 · MATLAB R2023b+ · 需要 Communications Toolbox

## 项目范围

AFDM 使用 DAFT 的 chirp 参数在双色散信道中获得分集增益。传统 EP 架构通过保护带隔离导频与数据，但会降低频谱效率。GI-Free 架构保留 DAFT 域单导频并移除保护带，使其余 $N-1$ 个子载波承载数据；代价是数据到导频干扰（ID2P）进入路径估计链路。

本工程的重点是把 GI-Free 从整数 Doppler、固定经验门限的早期设定推进到可复现实验主线。主线实现包含 OMP+CFAR 路径检测、分数 Doppler Dirichlet 核建模、动态导频功率跟踪、路径稳定度门控、置信门控清洗，以及最终路径诊断输出。

## 目录结构

`Core/` 保存可复用系统模块。GI-Free 主体位于 `Core/SystemFunctions/GiFree/`，EP 与 OFDM 基线分别位于 `Core/SystemFunctions/EmbeddedPilot/` 和 `Core/SystemFunctions/OfdmBaseline/`。通用物理时变信道函数为 `Core/buildTimeVaryingChannel.m`。

`Sim/MainSimulation.m` 是第四章主仿真入口，只负责公共参数、日志目录和 Part 1 至 Part 4 的分派。主干实验脚本位于 `Sim/experiments/`，共享辅助函数位于 `Sim/helpers/`。

`docs/memory/` 保存理论、公式、实现状态和第四章整合记录。涉及 AFDM/GI-Free 技术讨论前，应先读取该目录下的索引和核心文档。

## 主线实验

| Part | 脚本 | 作用 |
|------|------|------|
| Part 1 | `Sim/experiments/Exp1_AfdmVsOfdm.m` | AFDM-EP 与 CP-OFDM 在整数 Doppler 双色散信道下的 BER 基线对比 |
| Part 2 | `Sim/experiments/Exp2_CfarThreshold.m` | 固定导频阶段一检测器归因，比较 Zhou 固定门限、Zhou CFAR、OMP 固定门限和 OMP CFAR |
| Part 3 | `Sim/experiments/Exp3_DynamicPilotGain.m` | 固定导频与动态导频对比，输出 BER、NMSE 与路径检测概率 `Pd` |
| Part 4 | `Sim/experiments/Exp4_FractionalDoppler.m` | EP 与 GI-Free ProGuard-CLIP 在分数 Doppler 场景下的对比 |

## 核心实现

GI-Free 配置由 `GiFreeConfig` 管理，包含 Theorem-1 校验、动态导频、CFAR 目标虚警率和 R2 匹配统计量开关。`GiFreeEstimator` 提供 OMP+CFAR、固定门限 OMP 和完整复合响应归一化匹配统计量。`GiFreeReceiver` 实现 CLIP 接收流程，并输出 `rxDiag.finalEstimatedPaths` 与 `rxDiag.finalPathCount`，供 Exp3 统计 `Pd`。

R2 匹配统计量通过 `UseMatchedPilotMetric` 控制，默认值为 `false`。主干 CLIP 配置在 `createClipPilotConfigPair.m` 与 `Exp4_FractionalDoppler.m` 中显式开启该开关。

## 运行与验证

运行完整第四章主线：

```matlab
addpath(genpath(pwd));
MainSimulation;
```

运行当前轻量回归测试：

```matlab
addpath(genpath(pwd));
results = runtests('Tests/testSimulationRefactor.m');
assertSuccess(results);
```

## 环境要求

- MATLAB R2023b 或更高版本
- Communications Toolbox

