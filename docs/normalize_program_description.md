# T0 DSSD 归一化程序说明文档

## 1. 概述

`normalize.cpp` 是 brill2 数据分析框架中用于 T0 DSSD（双面硅条探测器）能量归一化（刻度）的核心程序。该程序的可执行文件名为 `./normalize`。

### 物理背景

DSSD 探测器有正面（front）和背面（back）两面的硅条，每个条（strip）的电子学增益可能不同。当同一个粒子穿过探测器时，其在正面和背面沉积的能量是相同的。归一化程序利用这一特性，以一面已归一化的条作为参考，对另一面的条进行能量刻度，使得所有条的能量响应一致。

### 归一化公式

```
E_normalized = p0 + p1 * E_raw
```

其中 `p0` 为偏移量（offset），`p1` 为比例系数（scale）。注意：代码中原本预留了二次项 `p2 * E_raw^2`，但当前被注释掉了，实际只使用线性归一化。

---

## 2. 命令行参数

| 参数 | 简写 | 类型 | 说明 |
|------|------|------|------|
| `--help` | `-h` | - | 打印帮助信息 |
| `--run` | `-r` | int | 起始 run 号（必填） |
| `--end-run` | `-e` | int | 结束 run 号（可选，默认等于 run） |
| `--trigger` | `-t` | string | 触发器类型 |
| `--config` | `-c` | string | 配置文件路径，默认 `config.toml` |
| `detector` | (位置参数) | vector\<string\> | 探测器名称列表，如 `t0d1 t0d2 t0d3 t0d4` |

### 使用示例

```bash
./normalize -r 100 -e 200 -t main -c config.toml t0d1 t0d2
```

---

## 3. 核心数据结构

### 3.1 DssdNormalizeParameters（归一化参数）

定义于 [dssd.h](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/t0/dssd.h#L13-L22)

```cpp
struct DssdNormalizeParameters {
    int front_strips = 0;          // 正面条数
    int back_strips = 0;           // 背面条数
    double front_p0[kMaxStrips];   // 正面各条 offset（默认 0.0）
    double front_p1[kMaxStrips];   // 正面各条 scale（默认 1.0）
    double front_p2[kMaxStrips];   // 正面各条二次项（默认 0.0，当前未使用）
    double back_p0[kMaxStrips];    // 背面各条 offset
    double back_p1[kMaxStrips];    // 背面各条 scale
    double back_p2[kMaxStrips];    // 背面各条二次项（当前未使用）
};
```

其中 `kMaxStrips = 128`。

### 3.2 NromalizeStripsConfig（归一化条配置）

定义于 [config.h](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/config.h#L28-L34)

```cpp
struct NromalizeStripsConfig {
    int norm_side = 0;             // 0=归一化正面, 1=归一化背面
    int ref[2];                    // 参考面条的范围 [起始条, 结束条]
    int norm[2];                   // 待归一化条的范围 [起始条, 结束条]
    double ref_energy[2];          // 参考面能量范围 [min, max]
    double norm_energy[2];         // 待归一化面能量范围 [min, max]
};
```

### 3.3 DssdEvent（DSSD 原始事件）

定义于 [dssd_event.h](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/event/ingot/dssd_event.h#L11-L19)

```cpp
struct DssdEvent {
    int front_num;                     // 正面 hits 数量
    int front_strip[kDssdMaxHits];     // 正面 hit 条号
    double front_energy[kDssdMaxHits]; // 正面 hit 能量
    double front_time[kDssdMaxHits];   // 正面 hit 时间
    int back_num;                      // 背面 hits 数量
    int back_strip[kDssdMaxHits];      // 背面 hit 条号
    double back_energy[kDssdMaxHits];  // 背面 hit 能量
    double back_time[kDssdMaxHits];    // 背面 hit 时间
};
```

其中 `kDssdMaxHits = 32`。

---

## 4. 核心算法流程

### 4.1 main() 函数流程

```
main()
  ├── 解析命令行参数
  ├── 验证探测器名称（只允许 t0d1, t0d2, t0d3, t0d4）
  ├── 加载配置文件 config.toml
  ├── 对每个探测器循环：
  │     ├── 构建 TChain，链入所有 run 的 ingot ROOT 文件
  │     ├── 初始化 DssdNormalizeParameters（p0=0, p1=1, p2=0）
  │     ├── 创建输出 ROOT 文件
  │     ├── 对每个 NromalizeStripsConfig 调用 NormalizeStrips()
  │     ├── 将归一化参数写入 front/back txt 文件
  │     └── 关闭 ROOT 文件
  └── 返回 0
```

### 4.2 NormalizeStrips() 函数详细流程

这是核心归一化函数，位于 [normalize.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/normalize.cpp#L41-L210)。

#### 第一阶段：填充 TGraph

```
对每条事件循环：
  ├── 跳过非单 hit 事件（front_num != 1 或 back_num != 1）
  ├── 特殊处理：t0d3 探测器中，跳过 back strip 18/19 且 front energy > 12000 的事件
  ├── 如果 norm_side == 0（归一化正面）：
  │     ├── 跳过：背面条不在参考范围 [ref[0], ref[1]] 内的事件
  │     ├── 跳过：正面条不在待归一化范围 [norm[0], norm[1]] 内的事件
  │     ├── 跳过：背/正面能量不在能量范围内的事件
  │     └── 填充 TGraph：ge[front_strip].AddPoint(front_raw_energy, back_normalized_energy)
  │         其中 back_normalized_energy = NormEnergy(parameters, side=1, back_strip, back_raw_energy)
  └── 如果 norm_side == 1（归一化背面）：
        ├── 跳过：正面条不在参考范围 [ref[0], ref[1]] 内的事件
        ├── 跳过：背面条不在待归一化范围 [norm[0], norm[1]] 内的事件
        ├── 跳过：正/背面能量不在能量范围内的事件
        └── 填充 TGraph：ge[back_strip].AddPoint(back_raw_energy, front_normalized_energy)
            其中 front_normalized_energy = NormEnergy(parameters, side=0, front_strip, front_raw_energy)
```

**关键理解**：
- 当 `norm_side=0`（归一化正面）时，假设**背面已经归一化好**，用背面归一化后的能量作为 Y 轴，正面原始能量作为 X 轴，拟合得到正面各条的归一化参数。
- 当 `norm_side=1`（归一化背面）时，反过来，用**正面归一化后的能量**作为参考。
- 因此，配置文件中必须先配置一面作为参考（已归一化），再配置另一面。通常先归一化一面（此时 p0=0, p1=1 即不做任何变换），然后再归一化另一面。

#### 第二阶段：拟合

```
对每个待归一化的条 i 循环：
  ├── 如果 TGraph 点数 > 10：
  │     ├── 创建 TF1("efit", "pol1", 0, 60000)
  │     ├── 初始值：p0=0.0, p1=1.0
  │     ├── 执行拟合：ge[i].Fit("efit", "QR+ ROB=0.8")
  │     │   - "Q" = 安静模式
  │     │   - "R" = 使用函数范围
  │     │   - "ROB=0.8" = 使用 Robust 拟合，排除 80% 外的异常点
  │     ├── 存储拟合结果到 parameters.front_p0[i]/front_p1[i] 或 back_p0[i]/back_p1[i]
  │     └── 打印拟合参数
  └── 如果点数 <= 10：打印警告，跳过该条
```

#### 第三阶段：残差分析

```
如果 norm 范围覆盖所有条（abs(norm[0]-norm[1]) == front_strips-1）：
  ├── 对每个条计算残差：res = NormEnergy(parameters, side, strip, raw_energy) - reference_energy
  ├── 填充残差直方图（每个条一个 + 总残差直方图）
  ├── 绘制 chi2/ndf vs strip 的 TGraph
  └── 写入 ROOT 文件
```

---

## 5. 输入输出文件

### 输入文件

| 文件 | 路径模式 | 说明 |
|------|----------|------|
| 原始数据 | `{workspace}/{paths.ingot}/{detector}_{trigger}{run:04d}.root` | DSSD 原始 ingot 数据 |
| 配置文件 | `config.toml` | 包含归一化配置 |

### 输出文件

**ROOT 文件**：`{workspace}/{paths.normalize}/{detector}_{trigger}{run:04d}_{end_run:04d}.root`

包含以下对象：
| 对象名 | 类型 | 说明 |
|--------|------|------|
| `g_f{strip}` | TGraph | 正面条 strip 的拟合散点图 |
| `g_b{strip}` | TGraph | 背面条 strip 的拟合散点图 |
| `res_f{strip}` | TH1D | 正面条 strip 的残差分布 |
| `res_b{strip}` | TH1D | 背面条 strip 的残差分布 |
| `h_total_res` | TH1D | 所有条的总残差分布 |
| `g_f_chi2ndf` | TGraph | 正面各条 chi2/ndf |
| `g_b_chi2ndf` | TGraph | 背面各条 chi2/ndf |

**TXT 参数文件**：
- `{workspace}/{paths.normalize}/{detector}_front_{trigger}{run:04d}.txt`
- `{workspace}/{paths.normalize}/{detector}_back_{trigger}{run:04d}.txt`

> **重要**：TXT 参数文件名中的 run 号**仅使用起始 run**（即 `-r` 参数），**不包含结束 run**。例如执行 `./normalize -r 60 -e 70 -t t1 t0d1`，生成的 TXT 文件名为 `t0d1_front_t1_0060.txt`，而非 `t0d1_front_t1_0060_0070.txt`。这与 ROOT 诊断文件的命名不同（ROOT 文件包含起始和结束 run）。

格式：
```
strip p0 p1 p2
0 0.5 1.02 0.0
1 -0.3 0.98 0.0
...
```

---

## 6. 配置文件示例（config.toml）

### 6.1 单段归一化示例

```toml
[normalize]

# 归一化参数文件与 run 号的对应关系
[[normalize.run]]
run = 0
use = 100   # run >= 0 时使用 run 100 的归一化参数

# t0d1 探测器归一化配置
[normalize.t0d1]
[[normalize.t0d1.strips]]
index = 0
norm_side = 0       # 归一化正面
ref = [0, 31]       # 参考面：背面条 0-31
norm = [0, 31]      # 归一化：正面条 0-31
ref_energy = [5000, 40000]   # 参考面能量范围
norm_energy = [5000, 40000]  # 归一化面能量范围

[[normalize.t0d1.strips]]
index = 1
norm_side = 1       # 归一化背面
ref = [0, 31]       # 参考面：正面条 0-31（已用上面结果归一化）
norm = [0, 31]      # 归一化：背面条 0-31
ref_energy = [5000, 40000]
norm_energy = [5000, 40000]
```

### 6.2 多段归一化示例（run 60-70 一段，run 70-80 一段）

假设探测器状态在 run 70 前后发生了变化，需要对两段 run 分别做归一化。执行命令：

```bash
# 第一段：run 60-70 归一化，TXT 文件名为 t0d1_front_t1_0060.txt
./normalize -r 60 -e 70 -t t1 t0d1

# 第二段：run 70-80 归一化，TXT 文件名为 t0d1_front_t1_0070.txt
./normalize -r 70 -e 80 -t t1 t0d1
```

生成的 TXT 文件：
```
normalize/t0d1_front_t1_0060.txt   ← 起始 run=60
normalize/t0d1_back_t1_0060.txt
normalize/t0d1_front_t1_0070.txt   ← 起始 run=70
normalize/t0d1_back_t1_0070.txt
```

对应的 `config.toml` 配置：

```toml
[normalize]

# run 0 ~ 69 使用 run 60 的归一化参数
[[normalize.run]]
run = 0
use = 60

# run 70 及以上使用 run 70 的归一化参数
[[normalize.run]]
run = 70
use = 70

# t0d1 探测器归一化配置（同上，此处省略 strips 配置）
# [normalize.t0d1]
# ...
```

**`normalize.runs` 映射逻辑**：
| 数据 run | 满足 `run >= start` 的最大 start | 使用的归一化参数 run |
|----------|----------------------------------|----------------------|
| 50 | 0 | 60 |
| 65 | 0 | 60 |
| 70 | 70 | 70 |
| 75 | 70 | 70 |
| 85 | 70 | 70 |

**注意**：由于代码中 `has_normalized` 的检查被注释掉了，当前配置中 strips 的顺序是：
1. 先 norm_side=0，此时背面参数是默认值（p0=0, p1=1），相当于用背面原始能量作为参考来归一化正面
2. 再 norm_side=1，此时正面参数已经被更新，用正面归一化后的能量作为参考来归一化背面

---

## 7. 已知问题与注意事项

### 7.1 has_normalized 检查被注释

代码中多处 `has_normalized` 的检查被注释掉了（[L90](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/normalize.cpp#L90), [L103](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/normalize.cpp#L103), [L126](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/normalize.cpp#L126)）。这意味着：
- 如果有多个 `NromalizeStripsConfig`，后一个配置会覆盖前一个配置的归一化参数
- 不会跳过已归一化的条

### 7.2 t0d3 特殊处理

```cpp
if (detector_name == "t0d3") {
    if (((bs == 18) || (bs == 19)) && (fe > 12000)) continue;
}
```

这是硬编码的探测器特殊处理，用于排除 t0d3 探测器背面 18、19 条上的异常高能事件。

### 7.3 二次项被注释

`p2` 参数（二次项系数）在拟合和计算中都被注释掉了，仅使用线性归一化。但参数文件仍会输出 `p2` 列（值为 0.0）。

### 7.4 硬编码数组大小

`TGraph ge[128]` 和 `TH1D res[128]` 硬编码为 128，与 `kMaxStrips` 一致，但如果探测器条数超过 128 会有问题。

### 7.5 拟合选项

使用 `ROB=0.8`（Robust 拟合），自动排除残差最大的 20% 数据点，这在有噪声条时能提高拟合稳定性。

### 7.6 残差分析仅在全条范围时触发

残差直方图只在 `abs(config.norm[0]-config.norm[1]) == parameters.front_strips-1` 时才生成，即只有当 norm 范围覆盖全部条时才会做残差分析。

---

## 8. 程序依赖关系

```
normalize.cpp
├── include/t0/dssd.h          → DssdNormalizeParameters, WriteDssdNormalizeParameters
├── include/config.h           → AppConfig, NromalizeStripsConfig, NormalizeConfig
├── include/event/ingot/dssd_event.h  → DssdEvent, SetupInput
├── include/utils.h            → JoinPath, TriggerInfix, IsJumpRun, LoadConfig
├── external/cxxopts.hpp       → 命令行参数解析
└── ROOT 库 (TFile, TTree, TChain, TGraph, TF1, TH1D, TString)
```

下游消费者：
- `estimate/estimate_normalize.cpp` → 读取归一化参数文件，应用归一化
- `match_dssd.cpp` → 读取归一化参数文件，进行 DSSD 前后条匹配（详见第 9 节）
- `src/t0/dssd.cpp` → `ApplyDssdNormalize()` 函数，在分析流程中应用归一化

---

## 9. match_dssd 中归一化参数的调用流程

`match_dssd.cpp` 是归一化参数的主要下游消费者，位于 [match_dssd.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/match_dssd.cpp)。该程序的目的是将归一化后的 DSSD 正反面 hits 进行匹配（即找出同一个粒子在正反两面产生的 hit 对），输出匹配后的事件。

### 9.1 归一化参数文件的查找逻辑

```cpp
int normalize_file_run = 0;
for (const auto &[start, use] : config.normalize.runs) {
    if (run >= start) normalize_file_run = use;
}
```

`match_dssd` **不直接使用当前 run 号的归一化文件**，而是通过 `config.normalize.runs` 配置表来查找应该使用哪个 run 的归一化参数。该表定义在 `config.toml` 中：

```toml
[[normalize.run]]
run = 0      # 从 run 0 开始
use = 100    # 使用 run 100 生成的归一化参数

[[normalize.run]]
run = 200    # 从 run 200 开始
use = 250    # 使用 run 250 生成的归一化参数
```

**查找规则**：遍历 `normalize.runs` 列表，对于当前数据 run 号，找到满足 `data_run >= start` 的最大 `start` 对应的 `use` 值。例如：
- 数据 run=50 → 使用 run 100 的归一化参数
- 数据 run=150 → 使用 run 100 的归一化参数  
- 数据 run=200 → 使用 run 250 的归一化参数
- 数据 run=300 → 使用 run 250 的归一化参数

**注意**：当前 `match_dssd.cpp` 中硬编码了 trigger 为 `t1`：
```cpp
TString front_path = TString::Format(
    "%s/%s_front_t1_%04d.txt",  // 硬编码 "t1"
    ...
);
```
这意味着 `match_dssd` 始终读取 `t1` trigger 的归一化参数，无论命令行传入的 `--trigger` 是什么。

### 9.2 归一化参数文件路径

读取的 TXT 文件路径为：
```
{workspace}/{paths.normalize}/{detector}_front_t1_{normalize_file_run:04d}.txt
{workspace}/{paths.normalize}/{detector}_back_t1_{normalize_file_run:04d}.txt
```

通过 `ReadDssdNormalizeParameters()` 函数（定义于 [dssd.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/t0/dssd.cpp#L134-L158)）读取。

### 9.3 match_dssd 完整处理流程

```
match_dssd main()
  ├── 解析命令行参数
  ├── 根据 config.normalize.runs 查找归一化参数文件 run 号
  ├── 对每个探测器循环：
  │     ├── 读取归一化参数 TXT 文件 → DssdNormalizeParameters
  │     ├── 打开原始 ingot 数据 ROOT 文件
  │     ├── 创建输出 ROOT 文件（match 目录）
  │     ├── 对每条事件循环：
  │     │     ├── ipt->GetEntry(entry)          // 读取原始事件
  │     │     ├── ApplyDssdNormalize()           // ★ 应用归一化参数
  │     │     ├── MatchDssdEvent()               // 正反面 hit 匹配
  │     │     └── opt.Fill()                     // 写入输出树
  │     ├── 写入 h_energy_diff 直方图
  │     └── 关闭文件
  └── 返回 0
```

### 9.4 ApplyDssdNormalize 的调用

在事件循环中，**每条事件都会调用** `ApplyDssdNormalize()`：

```cpp
for (long long entry = 0; entry < total; ++entry) {
    ipt->GetEntry(entry);
    brill::ApplyDssdNormalize(raw_event, parameters, normalized_event);
    brill::MatchDssdEvent(normalized_event, working_detector, match_event, &h_energy_diff);
    opt.Fill();
}
```

`ApplyDssdNormalize()` 函数（定义于 [dssd.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/t0/dssd.cpp#L160-L191)）对每个 hit 执行：

```cpp
output.front_energy[i] = NormalizeEnergy(
    input.front_energy[i],
    parameters.front_p0[strip],  // offset
    parameters.front_p1[strip],  // scale
    parameters.front_p2[strip]   // 二次项（当前为 0）
);
// 即: E_out = p0 + p1 * E_in + p2 * E_in²
```

### 9.5 MatchDssdEvent 的调用

归一化之后，`MatchDssdEvent()` 对正反面的 hits 进行匹配。匹配逻辑基于能量差：
- 对于归一化后的正面 hit 和背面 hit，计算 `|front_energy - back_energy|`
- 如果能量差小于 `match_tolerance`（可在 `config.toml` 的探测器配置中设置，或通过 `--window` 命令行参数覆盖），则认为它们来自同一个粒子
- 匹配成功的 hit 对写入 `DssdMatchEvent`

### 9.6 输出文件

| 输出 | 路径 | 说明 |
|------|------|------|
| 匹配事件 | `{workspace}/{paths.match}/{detector}_{trigger}{run:04d}.root` | 包含 `tree`（匹配后的事件）和 `h_energy_diff`（正反面能量差分布） |

### 9.7 典型工作流程

```
步骤1: ./normalize -r 100 -e 100 -t t1 t0d1
       → 生成 t0d1_front_t1_0100.txt 和 t0d1_back_t1_0100.txt

步骤2: 在 config.toml 中配置：
       [[normalize.run]]
       run = 0
       use = 100

步骤3: ./match_dssd -r 100 -t main t0d1
       → 读取 run 100 的归一化参数（根据 config 映射）
       → 读取 run 100 的原始数据
       → 应用归一化 → 匹配 → 输出
```

---

## 10. 归一化参数的使用方式

归一化参数被写入 TXT 文件后，由 `ReadDssdNormalizeParameters()` 读取，然后通过 `ApplyDssdNormalize()` 应用到原始数据：

```cpp
// 应用归一化（src/t0/dssd.cpp）
output.front_energy[i] = NormalizeEnergy(
    input.front_energy[i],
    parameters.front_p0[strip],  // offset
    parameters.front_p1[strip],  // scale
    parameters.front_p2[strip]   // 二次项（当前为 0）
);
// 即: E_out = p0 + p1 * E_in + p2 * E_in^2
```

在分析流程中，归一化后的数据用于后续的粒子鉴别（PID）、径迹重建等步骤。

---

*文档生成日期：2026-09-08*