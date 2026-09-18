# T0 能量刻度程序说明文档

## 1. 概述

brill2 数据分析框架提供三个 T0 探测器能量刻度程序：

| 程序 | 可执行文件 | 状态 |
|------|-----------|------|
| **calibrate_t0_v1** | `calibrate_t0_v1` | **标准程序** |
| calibrate_t0_v2 | `calibrate_t0_v2` | 实验性程序 |
| calibrate_t0 | `calibrate_t0` | 原始程序，已弃用 |

### 物理目标

T0 探测器系统由 5 层硅探测器（t0d1/t0d2/t0d3/t0d4/t0s）组成，每层探测器输出的原始信号为 ADC 道数。能量刻度的目标是为每层探测器确定一组线性刻度系数 `(p0, p1)`，使得：

```
E_calibrated[MeV] = p0 + p1 × E_raw[ADC]
```

### 刻度策略总览

所有版本均采用**基于 TCutG + 理论曲线拟合**方法：

1. 从 `pre_calibration` 程序生成的 TGraph 中读取各层对的原始 ADC 散点（`g_d1d2`、`g_d2d3`、`g_d3d4`、`g_d4s`）
2. 通过 `src/brill/Cut/` 目录下的 TCutG 截断文件筛选出目标粒子
3. 利用 catima 库计算粒子在多层硅探测器中的理论能损，建立理论 ΔE-E 曲线
4. 对实验数据点进行拟合以确定刻度系数

三个版本的区别在于**拟合策略**：原始版本一次性全局拟合 10 个参数；v1 将 d1 分离为两阶段（标准策略）；v2 采用四段式从深层到浅层逐步固定参数（实验性策略）。

---

## 2. calibrate_t0_v1（标准程序）

### 2.1 刻度策略：两阶段拟合

由于 d1 是 68 μm 薄探测器，其厚度不均匀性导致 d1-d2 全局拟合效果差。v1 将 d1-d2 分离出来独立处理：

```
Stage 1: d2-d3, d3-d4, d4-s1  联合拟合（8 参数: d2,d3,d4,s1）
                                      ↓ 固定 d2
Stage 2: d1-d2                  独立拟合（2 参数: d1，d2 固定）
```

**设计动机**：
- d2/d3/d4/s1 厚度（~1000 μm, ~1500 μm）相对均匀，联合拟合可互相约束获得稳定参数
- d1 厚度（68 μm）均匀性差，在 Stage 2 中用已固定的 d2 参数单独拟合，避免 d1 的不确定性污染其他层

### 2.2 命令行参数

| 参数 | 简写 | 类型 | 说明 |
|------|------|------|------|
| `--help` | `-h` | - | 打印帮助信息 |
| `--run` | `-r` | int | pre_calibration 的起始 run 号（必填） |
| `--trigger` | `-t` | string | 触发器类型 |
| `--config` | `-c` | string | 配置文件路径，默认 `config.toml` |

使用示例：
```bash
./calibrate_t0_v1 -r 57 -t t1
```

### 2.3 拟合参数初始值与限制

#### Stage 1（d2-d3 + d3-d4 + d4-s1 联合拟合，8 参数）

| 参数索引 | 对应探测器 | 初始值 | p0 限制 | p1 限制 |
|---------|-----------|--------|---------|---------|
| par[0], par[1] | t0d2 | (0.0, 0.006) | [-1.0, 1.0] | [0.0, 1.0] |
| par[2], par[3] | t0d3 | (0.0, 0.006) | [-1.0, 1.0] | [0.0, 1.0] |
| par[4], par[5] | t0d4 | (0.0, 0.003) | [-1.0, 1.0] | [0.0, 1.0] |
| par[6], par[7] | t0s | (0.0, 0.003) | [-1.0, 1.0] | [0.0, 1.0] |

TF1 范围：`[25000, 220000]`（跳过 d1-d2 的 0–25000 区间。layer=1,2,3 的 offset 起点分别是 26000, 102000, 180000）

#### Stage 2（d1-d2，d2 固定，2 参数）

| 参数索引 | 对应探测器 | 初始值 | p0 限制 | p1 限制 |
|---------|-----------|--------|---------|---------|
| par[0], par[1] | t0d1 | (0.0, 0.002) | [-1.0, 1.0] | [0.0, 1.0] |

TF1 范围：`[0, 25000]`（仅 d1-d2 对应区间）。固定参数 d2_p0, d2_p1 取自 Stage 1 结果。

### 2.4 拟合函数类

| 类名 | 用途 | 参数结构 |
|------|------|---------|
| `PidFitFuncStage1` | Stage 1：联合拟合 d2/d3/d4/s1 | 8 参数共享，`info.layer == 0` 跳过 |
| `PidFitFuncStage2` | Stage 2：固定 d2，拟合 d1 | 2 参数 + 固定 d2_p0/d2_p1 |

详细的拟合函数数学推导见第 6.4 节。

### 2.5 参数跨层对共享（全局拟合）

v1 的 Stage 1 中每个探测器（除端点外）同时出现在两个相邻层对中，使用同一组参数：

```
d2d3 层对: par[0],par[1] (d2) + par[2],par[3] (d3)
d3d4 层对: par[2],par[3] (d3) + par[4],par[5] (d4)
d4s  层对: par[4],par[5] (d4) + par[6],par[7] (t0s)
```

d2/d3/d4 的刻度系数同时受两个层对约束，属于**全局拟合**。

### 2.6 最终参数组装

```
final_parameters[10] = {
    stage2_pars[0], stage2_pars[1],  // d1 ← Stage 2
    stage1_pars[0], stage1_pars[1],  // d2 ← Stage 1
    stage1_pars[2], stage1_pars[3],  // d3 ← Stage 1
    stage1_pars[4], stage1_pars[5],  // d4 ← Stage 1
    stage1_pars[6], stage1_pars[7]   // s1 ← Stage 1
};
```

### 2.7 权重缩放配置

`kT0PidInfo` 数组第 7 个字段 `weight` 控制该 PID 在拟合中的缩放倍率。原理：

```
FillGcaliForLayers:   gcali.y  = weight × deep_ADC        （数据缩放）
PidFitFunc::operator: f(x)     = weight × prediction      （函数缩放）

χ² 贡献 = Σ(weight × y − weight × f)² = weight² × Σ(y − f)²
      → 有效权重 = weight²
```

即 `weight=2.0` 对应 4× 有效权重，`weight=5.0` 对应 25× 有效权重。

**当前配置**（[calibrate_t0_utils.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/t0/calibrate_t0_utils.cpp#L39-L53)）：

| 层对 | 粒子 | weight | 有效权重 (weight²) |
|------|------|--------|-------------------|
| d1d2 | ⁴He | 1.0 | 1× |
| d1d2 | ⁶Li | 1.0 | 1× |
| d2d3 | ⁴He | **5.0** | **25×** |
| d2d3 | ⁷Be | 2.0 | 4× |
| d2d3 | ¹²C | 1.0 | 1× |
| d3d4 | ¹H | **16.0** | **256×** |
| d3d4 | ⁴He | 4.0 | 16× |
| d3d4 | ⁷Be | 2.0 | 4× |
| d3d4 | ¹²C | 1.0 | 1× |
| d4s | ¹H | **8.0** | **64×** |
| d4s | ⁴He | 2.0 | 4× |
| d4s | ⁶Li | 1.0 | 1× |

**应用范围**：Stage 1（FillGcaliForLayers + PidFitFuncStage1）和 Stage 2（FillGcaliForLayers + PidFitFuncStage2）均已完整支持。

### 2.8 输出

详见第 9 节。

---

## 3. calibrate_t0_v2（实验性程序）

### 3.1 刻度策略：四段式逐步拟合

由于 d3 厚度可能存在测量误差（d3-d4 是拟合最差的层对），v2 采用从深层到浅层逐步固定参数的单向链式拟合：

```
Stage 1: d3-d4  BothFree（4 参数: d3, d4）
                    ↓ 固定 d4
Stage 2: d4-s1   FixLower（2 参数: s1，d4 固定）
                    ↓ 固定 d3（来自 Stage 1）
Stage 3: d2-d3   FixUpper（2 参数: d2，d3 固定）
                    ↓ 固定 d2
Stage 4: d1-d2   FixUpper（2 参数: d1，d2 固定）
```

**设计动机**：
- 先通过 d3-d4 联合拟合确定 d3 和 d4，让两者互相约束（而非让 d4 与 s1 联合拟合）
- 再用已固定的 d4 确定 s1，避免 s1 的不确定性反向传递到 d4
- 逐级向浅层推进，每一步只拟合一层

### 3.2 拟合函数类

v2 需要三个拟合类来处理"固定浅层"和"固定深层"两种不同情况：

| 类名 | 用途 | 固定/拟合关系 |
|------|------|-------------|
| `PidFitFuncBothFree` | Stage 1：同时拟合一对探测器的两层 | 无固定，4 参数 |
| `PidFitFuncFixLower` | Stage 2：固定浅层（d4），拟合深层（s1） | x=浅层ADC → Energy(layer, de) → (残余-p0_deep)/p1_deep |
| `PidFitFuncFixUpper` | Stage 3/4：固定深层（d3/d2），拟合浅层（d2/d1） | x=浅层ADC → Energy(layer, de) → (残余-p0_deep_fixed)/p1_deep_fixed |

`PidFitFuncFixLower` 和 `PidFitFuncFixUpper` 的核心区别在于哪一层的参数固定。对于 d4-s1，gcali 的 x 轴是 d4（浅层）ADC，拟合的是 s1（深层）——"固定浅层"，需要 `FixLower`。对于 d2-d3 和 d1-d2，gcali 的 x 轴是浅层 ADC，拟合的是浅层，"固定深层"，需要 `FixUpper`。

### 3.3 各阶段参数限制

| 阶段 | 层对 | 拟合层 | TF1 范围 | p0 限制 | p1 限制 |
|------|------|--------|----------|---------|---------|
| Stage 1 | d3-d4 | d3 | [100000, 180000] | [-1.0, 1.0] | [0.006, 0.008] |
| | | d4 | | [-1.0, 1.0] | [0.004, 0.006] |
| Stage 2 | d4-s1 | s1 | [180000, 220000] | [-1.0, 1.0] | [0.0, 0.1] |
| Stage 3 | d2-d3 | d2 | [25000, 102000] | [-1.0, 1.0] | [0.0, 0.1] |
| Stage 4 | d1-d2 | d1 | [0, 25000] | [-1.0, 1.0] | [0.0, 0.1] |

### 3.4 与 v1 的关键区别

| | v1 | v2 |
|---|---|---|
| 迭代方式 | d2-d4-s1 联合 + d1 独立 | d3-d4 → d4-s1 → d2-d3 → d1-d2 链式 |
| d3/d4 关系 | 与 d2、s1 联合同步约束 | d3-d4 先拟合，d4 固定后再拟合 s1 |
| 依赖方向 | Stage 1 中多探测器互相约束 | 单向链：d4→s1→d3→d2→d1 |
| d3/d4 限制 | 与其他层统一 [0,1] | 更严格独立限制（d3 p1 [0.006,0.008], d4 p1 [0.004,0.006]） |

---

## 4. calibrate_t0（原始程序，已弃用）

### 4.1 刻度策略：全局联合拟合

原始程序将全部 10 个参数（d1–s1）一次性全局拟合：

```
d1-d2, d2-d3, d3-d4, d4-s1  全局联合拟合（10 参数）
```

### 4.2 全局拟合参数表（10 参数）

| 参数索引 | 对应探测器 | 含义 |
|---------|-----------|------|
| par[0], par[1] | t0d1 | p0, p1 — d1 刻度系数 |
| par[2], par[3] | t0d2 | p0, p1 — d2 刻度系数 |
| par[4], par[5] | t0d3 | p0, p1 — d3 刻度系数 |
| par[6], par[7] | t0d4 | p0, p1 — d4 刻度系数 |
| par[8], par[9] | t0s | p0, p1 — t0s 刻度系数 |

### 4.3 初始刻度参数与限制

```cpp
double initial_calibration_parameters[10] = {
    0.0, 0.002,   // t0d1: p0=0, p1=0.002
    0.0, 0.006,   // t0d2: p0=0, p1=0.006
    0.0, 0.006,   // t0d3: p0=0, p1=0.006
    0.0, 0.003,   // t0d4: p0=0, p1=0.003
    0.0, 0.003    // t0s:  p0=0, p1=0.003
};
```

参数限制：p0 [0.0, 100.0]，p1 [0.0, 1.0]。

### 4.4 存在的问题

- d1 厚度（68 μm）不均匀导致 d1 参数不稳定，影响 d2/d3/d4/s1 的拟合
- 10 参数全局拟合自由度大，收敛性差
- 因此被 v1（标准程序）取代

---

## 5. 输入数据

### 5.1 输入数据：pre_calibration TGraph

程序读取 `estimate/pre_calibration_{trigger}{run}_*.root` 文件（如 `pre_calibration_t1057_0107.root`），该文件由 `pre_calibration` 程序生成。

文件包含 4 个 TGraph 对象，每个对应一个探测器层对：

| 对象名 | 层对 | x 轴 | y 轴 | 说明 |
|--------|------|------|------|------|
| `g_d1d2` | d1d2 | d2 原始 ADC | d1 原始 ADC | 符合 pre_calibration 条件的 d1-d2 事件 |
| `g_d2d3` | d2d3 | d3 原始 ADC | d2 原始 ADC | 符合 pre_calibration 条件的 d2-d3 事件 |
| `g_d3d4` | d3d4 | d4 原始 ADC | d3 原始 ADC | 符合 pre_calibration 条件的 d3-d4 事件 |
| `g_d4s` | d4s | s 原始 ADC | d4 原始 ADC | 严格 4-hit + s.valid 条件的 d4-s 事件 |

每个 TGraph 点的坐标惯例与 TCutG::IsInside(x,y) 一致，即 **x=深层探测器能量，y=浅层探测器能量**。

### 5.2 pre_calibration 程序的筛选条件

TGraph 中的数据点由 `pre_calibration` 程序根据以下条件筛选：

| 层对 | 数据点来源 | track 条件 |
|------|-----------|-----------|
| d1d2 | 取各自探测器的第一个 hit | `dx² + dy² ≤ 4 mm²`（圆形窗口） |
| d2d3 | 取各自探测器的第一个 hit | `dx² + dy² ≤ 4 mm²`（圆形窗口） |
| d3d4 | 取各自探测器的第一个 hit | `dx² + dy² ≤ 4 mm²`（圆形窗口） |
| d4s | 取 d4 第一个 hit + s 能量 | `d3.num==d4.num==1 && dx²+dy² ≤ 4 mm² && s.valid` |

> **距离阈值**由 `config.toml` 中 `[pre_calibration]` 段的 `max_distance_sq` 控制，默认值 4.0 mm²。

此筛选条件做了以下改进：
- **矩形 → 圆形**：`|dx|<2 && |dy|<2` 改为 `dx²+dy² ≤ 4`，物理上更合理（圆形区域内接正方形会多选 27% 的对角区域事件）
- **d4s 的 hit 条件**：从 `d1.num==d2.num==d3.num==d4.num==1`（全局 4-hit 一致）改为 `d3.num==d4.num==1`（仅约束相邻 d3-d4）

### 5.3 探测器配置

从 `config.toml` 读取探测器厚度参数（[config.h](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/config.h#L83-L85)）：

```toml
[t0]
silicon = ["t0d1", "t0d2", "t0d3", "t0d4", "t0s"]

[detectors.t0d1]
thickness_um = 68.0

[detectors.t0d2]
thickness_um = 1005.0

[detectors.t0d3]
thickness_um = 995.0

[detectors.t0d4]
thickness_um = 999.0

[detectors.t0s]
thickness_um = 1500.0
```

---

## 6. PID 截断（Cut）来源

### 6.1 TCutG 截断文件

在本版本中，PID 截断不再硬编码在程序中，而是通过独立的 TCutG 文件定义。

**文件存放位置**：`src/brill/Cut/`

**命名格式**：`{层对}_{粒子}_stop.C`

例如：
- `src/brill/Cut/d1d2_4He_stop.C`
- `src/brill/Cut/d2d3_4He_stop.C`
- `src/brill/Cut/d2d3_7Be_stop.C`
- `src/brill/Cut/d2d3_12C_stop.C`

**文件格式**（ROOT 宏，`.C` 后缀）：

```cpp
{
   std::vector<Double_t> cutg_vect0{
      x0, x1, x2, ...
   };
   std::vector<Double_t> cutg_vect1{
      y0, y1, y2, ...
   };
   TCutG *cutg = new TCutG("d1d2_4He_stop", N, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("e2");
   cutg->SetVarY("e1");
   cutg->SetFillStyle(1000);
   cutg->SetLineColor(2);
   cutg->SetLineWidth(2);
}
```

### 6.2 TCutG 加载机制

程序启动时，`LoadCuts()` 函数（[calibrate_t0_utils.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/t0/calibrate_t0_utils.cpp#L66-L118)）自动扫描 `src/brill/Cut/` 目录：

1. 对每个层对（`d1d2`、`d2d3`、`d3d4`、`d4s`），匹配以层对名开头的 `_stop.C` 文件
2. 从文件名中解析粒子名称（如 `d2d3_7Be_stop.C` → 粒子 `7Be`）
3. 解析粒子名为 (Z, A)（`4He`→{2,4}，`7Be`→{4,7}，`12C`→{6,12}）
4. 通过 `gROOT->ProcessLine(".x 文件名")` 执行宏，加载 TCutG
5. 通过 `gROOT->FindObject("d2d3_7Be_stop")` 获取 TCutG 指针并 Clone 保存

### 6.3 截断筛选逻辑（双重筛选：TCutG + left/right）

对每个加载的 cut，遍历对应层对的 TGraph 中所有数据点，同时使用 TCutG 和 pid_info 的 left/right 范围进行双重筛选：

```cpp
for (int pt = 0; pt < g->GetN(); ++pt) {
    double deep = g->GetPointX(pt);    // 深层探测器 ADC (TGraph x 轴)
    double shallow = g->GetPointY(pt);  // 浅层探测器 ADC (TGraph y 轴)
    if (cut.cut->IsInside(deep, shallow)   // TCutG 筛选：在 ΔE-E 截断曲线内
        && shallow > info->left            // left/right 筛选：浅层 ADC 在合理范围内
        && shallow < info->right) {
        gcali.AddPoint(shallow + info->offset, deep);
    }
}
```

双重筛选的意义：

| 筛选方式 | 作用 | 如果没有会怎样 |
|---------|------|---------------|
| `TCutG::IsInside` | ΔE-E 二维形状筛选，排除不在截断曲线内的点 | 无法区分不同粒子，不同粒子的散点混在一起 |
| `shallow > left && shallow < right` | 浅层能量一维范围筛选，确保数据点在拟合函数的有效区间内 | 数据点落入 PidFitFunc 的区间间隙（return 0），变成巨大 outlier 导致 chi2 飙升 |

**这与 reference 程序的行为一致**：reference 中同样同时检查 `event.layer`, `event.charge`, `event.mass`（粒子身份）**和** `event.energy[] > left && event.energy[] < right`（能量范围）。

命中的数据点被加入刻度用的 TGraph `gcali`，x 坐标加上层对 offset 用于后续 PidFitFunc 的区间路由。

### 6.4 pid_info 的双重角色

`pid_info` 表在本版本中承担**两个角色**：

| 角色 | 使用位置 | 作用 |
|------|---------|------|
| 散点筛选 | `main()` 中的 `gcali.AddPoint` 循环 | `left/right` 作为浅层 ADC 的额外约束，筛掉落在区间外的数据点 |
| 区间路由 | `PidFitFunc::operator()` | 将 gcali 的 x 坐标路由到正确粒子的理论曲线 |

两个角色共用同一组 `left/right/offset`，确保筛选条件和拟合路由完全一致。

### 6.5 各层对可用粒子

| 层对 | 可用粒子 | TCutG 文件 |
|------|---------|-----------|
| d1d2 | ⁴He | `d1d2_4He_stop.C` |
| d2d3 | ⁴He, ⁷Be, ¹²C | `d2d3_4He_stop.C`, `d2d3_7Be_stop.C`, `d2d3_12C_stop.C` |
| d3d4 | ⁴He, ⁷Be, ¹²C | `d3d4_4He_stop.C`, `d3d4_7Be_stop.C`, `d3d4_12C_stop.C` |
| d4s | ¹H, ⁴He, ⁶Li | `d4s_1H_stop.C`, `d4s_4He_stop.C`, `d4s_6Li_stop.C` |

> **注**：d1d2 的 ⁶Li 暂不参与拟合（`pid_info` 中已注释）。引入后拟合崩溃，需待 TCutG 优化后重新启用。

---

## 7. 理论曲线来源

### 7.1 能损计算链

理论曲线通过以下计算链得到：

```
catima 库（核物理能损计算）
    ↓
RangeEnergyCalculator（射程-能量关系，硅材料）
    ↓
DeltaEnergyCalculator（ΔE-E 理论曲线，多层硅探测器）
    ↓
PidFitFunc（多粒子 + 多层拟合函数）
```

### 7.2 RangeEnergyCalculator（射程-能量计算器）

定义于 [range_energy_calculator.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/energy_calculator/range_energy_calculator.cpp#L29-L47)，基于 [catima](https://github.com/hrosiak/catima) 库（C++ 版本的 ATIMA/STOPPING 能损计算库）。

**输入**：

| 输入 | 来源 | 说明 |
|------|------|------|
| `charge` (Z) | 自动从 cut 文件名推断 | 粒子电荷数，如 4He → Z=2 |
| `mass` (A) | 自动从 cut 文件名推断 | 粒子质量数，如 4He → A=4 |
| `material` | `SiliconMaterial()` | 硅材料：密度 2.329 g/cm³，单质 Si（Z=14, A=28） |
| `cache_path` | 由 `RangeCachePath()` 生成 | 缓存文件路径，见下方 |
| 最大能量 | `mass * 500.0` MeV | 粒子质量数 × 500 MeV，如 4He → 2000 MeV |

**输出（缓存文件）**：

```
<workspace>/energy_calculator/si_z{Z}_a{A}.root
```

例如 `/data/disk1/ribll2026_www_data/energy_calculator/si_z2_a4.root`（4He 在硅中的射程-能量曲线）。

文件中包含两个 TSpline3 对象（均为 1200 个点）：

| 对象名 | 方向 | x 轴 | y 轴 | 说明 |
|--------|------|------|------|------|
| `re` | 能量 → 射程 | 入射能量 (MeV) | 射程 (μm) | 给定能量，求粒子在硅中能走多远 |
| `er` | 射程 → 能量 | 射程 (μm) | 剩余能量 (MeV) | 给定剩余射程，求粒子还剩多少能量 |

**存储路径规则**（[range_energy_calculator.cpp:L54-L62](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/energy_calculator/range_energy_calculator.cpp#L54-L62)）：

```cpp
std::string RangeCachePath(const AppConfig &config, int charge, int mass) {
    return TString::Format(
        "%s/si_z%d_a%d.root",
        JoinPath(config.workspace, config.paths.energy_calculator).c_str(),
        charge, mass
    ).Data();
}
```

其中 `config.paths.energy_calculator` 默认为 `"energy_calculator"`（[config.h:L98](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/config.h#L98)），`config.workspace` 默认 `/data/disk1/ribll2026_www_data`。

**缓存机制**：构造函数中先调用 `LoadSplines()` 尝试从缓存文件读取，若文件不存在或读取失败则调用 `BuildSplines()` 重新生成并写入缓存。**首次运行 calibrate_t0 时自动生成，无需手动操作。**

**计算原理**（`BuildSplines`）：

1. 以 `(当前能量/最大能量)²` 的非线性步长在能量范围内取 1200 个点（低能区更密集）
2. 对每个能量点调用 `catima::range()` 计算射程，并将单位从 g/cm² 转换为 μm
3. 确保射程单调递增（若出现平台则加 ε=1e-6）
4. 构建两个方向的 TSpline3 三次样条插值

---

### 7.3 DeltaEnergyCalculator（ΔE-E 理论曲线计算器）

定义于 [delta_energy_calculator.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/energy_calculator/delta_energy_calculator.cpp#L62-L92)，计算粒子在相邻两层硅探测器中的 ΔE-E 理论关系。

**输入**：

| 输入 | 来源 | 说明 |
|------|------|------|
| `charge` (Z) | 调用方传入 | 粒子电荷数 |
| `mass` (A) | 调用方传入 | 粒子质量数 |
| `config` | `config.toml` | 读取 `t0.silicon` 探测器列表及各探测器厚度 |
| 探测器厚度 | `config.toml` → `detector.thickness_um` | 如 t0d1=68μm, t0d2=1005μm, t0d3=995μm, t0d4=999μm, t0s=1500μm |
| 射程-能量数据 | `RangeEnergyCalculator` 的缓存文件 | 即 `si_z{Z}_a{A}.root`（自动依赖） |

**输出（缓存文件）**：

```
<workspace>/energy_calculator/t0_delta_z{Z}_a{A}.root
```

例如 `/data/disk1/ribll2026_www_data/energy_calculator/t0_delta_z2_a4.root`（4He 在 T0 探测器中的 ΔE-E 理论曲线）。

文件中包含 4 个层对（d1d2, d2d3, d3d4, d4s），每个层对包含两个 TSpline3 对象：

| 对象名 | 方向 | x 轴 | y 轴 | 说明 |
|--------|------|------|------|------|
| `e_de_0` | ΔE → E | 在 d1 中沉积的能量 (MeV) | 穿透 d1 后的剩余能量 (MeV) | d1d2 层对 |
| `de_e_0` | E → ΔE | 穿透 d1 后的剩余能量 (MeV) | 在 d1 中沉积的能量 (MeV) | d1d2 层对（反向） |
| `e_de_1` | ΔE → E | 在 d2 中沉积的能量 (MeV) | 穿透 d2 后的剩余能量 (MeV) | d2d3 层对 |
| `de_e_1` | E → ΔE | 穿透 d2 后的剩余能量 (MeV) | 在 d2 中沉积的能量 (MeV) | d2d3 层对（反向） |
| `e_de_2` | ΔE → E | 在 d3 中沉积的能量 (MeV) | 穿透 d3 后的剩余能量 (MeV) | d3d4 层对 |
| `de_e_2` | E → ΔE | 穿透 d3 后的剩余能量 (MeV) | 在 d3 中沉积的能量 (MeV) | d3d4 层对（反向） |
| `e_de_3` | ΔE → E | 在 d4 中沉积的能量 (MeV) | 穿透 d4 后的剩余能量 (MeV) | d4s 层对 |
| `de_e_3` | E → ΔE | 穿透 d4 后的剩余能量 (MeV) | 在 d4 中沉积的能量 (MeV) | d4s 层对（反向） |

**缓存机制**：构造函数中先调用 `Load()` 尝试从缓存文件读取，若文件不存在或读取失败则调用 `Initialize()` 重新生成并写入缓存。**首次运行 calibrate_t0 时自动生成，无需手动操作。**

**计算原理**（`BuildSliceFunctions`）：

1. 确定粒子能量范围：从刚好能穿透第一层探测器（`thickness_first`）的最小能量，到刚好能穿透两层探测器总厚度（`thickness_first + thickness_second`）的最大能量
2. 以 **0.1 MeV** 步长遍历能量范围
3. 对每个入射能量 `E_total`：
   - `residual_range = Range(E_total) - thickness_first` — 穿透第一层后剩余的射程
   - `E_residual = Energy(residual_range)` — 剩余射程对应的能量
   - `ΔE = E_total - E_residual` — 在第一层中沉积的能量
4. 对生成的 (ΔE, E) 点对按 x 轴排序，确保单调性（解决 Bragg 峰导致的非单调问题）
5. 构建两个方向的 TSpline3 三次样条插值

**调用关系**：`DeltaEnergyCalculator` 内部依赖 `RangeEnergyCalculator`。在 `Initialize()` 中，先创建 `RangeEnergyCalculator` 对象（自动加载或生成 `si_z{Z}_a{A}.root`），然后利用其 `Range()` 和 `Energy()` 接口计算各层对的 ΔE-E 曲线。

### 7.4 PidFitFunc — 拟合函数数学推导

以 d2d3 层对为例，所有版本的拟合函数遵循相同的物理模型。拟合 TGraph 的坐标约定为：

```
gcali:  x = d2 的原始 ADC + offset   （浅层 = 先被粒子击中的探测器）
        y = d3 的原始 ADC            （深层 = 后被粒子击中的探测器）
```

拟合函数 `f(x)` 返回**预测的 d3 ADC 值**，Minuit 极小化 Σ(yᵢ − f(xᵢ))²：

```
步骤 1: 将 x 还原为浅层 ADC
        shallow_raw = x − offset                         (d2 的原始 ADC)

步骤 2: 对浅层 ADC 刻度 → 浅层沉积能量 ΔE
        ΔE = p0_d2 + p1_d2 × shallow_raw                (d2 中的沉积能量, MeV)

步骤 3: 通过理论曲线 ΔE → 深层剩余能量
        E_residual = Energy(layer=1, ΔE)                 (穿透 d2 后的剩余能量, MeV)

步骤 4: 将深层 MeV 转换回深层 ADC（预测值）
        deep_pred = (E_residual − p0_d3) / p1_d3         (预测的 d3 ADC)
```

ROOT 拟合器最小化 `(y_measured − deep_pred)²`。

#### 7.4.1 为什么参数 p1 的位置至关重要

上述公式中，**p1_shallow（p1_d2）出现在分子中**：

```
ΔE = p0_d2 + p1_d2 × shallow_raw
```

梯度：`∂(ΔE)/∂(p1_d2) = shallow_raw` ∝ ADC（典型值 ~10³–10⁴），梯度大小适中，Minuit 可以稳定探索。

而 **p1_deep（p1_d3）出现在分母中**：

```
deep_pred = (E_residual − p0_d3) / p1_d3
```

梯度：`∂(deep_pred)/∂(p1_d3) = −(E_residual − p0_d3) / p1_d3²`，**但 p1_d3 ~ 0.006**（数量级合理），且分子 `E_residual − p0_d3` 在拟合收敛时也较小（~0 附近），因此该梯度仍在可控范围。

---

#### 7.4.2 ★ 交换 x/y 轴导致拟合发散的根本原因

**如果将 TGraph 的 x/y 轴交换**（即 x = 深层 ADC + offset, y = 浅层 ADC）：

拟合函数变为 `f(x)` 返回**预测的浅层 ADC 值**：

```
步骤 1: 深层 ADC → 深层 MeV
        deep_mev = p0_deep + p1_deep × deep_raw

步骤 2: 深层 MeV → 浅层 ΔE
        ΔE = DeltaEnergy(layer, deep_mev)

步骤 3: 浅层 MeV → 浅层 ADC
        shallow_pred = (ΔE − p0_shallow) / p1_shallow
```

**问题出在第 3 步**：此时 **p1_shallow（如 p1_d1）出现在分母中**：

```
shallow_pred = (ΔE − p0_d1) / p1_d1
```

梯度：`∂(shallow_pred)/∂(p1_d1) = −(ΔE − p0_d1) / p1_d1²`

**致命问题**：p1_d1 的典型值约 **0.0007**（d1 是 68 μm 薄探测器，单位 ADC 对应能量极小），因此：

- `1 / p1_d1² ≈ 1 / (0.0007)² ≈ 2×10⁶`
- 当 ΔE 偏离 p0_d1 时，梯度被放大 200 万倍，Minuit 的梯度下降步长剧烈振荡
- 其他探测器（d2~d4）的 p1 也在 0.003–0.007 数量级，同样面临 `1/p1²` 放大问题
- Minuit 无法在这种病态梯度面上找到稳定的极小点 → **状态码 3（Abnormal termination），EDM 不收敛**

**对比总结**：

| 方案 | TGraph (x, y) | p1 在哪端 | 梯度数量级 | 稳定性 |
|------|---------------|-----------|-----------|--------|
| ✅ 当前版本 | (浅层ADC, 深层ADC) | p1_shallow 在分子 | ∝ ADC (~10⁴) | 稳定收敛 |
| ❌ 交换版本 | (深层ADC, 浅层ADC) | p1_shallow 在分母 | ∝ 1/p1² (~10⁶) | 发散 |

**本质上，哪个轴的 p1 出现在分母，哪个轴对应的探测器（通常更薄、p1 更小）就会放大梯度噪声。当前版本让较厚探测器（深层，p1 较大~0.006）的参数在分母，而薄探测器（浅层，p1 较小~0.0007）的参数在分子——这恰好是数值上更安全的方向。**

### 7.5 区间路由机制

`pid_info` 数组定义每个 (层对, 粒子) 组合对应的 x 区间和 offset：

```cpp
const std::vector<T0ParticlePidInfo> kT0PidInfo = {
    //layer charge mass left    right   offset  weight
    {0, 2,  4,  2250.0, 12000.0,       0, 1.0}, //d1d2 4He
    {0, 3,  6,  4355.0, 12000.0,   13000, 1.0}, //d1d2 6Li
    {1, 2,  4,  3890.0,  8250.0,   26000, 5.0}, //d2d3 4He
    {1, 4,  7, 10930.0, 17320.0,   36000, 2.0}, //d2d3 7Be
    {1, 6, 12, 22050.0, 45800.0,   55000, 1.0}, //d2d3 12C
    {2, 1,  1,   820.0,  1700.0,  102000, 16.0}, //d3d4 1H
    {2, 2,  4,  3400.0,  7100.0,  105000, 4.0}, //d3d4 4He
    {2, 4,  7,  9600.0, 18200.0,  113000, 2.0}, //d3d4 7Be
    {2, 6, 12, 22050.0, 45800.0,  133000, 1.0}, //d3d4 12C
    {3, 1,  1,   960.0,  2300.0,  183000, 8.0}, //d4s 1H
    {3, 2,  4,  3890.0,  9430.0,  187000, 2.0}, //d4s 4He
    {3, 3,  6,  7550.0, 16600.0,  203000, 1.0}, //d4s 6Li
};
```

当 `x` 落入某个 `[left+offset, right+offset]` 区间时，PidFitFunc 自动选择对应的 `(layer, charge, mass)` 进行拟合计算，并应用对应的 `weight` 缩放（参见 §2.7）。

**结构体字段说明**：

| 字段 | 类型 | 说明 |
|------|------|------|
| layer | int | 层对索引（0=d1d2, 1=d2d3, 2=d3d4, 3=d4s） |
| charge | int | 电荷数 Z |
| mass | int | 质量数 A |
| left | double | 浅层 ADC 下限 |
| right | double | 浅层 ADC 上限 |
| offset | double | 全局 x 区间偏移 |
| weight | double | 缩放倍率，有效权重 = weight² |

**offset 布局原则**：每个 (层对, 粒子) 组合分配一段独立的全局 x 区间，区间之间留有 ≥2000 ADC 的间隔（gap），确保各段互不重叠。所有 offset 单调递增。TF1 总范围 0–330000 覆盖所有区间。

### 7.6 ★ left/right 对理论曲线的截断效应

**这是本程序最容易被误解的设计点。** `left/right` 不仅用于筛选数据点和路由拟合区间，还对 PidFitFunc 内部使用的理论曲线产生了**截断效应**。

**拟合时（PidFitFunc）**：只使用理论曲线中落入 `[left+offset, right+offset]` 的一段。

```
PidFitFunc::operator() 的执行逻辑：
  if (x 不在任何 [left+offset, right+offset] 内)
      return 0.0;    ← 越界→残差 = (deep_ADC - 0)² = 极大！
  // 在区间内，正常计算：
  de = p0_shallow + p1_shallow × (x - offset);  → 浅层沉积能量
  e  = calculator→Energy(layer, de);            → 深层剩余能量（理论）
  return (e - p0_deep) / p1_deep;               → 预测深层 ADC
```

`DeltaEnergyCalculator::Energy(layer, de)` **内部使用的仍然是完整的理论曲线**，但调用入口 `de` 被 `left/right` 限制在了浅层 ADC 对应的能量范围内，因此实际上只查表了其中一段。

**TH2F 展示时（DrawTheoryCurves）**：绘制完整的理论曲线，不经过任何截断。

```
                 完整 de-e 理论曲线
                 ╱                ╲
                ╱                  ╲
       ΔE      ╱                    ╲                ← TH2F 绘制整条红线
      (浅层)  ╱                      ╲
             ╱                        ╲
            ╱    ┌──────────────┐      ╲
           ╱     │ left ~ right │       ╲              ← PidFitFunc 只取这一段
          ╱      │ PidFitFunc   │        ╲
         ╱       │ 拟合用区间    │         ╲
        ──────────────────────────────────
                 E (深层能量)
```

**后果对比**：

| | PidFitFunc（拟合） | TH2F 理论曲线（展示） |
|---|---|---|
| 曲线范围 | 被 `left/right` 截断 | 完整物理曲线 |
| 越界行为 | return 0 → 巨大残差 → 拟合失败 | 不存在越界 |
| 数据一致性 | 需要 left/right 覆盖所有数据点 | 总是覆盖整个图的能量范围 |

**关键约束**：`left/right` 必须**完整覆盖**所有通过 TCutG 筛选的数据点的浅层 ADC 范围。如果出现过小的 `left/right`（数据点落在区间外），`PidFitFunc` 返回 0，这些点会变成巨大的 outlier，导致 chi2 飙升、Minuit 无法收敛。

**合理设置 left/right 的原则**：
1. 对所有 TCutG 筛选后的数据点，统计其浅层 ADC 的最小值和最大值
2. `left` ≤ 最小浅层 ADC，`right` ≥ 最大浅层 ADC（留少许余量）
3. 相邻粒子的区间之间留 gap ≥ 2000 ADC，避免重叠
4. 同一层对的不同粒子使用相同的 p0/p1，共同约束该层对的刻度系数

---

## 8. 算法流程

### 8.1 v1 main() 函数流程

```
1. 解析命令行参数
2. 加载配置文件
3. 验证探测器配置存在
4. 提示重建理论曲线缓存（如有需要）
5. 在 estimate/ 目录中查找 pre_calibration 文件
6. 读取 g_d1d2/g_d2d3/g_d3d4/g_d4s 四个 TGraph
7. 扫描 src/brill/Cut/ 目录，加载所有符合命名格式的 TCutG
8. Stage 1：对 layer=1,2,3 的 cuts 执行 TCutG + left/right 双重筛选 → gcali_stage1
9. 构建 PidFitFuncStage1（8 参数），TF1 范围 [25000, 220000]
10. 拟合 Stage 1 → 获得 d2,d3,d4,s1 参数
11. Stage 2：对 layer=0 的 cuts 执行筛选 → gcali_stage2
12. 构建 PidFitFuncStage2（2 参数 + 固定 d2），TF1 范围 [0, 25000]
13. 拟合 Stage 2 → 获得 d1 参数
14. 组装 final_parameters[10]
15. 输出刻度参数到 calibration/t0_{run:04d}.txt
16. 校准 TH2F 并绘制叠加理论曲线的 Canvas
17. 保存 gcali 和拟合函数到 calibration/t0_{trigger}{run:04d}.root
```

### 8.2 v2 main() 函数流程

```
1-7. 同 v1
8.  Stage 1：对 layer=2 的 cuts 执行筛选 → gcali_stage1
9.  构建 PidFitFuncBothFree（4 参数 d3+d4），TF1 范围 [100000, 180000]
10. 拟合 Stage 1 → 获得 d3,d4 参数
9.  Stage 2：对 layer=3 的 cuts 执行筛选 → gcali_stage2
10. 构建 PidFitFuncFixLower（2 参数 s1 + 固定 d4），TF1 范围 [180000, 220000]
11. 拟合 Stage 2 → 获得 s1 参数
12. Stage 3：对 layer=1 的 cuts 执行筛选 → gcali_stage3
13. 构建 PidFitFuncFixUpper（2 参数 d2 + 固定 d3），TF1 范围 [25000, 102000]
14. 拟合 Stage 3 → 获得 d2 参数
15. Stage 4：对 layer=0 的 cuts 执行筛选 → gcali_stage4
16. 构建 PidFitFuncFixUpper（2 参数 d1 + 固定 d2），TF1 范围 [0, 25000]
17. 拟合 Stage 4 → 获得 d1 参数
18. 组装 final_parameters[10]
19-20. 同 v1
```

### 8.3 原始程序 main() 函数流程

```
1. 解析命令行参数
2. 加载配置文件
3. 在 estimate/ 目录中查找 pre_calibration 文件
4. 读取 g_d1d2/g_d2d3/g_d3d4/g_d4s 四个 TGraph
5. 扫描 src/brill/Cut/ 目录，加载所有 TCutG
6. 遍历每个 cut，对 TGraph 逐点执行双重筛选，命中点加入 gcali（所有 4 层对）
7. 从加载的 cut 中自动收集粒子列表，构建 PidFitFunc（10 参数）
8. 初始化刻度参数
9. 将 TGraph 与 TF1 进行全局拟合（10 参数）
10. 输出刻度参数到 calibration/t0_{run}.txt
11. 保存 gcali 到 calibration/t0_{trigger}{run}.root
```

---

## 9. 输出文件

### 9.1 calibration/t0_{run:04d}.txt

例如 `calibration/t0_0057.txt`（以起始 run 命名，4 位补齐）。

格式：
```
# layer p0 p1
0 <p0> <p1>
1 <p0> <p1>
2 <p0> <p1>
3 <p0> <p1>
4 <p0> <p1>
```

| layer 索引 | 对应探测器 |
|-----------|-----------|
| 0 | t0d1 |
| 1 | t0d2 |
| 2 | t0d3 |
| 3 | t0d4 |
| 4 | t0s |

### 9.2 calibration/t0_{trigger}{run:04d}.root

例如 `calibration/t0_t1057.root`。

v1 和 v2 输出文件包含：
- `gcali_stage1`–`gcali_stage4`：各阶段拟合用 TGraph（x 坐标已加 offset）
- `fcali_stage1`–`fcali_stage4`：拟合 TF1 对象
- 校准后 TH2F（`d1d2_cal`、`d2d3_cal`、`d3d4_cal`、`d4s_cal`）
- 叠加理论曲线的 Canvas（`c_d1d2_v1` 等）

---

## 10. 数据流

```
match/*.root + ingot/t0s_*.root
        ↓ pre_calibration
estimate/pre_calibration_{trigger}{run}_*.root
   （TGraph: g_d1d2, g_d2d3, g_d3d4, g_d4s）
        ↓ calibrate_t0_v1（TCutG 筛选 + 理论曲线两阶段拟合）
calibration/t0_{run:04d}.txt + t0_{trigger}{run:04d}.root
```

**数据流**：TCutG 替代了原来的 `track_t0` PID 步骤，直接从 pre_calibration 的原始 ADC 散点中筛选目标粒子进行刻度拟合。

---

## 11. 与 normalize 程序的关系

| 程序 | 刻度对象 | 刻度形式 | 输入 |
|------|---------|---------|------|
| `normalize` | 单个 DSSD 的各条 strip | per-strip: `E_norm = p0 + p1 * E_raw` | 原始 DSSD 数据 |
| `calibrate_t0_v1` | 整个探测器层的能量 | per-layer: `E_cal = p0 + p1 * E_raw` | pre_calibration TGraph + TCutG cuts |

**执行顺序**：`normalize`（条间归一化）→ `match`（正背面匹配）→ `pre_calibration`（生成 PID 图 + TGraph）→ （手动绘制 TCutG 截断）→ `calibrate_t0_v1`（层间能量刻度）。