# T0 能量刻度程序说明文档

## 1. 概述

`calibrate_t0.cpp` 是 brill2 数据分析框架中用于 T0 探测器系统能量刻度（ADC → MeV）的核心程序，可执行文件名为 `calibrate_t0`。

### 物理目标

T0 探测器系统由 5 层硅探测器（t0d1/t0d2/t0d3/t0d4/t0s）组成，每层探测器输出的原始信号为 ADC 道数。能量刻度的目标是为每层探测器确定一组线性刻度系数 `(p0, p1)`，使得：

```
E_calibrated[MeV] = p0 + p1 × E_raw[ADC]
```

### 刻度策略

本程序采用**基于 TCutG + 理论曲线拟合**方法：

1. 从 `pre_calibration` 程序生成的 TGraph 中读取各层对的原始 ADC 散点（`g_d1d2`、`g_d2d3`、`g_d3d4`、`g_d4s`）
2. 通过 `src/brill/Cut/` 目录下的 TCutG 截断文件筛选出目标粒子
3. 利用 catima 库计算粒子在多层硅探测器中的理论能损，建立理论 ΔE-E 曲线
4. 对实验数据点进行全局拟合以确定刻度系数

---

## 2. 命令行参数

| 参数 | 简写 | 类型 | 说明 |
|------|------|------|------|
| `--help` | `-h` | - | 打印帮助信息 |
| `--run` | `-r` | int | pre_calibration 的起始 run 号（必填） |
| `--trigger` | `-t` | string | 触发器类型 |
| `--config` | `-c` | string | 配置文件路径，默认 `config.toml` |

**注意**：不再需要 `-e`（end-run）参数，程序自动在 `estimate/` 目录中查找文件名以 `pre_calibration_{trigger}{run:04d}_` 开头的 ROOT 文件，从中读取 TGraph 数据点。

### 使用示例

```bash
./calibrate_t0 -r 57 -t t1
```

---

## 3. 数据来源

### 3.1 输入数据：pre_calibration TGraph

程序读取 `estimate/pre_calibration_{trigger}{run}_*.root` 文件（如 `pre_calibration_t1057_0107.root`），该文件由 `pre_calibration` 程序生成。

文件包含 4 个 TGraph 对象，每个对应一个探测器层对：

| 对象名 | 层对 | x 轴 | y 轴 | 说明 |
|--------|------|------|------|------|
| `g_d1d2` | d1d2 | d2 原始 ADC | d1 原始 ADC | 符合 pre_calibration 条件的 d1-d2 事件 |
| `g_d2d3` | d2d3 | d3 原始 ADC | d2 原始 ADC | 符合 pre_calibration 条件的 d2-d3 事件 |
| `g_d3d4` | d3d4 | d4 原始 ADC | d3 原始 ADC | 符合 pre_calibration 条件的 d3-d4 事件 |
| `g_d4s` | d4s | s 原始 ADC | d4 原始 ADC | 严格 4-hit + s.valid 条件的 d4-s 事件 |

每个 TGraph 点的坐标惯例与 TCutG::IsInside(x,y) 一致，即 **x=深层探测器能量，y=浅层探测器能量**。

### 3.2 pre_calibration 程序的筛选条件

TGraph 中的数据点由 `pre_calibration` 程序根据以下条件筛选：

| 层对 | 数据点来源 | track 条件 |
|------|-----------|-----------|
| d1d2 | 取各自探测器的第一个 hit | `|dx| < 2mm, |dy| < 2mm` |
| d2d3 | 取各自探测器的第一个 hit | `|dx| < 2mm, |dy| < 2mm` |
| d3d4 | 取各自探测器的第一个 hit | `|dx| < 2mm, |dy| < 2mm` |
| d4s | 取 d4 第一个 hit + s 能量 | `d1.num==d2.num==d3.num==d4.num==1 && s.valid` |

### 3.3 探测器配置

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

## 4. PID 截断（Cut）来源

### 4.1 TCutG 截断文件

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

### 4.2 TCutG 加载机制

程序启动时，`LoadCuts()` 函数（[calibrate_t0.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/calibrate_t0.cpp#L170-L214)）自动扫描 `src/brill/Cut/` 目录：

1. 对每个层对（`d1d2`、`d2d3`、`d3d4`、`d4s`），匹配以层对名开头的 `.C` 文件
2. 从文件名中解析粒子名称（如 `d2d3_7Be_stop.C` → 粒子 `7Be`）
3. 解析粒子名为 (Z, A)（`4He`→{2,4}，`7Be`→{4,7}，`12C`→{6,12}）
4. 通过 `gROOT->ProcessLine(".x 文件名")` 执行宏，加载 TCutG
5. 通过 `gROOT->FindObject("d2d3_7Be_stop")` 获取 TCutG 指针并 Clone 保存

### 4.3 截断筛选逻辑（双重筛选：TCutG + left/right）

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

### 4.4 pid_info 的双重角色

`pid_info` 表（[calibrate_t0.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/calibrate_t0.cpp#L51-L61)）在本版本中承担**两个角色**：

| 角色 | 使用位置 | 作用 |
|------|---------|------|
| 散点筛选 | `main()` 中的 `gcali.AddPoint` 循环 | `left/right` 作为浅层 ADC 的额外约束，筛掉落在区间外的数据点 |
| 区间路由 | `PidFitFunc::operator()` | 将 gcali 的 x 坐标路由到正确粒子的理论曲线 |

两个角色共用同一组 `left/right/offset`，确保筛选条件和拟合路由完全一致。

### 4.5 各层对可用粒子

| 层对 | 可用粒子 | TCutG 文件 |
|------|---------|-----------|
| d1d2 | ⁴He | `d1d2_4He_stop.C` |
| d2d3 | ⁴He, ⁷Be, ¹²C | `d2d3_4He_stop.C`, `d2d3_7Be_stop.C`, `d2d3_12C_stop.C` |
| d3d4 | ⁴He, ⁷Be, ¹²C | `d3d4_4He_stop.C`, `d3d4_7Be_stop.C`, `d3d4_12C_stop.C` |
| d4s | ¹H, ⁴He, ⁷Be | `d4s_1H_stop.C`, `d4s_4He_stop.C`, `d4s_7Be_stop.C` |

> **注**：d1d2 的 ⁶Li 暂不参与拟合（`pid_info` 中已注释）。引入后拟合崩溃，需待 TCutG 优化后重新启用。

---

## 5. 理论曲线来源

### 5.1 能损计算链

理论曲线通过以下计算链得到：

```
catima 库（核物理能损计算）
    ↓
RangeEnergyCalculator（射程-能量关系，硅材料）
    ↓
DeltaEnergyCalculator（ΔE-E 理论曲线，多层硅探测器）
    ↓
PidFitFunc（多粒子 + 多层全局拟合函数）
```

### 5.2 RangeEnergyCalculator（射程-能量计算器）

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

### 5.3 DeltaEnergyCalculator（ΔE-E 理论曲线计算器）

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
| `de_e_0` | ΔE → E | 在 d1 中沉积的能量 (MeV) | 穿透 d1 后的剩余能量 (MeV) | d1d2 层对 |
| `e_de_0` | E → ΔE | 穿透 d1 后的剩余能量 (MeV) | 在 d1 中沉积的能量 (MeV) | d1d2 层对（反向） |
| `de_e_1` | ΔE → E | 在 d2 中沉积的能量 (MeV) | 穿透 d2 后的剩余能量 (MeV) | d2d3 层对 |
| `e_de_1` | E → ΔE | 穿透 d2 后的剩余能量 (MeV) | 在 d2 中沉积的能量 (MeV) | d2d3 层对（反向） |
| `de_e_2` | ΔE → E | 在 d3 中沉积的能量 (MeV) | 穿透 d3 后的剩余能量 (MeV) | d3d4 层对 |
| `e_de_2` | E → ΔE | 穿透 d3 后的剩余能量 (MeV) | 在 d3 中沉积的能量 (MeV) | d3d4 层对（反向） |
| `de_e_3` | ΔE → E | 在 d4 中沉积的能量 (MeV) | 穿透 d4 后的剩余能量 (MeV) | d4s 层对 |
| `e_de_3` | E → ΔE | 穿透 d4 后的剩余能量 (MeV) | 在 d4 中沉积的能量 (MeV) | d4s 层对（反向） |

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

### 5.4 PidFitFunc — 多参数全局拟合函数

定义于 [calibrate_t0.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/calibrate_t0.cpp#L63-L101)，是传递给 ROOT TF1 的拟合函数，即 Minuit 极小化器的目标模型。

#### 5.4.1 拟合参数（共 10 个，每探测器 2 个）

| 参数索引 | 对应探测器 | 含义 |
|---------|-----------|------|
| par[0], par[1] | t0d1 | p0, p1 — d1 刻度系数 |
| par[2], par[3] | t0d2 | p0, p1 — d2 刻度系数 |
| par[4], par[5] | t0d3 | p0, p1 — d3 刻度系数 |
| par[6], par[7] | t0d4 | p0, p1 — d4 刻度系数 |
| par[8], par[9] | t0s | p0, p1 — t0s 刻度系数 |

#### 5.4.2 参数跨层对共享（全局拟合）

每个探测器（除 d1 和 t0s 两个端点外）同时出现在两个相邻层对中，使用同一组参数：

```
d1d2 层对: par[0],par[1] (d1) + par[2],par[3] (d2)
d2d3 层对: par[2],par[3] (d2) + par[4],par[5] (d3)
d3d4 层对: par[4],par[5] (d3) + par[6],par[7] (d4)
d4s  层对: par[6],par[7] (d4) + par[8],par[9] (t0s)
```

d2/d3/d4 的刻度系数同时受两个层对约束，属于**全局拟合**。

#### 5.4.3 拟合函数数学推导（当前正确版本）

以 d2d3 层对为例说明一号公式的物理含义。拟合 TGraph 的坐标约定为：

```
x = d2 的原始 ADC + offset
y = d3 的原始 ADC
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

#### 5.4.4 为什么参数 p1 的位置至关重要

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

#### 5.4.5 ★ 交换 x/y 轴导致拟合发散的根本原因

**如果将 TGraph 的 x/y 轴交换**（即 x = 深层 ADC + offset, y = 浅层 ADC）：

拟合函数变为 `f(x)` 返回**预测的浅层 ADC 值**：

```
步骤 1: 浅层改为深层，深层 ADC → 深层 MeV
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

#### 5.4.6 区间路由机制

`pid_info` 数组定义每个 (层对, 粒子) 组合对应的 x 区间和 offset：

```cpp
const std::vector<ParticlePidInfo> pid_info {
    {0, 2,  4,  2000.0, 11000.0,       0},  // d1d2: 4He  → 全局 x ∈ [    2000,    11000]
    //{0, 3,  6,  4500.0, 15000.0,   13000},  // d1d2: 6Li（暂注释，待优化cut）
    {1, 2,  4,  3900.0,  8000.0,   30000},  // d2d3: 4He  → 全局 x ∈ [   33900,    38000]
    {1, 4,  7, 11500.0, 19000.0,   39000},  // d2d3: 7Be  → 全局 x ∈ [   50500,    58000]
    {1, 6, 12, 22500.0, 46000.0,   59000},  // d2d3: 12C  → 全局 x ∈ [   81500,   105000]
    {2, 2,  4,  3400.0,  7000.0,  106000},  // d3d4: 4He  → 全局 x ∈ [  109400,   113000]
    {2, 4,  7,  9500.0, 19000.0,  114000},  // d3d4: 7Be  → 全局 x ∈ [  123500,   133000]
    {2, 6, 12, 22000.0, 40000.0,  134000},  // d3d4: 12C  → 全局 x ∈ [  156000,   174000]
    {3, 1,  1,  1000.0,  2300.0,  175000},  // d4s:  1H   → 全局 x ∈ [  176000,   177300]
    {3, 2,  4,  4000.0, 10000.0,  179000},  // d4s:  4He  → 全局 x ∈ [  183000,   189000]
    {3, 4,  7, 12000.0, 28000.0,  191000},  // d4s:  7Be  → 全局 x ∈ [  203000,   219000]
};
```

当 `x` 落入某个 `[left+offset, right+offset]` 区间时，PidFitFunc 自动选择对应的 `(layer, charge, mass)` 进行拟合计算。

**offset 布局原则**：每个 (层对, 粒子) 组合分配一段独立的全局 x 区间，区间之间留有 ≥2000 ADC 的间隔（gap），确保各段互不重叠。所有 offset 单调递增。TF1 总范围 0–330000 覆盖所有区间。

> **注**：d1d2 的 ⁶Li 目前被注释掉。经测试，引入 ⁶Li 后拟合完全崩溃（chi2 飙升 4 个数量级），即使 TCutG 不重叠且理论曲线准确也未能消除。根因尚在排查中，疑似与 ⁴He 和 ⁶Li 共享 d1/d2 参数时 Minuit 优化路径的数值不稳定性有关。待手动调整 cut 后重新启用。

#### 5.4.7 ★ left/right 对理论曲线的截断效应

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

**TH2D 展示时（GenerateTheoryCurve）**：绘制完整的理论曲线，不经过任何截断。

```
                 完整 de-e 理论曲线
                 ╱                ╲
                ╱                  ╲
       ΔE      ╱                    ╲                ← TH2D 绘制整条红线
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

| | PidFitFunc（拟合） | TH2D 理论曲线（展示） |
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

## 6. 算法流程

### 6.1 main() 函数流程

```
1. 解析命令行参数
2. 加载配置文件
3. 验证探测器配置存在
4. 在 estimate/ 目录中查找 pre_calibration 文件
5. 读取 g_d1d2/g_d2d3/g_d3d4/g_d4s 四个 TGraph
6. 扫描 src/brill/Cut/ 目录，加载所有符合命名格式的 TCutG
7. 遍历每个 cut，对 TGraph 逐点执行 **TCutG + left/right 双重筛选**，命中点加入 gcali
8. 从加载的 cut 中自动收集粒子列表，构建 PidFitFunc
9. 初始化刻度参数
10. 将 TGraph 与 TF1 进行全局拟合
11. 输出刻度参数到 calibration/t0_{run}.txt
12. 保存 gcali 到 calibration/t0_{trigger}{run}.root
```

### 6.2 初始刻度参数

```cpp
double initial_calibration_parameters[10] = {
    0.0, 0.002,   // t0d1: p0=0, p1=0.002
    0.0, 0.006,   // t0d2: p0=0, p1=0.006
    0.0, 0.006,   // t0d3: p0=0, p1=0.006
    0.0, 0.003,   // t0d4: p0=0, p1=0.003
    0.0, 0.003    // t0s:  p0=0, p1=0.003
};
```

所有参数均设置拟合边界：
- **p0（偏移量）**：`[0.0, 100.0]`（非负约束，避免限制拟合）  
- **p1（增益系数）**：`[0.0, 1.0]`（物理约束：能量/ADC 必须为正）

使用 `"R S"` 选项（TGraph 范围 + Strategy 2）进行拟合。

---

## 7. 输出文件

### 7.1 calibration/t0_{run}.txt

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

### 7.2 calibration/t0_{trigger}{run}.root

例如 `calibration/t0_t1057.root`。

包含拟合用的 TGraph `gcali`，保存所有用于拟合的数据点（x 坐标已加 offset），可用于可视化检查拟合质量。

---

## 8. 与 pre_calibration 的衔接

`calibrate_t0` 现在直接依赖 `pre_calibration` 的输出：

```
match/*.root + ingot/t0s_*.root
        ↓ pre_calibration
estimate/pre_calibration_t1{run}_*.root   (TGraph: g_d1d2, g_d2d3, g_d3d4, g_d4s)
        ↓ calibrate_t0 (TCutG 筛选 + 理论曲线拟合)
calibration/t0_{run}.txt
```

**数据流**：TCutG 替代了原来的 `track_t0` PID 步骤，直接从 pre_calibration 的原始 ADC 散点中筛选目标粒子进行刻度拟合。

---

## 9. 与 normalize 程序的关系

| 程序 | 刻度对象 | 刻度形式 | 输入 |
|------|---------|---------|------|
| `normalize` | 单个 DSSD 的各条 strip | per-strip: `E_norm = p0 + p1 * E_raw` | 原始 DSSD 数据 |
| `calibrate_t0` | 整个探测器层的能量 | per-layer: `E_cal = p0 + p1 * E_raw` | pre_calibration TGraph + TCutG cuts |

**执行顺序**：`normalize`（条间归一化）→ `match`（正背面匹配）→ `pre_calibration`（生成 PID 图 + TGraph）→ （手动绘制 TCutG 截断）→ `calibrate_t0`（层间能量刻度）。