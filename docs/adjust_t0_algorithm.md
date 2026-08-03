# T0 探测器位置矫正

## 概述

`adjust_t0_step1` 和 `adjust_t0_step2` 两个程序用于矫正 T0 探测器（D1-D4）的物理位置偏移。由于机械安装误差，探测器的实际位置与设计位置存在偏差，需要通过数据分析确定并修正。

### Step 1 vs Step 2

| | Step 1 | Step 2 |
|---|---|---|
| **目标** | D1, D2, D3 | D4 |
| **参考** | PPAC 径迹（束流） | D1-D2-D3 拟合直线 |
| **输入** | track 文件 + match 文件 | match 文件 |
| **适用** | 束流事件（main trigger） | 任意事件（反应事件可用） |
| **原因** | 束流方向不变，PPAC 径迹为准 | D4 无束流穿过，需用 D1-D3 外推 |

---

## Step 1: 束流矫正 D1/D2/D3

## 物理原理

### 束流方向不变假设

在束流实验中，束流方向几乎不会改变，所有束流粒子沿同一直线穿过探测器系统。PPAC 径迹重建得到的 `target_x`、`target_y`（z=0 处的位置）和 `dir_x`、`dir_y`（方向余弦）准确描述了束流粒子的运动方向。

### 残差计算

对于每个探测器 $i$（位于 $z_i$ 处），PPAC 径迹外推的理论位置为：

$$\text{predicted}_x = \text{target}_x + \text{dir}_x \cdot z_i$$
$$\text{predicted}_y = \text{target}_y + \text{dir}_y \cdot z_i$$

match 文件中的实际测量位置为 `match_events[i].x[0]`（加上 `config.toml` 中的 `x_offset_mm`）。残差定义为：

$$\text{residual}_x = \text{match}_x - \text{predicted}_x$$
$$\text{residual}_y = \text{match}_y - \text{predicted}_y$$

残差的平均值即为探测器相对于 PPAC 径迹标准的位置偏差，取反后得到应施加的 offset 矫正值。

## 算法流程

```
adjust_t0_step1
  │
  ├─ 1. 读取 PPAC track 文件 (ppac_XXXX.root)
  │     └─ [adjust_t0.cpp:L27-L34](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/event/t0/adjust_t0.cpp#L27-L34)
  │
  ├─ 2. 读取 T0D1/D2/D3 的 match 文件 (t0d1_XXXX.root, ...)
  │     └─ [adjust_t0.cpp:L40-L55](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/event/t0/adjust_t0.cpp#L40-L55)
  │
  ├─ 3. 逐事件筛选（总共 ~26M 事件）
  │     │
  │     ├─ 3a. PPAC track valid?    → ~90% 通过
  │     ├─ 3b. T0D1/D2/D3 同时 num==1?  → ~51% 通过
  │     └─ 3c. D2/D3 能量截断?     → ~20% 通过（取决于截断窗口）
  │     └─ [adjust_t0.cpp:L78-L114](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/event/t0/adjust_t0.cpp#L78-L114)
  │
  ├─ 4. 填充残差直方图（200 bins, -10~+10 mm）
  │     └─ [adjust_t0.cpp:L101-L112](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/event/t0/adjust_t0.cpp#L101-L112)
  │
  ├─ 5. 高斯拟合残差谱 → 峰值 = 偏移量
  │     └─ [adjust_t0.cpp:L147-L188](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/event/t0/adjust_t0.cpp#L147-L188)
  │
  └─ 6. 输出结果 & 保存 ROOT 文件
        └─ [adjust_t0.cpp:L190-L208](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/event/t0/adjust_t0.cpp#L190-L208)
```

## 数据结构

### 输入配置: T0AdjustConfig

[adjust_t0.h:L29-L37](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/event/t0/adjust_t0.h#L29-L37)

| 字段 | 类型 | 说明 |
|------|------|------|
| `track_path` | `string` | PPAC track 文件路径 |
| `match_dir` | `string` | match 文件目录 |
| `trigger` | `string` | 触发类型（`main` 或 `t1`） |
| `run` | `int` | Run 号 |
| `output_dir` | `string` | 输出目录（normalize/） |
| `detector_names` | `vector<string>` | 探测器名称列表 |
| `energy_cuts` | `vector<T0EnergyCut>` | 各探测器的能量截断 |

### 能量截断: T0EnergyCut

[adjust_t0.h:L23-L26](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/event/t0/adjust_t0.h#L23-L26)

| 字段 | 类型 | 说明 |
|------|------|------|
| `min` | `double` | 能量下限（ADC 通道） |
| `max` | `double` | 能量上限（ADC 通道），max=0 表示不截断 |

### 输出结果: T0AdjustResult

[adjust_t0.h:L12-L22](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/event/t0/adjust_t0.h#L12-L22)

| 字段 | 类型 | 说明 |
|------|------|------|
| `detector_name` | `string` | 探测器名称 |
| `dx`, `dy` | `double` | 偏移量（mm），`dx = -residual_mean` |
| `dsigma_x`, `dsigma_y` | `double` | 偏移量不确定度，<0 表示拟合失败 |
| `residual_x_mean` | `double` | 残差谱高斯拟合均值（mm） |
| `residual_x_sigma` | `double` | 残差谱高斯拟合宽度（mm） |
| `num_events` | `int` | 通过筛选的事件数 |

### 输入数据: PpacTrackEvent

[ppac_track_event.h:L12-L22](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/event/ppac/ppac_track_event.h#L12-L22)

| 字段 | 类型 | 说明 |
|------|------|------|
| `valid` | `int` | 径迹是否有效 |
| `target_x`, `target_y` | `double` | 束流在 z=0 处的位置（mm） |
| `dir_x`, `dir_y` | `double` | 束流方向（mm/mm） |
| `ppac_x[3]`, `ppac_y[3]` | `double` | 各 PPAC 上的位置 |

### 输入数据: DssdMatchEvent

[dssd_match_event.h:L8-L20](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/event/t0/dssd_match_event.h#L8-L20)

| 字段 | 类型 | 说明 |
|------|------|------|
| `num` | `int` | 匹配出的粒子数 |
| `energy[8]` | `double` | 粒子能量（ADC） |
| `x[8]`, `y[8]` | `double` | 粒子物理坐标（mm） |

## 事件筛选条件

### 条件 1: PPAC 径迹有效

`track_event.valid == 1`，即三个 PPAC 的 X 和 Y 方向均有有效 hit，且线性拟合成功。

### 条件 2: 三个探测器同时 num==1

`match_events[i].num == 1`（i = 0,1,2），即 T0D1、D2、D3 上均恰好匹配出一个粒子。这确保了束流粒子穿过所有三个探测器，且没有碎片或噪声干扰。

### 条件 3: D2/D3 能量截断

束流粒子在 D2 和 D3 中的能量沉积应在特定范围内。当前配置（[adjust_t0_step1.cpp:L80-L82](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/adjust_t0_step1.cpp#L80-L82)）：

| 探测器 | 能量范围 (ADC) | 说明 |
|--------|---------------|------|
| T0D1 | 无截断 | 薄探测器（68μm），能量沉积小 |
| T0D2 | 43000 - 47000 | 厚探测器（1005μm），束流粒子能量沉积集中 |
| T0D3 | 30000 - 35000 | 厚探测器（995μm） |

选择 D2/D3 的原因是它们厚度大（~1000μm），能量沉积足够大且稳定，适合用能量截断挑选束流粒子。D1 厚度仅 68μm，能量沉积太小，不做截断。

## 测试结果

| 探测器 | Run | 总事件 | 选中事件 | X 残差均值 | X σ | Y 残差均值 | Y σ | X offset | Y offset |
|--------|-----|--------|---------|-----------|-----|-----------|-----|----------|----------|
| T0D1 | 57 | 27.5M | 5.87M | -0.0945 | 0.711 | -1.961 | 0.627 | +0.095 | +1.961 |
| T0D1 | 60 | 26.4M | 2.77M | -0.104 | 0.713 | -1.947 | 0.625 | +0.104 | +1.947 |
| T0D1 | 64 | 15.3M | 0.93M | -0.099 | 0.711 | -1.977 | 0.632 | +0.099 | +1.977 |
| T0D2 | 57 | 27.5M | 5.87M | -1.003 | 0.476 | -2.406 | 0.370 | +1.003 | +2.406 |
| T0D2 | 60 | 26.4M | 2.77M | -1.007 | 0.482 | -2.395 | 0.369 | +1.007 | +2.395 |
| T0D2 | 64 | 15.3M | 0.93M | -1.004 | 0.476 | -2.416 | 0.373 | +1.004 | +2.416 |
| T0D3 | 57 | 27.5M | 5.87M | -0.103 | 0.728 | -1.918 | 0.622 | +0.103 | +1.918 |
| T0D3 | 60 | 26.4M | 2.77M | -0.109 | 0.733 | -1.929 | 0.623 | +0.109 | +1.929 |
| T0D3 | 64 | 15.3M | 0.93M | -0.108 | 0.731 | -1.930 | 0.627 | +0.108 | +1.930 |

（offset = -残差均值，单位为 mm）

### 结果解读

1. **Y 方向偏移显著**：所有探测器 Y 方向偏移约为 1.9-2.4 mm，说明探测器整体在 Y 方向有系统性的机械安装偏差。D2 的 Y 偏移最大（~2.4 mm），可能与 D2 的安装位置有关。

2. **X 方向偏移较小**：D1 和 D3 的 X 偏移约 0.1 mm，D2 的 X 偏移约 1.0 mm，说明 D2 在 X 方向也存在一定的偏移。

3. **跨 run 一致性极好**：run 57 和 run 60 的结果几乎一致，X 方向差异 < 0.01 mm，Y 方向差异 < 0.02 mm，表明结果具有极好的可重复性。

4. **D2 残差宽度最小**：D2 的残差 σ 约为 0.47-0.48 mm（X）和 0.37 mm（Y），明显小于 D1 和 D3（~0.7 mm）。这是因为 D2 厚度最大（1005μm），位置分辨率最高。

5. **能量截断效率**：run 57 约 5.9M 事件通过能量截断（~39% of num==1），说明能量窗口设置合理。之前的 388k 结果是因为能量截断窗口不同（更窄），新窗口（D2: 43000-47000, D3: 30000-35000）更合理。

## 程序用法

```bash
./build/bin/adjust_t0_step1 -r <run> [-t <trigger>] [-c <config>]
```

| 选项 | 说明 | 默认值 |
|------|------|--------|
| `-r, --run` | Run 号 | 必填 |
| `-t, --trigger` | 触发类型：`main`（束流）或 `t1`（反应） | `t1` |
| `-c, --config` | 配置文件路径 | `config.toml` |

**注意**：束流事件应使用 `main` trigger（非 `t1`），因为束流方向不变的前提只在主触发下成立。

输出文件保存至 `normalize/t0_offset_<trigger>_<run>.root`，包含：
- `h_res_x_<det>` / `h_res_y_<det>`：残差直方图
- `h_energy_<det>`：能量直方图（用于调试能量截断）

## 后续步骤

得到偏移量后，将其填入 `config.toml` 中对应探测器的 `x_offset_mm` / `y_offset_mm` 字段：

```toml
[detectors.t0d1]
x_offset_mm = 0.10   # 取 run57/60/64 平均值 (0.099)
y_offset_mm = 1.96   # 取 run57/60/64 平均值 (1.962)

[detectors.t0d2]
x_offset_mm = 1.00   # 取 run57/60/64 平均值 (1.005)
y_offset_mm = 2.41   # 取 run57/60/64 平均值 (2.406)

[detectors.t0d3]
x_offset_mm = 0.11   # 取 run57/60/64 平均值 (0.107)
y_offset_mm = 1.93   # 取 run57/60/64 平均值 (1.926)
```

这些 offset 会被所有使用 `SquareDetectorConfig` 的程序自动读入，在计算物理坐标时应用。

---

## Step 2: 径迹外推矫正 D4

### 物理原理

T0D4 位于探测器阵列最下游，束流事件中束流粒子通常不会穿过 D4（被 D3 后的阻挡物挡住），因此无法用 Step 1 的束流方法矫正 D4。

Step 2 利用反应事件（任意 trigger），对同时穿过 D1-D4 的事件：
1. 通过 D1/D2/D3 的位置拟合直线，作为粒子的真实径迹
2. 将径迹外推到 D4 的 z 位置，得到预测位置
3. 残差 = D4 实际位置 - 预测位置
4. 残差的平均值即为 D4 的偏移量

### 算法流程

```
adjust_t0_step2
  │
  ├─ 1. 读取 T0D1/D2/D3/D4 的 match 文件
  │     └─ [adjust_t0.cpp:L227-L244]
  │
  ├─ 2. 逐事件筛选，取前 10000 个 num==1 事件
  │     └─ [adjust_t0.cpp:L258-L264]
  │
  ├─ 3. 对每个事件，用 D1/D2/D3 拟合直线
  │     │  x = a + b·z  (最小二乘，3 点)
  │     │  y = a + b·z
  │     └─ [adjust_t0.cpp:L266-L288]
  │
  ├─ 4. 外推到 D4 的 z 位置 → 预测 (x, y)
  │     └─ [adjust_t0.cpp:L276, L288]
  │
  ├─ 5. 残差 = D4 位置 - 预测位置，填入直方图
  │     └─ [adjust_t0.cpp:L277, L289]
  │
  ├─ 6. 高斯拟合残差谱 → 峰值 = 偏移量
  │     └─ [adjust_t0.cpp:L303-L338]
  │
  └─ 7. 输出结果 & 保存 ROOT 文件
        └─ [adjust_t0.cpp:L342-L365]
```

### 线性拟合公式

对 D1/D2/D3 三点 $(z_1, x_1), (z_2, x_2), (z_3, x_3)$ 做最小二乘：

$$\bar{z} = \frac{z_1+z_2+z_3}{3}, \quad \bar{x} = \frac{x_1+x_2+x_3}{3}$$

$$b = \frac{\sum (z_i - \bar{z})(x_i - \bar{x})}{\sum (z_i - \bar{z})^2}, \quad a = \bar{x} - b\bar{z}$$

D4 预测位置：$x_4^{\text{pred}} = a + b \cdot z_4$

残差：$\Delta x = x_4^{\text{meas}} - x_4^{\text{pred}}$

### 程序用法

```bash
./build/bin/adjust_t0_step2 -r <run> [-t <trigger>] [-c <config>]
```

| 选项 | 说明 | 默认值 |
|------|------|--------|
| `-r, --run` | Run 号 | 必填 |
| `-t, --trigger` | 触发类型 | `t1` |
| `-c, --config` | 配置文件路径 | `config.toml` |

输出文件：`normalize/t0_offset_d4_<trigger>_<run>.root`，包含 `h_res_x_t0d4` / `h_res_y_t0d4` 直方图。

**注意**：Step 2 可以用于 `main` 或 `t1` trigger，因为 D4 上没有束流，只需 D1-D4 同时有 hit 即可。

---

## 源文件

| 文件 | 说明 |
|------|------|
| [adjust_t0_step1.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/adjust_t0_step1.cpp) | Step 1 主程序，参数解析与配置 |
| [adjust_t0_step2.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/bin/adjust_t0_step2.cpp) | Step 2 主程序，参数解析与配置 |
| [adjust_t0.h](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/include/event/t0/adjust_t0.h) | 数据结构定义与函数声明 |
| [adjust_t0.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/event/t0/adjust_t0.cpp) | Step 1 & Step 2 算法实现 |