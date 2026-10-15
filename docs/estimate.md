# Estimate 程序汇总

`estimate/` 目录下的程序用于估算、拟合和校准 T0 探测器的各项参数，是数据处理流程中的关键环节。

## 程序列表

| 程序 | 用途 | 输入源 | 输出 |
|------|------|--------|------|
| `estimate_t0_pid` | DSSD 层间 PID 谱 | `match/` (DSSD 匹配后) | 4 个 TH2F PID 直方图 |
| `estimate_t0_center` | 探测器中心偏移估算 | `match/` (DSSD 匹配后) | 控制台输出 offset 值 |
| `estimate_t0_straight` | Straight PID 拟合 | `match/` + `ingot/` (硅) | TH2F + TGraph + TF1 拟合参数 |
| `estimate_t0_csi_pid` | CsI PID 谱 | `ingot/` (硅 + CsI) | 36 个 TH2F PID 直方图 |
| `estimate_normalize` | DSSD 归一化能量 | `ingot/` (DSSD 原始) | 归一化后的 ROOT 文件 |
| `estimate_normalize_total_energy` | DSSD 归一化总能量 | `ingot/` (DSSD 原始) | 归一化总能量 ROOT 文件 |
| `check_beam_width` | 束流能散 + E-TOF/PPAC 关联 | `match/` + `ingot/` + `track/` | 分组 TH1D/TH2D + 拟合结果 |
| `check_tof_and_beam` | TOF-能量相对论校准 | `match/` + `ingot/` | TH1D/TH2D + 理论曲线 |

## 通用命令行参数

| 参数 | 简写 | 类型 | 必需 | 说明 |
|------|------|------|------|------|
| `--help` | `-h` | flag | 否 | 打印帮助信息 |
| `--run` | `-r` | int | 是 | 起始 run 号 |
| `--end-run` | `-e` | int | 否 | 终止 run 号（默认等于起始 run） |
| `--trigger` | `-t` | string | 否 | 触发类型（覆盖 config 中的设置） |
| `--config` | `-c` | string | 否 | 配置文件路径（默认 `config.toml`） |

---

## estimate_t0_pid

DSSD 匹配后估算 T0 PID。通过相邻探测器层之间的能量关联，生成二维 PID 谱。

### 输入

| 文件 | 路径 | 内容 |
|------|------|------|
| `t0d1_<trigger>_<run>.root` | `config.paths.match` | T0 D1 层匹配后数据 |
| `t0d2_<trigger>_<run>.root` | `config.paths.match` | T0 D2 层匹配后数据 |
| `t0d3_<trigger>_<run>.root` | `config.paths.match` | T0 D3 层匹配后数据 |
| `t0d4_<trigger>_<run>.root` | `config.paths.match` | T0 D4 层匹配后数据 |
| `t0s_<trigger>_<run>.root` | `config.paths.ingot` | T0 硅探测器数据 |

### 输出直方图

| 直方图 | 标题 | X 轴 | Y 轴 | 说明 |
|--------|------|------|------|------|
| `d1d2` | D1-D2 PID | D2 energy (0-60MeV) | D1 energy (0-30MeV) | D1-D2 层间能量关联 |
| `d2d3` | D2-D3 PID | D3 energy (0-50MeV) | D2 energy (0-80MeV) | D2-D3 层间能量关联 |
| `d3d4` | D3-D4 PID | D4 energy (0-45MeV) | D3 energy (0-45MeV) | D3-D4 层间能量关联 |
| `d4s` | D4-S PID | D4 energy (0-80MeV) | Si energy (0-45MeV) | D4-S 层间能量关联 |

输出路径：`<workspace>/<config.paths.estimate>/t0_pid_<trigger>_<run>_<end_run>.root`

### 筛选逻辑

- **D1-D2, D2-D3, D3-D4**：相邻层 hit 对的 x/y 差值在 track window 范围内才填入
- **D4-S**：仅当 D4 只有 1 个 hit 且 S 层 `valid == true` 时才填入

### 使用示例

```bash
./estimate_t0_pid -r 100 -t dt
./estimate_t0_pid -r 100 -e 200 -t dt
```

---

## estimate_t0_center

估算 T0 各探测器层的中心偏移量。通过相邻层之间 hit 位置的差值分布，拟合高斯函数得到 offset。

### 输入

4 层 DSSD 匹配后数据（`t0d1` ~ `t0d4`，来自 `config.paths.match`）。

### 工作原理

1. 对相邻层（D1-D2, D2-D3, D3-D4）的 hit 对计算位置差 `dx = right_x - left_x`, `dy = right_y - left_y`
2. 使用探测器 pitch 进行随机化处理（避免离散化效应）
3. 对差值直方图进行高斯拟合，输出 mean 和 sigma
4. 可选地，对 strip 级别（未随机化）也进行拟合

### 输出

控制台输出每对相邻层的 x/y 偏移量的高斯拟合结果（mean ± sigma）。

### 使用示例

```bash
./estimate_t0_center -r 100 -t dt
```

---

## estimate_t0_straight

估算 T0 Straight PID。利用 PID 拟合函数 `(0.5/A) * (sqrt(x² + 4A(Bx-C)²) - x)` 对相邻层能量关联数据进行拟合，得到各粒子（4He, 6He, 3H 等）的 A、B、C 参数。

### 输入

- DSSD 匹配后数据（`t0d1` ~ `t0d4`，来自 `config.paths.match`）
- 硅探测器数据（`t0s`，来自 `config.paths.ingot`）

### 工作流程

1. 对每个 slice（D1-D2, D2-D3, D3-D4, D4-S）生成 PID 和 Straight PID 直方图
2. 通过 track window 筛选相邻层 hit 对
3. 使用 `TCutG` 图形切割选择粒子区域
4. 对每个粒子的 `TGraph` 进行 PID 拟合，输出 A、B、C 参数
5. 使用 `4He` 的拟合结果作为默认参数

### 输出

- 控制台输出各 slice 各粒子的拟合参数
- ROOT 文件包含 TH2F、TGraph 和 TF1 对象

### 使用示例

```bash
./estimate_t0_straight -r 100 -t dt
```

---

## estimate_t0_csi_pid

估算 T0 CsI PID。生成硅探测器与 36 个 CsI 晶体之间的二维能量关联直方图。

### 输入

- 硅探测器数据（`t0s`，来自 `config.paths.ingot`）
- CsI 数据（`t0csi`，来自 `config.paths.ingot`）

### 输出

36 个 TH2F 直方图（`h0` ~ `h35`），每个对应一个 CsI 晶体，X 轴为 CsI 能量，Y 轴为硅能量。

输出路径：`<workspace>/<config.paths.estimate>/t0csi_pid_<trigger>_<run>_<end_run>.root`

### 筛选逻辑

仅当 `silicon.valid == true` 且对应 CsI 晶体 `valid[i] == true` 时才填入。

### 使用示例

```bash
./estimate_t0_csi_pid -r 100 -t dt
```

---

## estimate_normalize

对 DSSD 原始数据进行归一化处理。使用预先计算好的归一化参数，对 front/back 面各 strip 的能量进行校正。

### 输入

- DSSD 原始数据（来自 `config.paths.ingot`）
- 归一化参数文件（来自 `config.paths.normalize`）

### 位置参数

```bash
./estimate_normalize -r 100 -t dt t0d1 t0d2 t0d3 t0d4
```

支持同时处理多个探测器（`t0d1` `t0d2` `t0d3` `t0d4`）。

### 输出

归一化后的 ROOT 文件，包含 front/back 归一化能量，以及 front-back 能量差直方图。

输出路径：`<workspace>/<config.paths.estimate>/<detector>_normalize_<trigger>_<run>.root`

---

## estimate_normalize_total_energy

与 `estimate_normalize` 类似，但额外计算归一化后的 front/back 总能量。

### 输入

与 `estimate_normalize` 相同。

### 输出

归一化后的 ROOT 文件，包含每事件的前面总能量（`fe_t`）、背面总能量（`be_t`）、前/后面 hit 数（`num_f`/`num_b`），以及 front-back 能量差直方图。

输出路径：`<workspace>/<config.paths.estimate>/<detector>_normalize_total_energy_<trigger>_<run>.root`

---

## 数据处理流程

```
ingot/ (原始数据)
    │
    ├── estimate_normalize ──────────→ estimate/ (归一化能量)
    ├── estimate_normalize_total_energy → estimate/ (归一化总能量)
    │
match/ (DSSD匹配后)
    │
    ├── estimate_t0_pid ─────────────→ estimate/ (PID 直方图)
    ├── estimate_t0_center ──────────→ (控制台输出 offset)
    ├── estimate_t0_straight ────────→ estimate/ (Straight PID 拟合)
    │
ingot/ (硅 + CsI)
    │
    └── estimate_t0_csi_pid ─────────→ estimate/ (CsI PID 直方图)
```

---

## check_beam_width

束流能散检查和 E-TOF/PPAC2 关联分析程序。

### 功能概述

1. 对 d2/d3 能量进行刻度，筛选 14O 束流事件
2. 按每 `10^n` 个事件分组（最多 10 组），每组绘制一维束流能量谱
3. 每组内对能谱进行 Gaussian 拟合，获得峰值 E0 和 FWHM
4. **单 run 模式额外功能**：
   - 读取 ingot/beam 的 TOF 数据，绘制每组 E vs TOF 二维图
   - 对每组 E-vs-TOF 数据进行 robust 线性拟合 `tof = k·E + b`，在 ±1.5σ 能量范围内拟合
   - 读取 track/ppac 数据，绘制 PPAC2 的 x/y 位置、dir_x/dir_y 与 E 的二维关联图（仅 group 1）

### 命令行参数

| 参数 | 简写 | 类型 | 必需 | 说明 |
|------|------|------|------|------|
| `--help` | `-h` | flag | 否 | 打印帮助信息 |
| `--config` | `-c` | string | 否 | 配置文件路径（默认 `config.toml`） |
| `--run` | `-r` | int | 是 | run 号 |
| `--end-run` | `-e` | int | 否 | 终止 run 号，多 run 模式 |
| `--trigger` | `-t` | string | 否 | 触发类型 |
| `--group-power` | `-n` | int | 否 | 每组 `10^n` 事件（默认 6） |

### 输入

| 文件 | 路径 | 内容 |
|------|------|------|
| `t0d2_<trigger>_<run>.root` | `config.paths.match` | T0 D2 层匹配后数据 |
| `t0d3_<trigger>_<run>.root` | `config.paths.match` | T0 D3 层匹配后数据 |
| `beam_<trigger>_<run>.root` | `config.paths.ingot` | ingot 束流 TOF 数据 |
| `beam_<trigger>_<run>.root` | `config.paths.beam` | sorted 束流 PID 识别结果 |
| `ppac_<trigger>_<run>.root` | `config.paths.track` | PPAC track 数据（可选） |
| `t0_<calib_run>.txt` | `config.paths.calibration` | T0 能量刻度参数 |

### 筛选逻辑

- d2 和 d3 各只有 1 个 hit
- d2 和 d3 hit 位置差 `dx² + dy² < 4`
- 刻度后的 `(d3, d2)` 能量落入 14O 束流 cut（`cal_d2_d3_14O_beam_cut`）
- TOF/E-vs-TOF 填入需 `beam_event.valid && v_14O`

### 输出

输出路径：`<workspace>/<config.paths.beam>/check_beam_width_<trigger>_<run>.root`

| 对象 | 说明 |
|------|------|
| `h_beam_E_group{N}` | 每组的一维束流能量谱 |
| `h_E_vs_tof_group{N}` | 每组的 E vs TOF 二维图 |
| `h_ppac2_x_vs_E` | PPAC2 x 位置 vs E 二维图（group 1） |
| `h_ppac2_y_vs_E` | PPAC2 y 位置 vs E 二维图（group 1） |
| `h_dir_x_vs_E` | Track dir_x (斜率, 1000 bins, ±0.02 rad) vs E（group 1） |
| `h_dir_y_vs_E` | Track dir_y (斜率, 1000 bins, ±0.02 rad) vs E（group 1） |
| `c_beam_width` | 分组能谱叠画 TCanvas |
| `c_E_vs_tof_g{N}` | 每组 E-vs-TOF + 线性拟合红线 TCanvas |
| `g_peak` / `g_fwhm` | Gaussian 拟合峰值/FWHM 随组变化的 TGraph |

### 使用示例

```bash
# 单 run 模式
./check_beam_width -r 57 -t main

# 多 run 模式（不含 TOF/PPAC 功能）
./check_beam_width -r 57 -e 67 -t main

# 自定义每组 10^5 事件
./check_beam_width -r 57 -n 5
```

---

## check_tof_and_beam

束流 TOF 与能量的相对论运动学校准程序。

### 功能概述

1. 对 d1/d2/d3 能量进行刻度，筛选 14O 束流事件
2. 收集最多 **1e6** 个满足条件的事件
3. 对束流能量一维谱进行 Gaussian 拟合，得到中心值 **$E_0$**
4. 对 TOF 一维谱进行 Gaussian 拟合，得到中心值 **$t_0$**
5. 利用 14O 的相对论运动学关系，代入 $E_0$ 和 $t_0$ 解得 TOF 测量偏置 **$\Delta t$**
6. 绘制 E vs TOF 二维关联图，并将理论曲线 $E = f(t - \Delta t)$ 作为红线叠加

### 物理模型与计算方法

#### 已知常量

| 符号 | 值 | 含义 |
|------|-----|------|
| $m$ | $13044.84\ \mathrm{MeV}/c^2$ | 14O 核静止质量 |
| $L$ | $17000\ \mathrm{mm}$（可通过 `--distance` 指定） | TOF 飞行距离 |
| $c$ | $299.792458\ \mathrm{mm/ns}$ | 真空中光速 |

#### Step 1: 能量刻度与能损修正

读取 `<workspace>/<calibration>/t0_<calib_run>.txt` 中的刻度参数 $(p_0, p_1)$，对 d1 (layer=0)、d2 (layer=1)、d3 (layer=2) 的原始能量进行线性刻度：

$$
E_{\text{cal}, \ell} = p_0^{(\ell)} + p_1^{(\ell)} \cdot E_{\text{raw}, \ell}, \quad \ell \in \{0, 1, 2\}
$$

d1+d2+d3 沉积能量 $E_{\text{meas}} = E_{\text{cal},0} + E_{\text{cal},1} + E_{\text{cal},2}$。

束流实际路径为：**CD₂靶 (121.7μm)** → **d1 (68μm Si)** → d2 → d3。d1、d2、d3 三层均已实测能量，因此 $E_{\text{meas}}$ 近似为 d1 入射面（即 CD₂-Si 界面）处的动能。只需补全 CD₂ 靶中的能损，即可得到入射束流动能 $E$。

使用 catima 库构建 CD₂ 材料的 `RangeEnergyCalculator`（14O, Z=8, A=14）：

$$
R_{\text{CD}_2}(E_{\text{meas}}) \xrightarrow{+121.7\ \mu\text{m}} R_{\text{target}} \xrightarrow{\text{Energy}} E
$$

其中 CD₂ 材料参数：(C: Z=6, A=12, stoi=1) + (D: Z=1, A=2, stoi=2)，密度 $0.79867\ \mathrm{g/cm^3}$。

| 介质 | 厚度 (μm) | 材料 |
|------|-----------|------|
| CD₂ 靶 | 121.7 | $(^{12}\text{C}^2\text{H}_2)_n$ |

---

#### Step 2: 一维谱拟合

**束流能量**：先寻找 $h_{\text{beam\_E}}$ 的最大值 bin，设其中心为 $x_E^{\text{peak}}$，然后在 $[x_E^{\text{peak}}-6,\ x_E^{\text{peak}}+6]\ \mathrm{MeV}$ 范围内做 Gaussian 拟合：

$$
f_E(x) = A_E \exp\left(-\frac{(x - E_0)^2}{2\sigma_E^2}\right)
$$

得到 $E_0, \sigma_E$。

**TOF**：先寻找 $h_{\text{tof}}$ 的最大值 bin，设其中心为 $x_t^{\text{peak}}$，然后在 $[x_t^{\text{peak}}-1,\ x_t^{\text{peak}}+1]\ \mathrm{ns}$ 范围内做 Gaussian 拟合：

$$
f_t(x) = A_t \exp\left(-\frac{(x - t_0)^2}{2\sigma_t^2}\right)
$$

得到 $t_0, \sigma_t$。

---

#### Step 3: 相对论运动学

14O 动能 $E_k$ 与速度 $\beta = v/c$ 的关系（已确认非相对论效应不可忽略）：

$$
E_k = mc^2(\gamma - 1), \quad \gamma = \frac{1}{\sqrt{1 - \beta^2}}
$$

TOF 与真实速度的关系（$\Delta t$ 为测量偏置）：

$$
t_{\text{meas}} = t_{\text{true}} + \Delta t, \quad t_{\text{true}} = \frac{L}{\beta c}
$$

---

#### Step 4: 求解偏置 $\Delta t$

将 $E_0$ 代入运动学关系反推真实速度 $\beta_0$：

$$
\gamma_0 = 1 + \frac{E_0}{mc^2}
$$

$$
\beta_0 = \sqrt{1 - \frac{1}{\gamma_0^2}} = \sqrt{1 - \left(\frac{mc^2}{E_0 + mc^2}\right)^2}
$$

真实 TOF：

$$
t_{\text{true}} = \frac{L}{\beta_0 c}
$$

测量偏置：

$$
\Delta t = t_0 - t_{\text{true}} = t_0 - \frac{L}{c\beta_0}
$$

---

#### Step 5: 理论曲线

校准后的 TOF-能量关系（将 $\Delta t$ 代回）：

$$
E(t_{\text{meas}}) = mc^2\left(
\frac{1}{\sqrt{1 - \left(\dfrac{L}{(t_{\text{meas}} - \Delta t) \cdot c}\right)^2}} - 1
\right)
$$

该函数作为红色曲线绘制在 E vs TOF 二维图上。

---

#### Step 6: 门控分析

为研究 $\Delta E$ 和 $\Delta t$ 的来源，进行互斥门控：

**TOF-gated E**：以 $t_0 \pm 0.25\sigma_t$ 筛选事件，绘制 E 谱并 Gaussian 拟合，得到门控后 $E_0^{\text{gated}}$ 和 $\sigma_E^{\text{gated}}$。

**E-gated TOF**：以 $E_0 \pm 0.25\sigma_E$ 筛选事件，绘制 TOF 谱并 Gaussian 拟合，得到门控后 $\sigma_t^{\text{gated}}$。

---

#### Step 7: 参数汇总

| 参数 | 含义 | 来源 |
|------|------|------|
| $E_0$ | 束流动能中心值 | $h_{\text{beam\_E}}$ 的 Gaussian 拟合 |
| $\sigma_E$ | 能谱宽度 | 同上 |
| $t_0$ | TOF 测量中心值 | $h_{\text{tof}}$ 的 Gaussian 拟合 |
| $\sigma_t$ | TOF 宽度 | 同上 |
| $\beta_0$ | 由 $E_0$ 反推的真实速度 | $\beta_0 = \sqrt{1 - (mc^2/(E_0+mc^2))^2}$ |
| $t_{\text{true}}$ | 理论真实 TOF | $L / (\beta_0 c)$ |
| $\Delta t$ | TOF 测量偏置 | $t_0 - t_{\text{true}}$ |

### 命令行参数

| 参数 | 简写 | 类型 | 必需 | 说明 |
|------|------|------|------|------|
| `--help` | `-h` | flag | 否 | 打印帮助信息 |
| `--config` | `-c` | string | 否 | 配置文件路径（默认 `config.toml`） |
| `--run` | `-r` | int | 是 | run 号 |
| `--distance` | `-L` | double | 否 | TOF 飞行距离 (mm)，默认 17000 |

### 输入

| 文件 | 路径 | 内容 |
|------|------|------|
| `t0d1_<_run>.root` | `config.paths.match` | T0 D1 层匹配后数据 |
| `t0d2_<_run>.root` | `config.paths.match` | T0 D2 层匹配后数据 |
| `t0d3_<_run>.root` | `config.paths.match` | T0 D3 层匹配后数据 |
| `beam_<_run>.root` | `config.paths.ingot` | ingot 束流 TOF 数据 |
| `beam_<_run>.root` | `config.paths.beam` | sorted 束流 PID 识别（14O/13N/12C） |
| `t0_<calib_run>.txt` | `config.paths.calibration` | T0 能量刻度参数 |

### 筛选逻辑

- d1、d2、d3 各只有 1 个 hit（`num == 1`）
- d3 与 d1 hit 位置差满足 $dx^2 + dy^2 < 4$（欧氏距离 $< 2\ \mathrm{mm}$）
- 刻度后的 d2、d3 能量落入 `cal_d2_d3_14O_beam_cut` TCutG 切割区域内
- `beam_event.valid == true` 且粒子为 14O（`v_14O == true`）
- 最多收集 1e6 个事件

### 输出

输出路径：`<workspace>/<config.paths.beam>/check_tof_and_beam_<run>.root`

| 对象 | 说明 |
|------|------|
| `h_beam_E` | 束流能量一维谱 + Gaussian 拟合 |
| `h_tof` | TOF 一维谱 + Gaussian 拟合 |
| `h_E_vs_tof` | E vs TOF 二维关联图 + 理论红线 |
| `h_beam_E_tof_gated` | TOF 门控（±0.25σ）下的 E 一维谱 + Gaussian 拟合 |
| `h_tof_E_gated` | E 门控（±0.25σ）下的 TOF 一维谱 + Gaussian 拟合 |
| `c_E` | 能谱 TCanvas |
| `c_tof` | TOF TCanvas |
| `c_2d` | E-vs-TOF + 理论曲线 TCanvas，标注 $E_0$、$t_0$、$\Delta t$ |
| `c_E_tof_gated` | TOF 门控 E 谱 TCanvas |
| `c_tof_E_gated` | E 门控 TOF 谱 TCanvas |

输出参数文件：`<workspace>/<config.paths.beam>/check_tof_and_beam_<run>.txt`

包含五个空格分隔的浮点数（每行一份，首行为注释头）：

| 列 | 含义 |
|----|------|
| $E_0$ | 束流动能中心值 (MeV) |
| $\sigma_E$ | 能量宽度 (MeV) |
| $t_0$ | TOF 中心值 (ns) |
| $\sigma_t$ | TOF 宽度 (ns) |
| $\Delta t$ | TOF 测量偏置 (ns) |

可通过 `dt` 和实测 `tof` 反推能量：$E = mc^2(\gamma - 1)$，其中 $\beta = L / ((t_{\text{meas}} - \Delta t) \cdot c)$。

### 使用示例

```bash
# 使用默认飞行距离 17m
./check_tof_and_beam -r 57

# 自定义飞行距离
./check_tof_and_beam -r 57 -L 16000
```

---

## 编译

所有 estimate 程序已包含在 CMake 构建中：

```bash
cmake --build build --target estimate_t0_pid
cmake --build build --target estimate_t0_center
cmake --build build --target estimate_t0_straight
cmake --build build --target estimate_t0_csi_pid
cmake --build build --target estimate_normalize
cmake --build build --target estimate_normalize_total_energy
cmake --build build --target check_beam_width
cmake --build build --target check_tof_and_beam
```