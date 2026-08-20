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

## 编译

所有 estimate 程序已包含在 CMake 构建中：

```bash
cmake --build build --target estimate_t0_pid
cmake --build build --target estimate_t0_center
cmake --build build --target estimate_t0_straight
cmake --build build --target estimate_t0_csi_pid
cmake --build build --target estimate_normalize
cmake --build build --target estimate_normalize_total_energy
```