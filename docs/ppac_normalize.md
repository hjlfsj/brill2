# PPAC Normalize — PPAC 时间偏移校准

## 概述

`ppac_normalize` 用于校准 PPAC（Parallel Plate Avalanche Counter）探测器的时间偏移（offset），
计算延迟时间（delay time）和芯片峰位置与时间的线性刻度参数。适用于 **RIBLL2** 实验数据。

---

## 算法流程

### 1. 通道映射

每个 PPAC（共 3 个）有 5 个通道，定义在 `ppac_event.h` 的 `PpacChannelMap` 中：

| PPAC 序号 | X1 通道 | X2 通道 | Y1 通道 | Y2 通道 | 阳极通道 |
|-----------|---------|---------|---------|---------|----------|
| PPAC1     | 0       | 1       | 2       | 3       | 12       |
| PPAC2     | 4       | 5       | 6       | 7       | 13       |
| PPAC3     | 8       | 9       | 10      | 11      | 14       |

### 2. 两遍扫描（Two-Pass）

对每个 PPAC 独立执行以下流程：

#### Pass1：延迟时间谱

遍历所有事件，对每个事件检查 `valid` 标志和时间 > 0：

- **阳极效率统计**：阳极有效事件数 / 总事件数
- **阴极效率统计**：阳极有效且 X1/X2/Y1/Y2 均有效的事件数 / 阳极有效事件数
- 填充延迟时间直方图：

```
delay_x = time[X1] + time[X2] - 2 * time[anode]
delay_y = time[Y1] + time[Y2] - 2 * time[anode]
```

#### 延迟峰拟合

对 delay_x 和 delay_y 直方图分别做：
1. `TSpectrum::Search(sigma=0.7, threshold=0.9)` 寻峰
2. 以寻峰位置 ±1.2 ns 范围做高斯拟合 `gaus`
3. 记录 peak（均值）和 sigma（标准差）

#### Pass2：时间差谱

再次遍历所有事件，用 Pass1 得到的 3σ 截断：

```
|delay_x - peak_x| < 3 * sigma_x  →  fill dtx = time[X1] - time[X2]
|delay_y - peak_y| < 3 * sigma_y  →  fill dty = time[Y1] - time[Y2]
```

### 3. 芯片峰寻找与拟合

对 dtx/dty 直方图分别执行：

#### 寻峰

- X 方向：`TSpectrum::Search(sigma=2, threshold=0.2)`
- Y 方向：`TSpectrum::Search(sigma=3, threshold=0.001)`

#### 连丝检测（Wire Pair）

找到所有峰中**间距最小**的两个相邻峰，作为"连丝"：

```
gap = peak[i] - peak[i-1]
pair_idx = argmin(gap)
```

#### 8 个芯片峰选择

```
left_count = pair_idx - 1  （连丝左侧的峰数）

if left_count >= 10:
    start_idx = pair_idx - 9    → 取连丝左侧 8 个峰
else:
    start_idx = pair_idx + 1    → 取连丝右侧 8 个峰
```

对选出的 8 个峰，每个在 ±0.6 ns 范围内做高斯拟合，取拟合中心值作为芯片峰位置。

### 4. 位置刻度拟合

芯片的实际物理位置为固定值：

```
position = [-4, -3, -2, -1, 0, 1, 2, 3]
```

以 8 个芯片峰位置为 x，物理位置为 y，做线性拟合 `pol1`：

```
position = p0 + p1 * dt
```

输出拟合参数 p0、p1 和 chi2/ndf。

---

## 命令行用法

```bash
./ppac_normalize -r <run_number> [-t <trigger>] [-c <config_file>]
```

### 参数说明

| 参数 | 缩写 | 必需 | 默认值 | 说明 |
|------|------|------|--------|------|
| `--run` | `-r` | 是 | — | Run 号 |
| `--trigger` | `-t` | 否 | `main` | 触发类型（`main` 或 `t1`） |
| `--config` | `-c` | 否 | `config.toml` | 配置文件路径 |
| `--help` | `-h` | 否 | — | 打印帮助信息 |

### 使用示例

```bash
# main trigger
./ppac_normalize -r 510

# t1 trigger
./ppac_normalize -r 510 -t t1
```

---

## 输入

### 输入文件

```
{workspace}/ingot/ppac_{trigger}{run:04d}.root
```

例如：
- `ingot/ppac_0510.root`（main trigger）
- `ingot/ppac_t10510.root`（t1 trigger）

### 输入 TTree 格式

Tree 名：`tree`，包含 `PpacEvent` 结构：

```
flag:       int
valid[15]:  bool      — 通道有效性
time[15]:   double    — 时间（ns）
energy[15]: int       — 能量（ADC）
```

---

## 输出

### 输出 ROOT 文件

```
{workspace}/normalize/ppac_offset_{t1_}run{run:04d}.root
```

例如：`normalize/ppac_offset_0510.root`、`normalize/ppac_offset_t1_0510.root`

#### ROOT 文件内容（每个 PPAC 6 个对象）

| 对象名 | 类型 | 说明 |
|--------|------|------|
| `ppac1_Tx_Delay` | TH1D | X 方向延迟时间直方图 |
| `ppac1_Ty_Delay` | TH1D | Y 方向延迟时间直方图 |
| `ppac1_dtx12` | TH1D | X 方向时间差谱（dT = X1 - X2） |
| `ppac1_dty21` | TH1D | Y 方向时间差谱（dT = Y1 - Y2） |
| `ppac1_x_graph` | TGraph | X 方向位置刻度拟合图 |
| `ppac1_y_graph` | TGraph | Y 方向位置刻度拟合图 |

### 输出 TXT 文件

```
{workspace}/normalize/ppac_offset_{t1_}run{run:04d}.txt
```

#### TXT 文件格式（TSV，制表符分隔）

| 列名 | 说明 |
|------|------|
| `ppac_id` | PPAC 编号（1/2/3） |
| `delay_x_peak` | X 延迟时间峰位（ns） |
| `delay_x_sigma` | X 延迟时间峰宽（ns） |
| `delay_y_peak` | Y 延迟时间峰位（ns） |
| `delay_y_sigma` | Y 延迟时间峰宽（ns） |
| `position_x_p0` | X 位置刻度截距 |
| `position_x_p1` | X 位置刻度斜率 |
| `position_x_chi2/ndf` | X 拟合 chi2/ndf |
| `position_y_p0` | Y 位置刻度截距 |
| `position_y_p1` | Y 位置刻度斜率 |
| `position_y_chi2/ndf` | Y 拟合 chi2/ndf |
| `anode_eff(%)` | 阳极效率（%） |
| `x_cathode_eff(%)` | X 方向阴极效率（%） |
| `y_cathode_eff(%)` | Y 方向阴极效率（%） |

---

## 效率定义

| 效率 | 公式 |
|------|------|
| 阳极效率 | 阳极有效事件数 / 总事件数 × 100% |
| X 阴极效率 | 阳极有效且 X1/X2 有效事件数 / 阳极有效事件数 × 100% |
| Y 阴极效率 | 阳极有效且 Y1/Y2 有效事件数 / 阳极有效事件数 × 100% |

---

## 文件结构

```
src/brill/
├── bin/
│   ├── ppac_normalize.cpp          # 主程序
│   └── CMakeLists.txt              # 编译配置
├── include/
│   └── event/ingot/
│       └── ppac_event.h            # PpacEvent 结构体与 PpacChannelMap 定义
```

### 编译依赖

```
config ppac_event ROOT::Spectrum ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad
```