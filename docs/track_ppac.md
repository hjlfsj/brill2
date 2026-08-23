# track_ppac — PPAC 径迹重建

## 概述

`track_ppac` 从 ingot 原始 PPAC 数据中读取时间信息，应用 PPAC offset 参数将时间差转换为物理位置，并通过最小二乘拟合重建粒子径迹（靶点坐标和方向余弦）。

用法：`./track_ppac -r <run> -t <trigger> [--draw]`

---

## 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-h, --help` | 打印帮助信息 | |
| `-r, --run` | Run 号 | 必填 |
| `-t, --trigger` | 触发类型：`main` 或 `t1` | `main` |
| `-c, --config` | 配置文件路径 | `config.toml` |
| `--draw` | 启用径迹线绘制直方图 | 关闭 |

## 示例

```bash
./track_ppac -r 60 -t main
./track_ppac -r 60 -t main --draw
```

---

## 输入文件

| 文件 | 路径 | 内容 |
|------|------|------|
| `ppac_<trigger>_<run>.root` | `config.paths.ingot` | 原始 PPAC 时间数据 (`PpacEvent`) |
| `ppac_offset_<run>.txt` | `config.paths.normalize` | PPAC offset 参数（由 `ppac_normalize` 产生） |

### PPAC Offset 文件格式

第一行为 header，之后每行对应一个 PPAC（共 3 个）：

```
ppac delay_x_peak delay_x_sigma delay_y_peak delay_y_sigma pos_x_p0 pos_x_p1 chi2_x pos_y_p0 pos_y_p1 chi2_y anode_eff x_eff y_eff
1  60.90  0.87  72.30  11.57  0.5  2.0  1.2  0.3  1.8  1.1  0.95  0.90  0.88
2  ...
3  ...
```

---

## 输出文件

| 文件 | 路径 | 内容 |
|------|------|------|
| `ppac_<prefix><run>.root` | `config.paths.track` | PPAC 径迹重建结果 (`PpacTrackEvent`) |

prefix 规则：`main` trigger 时为空，`t1` trigger 时为 `t1_`。

### TTree 分支

| Branch | 类型 | 说明 |
|--------|------|------|
| `valid` | `int` | 径迹有效性（x 和 y 方向均成功拟合） |
| `x_used_ppac` | `unsigned int` | 位掩码，标记哪些 PPAC 的 x 位置被使用 |
| `y_used_ppac` | `unsigned int` | 位掩码，标记哪些 PPAC 的 y 位置被使用 |
| `target_x` | `double` | 径迹外推到 z=0 的 x 坐标 (mm) |
| `target_y` | `double` | 径迹外推到 z=0 的 y 坐标 (mm) |
| `dir_x` | `double` | x 方向余弦 (dx/dz) |
| `dir_y` | `double` | y 方向余弦 (dy/dz) |
| `ppac_x[3]` | `double` | 三个 PPAC 的 x 位置 (mm) |
| `ppac_y[3]` | `double` | 三个 PPAC 的 y 位置 (mm) |

### 可选直方图（`--draw` 时）

| 名称 | 内容 |
|------|------|
| `h_track_zx` | Z-X 径迹线投影，1500 bins (-1000, 500) × 100 bins (-50, 50) |
| `h_track_zy` | Z-Y 径迹线投影，1500 bins (-1000, 500) × 100 bins (-50, 50) |

---

## 算法流程

### 1. 读取 PPAC offset 参数

从 `ppac_offset_<run>.txt` 读取 3 个 PPAC 的 offset 参数（`PpacOffsetParams`），每个 PPAC 包含：

| 参数 | 符号 | 说明 |
|------|------|------|
| `delay_x_peak`, `delay_x_sigma` | μ, σ | x 方向延迟时间的高斯拟合参数 |
| `delay_y_peak`, `delay_y_sigma` | μ, σ | y 方向延迟时间的高斯拟合参数 |
| `position_x_p0`, `position_x_p1` | p0, p1 | x 位置刻度：`pos = p0 + p1 * dt` |
| `position_y_p0`, `position_y_p1` | p0, p1 | y 位置刻度：`pos = p0 + p1 * dt` |

### 2. 逐事件计算 PPAC 位置

对每个 PPAC（共 3 个），使用 5 个通道的时间信号：

| 通道 | 说明 |
|------|------|
| X1, X2 | x 方向两端时间 |
| Y1, Y2 | y 方向两端时间 |
| anode | 阳极时间 |

#### 有效性检查

- 阳极有效：`valid[anode] && time[anode] > 0`
- x 有效：阳极有效且 X1、X2 均有效
- y 有效：阳极有效且 Y1、Y2 均有效

#### 延迟时间筛选

```
delay_x = time[X1] + time[X2] - 2 * time[anode]
```

若 `|delay_x - delay_x_peak| < 3 * delay_x_sigma`，则计算位置：

```
dtx = time[X1] - time[X2]
ppac_x = position_x_p0 + position_x_p1 * dtx + x_offset_mm
```

y 方向同理。

### 3. 直线拟合

收集所有有效 PPAC 的 (z, x) 和 (z, y) 点，分别进行最小二乘直线拟合：

- **≥2 个有效点**：`pos = target + dir * z`，得到 `target_x/y` 和 `dir_x/y`
- **仅 1 个有效点**（且为 PPAC1）：`target = ppac_x[0]`，`dir = 0`
- **0 个有效点**：`valid = 0`

### 4. 可选径迹线绘制

当 `--draw` 启用时，对每个有效径迹，在 z ∈ [-1000, 500] mm 范围内按 `target_x + dir_x * z` 填充直方图，用于可视化检查径迹质量。

---

## PPAC offset 的作用阶段

PPAC offset 在 `track_ppac` 阶段被使用，具体流程：

```
ppac_normalize          →  计算 offset 参数（delay peak, position scale）
    │
    └→ normalize/ppac_offset_*.txt
            │
track_ppac              →  读取 offset，将时间差转换为物理位置
    │                       pos = p0 + p1 * dt + config_offset
    │
    └→ track/ppac_*.root
```

**PPAC offset 仅在 `track_ppac` 中使用**，后续程序（`adjust_t0`、`track_t0`、`GUI_track`）读取的都是 track 文件中的物理位置，不再接触原始时间数据。

---

## 编译

```cmake
add_executable(track_ppac track_ppac.cpp)
target_link_libraries(track_ppac
    PRIVATE
    config
    ppac_event
    ppac_track_event
    ppac_track
    ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf
)
```