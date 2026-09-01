# GUI_pid — DSSD 粒子鉴别可视化

## 功能概述

交互式 GUI 程序，用于 DSSD (t0d1-t0d4) 和 t0s 的粒子鉴别（PID）可视化。直接读取 match 文件和 ingot 文件，应用能量刻度后绘制二维关联图。

程序分为两个画布（均为 2×2 布局）：
- **主画布**（嵌入式）：4 张 hit0 能量关联直方图，展示 D1-D2、D2-D3、D3-D4、D4-T0S 的 hit0 能量关联
- **次画布**（独立 TCanvas）：4 张 hit1 能量关联直方图，展示 D1-D2、D2-D3、D3-D4、D4-T0S 的 hit1 关联。**子图按需绘制**：仅当左探测器存在 hit1 时才绘制对应子图，无 hit1 的探测器复用 hit0 能量

用法：`./GUI_pid -c <config>`

---

## 程序结构

| 路径 | 用途 |
|------|------|
| `src/brill/bin/GUI_pid.cpp` | 主程序 |
| `src/brill/include/d_6Li/extract.h` | D6LiCalibration 结构体、CalibrateD6LiEnergy、ReadD6LiCalibration |
| `src/brill/src/d_6Li/extract.cpp` | 刻度函数实现（编译为 `d_6Li_extract` 库） |
| `src/brill/include/event/t0/dssd_match_event.h` | DssdMatchEvent 结构体 |
| `src/brill/include/event/ingot/silicon_event.h` | SiliconEvent 结构体（t0s） |

---

## 运行方式

```bash
# 基本用法
./GUI_pid -c config.toml
```

---

## 界面说明

### 控制面板

| 控件 | 说明 |
|------|------|
| Trigger | 触发类型：`main` / `t1` |
| Run | 起始 run 号 ~ 结束 run 号 |
| d1_hit ~ d4_hit | 各探测器要求的 hit 数（上下箭头调节，范围 0-7） |
| Draw | 开始读取数据并绘图 |

**hit 数语义**：`d1_hit=2` 表示筛选 `d1.num == 2` 的事件（严格等于，非 ≥）。主画布使用 `energy[0]`（hit0），次画布使用 `energy[1]`（hit1）。

**约束**：前一层的 hit 数 ≥ 后一层（`hd1 >= hd2 >= hd3 >= hd4`），目前仅处理 hit=1 或 2。

### 主画布（嵌入式）

四张 hit0 二维直方图，标题标注 `(hit0)`：

| 位置 | 内容 | 坐标轴 |
|------|------|--------|
| 左上 | D1-D2 关联 | X: D2 Energy, Y: D1 Energy |
| 右上 | D2-D3 关联 | X: D3 Energy, Y: D2 Energy |
| 左下 | D3-D4 关联 | X: D4 Energy, Y: D3 Energy |
| 右下 | D4-T0S 关联 | X: T0S Energy, Y: D4 Energy |

### 次画布（独立 TCanvas）

四张 hit1 二维直方图，标题动态标注每个探测器使用的 hit。**子图按需绘制**：

| 条件 | 绘制 | 标题 |
|------|------|------|
| d1 有 hit1 | D1-D2 | `D1(hit1)-D2(hit1/hit0)` |
| d2 有 hit1 | D2-D3 | `D2(hit1)-D3(hit1/hit0)` |
| d3 有 hit1 | D3-D4 | `D3(hit1)-D4(hit1/hit0)` |
| d4 有 hit1 | D4-T0S | `D4(hit1)-T0S` |

右侧探测器若无 hit1，则复用其 hit0 能量。例如 `hit=(2,2,2,1)`：
- D1-D2: `D1(hit1)-D2(hit1)` ✓
- D2-D3: `D2(hit1)-D3(hit1)` ✓
- D3-D4: `D3(hit1)-D4(hit0)` ✓（d4 复用 hit0）
- D4-T0S: **不绘制**（d4 无 hit1）

### 进度显示

- `trigger=main` 时：每读完一个 run 刷新一次画布，状态栏实时更新
- `trigger=t1` 时：全部读完一次性绘制
- 终端打印每个 run 的百分比进度（10% 步进）

---

## 数据流

```
match/t0d1_*.root ─┐
match/t0d2_*.root  │
match/t0d3_*.root  ├──→ GUI_pid ──→ 主画布 + 次画布
match/t0d4_*.root  │
ingot/t0s_*.root  ─┘
        │
calibration/t0.txt ──→ 能量刻度 (p0 + p1 * E_raw)
```

**能量刻度**：复用 `d_6Li_extract` 库中的 `ReadD6LiCalibration` 和 `CalibrateD6LiEnergy`，刻度系数与 `extract_d_Li6` 完全一致。层映射：

| layer | 探测器 |
|-------|--------|
| 0 | t0d1 |
| 1 | t0d2 |
| 2 | t0d3 |
| 3 | t0d4 |
| 4 | t0s |

---

## 典型使用场景

```bash
# 查看 hit=2 的事件（每个探测器恰好 2 个 hit）
# 主画布显示 hit0 关联，次画布显示 hit1 关联
./GUI_pid -c config.toml
# 设置 d1_hit=2, d2_hit=2, d3_hit=2, d4_hit=2 → Draw

# 查看 hit=(2,2,2,1) 的事件
# 次画布 d4-t0s 不绘制，d3-d4 复用 d4 的 hit0
./GUI_pid -c config.toml
# 设置 d1_hit=2, d2_hit=2, d3_hit=2, d4_hit=1 → Draw
```