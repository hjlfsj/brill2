# GUI_track — PPAC + T0 径迹可视化

## 功能概述

交互式 GUI 程序，用于可视化 PPAC 径迹和 T0 DSSD 探测器 hit 的二维/三维分布。支持逐事件浏览、hit 数量过滤、多画布展示。T0 数据直接读取 DSSD match 文件（已废弃 T0 track 文件）。

用法：`./GUI_track -r <run> -t <trigger>`

---

## 数据来源

程序读取以下文件（路径均在 `data/` 目录下）：

| 文件 | 类型 | 说明 |
|------|------|------|
| `ppac_<trigger>_<run>.root` | `PpacTrackEvent` | PPAC 三层径迹，含靶点坐标和方向 |
| `t0d1_<trigger>_<run>.root` | `DssdMatchEvent` | D1 探测器正反面匹配结果 |
| `t0d2_<trigger>_<run>.root` | `DssdMatchEvent` | D2 探测器正反面匹配结果 |
| `t0d3_<trigger>_<run>.root` | `DssdMatchEvent` | D3 探测器正反面匹配结果 |
| `t0d4_<trigger>_<run>.root` | `DssdMatchEvent` | D4 探测器正反面匹配结果 |

每个 T0 match 文件由 `T0DetectorEntry` 结构体管理，包含 `file/tree/event` 指针以及探测器名称、z 位置、颜色、marker 样式。四个 T0 探测器存入 `t0_dets` 向量，索引 0/1/2/3 分别对应 D1/D2/D3/D4。

### PPAC 径迹 (PpacTrackEvent)

| Branch | 类型 | 说明 |
|--------|------|------|
| `valid` | int | 径迹有效性 |
| `target_x, target_y` | double | 径迹外推到靶点 (z=0) 的坐标 |
| `dir_x, dir_y` | double | 方向余弦 (dx/dz, dy/dz) |
| `ppac_x[3], ppac_y[3]` | double | 3 个 PPAC 的 hit 位置 |
| `x_used_ppac, y_used_ppac` | unsigned int | 位掩码，标记哪些 PPAC 被使用 |

### T0 match (DssdMatchEvent)

| Branch | 类型 | 说明 |
|--------|------|------|
| `num` | int | 该层匹配出的粒子数 |
| `energy[num]` | double | 粒子能量（原始值，未刻度） |
| `x[num], y[num], z[num]` | double | 粒子物理位置 (mm) |
| `front_strip[num], back_strip[num]` | double | 正反面加权条带位置 |
| `merge_tag[num]` | int | 配对类型 |
| `time[num]` | double | 时间 (ns) |

---

## 画布划分

程序创建 4 个独立画布：

| 画布 | 尺寸 | 内容 |
|------|------|------|
| canvas1 | 1600×800 | 主画布，z-x 和 z-y 径迹总览 |
| canvas2 | 600×600 | T0 原始事件详细信息 |
| canvas3 | 800×400 | 缩放视图，聚焦 D1-D4 区域 |
| canvas4 | 500×500 | 三维视图，靶点到 D2 hit 空间射线 |

### canvas1: 主画布

分为左右两个 pad，左侧为 z-x 投影，右侧为 z-y 投影。两者共享 z 轴范围 (-650, 150)，左右 pad 视觉上融合（移除相邻边框和标签）。左侧 pad 对应 z-x 投影，右侧 pad 对应 z-y 投影。两个 pad 的 x/y 轴范围固定为 (-35, 35)，使用 TH2F 框架控制坐标轴，不随数据变化。

### canvas2: T0 信息

左侧 pad 显示当前事件各探测器的原始数据，包括 D1-D4 每层每条 hit 的 strip 编号、能量、时间、物理坐标 (x, y, z)。右侧 pad 预留给拟合信息（当前为空白）。

### canvas3: 缩放视图

分为左右两个 pad，分别复制 z-x 和 z-y 投影，但绘图区域缩小到 z 方向 (80, 150)，聚焦 D1-D4 探测器区域。

### canvas4: 三维视图

TH3F 三维坐标系，坐标轴映射为：
- ROOT X 轴：z（束流方向，80-150 mm）
- ROOT Y 轴：x（-35, 35 mm）
- ROOT Z 轴：y（-35, 35 mm）

---

## 静态元素绘制

### 探测器位置线

在 z-x 和 z-y 画布中，用灰色虚线标注各探测器在 z 方向的位置，并绘制文字标签：

| 标注 | z 位置 (mm) | 颜色 |
|------|-------------|------|
| ppac1 | -625 | 灰色 |
| ppac2 | -531 | 灰色 |
| ppac3 | -389 | 灰色 |
| target | 0 | 灰色 |
| d1 | 98 | 灰色 |
| d2 | 118 | 灰色 |
| d3 | 128 | 灰色 |
| d4 | 138 | 灰色 |

### 靶点标记

在 z=0 处用黑色星形标记靶点位置（`target_x`, `target_y`）。

---

## 每事件动态绘制

### 主画布 (canvas1)

**z-x 投影**：
- PPAC 三层 hit 的 z-x 位置，蓝色粗线条
- T0 D1-D4 各层 hit 的 z-x 位置，每层用不同颜色圆点（D1=绿、D2=红、D3=蓝、D4=品红）
- 靶点位置，黑色星形
- 红色虚线：从靶点到每个 D2 hit 的连线，不同 hit 用不同颜色（红、绿、蓝、黄、品红、青...），绘制范围 z=(0, 150)

**z-y 投影**：与 z-x 对称，数据为 z-y 坐标。

### 缩放画布 (canvas3)

与主画布相同的绘制内容，但 z 范围缩小到 (80, 150)，x/y 范围保持 (-35, 35)。每条 D2 连线同样绘制。

### 三维画布 (canvas4)

使用 ROOT 3D 图形类绘制：

1. **TH3F 框架**：定义 3D 坐标空间，X=z(80-150), Y=x(-35,35), Z=y(-35,35)
2. **TPolyMarker3D 靶点**：黑色星形，位于 (target_x, target_y, 0)
3. **TPolyMarker3D D2 hit**：绿色圆点，位于各 D2 hit 的 (x, y, z) 位置
4. **TPolyLine3D 射线**：从靶点到每个 D2 hit 的虚线，射线终点外推到 z=140 平面，每条射线用不同颜色

### 绘制流程

`AddEventToCanvas` 函数对每个事件：
1. 将 PPAC hit 填充到 TGraph 中，加入 TMultiGraph
2. 遍历 `t0_dets` 向量中的每个 T0 探测器，从其 `DssdMatchEvent` 中读取 hit 的 x/y/z
3. 将所有 hit 填充到对应的 TGraph 中
4. 对 D2（`t0_dets[1]`）的 hit，绘制到各 D2 hit 的连线
5. 为 3D 画布创建 TPolyLine3D 和 TPolyMarker3D 对象

---

## 交互功能

### 按钮

| 按钮 | 功能 |
|------|------|
| Prev | 上一个事件 |
| Next | 下一个事件 |

按钮位于 canvas1 底部。

### 输入控件

| 控件 | 说明 |
|------|------|
| T0 hit_num | 过滤条件，使用 TGNumberEntry 实现 |

过滤逻辑：
- 值为 0：仅显示 T0 无 hit 的事件
- 值为 -1：显示所有事件，不过滤
- 其他值 N：仅显示 `num == N` 的事件

点击 Prev/Next 时，会跳过不满足过滤条件的事件。

---

## 事件导航流程

点击 Prev 或 Next 按钮时：

1. 根据当前 entry 索引和过滤条件，查找上一个/下一个满足条件的事件
2. 调用 `ClearCanvas` 清除所有画布上的动态图形对象（TGraph、TPolyLine3D、TPolyMarker3D 等）
3. 调用 `AddEventToCanvas` 绘制当前事件的动态数据
4. 更新所有画布，刷新显示
5. 更新 entry 编号显示

---

## PPAC offset 应用位置

PPAC 的 `x_offset_mm` 和 `y_offset_mm` 在 `track_ppac.cpp` 的 `TrackPpac` 函数中应用，而非 GUI_track 中。

### 应用时机

在 `TrackPpac` 函数中，对每个 PPAC（共 3 个），当 X 方向的 delay 信号有效时：

```cpp
track.ppac_x[ppac_idx] = offset[ppac_idx].position_x_p0
    + offset[ppac_idx].position_x_p1 * dtx
    + ppac_config.x_offset_mm[ppac_idx];
```

其中 `position_x_p0` 和 `position_x_p1` 是 delay → position 的刻度参数，`x_offset_mm[ppac_idx]` 是额外位置偏移。

### offset 来源

`config.toml` 中 `[detectors.ppac]` 节：

```toml
x_offset_mm = [0, -0.156, 0.078]
y_offset_mm = [0.0, 0.491, -0.246]
```

三个值分别对应 PPAC1、PPAC2、PPAC3。由 `adjust_ppac_position` 程序拟合，手动填入配置文件。

### 数据流

PPAC offset 在 `track_ppac.cpp` 中已加入，GUI_track 读取的 `PpacTrackEvent` 中 `ppac_x/ppac_y` 已经是修正后的值，GUI_track 不再做任何位置修正。

---

## T0 DSSD offset 应用位置

T0 DSSD 的 offset 在 `dssd_match.cpp` 的 `FillPhysicalPosition` 函数中应用，详见 `dssd_matching_algorithm.md`。GUI_track 读取的 `DssdMatchEvent` 中 `x/y/z` 已经是修正后的值，GUI_track 不再做任何位置修正。

### 与 PPAC offset 的对比

| 特性 | T0 DSSD offset | PPAC offset |
|------|---------------|-------------|
| 配置文件 | `[detectors.t0d1]` ~ `[detectors.t0d4]` | `[detectors.ppac]` |
| 数量 | 每个探测器独立 1 组 | 每个 PPAC 独立 1 组（共 3 组） |
| 应用位置 | `dssd_match.cpp` → `FillPhysicalPosition` | `ppac_track.cpp` → `TrackPpac` |
| 应用阶段 | DSSD 正反面匹配（strip→物理坐标） | PPAC 径迹重建（delay→position） |
| 计算程序 | `adjust_t0_step1`, `adjust_t0_step2` | `adjust_ppac_position` |

---

## 依赖

- ROOT: TGMainFrame, TCanvas, TPad, TH2F, TH3F, TMultiGraph, TGraph, TPolyLine3D, TPolyMarker3D, TGNumberEntry, TGTextButton, TPaveText
- glimmer 库: PpacTrackEvent, DssdMatchEvent
- `config.toml`: 探测器位置配置（z_mm, size_x_mm, size_y_mm 等）