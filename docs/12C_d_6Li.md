# 12C+d → 6Li+2α — 符合测量与 GUI 可视化

## 功能概述

12C+d → 6Li+2α 符合测量分析，研究 12C 与 d 靶反应产生 6Li 和两个 α 粒子的三体反应道。程序分为两个阶段：

- **extract_6Li_two4He**：从 match/beam/track/ingot 等原始数据中提取 6Li+2α 符合事件，经筛选、能量刻度和粒子分类后输出自包含的 ROOT 文件
- **GUI_12C_d_6Li**：交互式 GUI，读取 extract 阶段的输出文件，进行可视化筛选和物理分析

---

## 程序结构

| 路径 | 用途 |
|------|------|
| `src/brill/bin/C12_d_6Li/extract_6Li_two4He.cpp` | 数据提取主程序 |
| `src/brill/bin/GUI_12C_d_6Li.cpp` | 交互式可视化主程序 |
| `src/brill/include/C12_d_6Li/C12_d_6Li_event.h` | C12D6LiEvent 结构体定义 |
| `src/brill/src/C12_d_6Li/C12_d_6Li_event.cpp` | SetupInputC12D6Li / SetupOutputC12D6Li 实现 |
| `src/brill/include/C12_d_6Li/extract.h` | C12D6LiCalibration 结构体、刻度函数、筛选函数声明 |
| `src/brill/src/C12_d_6Li/extract.cpp` | PassC12D6LiCut / ClassifyTwo4He 实现 |
| `src/brill/src/C12_d_6Li/CMakeLists.txt` | C12D6Li_event 和 C12D6Li_extract 库构建 |
| `src/brill/include/physics/kinematics.h` | AngleBetween / AngleWithZ 角度计算函数 |
| `src/brill/Cut/cal_d1_d2_6Li_cut.C` | TCutG 截断，筛选 (E2_6Li, E1_6Li) 区域 |

---

## 探测器布局与能量定义

12C+d → 6Li+2α 反应中，6Li 和两个 α 粒子击中不同的探测器层：

| 探测器 | 碎片 | 能量变量 | 说明 |
|--------|------|----------|------|
| T0 D1 | 6Li + 4He1 + 4He2 | e1_6Li, e1_4He1, e1_4He2 | 第一层 DSSD，三个 hit 同时击中 |
| T0 D2 | 6Li + 4He1 + 4He2 | e2_6Li, e2_4He1, e2_4He2 | 第二层 DSSD，三个 hit |
| T0 D3 | 4He1 + 4He2 | e3_4He1, e3_4He2 | 第三层 DSSD，两个 hit（6Li 已停止） |
| T0 D4 | 4He1 + 4He2 | e4_4He1, e4_4He2 | 第四层 DSSD，两个 hit |
| T0 S1 | (CsI) | e5 | CsI 闪烁体，总残余能量 |

6Li 总动能：T_6Li = e1_6Li + e2_6Li
4He1 总动能：T_4He1 = e1_4He1 + e2_4He1 + e3_4He1 + e4_4He1
4He2 总动能：T_4He2 = e1_4He2 + e2_4He2 + e3_4He2 + e4_4He2

---

## extract_6Li_two4He — 数据提取

### 运行方式

| 参数 | 说明 |
|------|------|
| `-r, --run` | 起始 run 号（必填） |
| `-e, --end-run` | 结束 run 号（可选，默认等于起始 run） |
| `-t, --trigger` | 触发类型：`main` 或 `t1`（必填） |
| `-c, --config` | 配置文件路径（默认 `config.toml`） |
| `-h, --help` | 打印帮助信息 |

### 示例

```bash
./extract_6Li_two4He -r 57 -e 70 -t main
```

### 输入文件

每个 run 需要读取以下文件（从 workspace 对应目录）：

| 文件 | 来源 | 内容 |
|------|------|------|
| `match/t0d1_<trigger>_<run>.root` | match_dssd | D1 正反面匹配事件 |
| `match/t0d2_<trigger>_<run>.root` | match_dssd | D2 正反面匹配事件 |
| `match/t0d3_<trigger>_<run>.root` | match_dssd | D3 正反面匹配事件 |
| `match/t0d4_<trigger>_<run>.root` | match_dssd | D4 正反面匹配事件 |
| `ingot/t0s_<trigger>_<run>.root` | 原始数据 | S1 CsI 事件 |
| `beam/beam_<trigger>_<run>.root` | sort_beam | 束流粒子分类 |
| `track/ppac_<trigger>_<run>.root` | track_ppac | PPAC 径迹 |
| `calibration/t0.txt` | calibrate_t0 | T0 能量刻度参数 |

### 输出文件

输出到 `workspace/C12_d_6Li/` 目录（由 `config.toml` 中 `[paths]` 的 `c12_d_6Li` 配置，默认为 `C12_d_6Li`）：

```
extract_6Li_two4He_<trigger>_<run_start>_<run_end>.root
```

输出文件包含：
- **TTree `tree`**：C12D6LiEvent 结构体，含能量、角度、束流分类、PPAC 径迹等全部信息
- **8 个 TH2D**：6Li 和两个 α 粒子的各层能量关联图

### 事件筛选流程

```
读取各探测器数据
    │
    ├── PassC12D6LiCut 筛选：
    │   ├── d1.num == 3  （三个粒子击中 D1）
    │   ├── d2.num == 3  （三个粒子击中 D2）
    │   ├── d3.num == 2  （两个粒子击中 D3，6Li 已停止）
    │   ├── d4.num == 2  （两个粒子击中 D4）
    │   └── (E2_6Li, E1_6Li) 落在 cal_d1_d2_6Li_cut 内
    │
    ├── ClassifyTwo4He 粒子分类：
    │   ├── 6Li 识别：D1/D2 中能量最大的 hit 为 6Li
    │   ├── 4He1/4He2 识别：通过位置最近邻匹配将剩余 2 个 hit 区分为 4He1 和 4He2
    │   │   - D1 中除 6Li 外的两个 hit 按与 D2 中 4He1/4He2 位置最近原则分配
    │   │   - D3 和 D4 中的两个 hit 同样按最近邻分配
    │   └── 返回每个探测器上各粒子的 hit 索引
    │
    ├── 能量刻度：CalibrateC12D6LiEnergy(calib, layer, raw)
    │
    ├── 角度计算（依赖 PPAC 径迹）：
    │   ├── theta_beam, phi_beam：束流入射方向角
    │   ├── theta_6Li：6Li 相对于束流方向的散射角
    │   ├── theta_4He1：4He1 相对于束流方向的散射角
    │   ├── theta_4He2：4He2 相对于束流方向的散射角
    │   ├── opening_6Li_4He1：6Li 与 4He1 之间的张角
    │   ├── opening_6Li_4He2：6Li 与 4He2 之间的张角
    │   └── opening_4He1_4He2：4He1 与 4He2 之间的张角
    │
    └── 填充输出 TTree 和直方图
```

### 粒子分类算法

`ClassifyTwo4He` 函数实现了三体事件中 6Li 和两个 4He 的分类：

1. **6Li 识别**：在 D1 和 D2 中分别找出能量最大的 hit 作为 6Li。6Li 在探测器中的能量沉积通常比 α 大
2. **6Li Cut 验证**：刻度后的 (E2_6Li, E1_6Li) 必须落在 `cal_d1_d2_6Li_cut.C` 的 TCutG 区域内
3. **4He1/4He2 区分**：D2 中除 6Li 外的两个 hit 作为锚点（4He1, 4He2）。对 D1（剩余 2 个 hit）、D3、D4 分别独立处理：
   - 该层两个 hit H[0], H[1]，D2 两个锚点 A[0], A[1]
   - 组合 A：H[0]→A[0], H[1]→A[1]，sum_A = d²(H[0],A[0]) + d²(H[1],A[1])
   - 组合 B：H[0]→A[1], H[1]→A[0]，sum_B = d²(H[0],A[1]) + d²(H[1],A[0])
   - 取 sum 较小的组合，确定该层每个 hit 归属 4He1 还是 4He2
4. **相邻层位置一致性检查**：分类完成后，对 4He1 和 4He2 分别检查相邻探测器层之间的位置偏移：
   - D1↔D2：|dx| < 2mm 且 |dy| < 2mm
   - D2↔D3：|dx| < 2mm 且 |dy| < 2mm
   - D3↔D4：|dx| < 2mm 且 |dy| < 2mm
   - 任一粒子任一层不满足条件则丢弃该事件

---

## GUI_12C_d_6Li — 交互式可视化

### 运行方式

```bash
./GUI_12C_d_6Li -c config.toml
```

启动后通过 `File → Open` 选择 `C12_d_6Li/` 目录下的 `extract_6Li_two4He_*.root` 文件。

### 窗口布局

程序创建 2 个画布：

| 画布 | 标题 | 布局 | 说明 |
|------|------|------|------|
| 主画布（嵌入式） | `GUI_12C_d_6Li` | 2×3 | 6 张原始能量关联直方图，展示经束流筛选后的事件 |
| 次画布（独立） | `12C+d->6Li+2alpha Second` | 2×3 | 3 张能量关联直方图 + 2 张 6Li E-θ 图（全部/筛选）+ 1 个预留空位 |

### 主画布（嵌入式，2×3）

```
┌──────────────────────┬──────────────────────┬──────────────────────┐
│  E1_6Li vs E2_6Li    │  E1_4He1 vs E2_4He1  │  E1_4He2 vs E2_4He2  │
│  (E2_6Li, E1_6Li)    │  (E2_4He1, E1_4He1)  │  (E2_4He2, E1_4He2)  │
│  colz                │  colz                │  colz                │
├──────────────────────┼──────────────────────┼──────────────────────┤
│  E2_4He1 vs E3_4He1  │  E2_4He2 vs E3_4He2  │  E3_4He1 vs E4_4He1  │
│  (E3_4He1, E2_4He1)  │  (E3_4He2, E2_4He2)  │  (E4_4He1, E3_4He1)  │
│  colz                │  colz                │  colz                │
└──────────────────────┴──────────────────────┴──────────────────────┘
```

### 次画布（独立 TCanvas，2×3）

```
┌──────────────────────┬──────────────────────┬──────────────────────┐
│  E3_4He2 vs E4_4He2  │  E4sum vs E5          │  6Li E-θ (全部)      │
│  (E4_4He2, E3_4He2)  │  (E5, E4_4He1+        │  (θ_6Li, T_6Li)     │
│  colz                │   E4_4He2)            │  colz + 理论曲线     │
│                      │  colz                │                      │
├──────────────────────┼──────────────────────┼──────────────────────┤
│  (预留)               │  (预留)               │  6Li E-θ (4He cut)  │
│                      │                      │  (θ_6Li, T_6Li)     │
│                      │                      │  colz + 理论曲线     │
└──────────────────────┴──────────────────────┴──────────────────────┘
```

次画布中：
- `E4sum vs E5`：两个 α 粒子在 D4 上的能量之和与 CsI 中残余能量的关联图，用于评估粒子停止情况
- `6Li E-θ (全部)`：所有通过束流筛选的事件的 6Li 动能-散射角关联图，叠加三条理论曲线
- `6Li E-θ (4He cut)`：进一步要求 4He1 和 4He2 的 (e4, e3) 均落在 `cal_d3_d4_all_4He_cut` 内，同样叠加理论曲线

### 理论曲线

次画布的 6Li E-θ 图上叠加了三条 12C 激发态的理论运动学曲线，从 `assets/` 目录加载：

| 态 | 文件 | 颜色/线型 |
|----|------|----------|
| 0⁺（基态） | `12C_d_6Li_0+_e_theta.txt` | 蓝色实线 |
| 2⁺（4.44 MeV） | `12C_d_6Li_2+_e_theta.txt` | 绿色虚线 |
| 4⁺（14.08 MeV） | `12C_d_6Li_4+_e_theta.txt` | 品红点线 |

曲线通过 `brill::LoadEThetaCurve(path, "6Li")` 加载，提取 6Li 列的角度和每核子能量，换算为总能量绘制。

### 束流选择按钮

主窗口顶部提供 4 个束流选择按钮，**对所有画布同时生效**：

| 按钮 | 含义 | 影响 |
|------|------|------|
| `All` | 不过滤束流 | 显示全部事件 |
| `14O` | 仅 14O 束流 | 仅显示 is_14O=true 的事件 |
| `13N` | 仅 13N 束流 | 仅显示 is_13N=true 的事件 |
| `12C` | 仅 12C 束流 | 仅显示 is_12C=true 的事件 |

切换按钮时，主画布和次画布均自动重填，直方图标题中会显示当前束流标签（如 `(14O)`）。

### Run 号范围选择

主窗口顶部提供两个 `TGNumberEntry` 控件，用于指定 run number 范围：

- **Run Min**：run 号下限（闭区间），默认 `0` 表示不设下限
- **Run Max**：run 号上限（闭区间），默认 `0` 表示不设上限

当两个值均为 `0`（默认）时，加载文件中所有事件；否则仅加载 `run_number` 在 `[min, max]` 区间内的事件。

该筛选在 `File → Open` 时生效，加载时即过滤，所有画布均受影响。控制台会打印筛选后的事件数与跳过的数量。

### Draw 按钮

`Draw` 按钮位于顶部控制栏最右侧。按下后，使用当前 run 范围和束流选择重新加载当前文件并重绘所有画布，无需重新打开文件对话框。

---

## 刻度系数

在填充 TH2D 直方图前，程序会从 `<workspace>/<calibration>/t0.txt` 加载刻度系数，对每个探测器的能量应用刻度：

```
calibrated_energy = p0[layer] + p1[layer] × raw_energy
```

| layer 索引 | 对应探测器 |
|-----------|-----------|
| 0 | t0d1 |
| 1 | t0d2 |
| 2 | t0d3 |
| 3 | t0d4 |
| 4 | t0s1 (CsI) |

刻度文件格式：第一行为 header，之后每行 `index p0 p1`。

---

## C12D6LiEvent 数据结构

| Branch | 类型 | 说明 |
|--------|------|------|
| `run_number` | `int` | 事件所在 run 号 |
| `entry` | `Long64_t` | 事件在对应 run 中的 entry 号 |
| `e1_6Li` | `double` | 6Li 在 D1 的刻度后能量 (MeV) |
| `e2_6Li` | `double` | 6Li 在 D2 的刻度后能量 (MeV) |
| `e1_4He1` | `double` | 4He1 在 D1 的刻度后能量 (MeV) |
| `e2_4He1` | `double` | 4He1 在 D2 的刻度后能量 (MeV) |
| `e3_4He1` | `double` | 4He1 在 D3 的刻度后能量 (MeV) |
| `e4_4He1` | `double` | 4He1 在 D4 的刻度后能量 (MeV) |
| `e1_4He2` | `double` | 4He2 在 D1 的刻度后能量 (MeV) |
| `e2_4He2` | `double` | 4He2 在 D2 的刻度后能量 (MeV) |
| `e3_4He2` | `double` | 4He2 在 D3 的刻度后能量 (MeV) |
| `e4_4He2` | `double` | 4He2 在 D4 的刻度后能量 (MeV) |
| `e5` | `double` | CsI 闪烁体总残余能量 (MeV) |
| `is_14O` | `bool` | 束流分类：14O |
| `is_13N` | `bool` | 束流分类：13N |
| `is_12C` | `bool` | 束流分类：12C |
| `ppac_valid` | `bool` | PPAC 径迹是否有效 |
| `target_x` | `double` | PPAC 径迹在 target 处的 x 坐标 (mm) |
| `target_y` | `double` | PPAC 径迹在 target 处的 y 坐标 (mm) |
| `dir_x` | `double` | PPAC 径迹 x 方向的方向余弦 |
| `dir_y` | `double` | PPAC 径迹 y 方向的方向余弦 |
| `t0d2_6Li_x/y/z` | `double` | 6Li 在 D2 上的位置 (mm) |
| `t0d2_4He1_x/y/z` | `double` | 4He1 在 D2 上的位置 (mm) |
| `t0d2_4He2_x/y/z` | `double` | 4He2 在 D2 上的位置 (mm) |
| `theta_beam` | `double` | 束流与 z 轴夹角（度） |
| `phi_beam` | `double` | 束流方位角（度） |
| `theta_6Li` | `double` | 6Li 出射方向与束流方向夹角（度） |
| `theta_4He1` | `double` | 4He1 出射方向与束流方向夹角（度） |
| `theta_4He2` | `double` | 4He2 出射方向与束流方向夹角（度） |
| `opening_6Li_4He1` | `double` | 6Li 与 4He1 出射方向夹角（度） |
| `opening_6Li_4He2` | `double` | 6Li 与 4He2 出射方向夹角（度） |
| `opening_4He1_4He2` | `double` | 4He1 与 4He2 出射方向夹角（度） |

---

## 依赖关系

```
extract_6Li_two4He
    ├── C12D6Li_event (C12_d_6Li_event.cpp)
    ├── C12D6Li_extract (extract.cpp)
    ├── config
    ├── dssd_match_event
    ├── silicon_event
    ├── ppac_track_event
    ├── beam_sort
    ├── physics (kinematics.cpp)
    └── utils

GUI_12C_d_6Li
    ├── C12D6Li_event
    ├── config
    └── utils
```

---

## 数据处理流程

```
ingot/ (原始数据)
    │
    ├── match_dssd ──────────→ match/t0d{1,2,3,4}_*.root
    ├── track_ppac ──────────→ track/ppac_*.root
    ├── sort_beam ───────────→ beam/beam_*.root
    ├── calibrate_t0 ────────→ calibration/t0.txt
    │
    └── extract_6Li_two4He ─→ C12_d_6Li/extract_6Li_two4He_*.root
                                  │
                                  └── GUI_12C_d_6Li ──→ 可视化
```

---

## 配置文件

在 `config.toml` 的 `[paths]` 中添加：

```toml
[paths]
c12_d_6Li = "C12_d_6Li"
```

并在 workspace 下创建对应目录 `C12_d_6Li/`。