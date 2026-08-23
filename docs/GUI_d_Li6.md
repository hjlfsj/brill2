# GUI_d_Li6 — d+6Li 物理分析交互式可视化

## 功能概述

交互式 GUI 程序，用于 d+6Li 物理分析的可视化。直接读取 `extract_d_Li6` 输出的自包含 `D6LiEvent` ROOT 文件，无需回查 match/beam/calibration 等原始数据。

程序分为两个画布：
- **主界面**（嵌入式）：4 张原始能量相关二维直方图（2×2），展示 extract 阶段筛选出的全部事件
- **分析画布**（独立 TCanvas）：4 张物理分析图（2×2），包含 E-E 关联、E-θ 关联（含 Lise++ 参考线）和 θ-θ 关联

用法：`./GUI_d_Li6 [-c <config>]`

---

## 程序结构

| 路径 | 用途 |
|------|------|
| `src/brill/bin/GUI_d_Li6.cpp` | 主程序 |
| `src/brill/include/d_6Li/d_6Li_event.h` | D6LiEvent 结构体定义 |
| `src/brill/src/d_6Li/d_6Li_event.cpp` | SetupInput/Output 实现 |
| `src/brill/include/rebuild/rebuild_d_6Li.h` | 运动学计算函数声明 |
| `src/brill/src/rebuild/rebuild_d_6Li.cpp` | ComputeKinematics / LoadCutGFromFile 实现 |
| `src/brill/include/Lise++/e_theta.h` | Lise++ 曲线加载函数声明 |
| `src/brill/src/Lise++/e_theta.cpp` | LoadEThetaCurve / LoadThetaThetaCurve 实现 |
| `src/brill/include/physics/kinematics.h` | 角度计算通用函数 |

---

## 运行方式

### 参数

| 参数 | 说明 |
|------|------|
| `-c, --config` | 配置文件路径（默认 `config.toml`） |
| `-h, --help` | 打印帮助信息 |

### 示例

```bash
./GUI_d_Li6 -c config.toml
```

启动后通过 `File → Open` 选择 `d_Li6/` 目录下的 `extract_d_Li6_*.root` 文件。

---

## 窗口布局

程序创建 2 个窗口：

| 窗口 | 标题 | 尺寸 | 说明 |
|------|------|------|------|
| 主窗口 | `GUI_d_Li6` | 1200×900 | TGMainFrame，内含嵌入式 2×2 画布和菜单栏 |
| 分析画布 | `d+6Li Analysis` | 1200×900 | 独立 TCanvas，2×2 物理分析图 |

### 主界面（嵌入式画布）

```
┌────────────────────┬────────────────────┐
│  h_e1_10C_e2_10C   │  h_e1_6Li_e2_6Li   │
│  (E2_10C, E1_10C)  │  (E2_6Li, E1_6Li)  │
│  colz              │  colz              │
├────────────────────┼────────────────────┤
│  h_e2_10C_e3_10C   │  h_e3_10C_e4_10C   │
│  (E3_10C, E2_10C)  │  (E4_10C, E3_10C)  │
│  colz              │  colz              │
└────────────────────┴────────────────────┘
```

### 分析画布（独立 TCanvas）

```
┌────────────────────┬────────────────────┐
│  E_6Li vs E_10C    │  10C E-θ + ref    │
│  (E_10C, E_6Li)    │  (θ_10C, E_10C)    │
│  colz              │  + 红色 Lise++ 参考线│
├────────────────────┼────────────────────┤
│  6Li E-θ + ref     │  θ_10C vs θ_6Li   │
│  (θ_6Li, E_6Li)    │  (θ_10C, θ_6Li)    │
│  + 红色 Lise++ 参考线│  + 蓝色 Lise++ 参考线│
└────────────────────┴────────────────────┘
```

---

## 数据来源

程序**仅读取一个文件**：`extract_d_Li6` 的输出文件。该文件中的 `D6LiEvent` 已经是自包含的，所有能量已刻度、所有分类已完成。

| 文件 | 内容 |
|------|------|
| `d_Li6/extract_d_Li6_<trigger>_<run_start>_<run_end>.root` | TTree，含 31 个 branch 的 `D6LiEvent` |

额外读取的辅助文件（仅在打开文件时加载一次）：

| 文件 | 用途 |
|------|------|
| `src/brill/Cut/cal_d1_d2_6Li_cut.C` | TCutG，筛选 (E2_6Li, E1_6Li) 区域的 6Li 事件 |
| `assets/14O_d_6Li_0+_e_theta.txt` | Lise++ E-θ 参考线（10C 和 6Li） |
| `assets/14O_d_6Li_2+_theta_theta.txt` | Lise++ θ-θ 参考线 |

---

## 处理流程

### 1. 文件选择

通过 `File → Open` 打开 `d_Li6/` 目录下的 `extract_d_Li6_*.root` 文件，默认打开目录为 `paths.d_Li6`。

### 2. 加载事件列表

读取输入文件中的 `tree`，通过 `SetupInputD6Li` 将全部 `D6LiEvent` 加载到内存中的 `std::vector`。

### 3. 填充主界面直方图

遍历所有 `D6LiEvent`，填充 4 张原始能量相关图（不分束流类型）。

### 4. 加载分析辅助文件

- 加载 `cal_d1_d2_6Li_cut.C`（TCutG）
- 加载 Lise++ E-θ 参考线：10C 和 6Li 曲线
- 加载 Lise++ θ-θ 参考线

### 5. 填充分析画布直方图

二次遍历事件，同时满足以下条件的事件才进入分析画布：

- `ppac_valid == true`（PPAC 径迹有效）
- `(e2_6Li, e1_6Li)` 落在 `cal_d1_d2_6Li_cut` 区域内

对满足条件的事件计算运动学量：

- `E_6Li = e1_6Li + e2_6Li`（6Li 总能量）
- `E_10C = e1_10C + e2_10C + e3_10C + e4_10C`（10C 总能量）

### 6. 绘制

主界面 4 张图使用 `colz` 绘制，分析画布 4 张图同样使用 `colz` 绘制，并在 E-θ 和 θ-θ 图上叠加 Lise++ 参考线。

---

## 直方图

### 主界面（4 张）

| 名称 | x 轴 | y 轴 | 说明 |
|------|------|------|------|
| `h_e1_10C_e2_10C` | E2_10C (MeV) | E1_10C (MeV) | 10C 粒子 D1-D2 能量相关 |
| `h_e1_6Li_e2_6Li` | E2_6Li (MeV) | E1_6Li (MeV) | 6Li 粒子 D1-D2 能量相关 |
| `h_e2_10C_e3_10C` | E3_10C (MeV) | E2_10C (MeV) | 10C 粒子 D2-D3 能量相关 |
| `h_e3_10C_e4_10C` | E4_10C (MeV) | E3_10C (MeV) | 10C 粒子 D3-D4 能量相关 |

### 分析画布（4 张）

| 名称 | x 轴 | y 轴 | 说明 |
|------|------|------|------|
| `h_E_6Li_E_10C` | E_10C (MeV) | E_6Li (MeV) | 6Li 与 10C 总能量关联 |
| `h_10C_e_theta` | θ_10C (deg) | E_10C (MeV) | 10C E-θ 关联 + 红色 Lise++ 参考线 |
| `h_6Li_e_theta` | θ_6Li (deg) | E_6Li (MeV) | 6Li E-θ 关联 + 红色 Lise++ 参考线 |
| `h_theta_theta` | θ_10C (deg) | θ_6Li (deg) | θ-θ 关联 + 蓝色 Lise++ 参考线 |

---

## 筛选与分类

筛选和分类逻辑在 `extract_d_Li6` 中完成，详见 [d_6Li.md](d_6Li.md)。GUI 额外在分析画布中施加 `cal_d1_d2_6Li_cut` cut 和 `ppac_valid` 条件，用于选取高质量 d+6Li 反应事件。

---

## 进度输出

程序在终端输出详细运行进度：

```
Opening file: /data/disk1/ribll2026_www_data/d_Li6/extract_d_Li6_0057_0060.root
  Loading 43 events from d_Li6 file...
  Loaded 43 events
  Building main histograms...
  Main histograms filled: 43 events
  Loading cut file: src/brill/Cut/cal_d1_d2_6Li_cut.C
  Loading Lise++ reference curves...
  Building analysis histograms...
  Analysis: 43 events, 12 passed cut, 12 with ppac_valid
  Drawing...
  Done.
```

---

## 编译

在 `CMakeLists.txt` 中配置：

```cmake
add_executable(GUI_d_Li6 GUI_d_Li6.cpp)
target_link_libraries(
    GUI_d_Li6
    PRIVATE
    config
    d_6Li_event
    rebuild_physics
    physics
    lise_physics
    ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Gui
)
```

编译命令：

```bash
cmake --build build --target GUI_d_Li6 -j4
```

---

## 与旧版对比

| 方面 | 旧版 | 新版 |
|------|------|------|
| 输入文件 | 1 个 d_Li6 文件 + 逐 run 回查 match/beam/calibration/TCutG | 1 个 d_Li6 文件 |
| 依赖头文件 | 6 个 | 8 个（新增 rebuild_d_6Li、Lise++、physics） |
| I/O 次数 | 1 + 4×N_run + 1×N_run | 1 |
| 画布数量 | 4 个（All/14O/13N/12C） | 2 个（主界面 + 分析画布） |
| 直方图数量 | 16 张（4 束流 × 4 图） | 8 张（4 主界面 + 4 分析） |
| 分析功能 | 无 | E-E 关联、E-θ 关联、θ-θ 关联 + Lise++ 参考线 |
| 束流分类 | 4 组独立画布 | 不在主界面分束流，分析画布统一展示 |