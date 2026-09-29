# d+6Li 分析说明文档

## 1. 概述

本项目包含 d+6Li 物理分析的数据提取和 GUI 可视化程序：

- **extract_d_Li6**：从 match、beam、track 文件中提取并刻度数据，生成自包含的 `D6LiEvent` ROOT 文件
- **GUI_d_Li6**：交互式 GUI，读入 `D6LiEvent` 文件进行物理分析，包括束流筛选、粒子鉴别、能损修正和激发能谱重建

---

## extract_d_Li6 说明文档

## 1. 概述

`extract_d_Li6` 是 d+6Li 物理分析的数据提取程序，从 match、beam、track 文件中读取探测器数据和束流分类信息，在 extract 阶段一次性完成所有 I/O 和计算（刻度、筛选、分类），输出自包含的 `D6LiEvent` ROOT 文件，供 `GUI_d_Li6` 直接读取，无需回查原始数据。

## 2. 目录结构

| 路径 | 用途 |
|------|------|
| `src/brill/bin/d_6Li/extract_d_Li6.cpp` | 主程序 |
| `src/brill/include/d_6Li/extract.h` | 筛选条件函数声明 |
| `src/brill/src/d_6Li/extract.cpp` | 筛选条件函数实现 |
| `src/brill/include/d_6Li/d_6Li_event.h` | D6LiEvent 结构体定义 |
| `src/brill/src/d_6Li/d_6Li_event.cpp` | SetupInput/Output 实现 |
| `src/brill/include/physics/kinematics.h` | 角度计算通用函数 |
| `src/brill/src/physics/kinematics.cpp` | AngleBetween / AngleWithZ 实现 |

## 3. 运行方式

```bash
extract_d_Li6 [OPTION...]
```

### 参数

| 参数 | 说明 |
|------|------|
| `-h, --help` | 打印帮助信息 |
| `-r, --run` | 起始 run 号 |
| `-e, --end-run` | 结束 run 号（可选，默认等于 `-r`） |
| `-t, --trigger` | 触发类型 (如 `t1`, `main`) |
| `-c, --config` | 配置文件路径（默认 `config.toml`） |

### 示例

```bash
extract_d_Li6 -r 1000 -e 1010 -t t1
```

## 4. 输入文件

程序读取三类文件，以 `t0d1` 的 tree entry 作为主循环，按 entry 号同步读取所有文件：

| 文件 | 来源 | 内容 |
|------|------|------|
| `match/t0d1_<trigger>_<run>.root` | `match_dssd` | D1 `DssdMatchEvent` |
| `match/t0d2_<trigger>_<run>.root` | `match_dssd` | D2 `DssdMatchEvent` |
| `match/t0d3_<trigger>_<run>.root` | `match_dssd` | D3 `DssdMatchEvent` |
| `match/t0d4_<trigger>_<run>.root` | `match_dssd` | D4 `DssdMatchEvent` |
| `beam/beam_<trigger>_<run>.root` | `sort_beam` | `is_14O` / `is_13N` / `is_12C` |
| `track/ppac_<trigger>_<run>.root` | `track_ppac` | `PpacTrackEvent` |

## 5. 输出文件

- 文件名：`extract_d_Li6_<trigger>_<run_start>_<run_end>.root`
- 保存路径：`config.toml` 中 `paths.d_Li6` 指定的目录

### 输出内容

#### TTree：`tree` — 自包含的 `D6LiEvent`

所有能量已刻度（MeV），所有分类已完成，GUI 可直接读取无需回查原始数据。

| Branch | 类型 | 说明 |
|--------|------|------|
| `run_number` | `int` | 事件所在 run 号 |
| `entry` | `Long64_t` | 事件在对应 run 中的 entry 号 |
| `e1` | `double` | D1 第一个 hit 刻度后能量 (MeV) |
| `e2` | `double` | D2 第一个 hit 刻度后能量 (MeV) |
| `e3` | `double` | D3 第一个 hit 刻度后能量 (MeV) |
| `e4` | `double` | D4 第一个 hit 刻度后能量 (MeV) |
| `e1_10C` | `double` | D1 中 10C 分类 hit 的刻度后能量 |
| `e2_10C` | `double` | D2 中 10C 分类 hit 的刻度后能量 |
| `e3_10C` | `double` | D3 中 10C 分类 hit 的刻度后能量 |
| `e4_10C` | `double` | D4 中 10C 分类 hit 的刻度后能量 |
| `e1_6Li` | `double` | D1 中 6Li 分类 hit 的刻度后能量 |
| `e2_6Li` | `double` | D2 中 6Li 分类 hit 的刻度后能量 |
| `is_14O` | `bool` | 束流分类：14O（来自 `sort_beam`） |
| `is_13N` | `bool` | 束流分类：13N（来自 `sort_beam`） |
| `is_12C` | `bool` | 束流分类：12C（来自 `sort_beam`） |
| `ppac_valid` | `bool` | PPAC 径迹是否有效（来自 `track_ppac`） |
| `target_x` | `double` | PPAC 径迹在 target 处的 x 坐标 (mm) |
| `target_y` | `double` | PPAC 径迹在 target 处的 y 坐标 (mm) |
| `dir_x` | `double` | PPAC 径迹 x 方向的方向余弦 |
| `dir_y` | `double` | PPAC 径迹 y 方向的方向余弦 |
| `t0d2_10C_x` | `double` | 10C 粒子在 t0d2 上的 x 位置 (mm) |
| `t0d2_10C_y` | `double` | 10C 粒子在 t0d2 上的 y 位置 (mm) |
| `t0d2_10C_z` | `double` | 10C 粒子在 t0d2 上的 z 位置 (mm，即 DSSD2 的 z) |
| `t0d2_6Li_x` | `double` | 6Li 粒子在 t0d2 上的 x 位置 (mm) |
| `t0d2_6Li_y` | `double` | 6Li 粒子在 t0d2 上的 y 位置 (mm) |
| `t0d2_6Li_z` | `double` | 6Li 粒子在 t0d2 上的 z 位置 (mm，即 DSSD2 的 z) |
| `theta_beam` | `double` | 束流与 z 轴夹角（度，0-180），由 PPAC dir_x/dir_y 计算 |
| `phi_beam` | `double` | 束流方位角（度），atan2(dir_y, dir_x) |
| `theta_10C` | `double` | 10C 出射方向与束流方向夹角（度），ppac_valid=false 时置 0 |
| `theta_6Li` | `double` | 6Li 出射方向与束流方向夹角（度），ppac_valid=false 时置 0 |
| `opening_angle` | `double` | 10C 与 6Li 出射方向夹角（度），ppac_valid=false 时置 0 |

#### TH2D：能量相关图

| 名称 | x 轴 | y 轴 | 说明 |
|------|------|------|------|
| `h_d1d2` | D2 第一个 hit 能量 | D1 第一个 hit 能量 | D1-D2 能量相关 |
| `h_d2d3` | D3 第一个 hit 能量 | D2 第一个 hit 能量 | D2-D3 能量相关 |
| `h_d3d4` | D4 第一个 hit 能量 | D3 第一个 hit 能量 | D3-D4 能量相关 |
| `h_e1_10C_e2_10C` | E2_10C | E1_10C | 10C 粒子 D1-D2 能量相关 |
| `h_e1_6Li_e2_6Li` | E2_6Li | E1_6Li | 6Li 粒子 D1-D2 能量相关 |
| `h_e2_10C_e3_10C` | E3_10C | E2_10C | 10C 粒子 D2-D3 能量相关 |
| `h_e3_10C_e4_10C` | E4_10C | E3_10C | 10C 粒子 D3-D4 能量相关 |

## 6. 筛选条件

### 6.1 基础筛选 (`PassD6LiCut`)

定义在 [extract.cpp](file:///home/ribll2026/ribll2026_www/github_code/brill2/src/brill/src/d_6Li/extract.cpp)：

- D1 (`t0d1`)：`num == 2`
- D2 (`t0d2`)：`num == 2`
- D3 (`t0d3`)：`num == 1`
- D4 (`t0d4`)：`num == 1`
- 刻度后的 (E4, E3) 落在 `cal_d3_d4_10C_cut` 图形 cut 内

### 6.2 高级分类 (`ClassifyD6Li`)

在基础筛选通过后，进一步区分 10C 和 6Li 粒子：

1. 位置匹配：`abs(x_d3 - x_d4) < 2 mm` 且 `abs(y_d3 - y_d4) < 2 mm`
2. D2 分类：D2 的两个 hit 中有且仅有一个满足：
   - 刻度后的 (E3, E2) 落在 `cal_d2_d3_10C_cut` 图形 cut 内
   - `abs(x_d2 - x_d3) < 2 mm` 且 `abs(y_d2 - y_d3) < 2 mm`
   - 满足条件的为 `e2_10C`，另一个为 `e2_6Li`
3. D1 分类：计算 D1 两个 hit 到 `e2_10C` 位置的距离平方，距离较小的为 `e1_10C`，另一个为 `e1_6Li`

### 6.3 图形 Cut 文件

| 文件 | 用途 | x 轴 | y 轴 |
|------|------|------|------|
| `Cut/cal_d3_d4_10C_cut.C` | D3-D4 10C 粒子筛选 | E4 | E3 |
| `Cut/cal_d2_d3_10C_cut.C` | D2-D3 10C 粒子筛选 | E3 | E2 |

## 7. 刻度系数

在填充 TH2D 直方图前，程序会从 `<workspace>/<calibration>/t0.txt` 加载刻度系数，对每个探测器的能量应用刻度：

```
calibrated_energy = p0[layer] + p1[layer] × raw_energy
```

| layer 索引 | 对应探测器 |
|-----------|-----------|
| 1 | t0d1 |
| 2 | t0d2 |
| 3 | t0d3 |
| 4 | t0d4 |

刻度文件格式：第一行为 header，之后每行 `index p0 p1`。

## 8. 编译

在 `src/brill/bin/CMakeLists.txt` 中已添加编译目标，与其他 bin 程序（如 `GUI_track`）放在同一位置。

## 9. 数据流

```
track_ppac ──→ track/ppac_*.root ─────────────┐
match_dssd ──→ match/t0d{d}_*.root ───────────┤
sort_beam ───→ beam/beam_*.root ──────────────┤
calibrate_t0 → calibration/t0.txt ────────────┤
                                               │
                                               ▼
                                        extract_d_Li6
    ├── 逐 entry 同步读取 match + beam + PPAC track 文件
    ├── 加载刻度系数，对能量进行刻度
    ├── 应用 PassD6LiCut 筛选条件
    ├── 应用 ClassifyD6Li 分类（10C / 6Li）
    ├── 写入自包含 D6LiEvent（31 个 branch，含位置和角度）
    └── 填充 TH2D 能量相关图
                                               │
                                               ▼
                                    d_Li6/extract_d_Li6_*.root
                                               │
                                               ▼
                                         GUI_d_Li6
    (直接读取 D6LiEvent，无需回查任何原始文件)
```

---

## GUI_d_Li6 说明文档

### 概述

`GUI_d_Li6` 是 d+6Li 反应的交互式 GUI 物理分析程序，包含三个独立画布。程序通过 `File > Open` 载入 `extract_d_Li6` 生成的 ROOT 文件，支持束流选择、6Li 的 D1-D2 图形 cut 筛选，并在分析方法中集成了 CD₂ 靶能损修正和 ¹⁰C 激发能谱重建。

### 运行方式

```bash
GUI_d_Li6 -c config.toml
```

| 参数 | 说明 |
|------|------|
| `-c, --config` | 配置文件路径（默认 `config.toml`） |
| `-h, --help` | 打印帮助信息 |

### GUI 布局

- **顶部菜单栏**：File > Open（打开 ROOT 文件），File > Quit（退出）
- **束流选择按钮**：All / 14O / 13N / 12C，用于筛选束流成分
- **Run 范围筛选**：可选，限制加载的 run 号范围（默认加载所有）
- **三个独立画布**：主嵌入画布 + 两个弹出窗口

### 画布 1：粒子鉴别（2×2 布局）

嵌入在主窗口中的 `TRootEmbeddedCanvas`，展示四张 10C 和 6Li 的能量相关图：

| Pad | 直方图 | 内容 |
|-----|--------|------|
| 左上 | `h_e1_10C_e2_10C` | 10C 的 E1-E2 能量相关 |
| 右上 | `h_e1_6Li_e2_6Li` | 6Li 的 E1-E2 能量相关 |
| 左下 | `h_e2_10C_e3_10C` | 10C 的 E2-E3 能量相关 |
| 右下 | `h_e3_10C_e4_10C` | 10C 的 E3-E4 能量相关 |

右侧 6Li 的 E1-E2 图同时叠加 `Cut/cal_d1_d2_6Li_cut.C` 图形 cut。

### 画布 2：物理分析关联（2×2 布局）

弹出窗口 `canvas_analysis`，展示 E-E、E-θ、θ-θ 关联图，并叠加理论曲线：

| Pad | 直方图 | 内容 | 参考线 |
|-----|--------|------|--------|
| 左上 | `h_E_6Li_E_10C` | E₁₀C vs E₆Li 二维关联 | — |
| 右上 | `h_10C_e_theta` | ¹⁰C 能量-角度关联 | `14O_d_6Li_0+_e_theta.txt`（¹⁰C） |
| 左下 | `h_6Li_e_theta` | ⁶Li 能量-角度关联 | `14O_d_6Li_0+_e_theta.txt`（⁶Li） |
| 右下 | `h_theta_theta` | θ₁₀C vs θ₆Li | `14O_d_6Li_0+_theta_theta.txt` |

> **能损修正**：此画布中使用的 E₁₀C 和 E₆Li 已通过 CD₂ 靶能损函数补全。具体地，将探测器测得的动能 E 反推至穿出半靶厚（61 µm）之前的入射能量 E₀，考虑粒子出射角 θ 对有效路径长度的影响：
>
> ```
> E₀ = TargetEnergyLoss(Z, A).IncidentEnergy(E, θ, kCD2)
> E_eff = 61 / cos(θ)  [µm]
> ```

### 画布 3：¹⁰C 激发能谱（2×2 布局）

弹出窗口 `canvas_excitation`，2×2 布局，第一个 pad 展示 ¹⁰C 的激发能谱：

- **直方图**：`h_excitation`，范围为 -10 到 20 MeV（200 bins）
- **X 轴**：Eₓ(¹⁰C) (MeV)，**Y 轴**：Counts

其余 3 个 pad 预留待后续使用。

#### 重建方法

利用非相对论运动学和动量守恒重建入射 ¹⁴O 的动能，进而计算 ¹⁰C 的激发能。反应为：

```
¹⁴O + d → ¹⁰C + ⁶Li    (Q = -8.64 MeV, 吸热反应)
```

d 靶核在实验室系静止。步骤如下：

1. **动量计算**（非相对论，m = A × 931.5 MeV/c²）：

   ```
   p₁₀C = sqrt(2 · m₁₀C · E₁₀C_corr)
   p₆Li = sqrt(2 · m₆Li · E₆Li_corr)
   ```

2. **方向向量**（从靶点到 D2 位置）：

   ```
   d₁₀C = normalize(t0d2_10C_x - target_x, t0d2_10C_y - target_y, t0d2_10C_z - 0)
   d₆Li = normalize(t0d2_6Li_x - target_x, t0d2_6Li_y - target_y, t0d2_6Li_z - 0)
   ```

3. **动量守恒**（d 静止，p₄ = 0）：

   ```
   p₁₄O = p₁₀C · d₁₀C + p₆Li · d₆Li
   E₁₄O = |p₁₄O|² / (2 · m₁₄O)
   ```

4. **激发能**：

   ```
   Eₓ = E₁₄O - E₁₀C_corr - E₆Li_corr - Q      (Q = 8.64 MeV)
   ```

   其中 E₁₀C_corr 和 E₆Li_corr 为 CD₂ 靶修正后的能量。

### 筛选流程

1. **束流筛选**：根据按钮选择，筛选 `is_14O` / `is_13N` / `is_12C`
2. **PPAC 有效性**：要求 `ppac_valid == true`
3. **6Li 图形 cut**：要求 (E2₆Li, E1₆Li) 落在 `Cut/cal_d1_d2_6Li_cut.C` 内
4. **方向向量有效性**：要求 ¹⁰C 和 ⁶Li 的 D2 位置向量长度 > 0

### 编译

在 `src/brill/bin/CMakeLists.txt` 中已添加编译目标：

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

### 相关文件

| 文件 | 用途 |
|------|------|
| `src/brill/bin/GUI_d_Li6.cpp` | GUI 主程序 |
| `src/brill/include/Lise++/target_energy_loss.h` | CD₂/CH₂ 靶能损计算（TargetEnergyLoss 类） |
| `src/brill/src/Lise++/target_energy_loss.cpp` | 靶能损函数实现 |
| `src/brill/include/rebuild/rebuild_d_6Li.h` | ComputeKinematics / LoadCutGFromFile 声明 |
| `src/brill/src/rebuild/rebuild_d_6Li.cpp` | 运动学计算和 cut 加载实现 |
| `assets/14O_d_6Li_0+_e_theta.txt` | 基态反应 E-θ 理论曲线 |
| `assets/14O_d_6Li_0+_theta_theta.txt` | 基态反应 θ-θ 理论曲线 |
| `src/brill/Cut/cal_d1_d2_6Li_cut.C` | 6Li 的 D1-D2 图形 cut |