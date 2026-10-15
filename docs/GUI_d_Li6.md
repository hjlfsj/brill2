# GUI_d_Li6 — d+6Li 物理分析交互式可视化

## 功能概述

交互式 GUI 程序，用于 d+6Li 物理分析的可视化。直接读取 `extract_d_Li6` 输出的自包含 `D6LiEvent` ROOT 文件，无需回查 match/beam/calibration 等原始数据。

程序分为三个画布：
- **主界面**（嵌入式）：4 张原始能量相关二维直方图（2×2），展示 extract 阶段筛选出的全部事件
- **分析画布**（独立 TCanvas）：4 张物理分析图（2×2），包含 E-E 关联、E-θ 关联（含 Lise++ 参考线）和 θ-θ 关联
- **激发能画布**（独立 TCanvas）：3 张激发能谱（2×2），通过三种方法计算 ¹⁰C 激发能

用法：`./GUI_d_Li6 [-c <config>]`

---

## 程序结构

| 路径 | 用途 |
|------|------|
| `src/brill/bin/GUI_d_Li6.cpp` | 主程序（单文件，~820 行） |
| `src/brill/include/d_6Li/d_6Li_event.h` | D6LiEvent 结构体定义 |
| `src/brill/src/d_6Li/d_6Li_event.cpp` | SetupInput/Output 实现 |
| `src/brill/include/rebuild/rebuild_d_6Li.h` | 运动学计算函数声明 |
| `src/brill/src/rebuild/rebuild_d_6Li.cpp` | ComputeKinematics / LoadCutGFromFile 实现 |
| `src/brill/include/Lise++/e_theta.h` | Lise++ 曲线加载函数声明 |
| `src/brill/src/Lise++/e_theta.cpp` | LoadEThetaCurve / LoadThetaThetaCurve 实现 |
| `src/brill/include/physics/kinematics.h` | 角度计算通用函数 |
| `src/brill/include/energy_calculator/range_energy_calculator.h` | RangeEnergyCalculator（能损修正） |

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

程序创建 3 个窗口：

| 窗口 | 标题 | 尺寸 | 说明 |
|------|------|------|------|
| 主窗口 | `GUI_d_Li6` | 1200×900 | TGMainFrame，内含嵌入式 2×2 画布、菜单栏和控制面板 |
| 分析画布 | `d+6Li Analysis` | 1200×800 | 独立 TCanvas，2×2 物理分析图 |
| 激发能画布 | `10C Excitation Energy` | 800×600 | 独立 TCanvas，2×2 激发能谱（第 4 面板预留） |

### 控制面板

主窗口上方有一排控制组件，用于调整分析参数：

| 控件 | 类型 | 说明 |
|------|------|------|
| `All` / `14O` / `13N` / `12C` | TGTextButton | 束流类型筛选，切换后自动重分析已有数据 |
| `Run: [min] - [max]` | TGNumberEntry × 2 | 按 run 范围过滤事件（0 表示不限制） |
| `Bins: [N] /bin` | TGNumberEntry | 激发能直方图分 bin 数（5–500，默认 30） |
| `z(µm): [val]` | TGNumberEntry | 靶内反应深度偏移量（-61–61 µm，默认 0） |
| `Draw` | TGTextButton | 重新打开文件并全量重分析（适用修改参数后） |

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
│  colz              │  + Lise++ 参考线    │
├────────────────────┼────────────────────┤
│  6Li E-θ + ref     │  θ_10C vs θ_6Li   │
│  (θ_6Li, E_6Li)    │  (θ_10C, θ_6Li)    │
│  + Lise++ 参考线    │  + Lise++ 参考线    │
└────────────────────┴────────────────────┘
```

Lise++ 参考线：
- 绿色实线（`kGreen+2`）：0⁺ g.s. 参考线
- 红色虚线（`kRed, style=7`）：2⁺ 激发态参考线

### 激发能画布（独立 TCanvas）

```
┌──────────────────────┬──────────────────────┐
│  h_excitation_mm     │ h_excitation_invmass │
│  (missing mass, ¹⁴O) │  (invariant mass)    │
│  ¹⁰C Ex (MeV)        │  ¹⁰C Ex (MeV)         │
├──────────────────────┼──────────────────────┤
│ h_excitation_mm_fixed│       （预留）        │
│ (missing mass,       │                      │
│  E_b=35.6 MeV/u)     │                      │
└──────────────────────┴──────────────────────┘
```

---

## 数据来源

程序**仅读取一个数据文件**：`extract_d_Li6` 的输出文件。该文件中的 `D6LiEvent` 已经是自包含的，所有能量已刻度、所有分类已完成。

| 文件 | 内容 |
|------|------|
| `d_Li6/extract_d_Li6_<trigger>_<run_start>_<run_end>.root` | TTree，含 31 个 branch 的 `D6LiEvent` |

额外读取的辅助文件：

| 文件 | 用途 |
|------|------|
| `src/brill/Cut/cal_d1_d2_Li_large_cut.C` | TCutG，筛选 (E2_6Li, E1_6Li) 区域的 6Li 事件 |
| `assets/14O_d_6Li_0+_e_theta.txt` | Lise++ E-θ 参考线（0⁺，10C 和 6Li） |
| `assets/14O_d_6Li_2+_e_theta.txt` | Lise++ E-θ 参考线（2⁺，10C 和 6Li） |
| `assets/14O_d_6Li_0+_theta_theta.txt` | Lise++ θ-θ 参考线（0⁺） |
| `assets/14O_d_6Li_2+_theta_theta.txt` | Lise++ θ-θ 参考线（2⁺） |
| `beam/beam_<run>.root` | TOF 束流能量文件（仅 14O 模式，按 run 按需读取） |
| `<energy_calculator>/target_z6_a10_cd2.root` | ¹⁰C 在 CD₂ 靶中的 Range-Energy 缓存 |
| `<energy_calculator>/target_z3_a6_cd2.root` | ⁶Li 在 CD₂ 靶中的 Range-Energy 缓存 |
| `<energy_calculator>/target_z8_a14_cd2.root` | ¹⁴O 束流在 CD₂ 靶中的 Range-Energy 缓存（方法一和方法三共用） |

---

## 处理流程

### 1. 文件选择

通过 `File → Open` 打开 `d_Li6/` 目录下的 `extract_d_Li6_*.root` 文件，默认打开目录为 `paths.d_Li6`。

### 2. 加载事件列表

读取输入文件中的 `tree`，通过 `SetupInputD6Li` 将全部 `D6LiEvent` 加载到内存中的 `std::vector`。支持按 run 范围过滤（`Run` 控件）。

### 3. 填充主界面直方图

根据当前束流类型（All/14O/13N/12C）过滤并填充 4 张原始能量相关二维图。

### 4. 加载分析辅助文件

- 加载 `cal_d1_d2_Li_large_cut.C`（TCutG）
- 加载 Lise++ E-θ 参考线（0⁺ 和 2⁺，各含 10C 和 6Li 曲线）
- 加载 Lise++ θ-θ 参考线（0⁺ 和 2⁺）
- 初始化 RangeEnergyCalculator（¹⁰C、⁶Li 在 CD₂ 靶中的能损修正）
- 初始化 ¹⁴O 束流的 RangeEnergyCalculator 缓存（所有模式均尝试加载，用于方法三；方法一仅 14O 模式执行）

### 5. 填充分析画布直方图

遍历事件，同时满足以下条件的事件才进入分析：

- 通过束流筛选（`PassBeamFilter`）
- `ppac_valid == true`（PPAC 径迹有效）
- `(e2_6Li, e1_6Li)` 落在 `cal_d1_d2_Li_large_cut` 区域内

对满足条件的事件进行能损修正和运动学计算。

---

## 能损修正

由于产物粒子穿过剩余 CD₂ 靶物质，需进行能损修正以还原反应点的能量。

| 粒子 | RangeEnergyCalculator | 缓存文件 |
|------|----------------------|----------|
| ¹⁰C (Z=6, A=10) | `calc_10C` | `target_z6_a10_cd2.root` |
| ⁶Li (Z=3, A=6) | `calc_6Li` | `target_z3_a6_cd2.root` |
| ¹⁴O (Z=8, A=14) | `s_beam_calc`（所有模式） | `target_z8_a14_cd2.root` |

修正公式：

- `exit_thick = half_thick_um / 2 - z`（穿过剩余靶的半厚度，z 为反应深度偏移）
- `range_rxn = range(E_measured) + exit_thick / cos(θ)`（还原到反应点的 range）
- `E_corrected = Energy(range_rxn)`（由 range 反算能量）

其中 `half_thick_um = 122.0 / 2 = 61 µm`，`z` 由控件 `z(µm)` 设置。

---

## 激发能计算

激发能通过三种方法计算：

### 方法一：Missing Mass — TOF 束流能量（仅 14O 模式）

仅当 `g_beam_filter == 1`（14O 模式）且束流 RangeEnergyCalculator 缓存存在时执行。

1. 通过 TOF 读取束流能量 `E_beam_raw`（按 run、entry 查询 `beam/beam_<run>.root`）
2. 束流能损修正：`beam_range_rxn = Range(E_beam_raw) - beam_depth`，`E_beam = Energy(beam_range_rxn)`
3. 构造四动量：
   - `P_A = (p_A * dir_x/norm, p_A * dir_y/norm, p_A * 1/norm, E_A)`（¹⁴O 束流）
   - `P_B = (0, 0, 0, m_d)`（静止氘靶）
   - `P_6Li`（测量的 ⁶Li）
4. `P_10C_star = P_A + P_B - P_6Li`（反冲 ¹⁰C 的四动量）
5. `E_x_mm = M(P_10C_star) - m_10C`（激发能 = 不变量质量 - 基态质量）

### 方法二：Invariant Mass（所有束流类型）

对所有束流类型均执行：

1. 由测量的 ¹⁰C 和 ⁶Li 动量构造四动量 `P_10C`、`P_6Li`
2. 反推 ¹⁴O 四动量：`P_14O_reco = P_10C + P_6Li`（动量守恒），能量由质量壳关系计算
3. `P_10C_star_inv = P_14O_reco + P_d - P_6Li`
4. `E_x_inv = M(P_10C_star_inv) - m_10C`

### 方法三：Missing Mass — 固定束流能量（所有束流类型）

对所有束流类型均执行（不依赖 `g_beam_filter`），使用固定束流能量 **35.6 MeV/u**（对应 ¹⁴O 名义束流能量），无需读取 beam 文件。

1. `E_beam_fixed = 35.6 × 14 = 498.4 MeV`（动能）
2. 束流能损修正：`beam_range_rxn = Range(E_beam_fixed) - beam_depth`，`E_beam_at_rxn = Energy(beam_range_rxn)`
3. 构造四动量（与本方法一相同）：
   - `P_A`：固定能量的 ¹⁴O 束流
   - `P_B = (0, 0, 0, m_d)`：静止氘靶
   - `P_6Li`：测量的 ⁶Li
4. `P_10C_star = P_A + P_B - P_6Li`
5. `E_x_mm_fixed = M(P_10C_star) - m_10C`

> **注意**：方法三在所有束流类型下均用 ¹⁴O 的 35.6 MeV/u 作为近似束流能量，对于 13N、12C 数据属于交叉检验，不作为精确测量。

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
| `h_E_6Li_E_10C` | E_10C (MeV) | E_6Li (MeV) | 6Li 与 10C 总能损修正后能量关联 |
| `h_10C_e_theta` | θ_10C (deg) | E_10C (MeV) | 10C E-θ 关联 + Lise++ 参考线 |
| `h_6Li_e_theta` | θ_6Li (deg) | E_6Li (MeV) | 6Li E-θ 关联 + Lise++ 参考线 |
| `h_theta_theta` | θ_10C (deg) | θ_6Li (deg) | θ-θ 关联 + Lise++ 参考线 |

### 激发能画布（3 张，2×2 布局，第 4 面板预留）

| 名称 | 说明 | 填充条件 |
|------|------|----------|
| `h_excitation_mm` | ¹⁰C 激发能（方法一：TOF missing mass） | 仅 14O 模式 + 束流缓存存在 |
| `h_excitation_invmass` | ¹⁰C 激发能（方法二：invariant mass） | 所有束流类型 |
| `h_excitation_mm_fixed` | ¹⁰C 激发能（方法三：固定 E_b=35.6 MeV/u missing mass） | 所有束流类型 + 束流缓存存在 |

激发能直方图范围：-10 ~ 20 MeV，bin 数由 `Bins` 控件设置（默认 30，即 ~1 MeV/bin）。

---

## 束流能量读取

`ReadBeamEnergy(run, entry, beam_dir)` 按需从 `beam/beam_<run>.root` 读取 TOF 束流能量，采用带缓存的惰性加载策略：

- 首次访问某个 run 时，打开对应 `beam_<run>.root` 文件并保持打开状态
- 后续访问同一 run 时直接使用已缓存的 `TTree*`
- `BeamFileCache` 结构体使用移动语义管理 `TFile*` 生命周期，确保 map 中的缓存独占文件所有权

### BeamFileCache 设计

```cpp
struct BeamFileCache {
    TFile *file = nullptr;
    TTree *tree = nullptr;
    double E_beam = 0.0;
    int run = -1;
    ~BeamFileCache() { if (file) { file->Close(); delete file; } }
    // 移动语义：转移所有权，源对象置空
    BeamFileCache(BeamFileCache&& other) noexcept;
    BeamFileCache& operator=(BeamFileCache&& other) noexcept;
    // 禁止拷贝：防止多个对象共享同一 TFile* 导致 double-free
    BeamFileCache(const BeamFileCache&) = delete;
    BeamFileCache& operator=(const BeamFileCache&) = delete;
};
```

**关键**：移动语义确保 `g_beam_cache[run] = std::move(cache)` 后，局部 `cache` 的 `file` 指针被置空，析构时无操作；map 中的条目独占文件所有权。若无移动语义，自定义析构函数会阻止编译器生成移动构造函数，导致退化为拷贝——两个对象共享同一 `TFile*`，局部对象析构时释放文件，造成 map 中悬垂指针的段错误。

---

## 按钮交互逻辑

| 操作 | 触发函数 | 行为 |
|------|---------|------|
| `File → Open` | `OnFileOpen()` | 重新读取文件 → 加载事件 → 填充主界面 → 填充分析 → 绘制 |
| 束流按钮（All/14O/13N/12C） | `OnBeamChanged()` | 仅重新筛选已有事件并重分析（不重读文件） |
| `Draw` 按钮 | `OnFileOpen()` | 与 File→Open 相同，重读文件并全量重分析 |
| 修改 Run/Bins/z 参数后点 `Draw` | `OnFileOpen()` | 用新参数重新加载文件并分析 |

**注意**：束流按钮切换时不重读文件——若修改了 Run/Bins/z 参数，需先点击 `Draw` 使其生效，再切换束流。

---

## 筛选与分类

筛选和分类逻辑在 `extract_d_Li6` 中完成，详见 [d_6Li.md](d_6Li.md)。GUI 额外施加：

| 条件 | 适用画布 |
|------|---------|
| `PassBeamFilter`（is_14O / is_13N / is_12C） | 主界面、分析画布、激发能画布 |
| `ppac_valid == true` | 分析画布、激发能画布 |
| `(e2_6Li, e1_6Li)` 在 `cal_d1_d2_Li_large_cut` 内 | 分析画布、激发能画布 |

---

## 进度输出

程序在终端输出详细运行进度：

```
Opening file: /data/disk1/ribll2026_www_data/d_Li6/extract_d_Li6_0057_0060.root
  Loading 43 events from d_Li6 file...
  Loaded 43 events
  Done: 43 events, passed=43
  Loading cut: src/brill/Cut/cal_d1_d2_Li_large_cut.C
  --- Selected events (run, entry) ---
  Analysis done: 43 events, passed=12 (27%)
  Done.
```

---

## 常见问题

### 14O 模式崩溃

如曾遇到过切换到 14O 时崩溃，原因是 `BeamFileCache` 自定义析构函数阻止了移动构造函数的隐式生成，导致 `std::move` 退化为拷贝，文件被提前释放。已通过显式定义移动语义修复。

### 束流缓存文件不存在

若 `target_z8_a14_cd2.root` 不存在，程序将打印 WARNING 并跳过 missing mass 方法（方法一）。该文件由 `GUI_10C_4He` 在 14O 模式下首次运行时自动生成。请确保先运行 `GUI_10C_4He` 的 14O 模式。

### Draw 按钮行为

点击 `Draw` 会重新打开文件（`OnFileOpen()`），而非仅重新分析已有事件。若只需切换束流筛选，使用束流按钮（All/14O/13N/12C）即可，无需重读文件。

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
    energy_calculator
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
| 依赖头文件 | 6 个 | 9 个（新增 rebuild_d_6Li、Lise++、physics、energy_calculator） |
| I/O 次数 | 1 + 4×N_run + 1×N_run | 1（+ beam 文件按 run 惰性加载） |
| 画布数量 | 4 个（All/14O/13N/12C） | 3 个（主界面 + 分析 + 激发能） |
| 直方图数量 | 16 张（4 束流 × 4 图） | 10 张（4 主界面 + 4 分析 + 2 激发能） |
| 分析功能 | 无 | E-E 关联、E-θ 关联、θ-θ 关联 + Lise++ 参考线（0⁺+2⁺） |
| 激发能 | 无 | Missing mass + invariant mass 双方法 |
| 能损修正 | 无 | CD₂ 靶中 ¹⁰C/⁶Li/¹⁴O 能损修正 |
| 束流分类 | 4 组独立画布 | 束流按钮动态切换，无需重读文件 |
| 交互控件 | 无 | Run 范围、Bins、z 深度 + Draw 重分析 |