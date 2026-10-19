









## 三Alpha事件提取与分析方法

### 概述

本项目用于提取和分析 T0 探测器中的三alpha（⁴He + ⁴He + ⁴He）事件，目标物理过程为 ¹²C → ³α 衰变。程序分为两个阶段：

- **extract_three_alpha**：从 match、beam、track 文件中读取探测器数据，进行径迹匹配、PID cut 验证、能量修正，输出自包含的 ROOT 文件
- **GUI_three_alpha**：交互式 GUI，读取 extract 阶段的输出文件，进行 PID 分布可视化和 ⁸Be/¹²C 激发能谱（不变质量法）分析

---

### 目录结构

| 路径 | 用途 |
|------|------|
| `src/brill/bin/three_alpha/extract_three_alpha.cpp` | 数据提取主程序 |
| `src/brill/bin/three_alpha/GUI_three_alpha.cpp` | 交互式可视化主程序 |
| `src/brill/include/three_alpha/three_alpha_event.h` | `ThreeAlphaEvent` 结构体定义（输出事件格式） |
| `src/brill/src/three_alpha/three_alpha_event.cpp` | `SetupInputThreeAlpha` / `SetupOutputThreeAlpha` 实现 |
| `src/brill/include/three_alpha/extract.h` | 校准、筛选、分类、能量反推函数声明 |
| `src/brill/src/three_alpha/extract.cpp` | 三粒子径迹匹配、PID cut 验证、二分法穿透能量反推等核心逻辑 |
| `src/brill/src/three_alpha/CMakeLists.txt` | `three_alpha_event` 和 `three_alpha_extract` 库构建 |
| `src/brill/include/energy_calculator/lost_energy_calculator.h` | `LostEnergyCalculator`：穿过给定厚度后的能量损失计算 |
| `src/brill/include/energy_calculator/range_energy_calculator.h` | `RangeEnergyCalculator`：range-energy 关系（catima） |
| `src/brill/include/rebuild/nuclear_data.h` | `GetMass(z, a)`：从 AME2020 表获取原子核静止质量 |
| `src/brill/include/physics/kinematics.h` | `AngleBetween` / `AngleWithZ`：角度计算通用函数 |

---

### 依赖的 config.toml 配置

```toml
[paths]
three_alpha = "three_alpha"
```

程序输出到 `<workspace>/three_alpha/` 目录。

---

### 依赖的 TCutG 文件（需用户自行创建）

文件位于 `src/brill/Cut/` 下，通过 ROOT macro (.C) 格式加载，在 `extract_three_alpha` 启动时一次性读入内存：

| 文件 | x 轴 | y 轴 | 用途 |
|------|------|------|------|
| `cal_d1_d2_tran_He4_cut.C` | D2 能量 (MeV) | D1 能量 (MeV) | D1-D2 层面 ⁴He 透射 PID 选区和α鉴别 |
| `cal_d2_d3_tran_He4_cut.C` | D3 能量 (MeV) | D2 能量 (MeV) | D2-D3 层面 ⁴He 透射 PID 选区 |
| `cal_d3_d4_stop_He4_cut.C` | D4 能量 (MeV) | D3 能量 (MeV) | D3-D4 层面 ⁴He **阻停** PID 选区 |
| `cal_d3_d4_tran_He4_cut.C` | D4 能量 (MeV) | D3 能量 (MeV) | D3-D4 层面 ⁴He **穿透** PID 选区 |

> **注意**：d3-d4 有两种 cut（stop / tran），因为粒子可能在 d4 层阻停或穿透。其他层面（d1-d2、d2-d3）只有一种 tran cut，因为它们必须穿透当前层进入下一层。

---

## extract_three_alpha — 数据提取

### 运行方式

| 参数 | 说明 |
|------|------|
| `-r, --run` | 起始 run 号（必填） |
| `-e, --end-run` | 结束 run 号（可选，默认等于 `-r`） |
| `-t, --trigger` | 触发类型：`main` 或 `t1`（必填） |
| `-c, --config` | 配置文件路径（默认 `config.toml`） |
| `-h, --help` | 打印帮助信息 |

### 示例

```bash
# 单 run 提取
./extract_three_alpha -r 1234 -t main

# run 区间提取
./extract_three_alpha -r 1234 -e 1240 -t t1
```

### 输入文件

程序以 `t0d1` 的 match tree entry 作为主循环，按 entry 号同步读取所有文件：

| 文件 | 来源程序 | 内容 |
|------|----------|------|
| `match/t0d1_<trigger>_<run>.root` | `match_dssd` | D1 层正反面匹配事件 (`DssdMatchEvent`) |
| `match/t0d2_<trigger>_<run>.root` | `match_dssd` | D2 层正反面匹配事件 |
| `match/t0d3_<trigger>_<run>.root` | `match_dssd` | D3 层正反面匹配事件 |
| `match/t0d4_<trigger>_<run>.root` | `match_dssd` | D4 层正反面匹配事件 |
| `beam/beam_<trigger>_<run>.root` | `sort_beam` | 束流粒子分类：`is_14O` / `is_13N` / `is_12C` |
| `track/ppac_<trigger>_<run>.root` | `track_ppac` | PPAC 径迹：靶点坐标 (`target_x`, `target_y`)、方向余弦 (`dir_x`, `dir_y`) |
| `calibration/t0_<calib_run>.txt` | `calibrate_t0` | T0 探测器四层（d1-d4）的线性刻度参数：每行 `index p0 p1` |

### 输出文件

- 文件名：`extract_three_alpha_<trigger>_<run_start>_<run_end>.root`
- 保存路径：`<workspace>/three_alpha/`

#### TTree `tree` 的 branch 列表

所有能量已刻度（MeV），三个 alpha 粒子按 `et`（修正总能量）从大到小排序后分别标记为 He1、He2、He3。

| Branch | 类型 | 说明 |
|--------|------|------|
| `run_number` | `int` | 事件所在 run 号 |
| `entry` | `Long64_t` | 事件在对应 run 中的 entry 号 |
| **能量信息** |||
| `e1_He1 / e1_He2 / e1_He3` | `double` | 三个 alpha 在 D1 层的沉积能量 (MeV) |
| `e2_He1 / e2_He2 / e2_He3` | `double` | 三个 alpha 在 D2 层的沉积能量 (MeV) |
| `e3_He1 / e3_He2 / e3_He3` | `double` | 三个 alpha 在 D3 层的沉积能量 (MeV) |
| `e4_He1 / e4_He2 / e4_He3` | `double` | 三个 alpha 在 D4 层的沉积能量 (MeV) |
| `et_He1 / et_He2 / et_He3` | `double` | 三个 alpha 的修正总动能 (MeV)（阻停=各层之和，穿透=入射能量反推值） |
| `e_sum_He1 / e_sum_He2 / e_sum_He3` | `double` | 三个 alpha 四层沉积能量之和 (MeV)，即 e1+e2+e3+e4 |
| `stop_He1 / stop_He2 / stop_He3` | `bool` | 三个 alpha 是否在 D3-D4 stop cut 内（true=阻停，false=穿透） |
| **束流信息** |||
| `is_14O / is_13N / is_12C` | `bool` | 束流粒子分类（来自 `sort_beam`） |
| **PPAC 径迹信息** |||
| `ppac_valid` | `bool` | PPAC 径迹是否有效 |
| `target_x / target_y` | `double` | PPAC 径迹在靶点的坐标 (mm) |
| `dir_x / dir_y` | `double` | PPAC 径迹的方向余弦 |
| **D2 位置信息** |||
| `t0d2_He1_x / t0d2_He1_y / t0d2_He1_z` | `double` | He1 在 D2 上的击中位置 (mm) |
| `t0d2_He2_x / t0d2_He2_y / t0d2_He2_z` | `double` | He2 在 D2 上的击中位置 (mm) |
| `t0d2_He3_x / t0d2_He3_y / t0d2_He3_z` | `double` | He3 在 D2 上的击中位置 (mm) |
| **角度信息** |||
| `theta_beam` | `double` | 束流与 z 轴夹角（度），由 PPAC dir_x/dir_y 计算 |
| `phi_beam` | `double` | 束流方位角（度），atan2(dir_y, dir_x) |
| `theta_He1 / theta_He2 / theta_He3` | `double` | 三个 alpha 出射方向与束流方向的夹角（度） |
| `opening_He1_He2 / opening_He1_He3 / opening_He2_He3` | `double` | 三个 alpha 粒子两两之间的张角（度） |

#### 输出的 TH2D

| 名称 | x 轴 | y 轴 | 说明 |
|------|------|------|------|
| `h_e1_He1_e2_He1` | E2_He1 | E1_He1 | He1 D1-D2 能量关联 |
| `h_e2_He1_e3_He1` | E3_He1 | E2_He1 | He1 D2-D3 能量关联 |
| `h_e3_He1_e4_He1` | E4_He1 | E3_He1 | He1 D3-D4 能量关联 |
| `h_e1_He2_e2_He2` ~ `h_e3_He2_e4_He2` | 同上 | 同上 | He2 三层能量关联 |
| `h_e1_He3_e2_He3` ~ `h_e3_He3_e4_He3` | 同上 | 同上 | He3 三层能量关联 |

---

### 事件筛选流程（详细）

```
程序启动：
  │
  ├── 从 config.toml 加载配置，读取 D1-D4 厚度（um）
  ├── 加载 4 个 TCutG：d12_tran, d23_tran, d34_stop, d34_tran
  │
  ├── 打印硅探测器总厚度和各层厚度
  │
  ├── [验证] 计算 50 MeV 沉积能量对应的入射能量，打印结果
  │     （用于用户验证穿透能量反推的 catima 计算是否合理）
  │
  └── 遍历 run_start 到 run_end
        │
        ├── 跳过 jump_run（IsJumpRun）
        │
        ├── 打开 match/t0d1~t0d4 的 root 文件 → 读取 DssdMatchEvent
        ├── 打开 beam 文件 → 读取 is_14O / is_13N / is_12C
        ├── 打开 ppac track 文件 → 读取 PpacTrackEvent
        ├── 打开 calibration/t0_<calib_run>.txt → ReadThreeAlphaCalibration
        │
        └── 以 d1 entry 为主循环：
              │
              ├── 1. ppac_valid 检查
              │     → ppac_event.valid == 0 → skip
              │
              ├── 2. PassThreeAlphaCut 基础筛选
              │     → d1.num != 3 → skip
              │     → d2.num != 3 → skip
              │     → d3.num != 3 → skip
              │     → d4.num != 3 → skip (必须四层均等于3)
              │
              ├── 3. ClassifyThreeAlpha 径迹匹配 + PID cut 验证
              │     │
              │     ├── 3.1 D2 排序：将 d2 的三个 hit 按能量从大到小排序
              │     │       得到 d2_He[0], d2_He[1], d2_He[2]
              │     │
              │     ├── 3.2 跨层匹配（AssignThree，对 d1/d3/d4 分别执行）
              │     │       对于 d1 的三个 hit {0,1,2} 和 d2 的三个 hit：
              │     │         遍历 3! = 6 种排列，计算距离平方和：
              │     │           Σ_{i=0,1,2} [ (x_l - x_{d2_i})² + (y_l - y_{d2_i})² ]
              │     │         取距离平方和最小的排列作为匹配方案
              │     │       这样 d1/d3/d4 的三个 hit 分别与 d2 的三个 hit 一一对应
              │     │
              │     └── 3.3 PID cut 验证（对三个粒子逐个检查，必须全部通过）
              │           ├── cal_d1_d2_tran_cut → (E2_cal, E1_cal) 必须在 cut 内
              │           ├── cal_d2_d3_tran_cut → (E3_cal, E2_cal) 必须在 cut 内
              │           └── cal_d3_d4_stop/tran_cut →
              │                 (E4_cal, E3_cal) 必须在 stop cut 或 tran cut 内（二者至少满足其一）
              │          任一粒子任一 cut 不通过 → 整个事件舍弃
              │
              ├── 4. 能量计算
              │     │
              │     ├── 对每个粒子：e_sum = e1 + e2 + e3 + e4
              │     │
              │     ├── 判定 stop/tran：
              │     │     (E4, E3) 在 d34_stop_cut 内 → stop 粒子
              │     │       → et = e_sum（直接相加）
              │     │     (E4, E3) 在 d34_tran_cut 内 → tran 粒子
              │     │       → et = FindIncidentEnergy(e_sum, total_thickness_um)
              │     │         （通过 LostEnergyCalculator + 二分法，在 [e_sum, e_sum×5,...]
              │     │           区间内 80 次迭代，找到入射能量使穿透四层硅的损失恰好为 e_sum）
              │     │
              │     └── 最终排序：将三个 alpha 按 et 从大到小重新排序
              │           He1 = et 最大，He2 = et 中等，He3 = et 最小
              │
              ├── 5. 角度计算（基于 D2 位置和 PPAC 靶点）
              │     ├── beam_dir = (dir_x, dir_y, 1.0) → theta_beam, phi_beam
              │     ├── dir_He = (t0d2_x - target_x, t0d2_y - target_y, t0d2_z - 0)
              │     └── theta_He = AngleBetween(dir_He, beam_dir)    (散射角)
              │         opening = AngleBetween(dir_HeA, dir_HeB)     (张角)
              │
              └── 6. 填充输出 TTree 和直方图
```

---

### 核心算法说明

#### 三粒子跨层匹配（AssignThree）

由于 T0 DSSD 每层有三个 hit，但无法直接判断哪个 hit 属于同一个粒子，需要跨层进行径迹关联。算法流程：

1. 将 D2 的三个 hit 按能量从大到小排序，得到基准顺序 d2_He[0], d2_He[1], d2_He[2]
2. 对 D1（同理对 D3、D4），其有三个 hit {0, 1, 2}，有 3! = 6 种排列方式
3. 对每种排列 P，计算与 D2 基准的总距离平方和：Σ (x 差² + y 差²)
4. 选取总距离最小（即 x-y 位置最接近）的排列
5. 最终每层输出三个 hit 的索引，分别匹配到 d2 的三个 hit

这种方法基于物理直觉：同一个粒子穿过各层 DSSD 时，其 x-y 坐标应最接近，而不是通过靶点投影来确定。

#### 穿透能量反推（FindIncidentEnergy）

当粒子穿透 D4（即 d3-d4 在 tran cut 区域），无法直接测量总动能。利用 catima 能量损失计算进行反推：

- **输入**：四层总沉积能量 `e_total`（已刻度，MeV）
- **目标**：找到入射动能 E_true，使粒子穿透 D1-D4 四层硅时沉积的能量恰好等于 e_total
- **方法**：
  1. 创建 `LostEnergyCalculator(2, 4, Silicon, total_thickness_um)`，计算给定入射动能穿过给定厚度的能量损失
  2. 二分法搜索：下界 lo = e_total，上界 hi = e_total × 5（逐步倍增直到 E_lost(hi) < e_total）
  3. 在 [lo, hi] 区间内 80 次二分迭代，收敛误差 ≪ 0.1 MeV

---

## GUI_three_alpha — 交互式可视化

### 运行方式

```bash
./GUI_three_alpha -c config.toml
```

### 窗口总览

启动后共 **5 个窗口**（1 个主窗口 + 4 个独立弹出窗口），每个窗口尺寸均为 1700×900，排列如下：

```
(50,50)         (1770,50)        (3490,50)
┌──────────┐    ┌──────────┐     ┌──────────┐
│  Main    │    │  He1 PID │     │  He3 PID │
│ 2×2 分析 │    │  2×2     │     │  2×2     │
└──────────┘    └──────────┘     └──────────┘
                ┌──────────┐     ┌──────────┐
                │  He2 PID │     │ 8BePairs │
                │  2×2     │     │  2×2     │
                └──────────┘     └──────────┘
                (1770,1010)      (3490,1010)
```

#### 窗口 1 — 主窗口 (Main)

**控制栏**（顶部水平条）：

| 控件 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `[All] [14O] [13N] [12C]` | 互斥按钮 | All | 束流过滤，选中后自动重绘 |
| `Run:` | 标签 | — | — |
| `[0] - [100]` | 整数输入 | 0, 100 | run 号范围 |
| `[Draw]` | 按钮 | — | 加载数据并绘图（首次）或重绘 |
| `n1:` | 标签 | — | — |
| `[300]` | 整数输入 | 300 | ⁸Be 一维谱分 bin 数（范围 10-10000） |
| `n2:` | 标签 | — | — |
| `[300]` | 整数输入 | 300 | ¹²C 一维谱分 bin 数（范围 10-10000） |
| `x1:` | 实数输入 | -0.5 | ⁸Be 门控下限（范围 -30~30，可负数） |
| `x2:` | 实数输入 | 0.5 | ⁸Be 门控上限（范围 -30~30，可负数） |

> 修改 n1/n2/x1/x2 后按 Enter 或点击焦点外自动触发重绘。

**分析面板**（占据窗口主体，内嵌 2×2 布局）：

```
┌─────────────────────────┬──────────────────────────────┐
│ pad (1,1)               │ pad (1,2)                    │
│ ⁸Be Excitation (1D)     │ ¹²C vs ⁸Be Excitation (2D)   │
│ n1 bins, -1 ~ 29 MeV    │ x = Ex(¹²C), n2, -1~29      │
│ 两条红色竖线：x=x1, x=x2│ y = Ex(⁸Be), n1, -1~29      │
├─────────────────────────┼──────────────────────────────┤
│ pad (2,1)               │ pad (2,2)                    │
│ ¹²C Excitation (1D)     │ ¹²C Excitation (gated)       │
│ n2 bins, -1 ~ 29 MeV    │ 仅 ⁸Be ∈ [x1, x2] 的事例    │
│                         │ n2 bins, -1 ~ 29 MeV         │
└─────────────────────────┴──────────────────────────────┘
```

- **pad (1,1) — ⁸Be 激发能谱**：选取开角最小的一对 αα 重建 ⁸Be，纵轴自动范围。叠加两条红色竖直线 x=x1 和 x=x2，高度为直方图最大值的 1.15 倍，用于可视化门控范围。
- **pad (1,2) — 二维关联图**：x 轴为 Ex(¹²C)，y 轴为 Ex(⁸Be)。每格像素尺寸 = (30+1)/n2 × (30+1)/n1。
- **pad (2,1) — ¹²C 激发能谱**：三个 α 四动量之和重建 ¹²C。
- **pad (2,2) — gated ¹²C 谱**：仅当最小开角对的 ⁸Be 激发能落在 [x1, x2] 内时才填充。

**状态栏**（底部）：显示当前进度和状态信息。

#### 窗口 2-4 — He1/He2/He3 PID 窗口

三个独立弹出窗口，窗口标题分别为 `He1 PID`、`He2 PID`、`He3 PID`，每个内部 2×2 布局：

```
┌──────────────────┬──────────────────┐
│ E1 vs E2 (colz)  │ E2 vs E3 (colz)  │
│ 1000×1000 bins   │ 1000×1000 bins   │
│ x: 0~400 MeV     │ x: 0~350 MeV     │
│ y: 0~25 MeV      │ y: 0~400 MeV     │
├──────────────────┼──────────────────┤
│ E3 vs E4 (colz)  │ E_true vs θ      │
│ 1000×1000 bins   │ 50×100 bins      │
│ x: 0~250 MeV     │ x: 0~50 deg      │
│ y: 0~300 MeV     │ y: 0~200 MeV     │
└──────────────────┴──────────────────┘
```

每个窗口中：
- **E1-E2 / E2-E3 / E3-E4**：相邻硅层的能量关联图，用 colz 绘制。各层可区分阻停（斜线带）和穿透（横向聚集）粒子。
- **E_true vs θ**：修正总动能与束流散射角的二维分布。

#### 窗口 5 — 8Be Pairs 诊断窗口

独立弹出窗口，标题 `8Be Pairs`，2×2 布局，专门显示三种 αα 组合的 ⁸Be 激发能谱。分 bin 和范围跟随主界面的 **n1** 参数（默认 300, -1~29 MeV）。

```
┌─────────────────────────┬─────────────────────────┐
│ ⁸Be Exc: α₁α₂ pair     │ ⁸Be Exc: α₁α₃ pair     │
│ n1 bins, -1 ~ 29 MeV   │ n1 bins, -1 ~ 29 MeV   │
├─────────────────────────┼─────────────────────────┤
│ ⁸Be Exc: α₂α₃ pair     │        (empty)          │
│ n1 bins, -1 ~ 29 MeV   │                         │
└─────────────────────────┴─────────────────────────┘
```

- **α₁α₂**：He1 + He2 的四动量之和重建 ⁸Be
- **α₁α₃**：He1 + He3 的四动量之和重建 ⁸Be
- **α₂α₃**：He2 + He3 的四动量之和重建 ⁸Be
- 三种组合各自独立计算，不受主窗口"最小开角"选择策略的影响。该窗口用于诊断不同 αα 配对下 ⁸Be 激发能峰是否一致。
- 右下角 pad 预留为空。

### 束流过滤逻辑

四个按钮为互斥切换模式，再次点击同一按钮取消选中：
- `All`：显示所有事件（默认）
- `14O`：仅显示 `is_14O == true` 的事件
- `13N`：仅显示 `is_13N == true` 的事件
- `12C`：仅显示 `is_12C == true` 的事件

切换束流过滤后自动重绘所有窗口。

### 激发能不变质量计算

不需要束流能量，仅需三个 alpha 的四动量：

1. 从输出 tree 读取各 alpha 的修正总动能 `et` 和 D2 位置/靶点坐标
2. 构建四动量：
   - 方向：从靶点指向 D2 击中位置（归一化）
   - 能量：E_tot = et + m(⁴He)
   - 动量模值：√(E_tot² − m(⁴He)²)
3. **⁸Be 重建**：选取三个 alpha 中两两开角最小的一对，P_8Be = P_α[a] + P_α[b]，Ex = M_inv(8Be) − m(⁸Be_gs)
4. **¹²C 重建**：P_12C = P_α[0] + P_α[1] + P_α[2]，Ex = M_inv(12C) − m(¹²C_gs)
5. 静止质量来源于 AME2020 数据表（`nuclear_data.h` 中的 `GetMass` 函数）

------------此行以上为说明文档--------------------------------------------------

------------此行以下为提示词--------------------------------------------------
目前需要构建一个新的程序，用于提取T0探测器的三alpha事件。程序的总体和10C_4He，d_6Li程序有相同的架构，只是在提取事件时，需要根据三alpha事件的特征进行提取。你首先应该读取，这两个程序的源文件和说明文档，了解架构。目前我已经创建了相关的文件夹，具体文件则需要你根据需要进行创建。
- src/brill/bin/three_alpha/：用于创建extract_three_alpha.cpp，GUI_three_alpha.cpp这两个主程序，extract_three_alpha.cpp用于提取三alpha事件，GUI_three_alpha.cpp用于创建GUI分析界面
- src/brill/include/three_alpha/: 用于创建所需要的头文件，用于在extract_three_alpha.cpp和GUI_three_alpha.cpp中进行调用
- src/brill/src/three_alpha/: 用于创建所需要的源文件，用于在extract_three_alpha.cpp和GUI_three_alpha.cpp中进行调用

数据筛选条件：

1. 我们需要筛选三个alpha的事件，同时只考虑d1,,d2,d3,d4四层硅条，因为s1无法分辨多重性事件
2. 事件类型 d1.num = d2.num = d3.num = 3, d4.num = k, 这里不对k进行要求。另一个条件要求 ppac_valid，即有靶点
3. 然后对粒子径迹进行分裂，将d2上的三个粒子的位置和靶点进行连线，然后预测d1,d3,d4上的理论位置，然后把d1,d3,d4上的粒子的实际位置，按照
最小距离进行匹配。即d2_p1,和d1_p2,d3_4的距离最近，就把这三个粒子匹配到一起，进行分类为三个可能的粒子。
4. 然后进行判断，(d1_He1_e, d2He1_e)是否位于cal_d1_d2_tran_He4_cut.C中
                (d2He1_e, d3He1_e)是否位于cal_d2_d3_tran_He4_cut.C中
                (d3He1_e, d4He1_e)是否位于cal_d3_d4_tran_He4_cut.C/cal_d3_d4_stop_He4_cut.C中（cut我自己创建，和其他程序的cut在同一个文件夹下）

5. 如果判断为三alpha事件，就提取下这个事件的特征，包括事件时间，事件位置，事件能量等。同时根据该alppha粒子
    是否位于d3-d4的阻停cut还是穿透区域，阻停区域和穿透区域的能量计算方法不同。如果阻停直接相加即可，如果穿透，则需要使用理论曲线进行计算。

6. 仍然配置多个run ，main/t1运行方式


以上是extract_three_alpha.cpp的执行计划。

对于GUI_three_alpha.cpp，我们需要创建一个GUI界面，用于分析三alpha事件的特征。

1. GUI界面主窗口和d_Li6类似，配置all/14O/12C/13N的束流情况，run区间，draw按钮。这几个就够了
2. 主窗口绘制 He1的 d1e-d2e，d2e-d3e，d3e-d4e的pid分布直方图。和E-theta分布
3. 第二个窗口绘制 He2的 d1e-d2e，d2e-d3e，d3e-d4e的pid分布直方图。和E-theta分布
4. 第三个窗口绘制 He3的 d1e-d2e，d2e-d3e，d3e-d4e的pid分布直方图。和E-theta分布
5. 第四个窗口为2*2模式，第一张图绘制8Be的激发能谱（不变质量法），选取夹角最小的两个alpha事件进行绘制
第二张图绘制刚才的Ex_8Be(夹角最小的两个alpha)：Ex_12C（三个alpha的不变质量谱）。剩下两个先不绘制



我推翻刚才的一个大前提，我现在要求d4_num==3
建议：
1. 同意
2.按照d4_num==3,直接处理
3. 穿透区域的能量计算你仍然需要得到E_total，然后E_total为穿透四层硅沉积的能量，然后你需要计算E_true,即该粒子穿透四层硅沉积的能量为E_total
4. 对

疑问：
4.不需要
5. 12C衰变为3alpha，的不变质量法计算不需要束流能量，4He没有激发态




建议：

1.可以2.可以3.同意
4. 同意，我们的角度计算和之前一样仍然是基于d2的x,y位置

疑问：

1.精确。没错运行extract_three_alpha.cpp时候第一行打印一个 在d1,d2,d3,d4上沉积总能量为50MeV的alpha粒子的实际能量为多少的结果，以便于我验证
2.先按d2能量，最后再按总能量重新排一遍，从大到小
3.对，需要这些判定均通过
4.必须全部通过



---以下部分为你的执行计划--------------------------------------------------