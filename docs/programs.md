# 可执行程序一览

## 主程序 (bin/)

| 程序 | 功能 | 用法示例 |
|------|------|---------|
| `normalize` | DSSD 归一化 | `./normalize -r 60 -e 70 -t main` |
| `ppac_normalize` | PPAC offset 标定 | `./ppac_normalize -r 60 -t main` |
| `track_ppac` | PPAC 径迹重建 | `./track_ppac -r 60 -t main` |
| `match_dssd` | DSSD 正反面匹配 | `./match_dssd -r 60 -t main` |
| `check_match` | DSSD 匹配检查器 (GUI) | `./check_match -r 60 -t main -d t0d1` |
| `check_ppac_track` | PPAC 径迹检查器 (GUI) | `./check_ppac_track -r 60 -t main` |
| `GUI_track` | PPAC + T0 径迹查看器 (GUI) | `./GUI_track -r 60 -t main` |
| `adjust_ppac_position` | PPAC 位置偏移矫正 | `./adjust_ppac_position -r 60 -t main -n 5` |
| `adjust_t0_step1` | T0 D1/D2/D3 位置矫正（束流） | `./adjust_t0_step1 -r 60 -t main` |
| `adjust_t0_step2` | T0 D4 位置矫正（径迹外推） | `./adjust_t0_step2 -r 60 -t main` |
| `track_t0` | T0 粒子径迹与 PID | `./track_t0 -r 60 -t main` |
| `rebuild_t0` | T0 粒子重建 | `./rebuild_t0 -r 60 -t main` |
| `rebuild_2alpha_2p` | 2α+2p 反应重建 | `./rebuild_2alpha_2p -r 60 -e 70 -t main` |
| `calibrate_t0` | T0 能量刻度 | `./calibrate_t0 -r 60 -e 70 -t main` |
| `sort_beam` | 束流粒子分类 | `./sort_beam -r 57 -t main` |
| `extract_d_Li6` | d+6Li 事件提取 | `./extract_d_Li6 -r 57 -e 60 -t main` |
| `GUI_d_Li6` | d+6Li 交互式可视化 | `./GUI_d_Li6 -c config.toml` |
| `extract_10C_4He` | 10C+4He 事件提取 | `./extract_10C_4He -r 57 -e 70 -t main` |
| `GUI_10C_4He` | 10C+4He 交互式可视化 | `./GUI_10C_4He -c config.toml` |
| `GUI_pid` | DSSD PID 可视化 (GUI) | `./GUI_pid -c config.toml` |

## d_6Li/

| 程序 | 功能 | 用法示例 |
|------|------|---------|
| `extract_d_Li6` | d+6Li 数据提取 | `./extract_d_Li6 -r 57 -e 60 -t main` |

## 10C+4He/

| 程序 | 功能 | 用法示例 |
|------|------|---------|
| `extract_10C_4He` | 10C+4He 数据提取 | `./extract_10C_4He -r 57 -e 70 -t main [--d1_hit 1]` |
| `GUI_10C_4He` | 10C+4He 交互式可视化 | `./GUI_10C_4He -c config.toml` |

## estimate/

| 程序 | 功能 | 用法示例 |
|------|------|---------|
| `estimate_normalize` | 估计归一化参数 | `./estimate_normalize -c config.toml` |
| `estimate_normalize_total_energy` | 估计总能量归一化 | `./estimate_normalize_total_energy -c config.toml` |
| `estimate_t0_center` | 估计 T0 探测器中心 | `./estimate_t0_center -r 60 -e 70 -t main` |
| `estimate_t0_pid` | 估计 T0 PID 参数 | `./estimate_t0_pid -r 60 -e 70 -t main` |
| `estimate_t0_straight` | 估计 T0 直线刻度 | `./estimate_t0_straight -r 60 -e 70 -t main` |
| `estimate_t0_csi_pid` | 估计 T0 CsI PID | `./estimate_t0_csi_pid -r 60 -e 70 -t main` |

## rebuild/

| 程序 | 功能 | 用法示例 |
|------|------|---------|
| `rebuild_2alpha` | 2α 反应重建 | `./rebuild_2alpha -r 60 -e 70 -t main` |
| `rebuild_3alpha` | 3α 反应重建 | `./rebuild_3alpha -r 60 -e 70 -t main` |

## 通用选项

| 选项 | 说明 |
|------|------|
| `-r, --run` | 起始 run 号（必填） |
| `-e, --end-run` | 结束 run 号（部分程序支持） |
| `-t, --trigger` | 触发类型：`main` 或 `t1` |
| `-c, --config` | 配置文件路径，默认 `config.toml` |
| `-h, --help` | 查看帮助 |

---

## 数据处理流程

下图展示完整的数据处理管线，从原始数据到物理分析结果。特别注意三种校准参数（PPAC offset、T0 offset、T0 刻度）分别在什么阶段被使用。

### 流程图

```
ingot/ (原始数据)
  │
  ├─── [1] normalize ─────────────────────────→ normalize/ (DSSD 条带归一化参数)
  │     每 strip 独立拟合 front/back 能量关系
  │
  ├─── [2] ppac_normalize ────────────────────→ normalize/ (PPAC offset 参数)
  │     计算 delay_x/y 峰位和 position 刻度系数
  │
  └─── [3] track_ppac ────────────────────────→ track/ (PPAC 径迹)
        输入: ingot/ppac_*.root + normalize/ppac_offset_*.txt
        ★ PPAC offset 在此应用：将时间差转换为物理位置
        pos = p0 + p1 * dt + x_offset_mm

  ┌─── [4] match_dssd ────────────────────────→ match/ (DSSD 正反面匹配)
  │     输入: normalize/ (归一化 DSSD)
  │     ★ T0 位置 offset 在此应用：x_offset_mm / y_offset_mm
  │     (从 config.toml 读取，由 adjust_t0_step1/2 标定)
  │
  ├─── [5] adjust_t0_step1 ───────────────────→ normalize/ (T0 D1/D2/D3 位置 offset)
  │     输入: track/ppac_*.root + match/t0d{1,2,3}_*.root
  │     用 PPAC 径迹外推残差标定 D1/D2/D3 位置偏移
  │
  └─── [6] adjust_t0_step2 ───────────────────→ normalize/ (T0 D4 位置 offset)
        输入: match/t0d{1,2,3,4}_*.root
        用 D1-D2-D3 拟合外推残差标定 D4 位置偏移

  ┌─── [7] track_t0 ──────────────────────────→ track/ (T0 粒子径迹)
  │     输入: track/ppac_*.root + match/t0d{1,2,3,4}_*.root
  │     从 D1 播种，按 track_window 匹配 D2/D3/D4
  │
  ├─── [8] rebuild_t0 ─────────────────────────→ particle/ (T0 粒子重建)
  │     输入: track/t0_*.root + calibration/t0.txt
  │     ★ T0 能量刻度在此应用：E_cal = p0 + p1 * E_raw
  │
  └─── [9] calibrate_t0 ───────────────────────→ calibration/t0.txt
        输入: particle/t0_*.root
        基于 PID 拟合产生能量刻度系数

  ┌─── [10] sort_beam ─────────────────────────→ beam/ (束流粒子分类)
  │     输入: ingot/beam_<trigger>_*.root
  │     对 TOF 直方图寻峰、高斯拟合，标记 14O/13N/12C
  │
  └─── [11] extract_d_Li6 ─────────────────────→ d_Li6/ (自包含 D6LiEvent)
        输入: match/t0d{1,2,3,4}_*.root + beam/beam_*.root + track/ppac_*.root + calibration/t0.txt
        ★ T0 能量刻度在此应用：E_cal = p0 + p1 * E_raw
        ★ 读取 PPAC 径迹和束流分类，一并写入 D6LiEvent
        输出自包含文件，GUI 无需回查原始数据

  ┌─── [12] GUI_d_Li6 ───────────────────────── (交互式物理分析可视化)
        ★ 输入: d_Li6/extract_d_Li6_*.root + Cut/cal_d1_d2_6Li_cut.C + assets/ Lise++ 参考线
        ★ 直接读取自包含 D6LiEvent，无需回查 match/beam/calibration
        主界面 4 张原始能量相关图 + 分析画布 4 张物理分析图（E-E、E-θ、θ-θ + Lise++ 参考线）

  ┌─── [13] GUI_pid ─────────────────────────── (DSSD PID 可视化)
        ★ 输入: match/t0d{1,2,3,4}_*.root + ingot/t0s_*.root + calibration/t0.txt
        ★ 直接读取 match 文件，应用能量刻度后绘制
        主画布 4 张 hit0 关联图 + 次画布 4 张 hit1 关联图（按需绘制）
        hd 字段表示严格 hit 数（==），主画布用 hit0，次画布用 hit1
```

### 三种校准参数的使用阶段

#### PPAC offset（延迟时间 → 物理位置）

| 阶段 | 程序 | 说明 |
|------|------|------|
| **标定** | `ppac_normalize` | 计算 delay_x/y 峰位 (μ, σ) 和 position 刻度 (p0, p1)，输出 `normalize/ppac_offset_*.txt` |
| **使用** | `track_ppac` | 将 PPAC 时间差转换为物理位置：`pos = p0 + p1 * dt + x_offset_mm`，输出 `track/ppac_*.root` |

PPAC offset 是**一次性**的：在 `track_ppac` 中转换后，后续所有程序只使用 track 文件中的物理坐标，不再接触原始时间。

#### T0 位置 offset（探测器物理位置矫正）

| 阶段 | 程序 | 说明 |
|------|------|------|
| **标定** | `adjust_t0_step1` | 用 PPAC 束流径迹外推，标定 D1/D2/D3 的 x/y 偏移，输出 `normalize/t0_offset_*.root` |
| **标定** | `adjust_t0_step2` | 用 D1-D2-D3 拟合外推，标定 D4 的 x/y 偏移，输出 `normalize/t0_offset_d4_*.root` |
| **使用** | `match_dssd` | 从 `config.toml` 读取 `x_offset_mm / y_offset_mm`，在计算物理位置时应用 |

T0 位置 offset 是**迭代**的：标定结果需手动更新到 `config.toml` 后重新运行 `match_dssd`，再重新标定，直到收敛。

#### T0 能量刻度（ADC → MeV）

| 阶段 | 程序 | 说明 |
|------|------|------|
| **标定** | `calibrate_t0` | 基于 PID 的粒子鉴别，拟合每层探测器的能量刻度系数 (p0, p1)，输出 `calibration/t0.txt` |
| **使用** | `rebuild_t0` | 在粒子重建阶段应用：`E_cal = p0 + p1 * E_raw` |
| **使用** | `extract_d_Li6` | 在填充直方图前对每个 hit 的能量进行刻度，写入 D6LiEvent |
| **使用** | `GUI_d_Li6` | 直接读取 D6LiEvent 中已刻度的能量，分析画布额外施加 ppac_valid 和 D1-D2 6Li cut |
| **使用** | `GUI_pid` | 读取 match 文件，应用 T0 刻度后绘制 PID 关联图 |

T0 能量刻度是**迭代**的：`calibrate_t0` → `rebuild_t0` → `calibrate_t0` 循环，直到刻度系数收敛。`extract_d_Li6` 使用收敛后的刻度系数将能量写入 D6LiEvent，`GUI_d_Li6` 直接读取已刻度的值。

### 数据目录一览

| 目录 | 内容 | 关键文件 |
|------|------|---------|
| `ingot/` | 原始数据 | `ppac_*.root`, `beam_*.root`, `t0d1_*.root` ... |
| `normalize/` | 归一化参数 | `ppac_offset_*.txt`, `t0_offset_*.root`, `*_front.txt`, `*_back.txt` |
| `track/` | 径迹重建 | `ppac_*.root`, `t0_*.root` |
| `match/` | DSSD 匹配 | `t0d1_*.root`, `t0d2_*.root`, `t0d3_*.root`, `t0d4_*.root` |
| `particle/` | 粒子重建 | `t0_*.root` |
| `calibration/` | 能量刻度 | `t0.txt` |
| `beam/` | 束流分类 | `beam_*.root` |
| `d_Li6/` | d+6Li 分析 | `extract_d_Li6_*.root`（自包含 D6LiEvent，31 branch） |

### 典型运行顺序

```bash
# 1. DSSD 归一化
./normalize -r 57 -e 70 -t main t0d1 t0d2 t0d3 t0d4

# 2. PPAC offset 标定
./ppac_normalize -r 57 -t main

# 3. PPAC 径迹重建
./track_ppac -r 57 -t main

# 4. DSSD 匹配（首次使用 config 中的初始 offset）
./match_dssd -r 57 -t main

# 5. T0 位置标定（需用 main trigger 束流数据）
./adjust_t0_step1 -r 57 -t main
./adjust_t0_step2 -r 57 -t main
# → 将结果更新到 config.toml → 重新运行 match_dssd

# 6. 束流粒子分类
./sort_beam -r 57 -t main

# 7. T0 能量刻度
./calibrate_t0 -r 57 -e 70 -t main

# 8. d+6Li 事件提取
./extract_d_Li6 -r 57 -e 60 -t main

# 9. GUI 可视化（直接读取自包含 D6LiEvent，无需重新运行 extract）
./GUI_d_Li6 -c config.toml

# 10. DSSD PID 可视化（直接读取 match 文件，应用能量刻度）
./GUI_pid -c config.toml
```