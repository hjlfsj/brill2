# extract_d_Li6 说明文档

## 1. 概述

`extract_d_Li6` 是 d+6Li 物理分析的数据提取程序，从 track 和 match 文件中读取探测器数据，按筛选条件提取事件，输出包含命中事件和 2D 能量相关图的 ROOT 文件。

## 2. 目录结构

| 路径 | 用途 |
|------|------|
| `src/brill/bin/d_6Li/extract_d_Li6.cpp` | 主程序 |
| `src/brill/include/d_6Li/extract.h` | 筛选条件函数声明 |
| `src/brill/src/d_6Li/extract.cpp` | 筛选条件函数实现 |

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

与 `GUI_track` 使用相同的输入文件：

| 文件 | 来源 | 内容 |
|------|------|------|
| `track/ppac_<trigger>_<run>.root` | `track_ppac` | `PpacTrackEvent` |
| `match/t0d1_<trigger>_<run>.root` | `match_dssd` | `DssdMatchEvent` |
| `match/t0d2_<trigger>_<run>.root` | `match_dssd` | `DssdMatchEvent` |
| `match/t0d3_<trigger>_<run>.root` | `match_dssd` | `DssdMatchEvent` |
| `match/t0d4_<trigger>_<run>.root` | `match_dssd` | `DssdMatchEvent` |

## 5. 输出文件

- 文件名：`extract_d_Li6_<trigger>_<run_start>_<run_end>.root`
- 保存路径：`config.toml` 中 `paths.d_Li6` 指定的目录

### 输出内容

#### TTree：`tree`

| Branch | 类型 | 说明 |
|--------|------|------|
| `run_number` | `int` | 事件所在 run 号 |
| `entry` | `Long64_t` | 事件在对应 run 中的 entry 号 |

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
track_ppac → track/ppac_*.root
match_dssd → match/t0d{d}_*.root
    │
    ▼
extract_d_Li6
    ├── 读取各 run 的 PPAC track + DSSD match 文件
    ├── 逐 event 应用 PassD6LiCut 筛选条件
    ├── 通过的事件写入 output_tree (run_number, entry)
    └── 填充 TH2D 能量相关图
    │
    ▼
d_Li6/extract_d_Li6_*.root
```