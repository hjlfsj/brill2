# sort_beam 说明文档

## 1. 程序功能

从 `ingot` 目录读取 beam 文件，按 TOF 分布寻峰（14O、13N、12C），高斯拟合后标记每个事件属于哪个束流粒子，输出到 `beam` 文件夹。

## 2. 命令行参数

| 参数 | 说明 |
|------|------|
| `-r, --run` | run 号（必填，单 run 运行） |
| `-t, --trigger` | 触发类型（必填） |
| `-c, --config` | 配置文件路径（默认 `config.toml`） |
| `-h, --help` | 打印帮助信息 |

## 3. 运行示例

```bash
sort_beam -r 57 -t t1
```

## 4. 工作流

`sort_beam` 与 `check_tof_and_beam` 存在相互依赖关系：

- `sort_beam` 负责 TOF 寻峰和粒子 PID 分类（14O/13N/12C），输出 sorted beam 文件
- `check_tof_and_beam` 读取 sorted beam 文件中的 `v_14O` 标记筛选 14O 事件，通过相对论运动学计算 TOF 测量偏置 $\Delta t$，输出 `check_tof_and_beam_<run>.txt`
- `sort_beam` 再次运行时读取 $\Delta t$，为 14O 事件计算准确的束流能量 $E$

**执行顺序**：先运行 `check_tof_and_beam`，再运行 `sort_beam` 即可。

对于**首次处理**一个 run（尚无 sorted beam 文件）：
1. 先跑一次 `sort_beam`（此时无 $\Delta t$，仅产出 PID 分类，$E = -1$）
2. 跑 `check_tof_and_beam`（利用 PID 结果做相对论校准，得到 $\Delta t$）
3. 再跑一次 `sort_beam`（此时有 $\Delta t$，正确计算 14O 能量）

## 5. 输入文件

| 文件 | 来源 | 内容 |
|------|------|------|
| `ingot/beam_<trigger>_<run>.root` | forge beam | `BeamEvent` |

## 6. 输出文件

- 文件名：`beam_<trigger>_<run>.root`
- 保存路径：`beam/` 目录

## 7. 输出内容

### 7.1 TH1D

| 名称 | 内容 |
|------|------|
| `h_tof` | TOF 直方图（5000 bins, -500~500 ns），含高斯拟合曲线和标注 |

### 7.2 TCanvas

| 名称 | 内容 |
|------|------|
| `c1` | TOF 直方图 + 峰区竖线 + 粒子标注，用于可视化检查 |

### 7.3 TTree

| Branch | 类型 | 含义 |
|--------|------|------|
| `14O_valid` | `Bool_t` | TOF 落在 14O 峰的 5σ 范围内 |
| `13N_valid` | `Bool_t` | TOF 落在 13N 峰的 5σ 范围内 |
| `12C_valid` | `Bool_t` | TOF 落在 12C 峰的 5σ 范围内 |
| `tof` | `Double_t` | 原始 TOF 值 (ns) |
| `E` | `Double_t` | 14O 束流能量 (MeV)，仅当 `14O_valid==true` 且 `check_tof_and_beam_<run>.txt` 存在时有效，否则为 -1 |

E 的计算依赖 `beam/check_tof_and_beam_<run>.txt` 中的 $\Delta t$ 参数：

$$
E = mc^2\left(\frac{1}{\sqrt{1-\beta^2}} - 1\right), \quad
\beta = \frac{L}{(t_{\text{tof}} - \Delta t) \cdot c}
$$

其中 $m = 13044.84\ \mathrm{MeV}/c^2$, $L = 17000\ \mathrm{mm}$, $c = 299.792458\ \mathrm{mm/ns}$。若 txt 文件不存在，打印警告且所有事件的 `E = -1`。

## 8. 算法流程

1. 读取单个 run 的 `beam_<trigger>_<run>.root`，用 `tree->Draw("tof>>h", "valid")` 生成 TOF 直方图
2. 用 `TSpectrum` 寻峰，取高度最高的 3 个峰，按 x 位置从左到右排序
3. 峰位分配：最高峰 = 14O，右邻 = 13N，最右 = 12C
4. 对每个峰在 ±3 ns 范围内做高斯拟合，得 mean 和 sigma
5. 计算 5σ 区间，检查相邻峰是否重叠：若重叠则取峰位中点作为分割边界，确保三个 bool 值互斥（最多一个为真）
6. 在直方图上绘制：
   - 红色虚线：14O 的 5σ 边界
   - 蓝色虚线：13N 的 5σ 边界
   - 绿色虚线：12C 的 5σ 边界
   - 黑色虚线：峰间分割线
   - 文字标注：14O、13N、12C
7. 遍历同一 run 的所有事件，对每个事件判断 TOF 落在哪个区间，填入对应的 bool 变量

## 9. 库文件

束流排序的核心逻辑已提取为独立库，供后续程序调用。

### 9.1 文件结构

| 文件 | 用途 |
|------|------|
| `include/event/beam/beam_sort.h` | 头文件：`BeamSortResult` 结构体、函数声明 |
| `src/event/beam/beam_sort.cpp` | 源文件：函数实现 |

### 9.2 数据结构

```cpp
struct BeamSortResult {
    double mean_14O, sigma_14O, x_low_14O, x_high_14O;
    double mean_13N, sigma_13N, x_low_13N, x_high_13N;
    double mean_12C, sigma_12C, x_low_12C, x_high_12C;
};
```

### 9.3 API

| 函数 | 返回值 | 说明 |
|------|--------|------|
| `SortBeamTOF(TH1D *h, BeamSortResult &result)` | `int` | 寻峰→高斯拟合→重叠检测，返回 0 成功 |
| `SetupInputSortBeamTree(TTree *tree, bool &v_14O, bool &v_13N, bool &v_12C)` | `void` | 从已有 tree 读取三个 branch |
| `SetupOutputSortBeamTree(TTree *tree, bool &v_14O, bool &v_13N, bool &v_12C)` | `void` | 创建输出 TTree 的三个 branch |
| `ClassifyBeam(double tof, bool valid, const BeamSortResult &result, bool &v_14O, bool &v_13N, bool &v_12C)` | `void` | 分类单个事件，三值互斥 |

### 9.4 调用示例

```cpp
#include "include/event/beam/beam_sort.h"

// 已有 TOF 直方图 h_tof
brill::BeamSortResult sort_result;
if (brill::SortBeamTOF(h_tof, sort_result) == 0) {
    // 创建输出 tree
    bool v_14O, v_13N, v_12C;
    brill::SetupOutputSortBeamTree(output_tree, v_14O, v_13N, v_12C);

    // 遍历事件分类
    for (Long64_t i = 0; i < n_entries; ++i) {
        beam_tree->GetEntry(i);
        brill::ClassifyBeam(beam_event.tof, beam_event.valid,
                            sort_result, v_14O, v_13N, v_12C);
        output_tree->Fill();
    }
}
```