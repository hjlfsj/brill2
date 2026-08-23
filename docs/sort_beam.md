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

## 4. 输入文件

| 文件 | 来源 | 内容 |
|------|------|------|
| `ingot/beam_<trigger>_<run>.root` | forge beam | `BeamEvent` |

## 5. 输出文件

- 文件名：`beam_<trigger>_<run>.root`
- 保存路径：`beam/` 目录

## 6. 输出内容

### 6.1 TH1D

| 名称 | 内容 |
|------|------|
| `h_tof` | TOF 直方图（5000 bins, -500~500 ns），含高斯拟合曲线和标注 |

### 6.2 TCanvas

| 名称 | 内容 |
|------|------|
| `c1` | TOF 直方图 + 峰区竖线 + 粒子标注，用于可视化检查 |

### 6.3 TTree

| Branch | 类型 | 含义 |
|--------|------|------|
| `14O_valid` | `Bool_t` | TOF 落在 14O 峰的 5σ 范围内 |
| `13N_valid` | `Bool_t` | TOF 落在 13N 峰的 5σ 范围内 |
| `12C_valid` | `Bool_t` | TOF 落在 12C 峰的 5σ 范围内 |

## 7. 算法流程

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

## 8. 库文件

束流排序的核心逻辑已提取为独立库，供后续程序调用。

### 8.1 文件结构

| 文件 | 用途 |
|------|------|
| `include/event/beam/beam_sort.h` | 头文件：`BeamSortResult` 结构体、函数声明 |
| `src/event/beam/beam_sort.cpp` | 源文件：函数实现 |

### 8.2 数据结构

```cpp
struct BeamSortResult {
    double mean_14O, sigma_14O, x_low_14O, x_high_14O;
    double mean_13N, sigma_13N, x_low_13N, x_high_13N;
    double mean_12C, sigma_12C, x_low_12C, x_high_12C;
};
```

### 8.3 API

| 函数 | 返回值 | 说明 |
|------|--------|------|
| `SortBeamTOF(TH1D *h, BeamSortResult &result)` | `int` | 寻峰→高斯拟合→重叠检测，返回 0 成功 |
| `SetupInputSortBeamTree(TTree *tree, bool &v_14O, bool &v_13N, bool &v_12C)` | `void` | 从已有 tree 读取三个 branch |
| `SetupOutputSortBeamTree(TTree *tree, bool &v_14O, bool &v_13N, bool &v_12C)` | `void` | 创建输出 TTree 的三个 branch |
| `ClassifyBeam(double tof, bool valid, const BeamSortResult &result, bool &v_14O, bool &v_13N, bool &v_12C)` | `void` | 分类单个事件，三值互斥 |

### 8.4 调用示例

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