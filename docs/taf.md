# TAFD — 扇形硅探测器（T0 Array of Forward Detectors）

## 功能概述

TAFD 是 ingot 实验中位于束流下游的前角扇形硅条探测器阵列。由 6 个独立的扇形硅探测器（TAFD0–TAFD5）组成，围绕束流轴线（z 轴）环形分布。每片探测器为双面硅条结构，正面 16 条环形（ring），背面 8 条扇形（pie），用于探测前角散射粒子的能量、位置和时间信息。

TAFD 同时承担触发探测器功能，`t1` 触发即为 TAFD 给出的触发信号。

---

## 探测器几何结构

### 整体参数

| 参数 | 值 | 说明 |
|------|-----|------|
| z 位置 | 135.0 mm | 探测器平面与 z 轴交点 |
| 平面方向 | 垂直于 z 轴 | 位于 x-y 平面 |
| 探测器圆心距 z 轴 | 34.4 mm | 各探测器圆心不在 z 轴上 |
| 探测器内径 | 32.6 mm | 以探测器自身圆心为原点 |
| 探测器外径 | 135.1 mm | 以探测器自身圆心为原点 |
| 最内侧条距 z 轴 | 67.0 mm | = 34.4 + 32.6 |
| 最外侧条距 z 轴 | 169.5 mm | = 34.4 + 135.1 |
| 距 z 轴 135 mm 处张角 | 26.4° – 51.5° | |

### 正面 Ring 条（径向）

| 参数 | 值 |
|------|-----|
| 条数 | 16（ring 0 – ring 15） |
| 单条径向宽度 | 6.40625 mm `= (135.1 - 32.6) / 16` |
| ring i 中心半径 | `32.6 + (i + 0.5) × 6.40625` mm |
| 方向 | ring 0 半径最小，ring 15 半径最大 |

### 背面 Pie 条（角度）

| 参数 | 值 |
|------|-----|
| 条数 | 8（pie 0 – pie 7） |
| 扇形总张角 | 54.4° |
| 单条张角 | 6.8° |
| 方向 | pie 0→7 顺时针（角度减小方向） |
| 对称轴 | 与探测器圆心 → z 轴连线方向一致 |

### 各探测器角度位置

| 探测器 | 圆心角度 | 位置描述 |
|--------|----------|----------|
| TAFD0 | 270° | -y 轴方向 |
| TAFD1 | 330° | 第四象限 |
| TAFD2 | 30° | 第一象限 |
| TAFD3 | 90° | +y 轴方向 |
| TAFD4 | 150° | 第二象限 |
| TAFD5 | 210° | 第三象限 |

```
                TAFD3 (90°)
                    ▲
                    │
    TAFD4 (150°)    │    TAFD2 (30°)
            ╲       │       ╱
             ╲      │      ╱
              ╲     │     ╱
               ╲    │    ╱
                ╲   │   ╱
                 ╲  │  ╱
                  ╲ │ ╱
    ────────────────●────────────────   z 轴（原点）
                   ╱│╲
                  ╱ │ ╲
                 ╱  │  ╲
                ╱   │   ╲
               ╱    │    ╲
              ╱     │     ╲
             ╱      │      ╲
    TAFD5 (210°)    │    TAFD1 (330°)
                    │
                    ▼
                TAFD0 (270°)
```

---

## 数据源和结构

### 数据文件

| 文件路径 | 说明 |
|---------|------|
| `ingot/tafd_t1_<run>.root` | t1 触发的 TAFD 原始数据 |
| `ingot/tafd_main_<run>.root` | main 触发的 TAFD 数据 |

### TTree 结构

ROOT 文件中 TTree 名为 `tree`，其 branch 对应 `brill::TafdEvent` 结构体：

| Branch | 类型 | 数组长度 | 含义 |
|--------|------|----------|------|
| `flag/I` | int | 1 | 事件标记 |
| `v[6]/O` | bool | 6 | 各 TAFD 通道是否有效 |
| `fe[6]/D` | double | 6 | 正面（ring）原始 ADC 通道（需刻度） |
| `be[6]/D` | double | 6 | 背面（pie）原始 ADC 通道（需刻度） |
| `fs[6]/I` | int | 6 | 正面条带号 (ring, 0–15) |
| `bs[6]/I` | int | 6 | 背面条带号 (pie, 0–7) |
| `ft[6]/D` | double | 6 | 正面时间 (ns) |

> TafdEvent 定义文件：[src/brill/include/event/ingot/tafd_event.h](../src/brill/include/event/ingot/tafd_event.h)

---

## 程序结构

| 路径 | 用途 |
|------|------|
| `src/brill/include/config.h` | TAFDConfig 结构体定义 |
| `src/brill/src/config.cpp` | LoadTAFD() — 从 config.toml 加载 TAFD 配置 |
| `src/brill/include/event/ingot/tafd_event.h` | TafdEvent 结构体，raw 数据读取 |
| `src/brill/src/event/ingot/tafd_event.cpp` | SetupInput/SetupOutput() — TTree branch 绑定 |
| `src/brill/include/10C+4He/tafd.h` | 位置计算 + 能量刻度 + 目标筛选 全部 API 声明 |
| `src/brill/src/10C+4He/tafd.cpp` | 上述函数的完整实现 |
| `src/brill/src/10C+4He/CMakeLists.txt` | C10He4_tafd 库构建规则 |
| `calib_vme_alpha_0255.txt` | α 源能量刻度参数 (k, b) |

---

## 位置计算 API

### 配置文件

TAFD 几何参数配置于 `config.toml` 的 `[tafd]` 段：

```toml
[tafd]
z_mm = 135.0
num_rings = 16
num_pies = 8
inner_radius_mm = 32.6
outer_radius_mm = 135.1
center_distance_mm = 34.4
fan_angle_deg = 54.4
angles_deg = [270.0, 330.0, 30.0, 90.0, 150.0, 210.0]
```

### 函数接口

```cpp
namespace brill {

// 获取 ring 条中心到探测器圆心的径向距离 (mm)
// ring: 0–15
double TAFDGetRingRadius(const TAFDConfig &config, int ring);

// 获取 pie 条中心相对于探测器对称轴的角度 (°)
// pie: 0–7，顺时针为正方向对应角度减小
// 返回正值表示在对称轴的逆时针侧
double TAFDGetPieAngle(const TAFDConfig &config, int pie);

// 获取 TAFD 条击中点在实验室坐标系中的 (x, y, z) 坐标 (mm)
// tafd_index: 0–5 (对应 TAFD0–TAFD5)
// ring: 0–15 (正面条带)
// pie:  0–7  (背面条带)
// x, y, z: 输出坐标 (mm)，z 始终等于 z_mm
void TAFDGetPosition(
    const TAFDConfig &config,
    int tafd_index,
    int ring,
    int pie,
    double &x,
    double &y,
    double &z
);

}  // namespace brill
```

### 几何计算逻辑

```
实验室坐标系：+z = 束流方向，x-y = 探测器平面

对于 TAFD[i]（圆心角度 θ_i = angles_deg[i]）：

  1. 探测器圆心坐标：
     cx = center_distance_mm × cos(θ_i)
     cy = center_distance_mm × sin(θ_i)

  2. Ring 条 r 的径向距离（从探测器圆心算起）：
     pitch = (outer_radius_mm - inner_radius_mm) / num_rings
     r = inner_radius_mm + (ring + 0.5) × pitch

  3. Pie 条 j 相对对称轴的角度：
     φ_rel = fan_angle/2 − pie_angle/2 − j × pie_angle
           = 27.2° − 3.4° − j × 6.8°
           = 23.8° − j × 6.8°

  4. 全局坐标：
     φ_global = θ_i + φ_rel
     x = cx + r × cos(φ_global)
     y = cy + r × sin(φ_global)
     z = z_mm
```

### 使用示例

```cpp
#include "include/10C+4He/tafd.h"
#include "include/config.h"

brill::AppConfig config;
brill::LoadConfig("config.toml", config);

int tafd_index = 0;   // TAFD0
int ring = 5;         // ring 5
int pie = 3;          // pie 3
double x, y, z;
brill::TAFDGetPosition(config.tafd, tafd_index, ring, pie, x, y, z);

double r = brill::TAFDGetRingRadius(config.tafd, ring);
double phi_rel = brill::TAFDGetPieAngle(config.tafd, pie);
```

### CMake 链接

```cmake
target_link_libraries(your_target PUBLIC C10He4_tafd config)
```

---

## 能量刻度 API

### 刻度文件格式

`calib_vme_alpha_0255.txt` 由 α 源（5.486 MeV, 5.805 MeV）刻度得到，公式 `E = k × ADC + b` (MeV)，按 `detector/face/strip` 三级索引：

```
# 列: detector face strip peak_low_ch peak_low_e peak_high_ch peak_high_e k b valid entries
tafd0  front  0  922.5  0.1  969.7  0.2  0.006759  -0.7496  yes  5189
tafd0  back   0  844.7  0.2  891.0  0.3  0.006896  -0.3392  yes  14851
...
tafd5  front  1  -  -  -  -  -  -  no  3758     ← 无效条
t1du   front  0  459.9  0.2  483.4  0.2  0.013598  -0.7677  yes  1909   ← 非 TAFD，自动跳过
```

### TafdCalibration 结构体

```cpp
struct TafdCalibration {
    double front_k[6][16] = {};    // k 参数，[tafd][strip]
    double front_b[6][16] = {};    // b 参数
    double back_k[6][8] = {};
    double back_b[6][8] = {};
    bool front_valid[6][16] = {};  // 有效性（排除 tafd5 strip1 等坏条）
    bool back_valid[6][8] = {};
};
```

### ReadTafdCalibration

```cpp
// 解析 calib_vme_alpha_0255.txt，自动跳过 # 注释行和非 tafd 探测器
// 返回 0 成功，-1 失败
int ReadTafdCalibration(const std::string &path, TafdCalibration &calib);
```

---

## 目标事件筛选 API

### 处理逻辑

```
对每个 TTree entry：

  for i in 0..5:                           // 遍历 6 片 TAFD
    if !valid[i]                          → 跳过
    if strip 越界                          → 跳过
    if calib 无效 (如 tafd5 strip1)        → 跳过
    fe_cal = k_front × ADC + b_front       // 正面刻度
    be_cal = k_back  × ADC + b_back        // 背面刻度
    if |fe_cal − be_cal| ≥ 0.2 MeV        → 丢弃（噪声）
    加入候选列表

  从候选中选 fe_cal 最大的 → 目标事件
  调用 TAFDGetPosition → 计算 (x, y, z)
  返回 {能量, x, y, z, tafd_index, ring, pie}
```

> 设计原则：正面和背面是同一粒子在同一条带上沉积的能量，两者应一致。差值 < 0.2 MeV 作为一致性判据；其他能量小的候选大概率是噪声，取最大能量作为目标事件。

### TafdTarget 结构体

```cpp
struct TafdTarget {
    bool found = false;       // 是否找到有效目标
    double energy = 0.0;      // 正面刻度能量 (MeV)
    double x = 0.0;           // 实验室系 x (mm)
    double y = 0.0;           // 实验室系 y (mm)
    double z = 0.0;           // 实验室系 z (mm)
    int tafd_index = -1;      // 0–5，哪个 TAFD
    int ring = -1;            // 0–15，正面 strip
    int pie = -1;             // 0–7，背面 strip
};
```

### FindTafdTarget

```cpp
// 对一个 TTree entry 执行刻度、正反面一致性筛选、选最大能量
// 若找不到有效目标，返回 target.found == false
TafdTarget FindTafdTarget(
    const TafdEvent &event,
    const TafdCalibration &calib,
    const TAFDConfig &config
);
```

---

## 完整使用示例

```cpp
#include "include/10C+4He/tafd.h"
#include "include/config.h"
#include "include/event/ingot/tafd_event.h"

#include <TFile.h>
#include <TTree.h>

int main() {
    // 1. 加载配置
    brill::AppConfig config;
    brill::LoadConfig("config.toml", config);

    // 2. 读取刻度文件
    brill::TafdCalibration calib;
    if (brill::ReadTafdCalibration("calib_vme_alpha_0255.txt", calib) != 0) {
        return -1;
    }

    // 3. 打开数据文件
    TFile f("ingot/tafd_t1_0057.root");
    TTree *tree = (TTree*)f.Get("tree");
    brill::TafdEvent event;
    brill::SetupInput(tree, event);

    // 4. 逐 entry 处理
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        auto target = brill::FindTafdTarget(event, calib, config.tafd);
        if (target.found) {
            printf("entry=%lld TAFD%d ring=%d pie=%d E=%.3f MeV pos=(%.1f, %.1f, %.1f)\n",
                   i, target.tafd_index, target.ring, target.pie,
                   target.energy, target.x, target.y, target.z);
        }
    }
    return 0;
}
```

---

## 触发与物理含义

### 触发类型

| 触发 | 含义 | TAFD 角色 |
|------|------|-----------|
| `main` | 主触发（T0 DSSD 等） | TAFD 不做触发，被动记录 |
| `t1` | TAFD 触发 | TAFD 给出触发信号——只有 TAFD 探测到粒子时才记录事件 |

### 物理意义

- **前角粒子探测**：TAFD 位于束流下游（z=135mm），覆盖实验室系极角 26.4°–51.5° 范围
- **粒子鉴别**：正面 (ring) 和背面 (pie) 能量信息 + 位置信息可用于粒子种类鉴别
- **触发系统**：t1 触发模式确保束流粒子或反应产物击中 TAFD 时才记录，适合研究前角出射道
- **位置重建**：ring（径向）+ pie（角度）交叉确定粒子击中位置，结合 z 位置可用于角度重建

---

## TAFCSI — CsI 闪烁体探测器

每块 TAFD 硅探测器背后对应两块 CsI 闪烁体，用于测量穿透 TAFD 后粒子的残余能量，提供粒子鉴别能力。

### 几何对应

| TAFD | CsI 编号 |
|------|----------|
| TAFD0 | CsI 0, 1 |
| TAFD1 | CsI 2, 3 |
| TAFD2 | CsI 4, 5 |
| TAFD3 | CsI 6, 7 |
| TAFD4 | CsI 8, 9 |
| TAFD5 | CsI 10, 11 |

### 数据文件

| 文件路径 | 说明 |
|---------|------|
| `ingot/tafcsi_t1_<run>.root` | t1 触发的 TAFCSI 原始数据 |

### TTree 结构

ROOT 文件中 TTree 名为 `tree`，其 branch 对应 `brill::CsiEvent<12>` 模板：

| Branch | 类型 | 数组长度 | 含义 |
|--------|------|----------|------|
| `flag/l` | long | 1 | 事件标记 |
| `v[12]/O` | bool | 12 | 各 CsI 通道是否有效 |
| `t[12]/D` | double | 12 | 时间 (ns) |
| `e[12]/I` | int | 12 | 能量 (ADC channel，不刻度) |

> entry 号与 `tafd_t1_<run>.root` 一一对应。

### 数据读取 API

```cpp
namespace brill {

constexpr int kMaxTafcsi = 12;

struct TafcsiTarget {
    bool found = false;
    int energy = 0;    // raw ADC，取两块 CsI 中较大的信号
    int channel = -1;  // 被选中的 CsI 通道号 (0-11)
};

// 对给定 tafd_index (0–5)，从 CsI 数据中提取对应信号
TafcsiTarget FindTafcsiTarget(const CsiEvent<kMaxTafcsi> &event, int tafd_index);

}  // namespace brill
```

**处理逻辑**：根据 `tafd_index` 定位 `CsI[2*idx]` 和 `CsI[2*idx+1]`，两者均有效取较大值，仅一个有效取该值，均无效返回 `found=false`。

### 使用示例

```cpp
#include "include/10C+4He/tafd.h"
#include "include/10C+4He/tafcsi.h"

TFile f_tafd("tafd_t1_0057.root");
TFile f_csi("tafcsi_t1_0057.root");
TTree *t_tafd = (TTree*)f_tafd.Get("tree");
TTree *t_csi  = (TTree*)f_csi.Get("tree");

TafdEvent tafd_ev;
CsiEvent<12> csi_ev;
SetupInput(t_tafd, tafd_ev, "");
SetupInput(t_csi, csi_ev, "");

for (Long64_t i = 0; i < t_tafd->GetEntries(); ++i) {
    t_tafd->GetEntry(i);
    auto target = FindTafdTarget(tafd_ev, calib, config.tafd);
    if (!target.found) continue;

    t_csi->GetEntry(i);
    auto csi = FindTafcsiTarget(csi_ev, target.tafd_index);
    if (csi.found) {
        printf("TAFD%d: E_Si=%.2f MeV, CsI%02d=%d ch\n",
               target.tafd_index, target.energy, csi.channel, csi.energy);
    }
}
```

## 探测器版本差异

| 探测器 | 厚度 | 适用场景 |
|--------|------|----------|
| TAFD0、TAFD1 | 薄 (~68 μm) | 适合低能粒子、高能量分辨 |
| TAFD2、TAFD3 | 厚 (~300 μm) | 适合穿透性较强的粒子 |
| TAFD4、TAFD5 | — | 刻度文件中 tafd5 ch1 标记为 bad |

> 厚度信息来源于 `calib_vme_alpha_0255.txt` 的 α 源刻度数据，实际厚度以实验记录为准。