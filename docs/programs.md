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
| `adjust_ppac_position` | PPAC 位置偏移矫正 | `./adjust_ppac_position -r 60 -t main -n 5` |
| `track_t0` | T0 粒子径迹与 PID | `./track_t0 -r 60 -t main` |
| `rebuild_t0` | T0 粒子重建 | `./rebuild_t0 -r 60 -t main` |
| `rebuild_2alpha_2p` | 2α+2p 反应重建 | `./rebuild_2alpha_2p -r 60 -e 70 -t main` |
| `calibrate_t0` | T0 能量刻度 | `./calibrate_t0 -r 60 -e 70 -t main` |

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