


vibe coding order


1. d_6Li

在brill2/src/brill/bin/ 目录下,创建d_6Li目录，用于存在d_6Li物理分析主程序相关的代码
在brill2/src/brill/include/ 目录下,创建d_6Li目录，用于存在d_6Li物理分析相关的头文件
在brill2/src/brill/src/ 目录下,创建d_6Li目录，用于存在d_6Li物理分析相关的源文件
- 在brill2/src/brill/bin/d_6Li目录下,创建extract_d_Li6.cpp文件，用于提取d_6Li物理分析相关的数据。该文件所依赖的
头文件和源文件在src和include目录下自行创建。extract_d_Li6.cpp文件作为主程序，运行方式为：
Usage:
extract_d_Li6 [OPTION...]

  -h, --help             Print help information.
  -r, --run run          Start run number.
  -e, --end-run run      End run number.
  -t, --trigger trigger  Trigger type.
  -c, --config file      Config file path. (default: config.toml)

生成的文件命名为extract_d_Li6_trigger_runstart_endrun.root, 保存到配置文件中的 d_Li6 目录下

- 生成的root文件中包含两个branch，一个为筛选出的事件的run_number,另一个为该事件在它的run中的entry号。

- 输入文件和GUI_track中的文件一致,包含ppac的track和dssd的match。需要你灵活调用已有的函数。

- 输出文件中存几张TH2D：1.d1的第一个hit的能量（y轴）：d2的第一个hit的能量（x轴）
                    2. d2的第一个hit的能量（y轴）：d3的第一个hit的能量（x轴）
                    3. d3的第一个hit的能量（y轴）：d4的第一个hit的能量（x轴）

- 目前设置简单的筛选条件：d1,d2的hitnum==2，d3,d4的hitnum==1。单独使用一个文件，函数用来设置条件，目前的条件只是用来调试，之后会使用非常复杂的条件！


1. d_6Li_cut

现在需要加入另一个cut:cal_d2_d3_10C_cut.C.
在目前(e3,e4)已经在d3d4_cut的情况下，继续加入条件：
- abs (x_e3 - x_e4) < 2 && abs (y_e3 - y_e4) < 2 ，即x,y坐标在2个单位内。
- （e2_1,e3）和（e2_2,e3）有且只有一个在cal_d2_d3_10C_cut中,并且在cut中的事件满足x,y坐标在2个单位内。
- 如果满足以上条件则满足条件的视为 e2_10C，另一个为e2_6Li
- 此时判断e1上的两个hit，选择和e2_10C的x,y距离的平方较小的那个hit，作为e1_10C，另一个作为e1_6Li

- 最后创建TH2D：1. e1_10C的能量（y轴）：e2_10C的能量（x轴）
                    2. e1_6Li的能量（y轴）：e2_6Li的能量（x轴）
                    3. e2_10C的能量（y轴）：e3_10C的能量（x轴）
                    4. e3_10C的能量（y轴）：e4_10C的能量（x轴）

- 尽量加入可扩展的函数，之后还会加入更复杂的逻辑条件。


