


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


2. beam


从工作区ingot读入 beam文件，依旧按照trigger分类，输出文件放到beam文件夹中。可执行程序的结构和之前的程序保持一致，主程序放在bin文件夹下，命名为sort_beam.cpp 
 该程序读入ingot中的beam文件信息，首先绘制一个直方图tree->Draw("tof>>(5000,-500,500)","valid")，然后在直方图上寻三个峰，其中最小的峰高度至少要有最高峰的10%。然后最高峰为 14O，14O右侧为13N，13N右侧为12C（12C的峰位最低）. 
 然后在直方图上进行高斯拟合，在对应的峰位标注14O, 13N，12C以便之后检查。然后直方图存储到生成的root文件中。最后给新生成的root 文件中 设置三个Branch的bool变量，分别为 14O_valid，13N_valid，12C_valid，表示这个事件是否落在beam的峰内(以拟合结果的5sigma为标准，注意不要重合)，在直方图上用竖线表示一下分割区域。 





 Ok，我们现在需要大规模修改GUI_d_Li6进行真正的物理分析了。
首先，目前的主界面的四张图不动，其他的三个画布全部去除
然后，在rebuild目录下创建新文件rebuild_d_6Li，用于写入之后的计算函数
然后，根据extract_d_Li6_0057_0092.root文件中的信息，进行计算。首先读取cal_d1_d2_6Li_cut.C 这个cut文件，判断(e2_6Li, e1_6Li)是否在cut区域中，如果在，则对事件进行计算，首先计算6Li和10C粒子的总能量 E_6Li，E_10C，然后





3. 10C+4He
和 d_Li6一致，在bin/10C+4He目录下创建extract_10C_4He.cpp文件，用于提取10C+4He物理分析相关的数据。关于数据结构和筛选条件的源文件和头文件
放置到src/brill/src/10C+4He和src/brill/include/10C+4He目录下。整体结构和d_6Li一致。
但是需要注意的是，10C+4He的筛选条件和d_6Li不同，需要根据实际情况进行调整: 10C阻停在d3中，所以只有e1,e2,e3.4He阻停在d4或者s1中。
因此相比d_Li6，我们需要额外在ingot目录下读取t0s的文件。4He则有e1,e2,e3,e4,e5。e5即为t0s的能量。这里的能量都采用刻度之后的能量。
和d_Li6一样，我们也使用d2的位置信息计算10C和4He的角度信息。

同时在源文件中命名Pass10C_d3_4He_s1Cut()，用于筛选10C+4He的事件,该cut用于选择10C阻停在d3中，4He阻停在s1中的事件，
然后extract_10C_4He.cpp文件中目前选择该函数，之后可以扩展为其他函数。

Pass10C_d3_4He_s1Cut()的筛选条件：

- d1.num==2&&d2.num==2&&d3.num==2&&d4.num==1
- 我们知道match中的每一层hit都已经按照能量大小进行了排序(检查一下)，判断(d3hit[0],d2hit[0])是否在cal_d2_d3_stop_10C_cut(cut文件的存储路径和
d_Li6一致)。
- 判断t0s 是否有响应 (t0s_valid)
- 满足以上三个条件则为目标事件，其中d1,d2,d3,的hit[0](即能量较大的hit)为10C的e1,e2,e3.其余的为4He的e1,e2,e3,e4,e5.
- 然后即可进行角度，ppac,tof的相关计算
