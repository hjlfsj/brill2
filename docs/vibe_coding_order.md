


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


4. C12_d_6Li
接下来进行 该反应道的分析，该反应道的末态产物是6Li和两个能量和角度都很接近的 4He.
我们需要参考10C+4He， d_Li6的分析方法
/home/ribll2026/ribll2026_www/github_code/brill2/docs/10C_4He.md
/home/ribll2026/ribll2026_www/github_code/brill2/docs/d_6Li.md
然后在bin/C12_d_6Li目录下创建extract_6Li_two4He.cpp文件作为主程序，用于提取C12_d_6Li物理分析相关的数据。
其他配套的程序在src/brill/src/C12_d_6Li和src/brill/include/C12_d_6Li目录下。
extract_6Li_two4He.cpp 中的筛选条件：

- d1.num==3&&d2.num==3&&d3.num==2&&d4.num==2
- 选取d1,d2上能量最大的hit，该粒子的pid需要满足 6Li的cut(在brill2/src/brill/Cut/cal_d1_d2_6Li_cut.C)

满足以上条件的事件，即可进行角度，ppac,tof的相关计算。
这时候6Li已经确定了，因为6Li阻停在了d2中。
同时，需要对d1,d2,d3,d4上剩下的两个hit进行分类，分别组成两个粒子.此时使用距离最近的原则
d2上的剩下两个hit分别为He1,He2,然后d1，d3,d4上的两个hit分别和d2上的hit进行距离平方计算，一共两种情况，选择平方和较小的组合,并入
这两个He的信息中。同时s1中的能量为两个alpha的能量和，也应该记录下来(因为s1只有一个通道，没有分条)。
角度计算需要分别使用d2的位置信息，计算三个粒子和束流之间的夹角。同时计算两两之间的张角。
保存8张二维图：
- 1. e1_6Li的能量（y轴）：e2_6Li的能量（x轴）
- 2. e1_4He1的能量（y轴）：e2_4He1的能量（x轴）
- 3. e1_4He2的能量（y轴）：e2_4He2的能量（x轴）
- 4. e2_4He1的能量（y轴）：e3_4He1的能量（x轴）
- 5. e2_4He2的能量（y轴）：e3_4He2的能量（x轴）
- 6. e3_4He1的能量（y轴）：e4_4He1的能量（x轴）
- 7. e3_4He2的能量（y轴）：e4_4He2的能量（x轴）
- 8. e4_4He1+e4_4He2的能量（y轴）：e5的能量（x轴）



5. 新的刻度输入文件
我们需要通过match得到的文件，生成新的刻度输入文件。步骤分为两步，首先根据一定的条件绘制二维pid图
以便于手动检验cut的选择是否正确，其次是把相同条件下筛选出来的点填入TGraph中,作为后续的拟合数据。

- 首先是d1-d2的pid图，选取两层上的第一个hit，然后要求dx<2 mm, dy<2 mm,将满足条件的点填入（d1_e[0],d2_e[0]）分别填入
TH2D 和 TGraph中
- 然后是d2-d3的pid图，选取两层上的第一个hit，然后要求dx<2 mm, dy<2 mm,将满足条件的点填入（d2_e[0],d3_e[0]）分别填入
TH2D 和 TGraph中
- 最后是d3-d4的pid图，选取两层上的第一个hit，然后要求dx<2 mm, dy<2 mm,将满足条件的点填入（d3_e[0],d4_e[0]）分别填入
TH2D 和 TGraph中
- 最后是d4-s1的pid图，此时要求d1_hit=d2_hit=d3_hit=d4_hit=1,将满足条件的点填入（d4_e[0],s1_e[0]）分别填入
TH2D 和 TGraph中

其他配置和运行方式和normalize,match等程序一致，把源文件和生成的文件放进estimate目录下，命名为pre_calibration。
先陈述行动方案，并提出建议和疑问



6. 新的刻度 calibration_t0

关于程序的整体输入输出框架在/home/ribll2026/ribll2026_www/github_code/brill2/docs/calibrate_t0.md 文件中已经写明
关于算法/home/ribll2026/ribll2026_www/github_code/brill2/reference/calibrate_t0.cpp， 可以以这个旧的程序作为参考

我这里要重申一下拟合算法：

- 输入的拟合数据为TGraph，其中x轴为下一层的adc,y轴为上一层的adc。例如(d1_e,d2_e)
- 生成的delta_e的理论曲线也应该为相应粒子的 （下一层沉积能量， 上一层沉积能量）
- 然后将理论曲线拼接一个多参数的TF1，和 拼接的 TGraph中的散点进行拟合

在得到系数之后，使用系数重新计算输入文件中的四个TH2D中的点,进行刻度，并绘制常见粒子的理论曲线进行对照


7. calibration_t0_v1
目前的刻度效果并不理想，尤其是对d1的效果很差，原因很可能是d1的厚度太小，厚度存在不均匀性，并且
de-e的pid分辨较差，因此我们需要新写一个calibration_t0_v1.cpp程序，以原先的程序为基础，进行改进，升级
t0的刻度方案。有以下几点：

- 首先使用八个共同参数,利用d2-d3,d3-d4,d4-s1的pid图，进行拟合,得到除了d1的其他参数的系数。
- 在得到准确的d2刻度系数之后，再使用d1-d2的pid图，进行拟合，得到d1的刻度系数。便可以规避d1的不准确性对其他层
的影响。本次拟合只拟合d1的 p0,p1
- 以上两次拟合便得到完整的t0刻度系数
- 另外，我们目前采用的各个粒子的delta_e的理论曲线是根据config文件中硅的理论厚度得到，但是
我们可能有时候会改变硅的厚度，因此需要根据实际情况，重新计算理论曲线。因此在刻度开始前，可以对话
询问用户是否需要根据config文件中的厚度，重新计算de-e的理论曲线.


架构建议：
1. 同意，把独立的头文件和源文件放到brill2/src/brill/include/t0和brill2/src/brill/src/t0目录下，至于
calibrate_t0.cpp就维持原状吧。复用函数都只给v1使用
2. 正确
3. 没错，但是可以完全沿用我之前给你的offset表，无非就是d1d2的只在第二阶段用就行
4. 可以。但是si_z*_a*.root 文件是和厚度无关的，删除重建的应该是t0_delta_z*_a*.root文件
算法建议：
5. 同意
6. 沿用之前的限制即可，并且这两者粒子的pid没有重叠区域

疑问：
1. 使用所有除了d1-d2以外的所有cut，同时要留有之后增加新粒子，新cut的空间。
2. 两种都使用，但是要留有之后增加新粒子，新cut的空间。
3. 覆盖v0的输出文件
4. 不变
5. 并存，可以执行程序命名为calibration_t0_v1
