# WCP参数计算说明

​		首先说明一下, WCP全称Warren-Cowley parameters, 用于描述合金模型的短程有序, 对于在第一近邻处的WCP而言, 其计算公式如下:
$$
WCP_{mn} = 1 - \frac{Z_{mn}}{\chi_nZ_m}
$$
​		其中$WCP_{mn}$为第m种原子和第n种原子之间的WCP系数, 而$Z_{mn}$为m类型的原子的第一近邻中所有的n类型的原子的个数, $\chi_n$为n类型原子的在模型中的浓度, $Z_m$为m类型原子中的第一近邻中所有的m类原子,  其取值范围在(-2.0, 1.0)之间, 如果$WCP_{mn}$越接近1.0, 说明m和n元素越排斥, $WCP_{mn}$值越接近-2.0, 说明m和n元素越吸引, 如果$WCP_{mn}$接近于0, 说明mn原子对相对来说是随机分布的. 

​		这里也有必要解释一下第一近邻的意思, 看着很高大上, 说白了就是RDF曲线的第一个波峰的横坐标, 而另一个我们熟知的算法也利用到了第一近邻这个数据, 即CNA(Common neighbor analysis), 而里面也提及了第一近邻的计算[1](https://linkinghub.elsevier.com/retrieve/pii/0927025694901090). 那么我们岂不是每一次计算都需要先算一次RDF再去算WCP? 其实在lammps手册中, 提及到对于金属晶体而言, 波峰的横坐标与其晶格常数有一个**较为精确**的数学关系.[2](https://docs.lammps.org/compute_cna_atom.html)
$$
r_c^{fcc} = \frac{1}{2}(\frac{\sqrt2}{2} + 1 )a \approx 0.8536a \tag{1}
$$

$$
r_c^{bcc} = \frac{1}{2}(\sqrt2 + 1 )a \approx 1.207a \tag{2}
$$

$$
r_c^{hcp} = \frac{1}{2}\left(1 + \sqrt{\frac{4 + 2x^2}{3}} \right) a \tag{3}
$$

​		在(3)式中, x表示$\frac{c}{a}$, 对于大多数的HCP晶体来说, 1.6333是一个比较常用的值, 而在lammps的源代码中, 就默认了x均为1.63333, 但是到了具体问题上, 最好去查询资料以确定相关的参数. 

​		那么现在就很明了了, 我们只需要模型文件, 晶格常数和晶格类型, 理论上就可以计算出所有合金晶体模型的WCP系数矩阵了. 而在本程序中, 目前只支持**FCC**和**BCC**的合金的模型的计算, HCP晶格类型的合金, 由于印象中在文献中看到的比较少, 所以暂时没有加入计算中, 同时默认盒子的3个方向均为周期性边界条件, 相关的系数均按照Fortran的namelist的格式进行输入. 其格式如下:

```
&inputparameter
filename =    'CoCuFeNiPd-2m.txt'
ltype    =    'fcc'
lattice  =    3.6d0
/

```

​		你需要以& + namelist的名字开头, 并以 / 号结尾, 之后一一为对应的变量赋值, 就像在Fortran中那样即可, 注意字符串变量可以用单引号或者双引号括起来, 3.6d0表明这是一个双精度的实数, namelist的变量赋值的语法和在Fortran中编程基本一致. 另一个很重要的点就是最后面还需要额外的空一行, 否则程序读取的时候会报错.



# 编译和运行方法

​		本程序只有2个文件, 其实手动敲命令输入也可以, 但是我还是写了Makefile, 准备了2种不同的编译方式, 分别是串行版本和并行版本, 当然这里推荐并行编译, 因为现在openMP已经是gfortran的标准库了, 基本不存在环境配置的问题, 在Linux的**Ubuntu**版本上, 可以使用:

```shell
sudo apt install gfortran make
```

​		安装gfortran和make, 并使用:

```shell
which gfortran
```

​		查看gfortran是否安装成功. 软件安装准备完毕之后就是开始安装软件了.

## 快速开始

​		这一节是专门针对Linux和Fortran老手的快速开始教程, 直接输入:

```
gfortran mod_type.f90 wcsro.f90 -o wcs -fopenmp -Wall -Wno-unused-function -Wno-maybe-uninitialized -Wno-unused-variable
```

​		之后就可以获得可执行文件wcs了, 使用方法就不再赘述.

## 并行版本编译和使用

​		其实并行版本是**默认的**编译选项, 只需要输入: 

```shell
make
```

​		就可以编译成功, 同时文件夹下面会出现一个名为wcp_op的可执行文件, 假设我们的namelist的文件内容如下, 且文件名为input.txt, 需要用户输入文件名称(**filename**), 晶格类型(**ltype**), 以及晶格常数(**lattice**):

```
&inputparameter
filename =    'CoCuFeNiPd-2m.txt'
ltype    =    'fcc'
lattice  =    3.6d0
/
```

​		openMP并行, 可以通过以下命令设置并行核数, 这里我们假设调用4核并行, 

```shell
export OMP_NUM_THREADS=4
```

​		输入如下命令运行:

```shell
./wcp_op input.txt
```

​		运行结果:

```
    5    0.70    0.05   -0.95    0.78   -0.58
    4   -0.31    0.51    0.59   -1.57    0.78
    3    0.45    0.36   -0.45    0.59   -0.95
    2    0.51   -1.43    0.36    0.51    0.05
    1   -1.35    0.51    0.45   -0.31    0.70
            1       2       3       4       5
```

​		当然这里说明一下, 如果不预先设定openMP的并行核数, 其实Fortran默认直接使用能调用的最大核数进行并行. 所以如果不是那种数十个乃至100多个核的服务器, 没必要提前设定并行核数.

## 串行版本编译和使用

​		其实不推荐串行版本, 因为速度较慢, 不过说是说速度慢, 其实就本文件夹中的3个模型文件而言, 和并行比起来就是一个0.4s和一个2s的区别, 基本上都是一瞬间完成, 串行版本需要额外的指令输入:

```shell
make wcp_se
```

​		之后发现文件夹下面多出了一个名为wcp_se的可执行文件, 运行方式与并行版完全相同, 假设文件名和文件内容与并行版本软件完全相同, 则不需要设置openMP的核数, 命令同样是:

```
./wcp_se input.txt
```

​		运行结果如下:

```
    5    0.70    0.05   -0.95    0.78   -0.58
    4   -0.31    0.51    0.59   -1.57    0.78
    3    0.45    0.36   -0.45    0.59   -0.95
    2    0.51   -1.43    0.36    0.51    0.05
    1   -1.35    0.51    0.45   -0.31    0.70
            1       2       3       4       5
```

​		毫不意外, 结果与并行相同但是速度差了不少.

# 补充说明

​		附带的文件夹里面提供了一篇研究CoCuFeNiPd短程有序的NC顶刊和文章里面的3个模型, 用户可以试着修改输入文件计算3个模型的WCP矩阵和文章中的WCP矩阵对比, 这里由于文章中提供的3个模型没有给出1-5号原子对应哪些元素, 所以只能以数字1, 2, 3表示不同类型的原子, **同时用户如果想要计算其他的合金模型, 模型文件也必须严格按照文件夹中提供的3个模型文件的格式, 否则计算就会失败**. 同时对比了我的计算结果和文章中的计算结果, 数值都完全可以对上, 误差在$\pm 0.01$以内.

​		我的计算结果:

​		模型文件: CoCuFeNiPd-0.txt

```
    5   -0.00   -0.01   -0.00    0.00    0.01
    4    0.01   -0.01    0.02   -0.02    0.00
    3    0.01    0.01   -0.04    0.02   -0.00
    2    0.00    0.00    0.01   -0.01   -0.01
    1   -0.02    0.00    0.01    0.01   -0.00
            1       2       3       4       5
```

​		文献结果:

​                                                                		<img src=".\picture\WCP_RM.png" style="zoom: 33%;" />

​		模型文件: CoCuFeNiPd-2m.txt

```
    5    0.70    0.05   -0.95    0.78   -0.58
    4   -0.31    0.51    0.59   -1.57    0.78
    3    0.45    0.36   -0.45    0.59   -0.95
    2    0.51   -1.43    0.36    0.51    0.05
    1   -1.35    0.51    0.45   -0.31    0.70
            1       2       3       4       5
```

文献结果:

<img src=".\picture\WCP_2M.png" style="zoom: 33%;" />

​		模型文件: CoCuFeNiPd-4m.txt

	    5    0.66    0.20   -1.04    0.85   -0.67
	    4   -0.30    0.51    0.63   -1.69    0.85
	    3    0.39    0.58   -0.57    0.63   -1.04
	    2    0.64   -1.94    0.58    0.51    0.20
	    1   -1.39    0.64    0.39   -0.30    0.66
	            1       2       3       4       5

​		文献结果:

<img src=".\picture\WCP_4M.png" style="zoom: 33%;" />