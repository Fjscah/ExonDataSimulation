---
title: 外显子数据生成模拟
---
# 外显子数据生成模拟
## 介绍

    这是一个模拟多个基因型以某种比列混杂,生成外显子测序数据的脚本.主要用于癌症比列分析软件的准确性检测
    
    癌症的变异类型多样,有的程度轻微,有的程度巨大,而且一个肿瘤组织中混杂着正常和多种变异类型的细胞.因而出现了很多分析癌变细胞比例的软件.但是这些软件分析结果差异很大,与人工分析细胞形态得出的分析结果也相差甚远.因为实验细胞样的真正比例无从事先知道,因而需从设置变异参数出发,确定变异比例,模拟生成外显子测序数据,来检测分析软件的准确性

<<<<<<< HEAD
    这个脚本的过程大致是:根据染色体碱基序列和外显子注释获得外显子序列->由原始的fastq文件获得碱基质量分布->根据外显子序列设置突变情况(各个变异型细胞比例,它们的snp,cnv突变情况)->最后根据突变文件,碱基质量分布,外显子序列生成相应fastq文件.而后用于检测分析软件的准确性

* 程序分为两部分:
  
    1. 初始化文件 exonlist.py :
   
        用于获得外显子序列和碱基质量分布
    2. fastq生成文件 readout.py:
        
        外显子序列和碱基质量分布
=======
文件已准备了2.3.步骤生成的数据,可以跳过2.3.步,直接到第4.步.

如果你要更改设置 setting.py, 2.3.步骤不可缺少
## 2. 数据初始化
>>>>>>> b4d11df07c2e3febb9d1bdde95e67504ebaeaa20

    3. 突变文件目前需要手动设置,没有相关程序
   
* 数据生成需要一些文件:

    共需要准备4个文件:

        基因组文件,基因注释文件,fastq参考文件,突变设置文件

    产生一些中间文件:
        
        染色体基因组文件,外显子基因注释文件,碱基质量分布文件,外显子初始化文件
    
    最后得到PE的fastq文件
    
## 环境准备

    python3 

## 数据初始化


这个步骤是通过 基因组文件 , 基因注释文件 和 fastq参考文件 , 提取外显子信息和碱基质量分布, 供生成测序下机数据使用

前面所提到的中间文件都是在这里生成


1. 基因组文件

    fna格式, 可以在 UCSC ( http://genome.ucsc.edu/ ), NCBI ,下载 hg37(19)/hg38 文件 . 请确保染色体文件在同一文件中,并且按照从前往后的顺序. UCSC是每个染色体一个文件,需要将他们合并,合并完成后就是由完整染色体构成的基因组文件
    
2. 基因组注释文件

    gff格式, 可以在 Gencode( https://www.gencodegenes.org/human/ ) , NCBI下载hg37(19)/38 的基因文件 . 注释文件的染色体顺序也是从前往后,即与基因组文件保持一致

3. 碱基质量分布参考文件
   
   因为从现有数据分析上(还没有查资料,只是根据分布图得出),还没有找到碱基质量的匹配的分布,因而目前只能先按照碱基的频率分布获得碱基质量的分布. 这就需要有参考的fastq文件来得到碱基质量的分布情况

4. 运行exonlist.py
   
    执行 exonlist.py
    ```
    python3 exonlist.py
    ```
    会给出操作说明:
    ```
    init -seq -filex -file1 -ver    : get chromosome sequence from filex, generate file1
    init -ano -filey -file2 -ver    : get exon annotation from filey, generate file2
    init -list -file1 -file2 -file3 -ver :get exonlist from file1,file2, generate file3
    init -all                       :excute above all according setting,py
    phred -file -XX                 : get qphred frequencies from file,generate 'phred.json'
    view -file3                     : view exonlist's exon sequence
    ```
    首先将基因组文件过滤,只获得染色体的基因(如果你的文件已经是完整的染色体文件,没有参杂其他片段,可以跳过这一步).
    
    比如从NCBI下载的fna文件(NCBI_gh38.fna),既有染色体也有其他序列,要获得过滤的文件( NCBI_hg38ref.fna ), 则按如下执行,后面的'NCBI'表示是针对NCBI的格式
    ```
    >init -seq -NCBI_gh38.fna -NCBI_hg38ref.fna -NCBI
    ```
    接着将基因组注释文件过滤,获得外显子的注释文件(如果你的注释文件就是的每行都是你的目的序列,即外显子的标注,可以跳过这一步)

    比如从NCBI下载的gff文件(NCBI_hg38.gff),只要其中外显子的注释,生成NCBI_hg38exon.gff , 则按如下输入, 'NCBI'表示这是NCBI常用的注释格式
    ```
    >init -ano -NCBI_hg38.gff -NCBI_hg38exon.gff -NCBI
    ```
    最后根据前面生成好的两个文件,就可以获得初始化的外显子序列

    就接上面示例操作,要根据NCBI_hg38ref.fna和NCBI_hg38exon.gff获得外显子序列NCBI_hg38list.txt, 则如下输入
    ```
    >init -list -NCBI_hg38ref.fna -NCBI_hg38exon.gff -NCBI_hg38list.txt -NCBI
    ```
    对于上面三部如果你已经在 setting.py 中设置(详细参考请见:),就可以直接执行
    ```
    >init -all
    ```
    接下来还要获得碱基的质量分布

    比如准备好了source100.fastq文件 ,这个文件的碱基质量是 ASCII+33 转换的(即碱基质量+33=它对应在fastq中的质量符号的值), 按照如下输入 , 如果文件是ASCII+64,则改为+64即可
    ```
    >phred -source100.fastq -33
    ```
    还剩下突变设置文件要准备,为辅助突变文件的书写,可用'view'查看外显子序列. (突变格式要求和view查看详细请见:)

    比如查看上面获得NCBI_hg38list.txt外显子序列,如下输入,即可转到view查看
    ```
    >view -NCBI_hg38list.txt
    ```
    完成后, 'exit' 退出
    ```
    >exit
    ```


## 输出fastq文件

需要 外显子序列文件,碱基质量分布文件和突变设置 3个文件

1. 运行readout.py
   
    执行 readout.py 
    ```
    python3 readout.py
    ```
    会询问需要通过哪种方式创造突变
    ```
    which mutation creating way do you want?
        1. set all typess of mutation's total number
        2. set every mutatition manually
        3. import mutation file
        >
    ```
<<<<<<< HEAD
    目前只完成了第三种,请输入3
=======
    手动输入或编写的test.txt的突变参数参考下一步(2. 设置突变参数) 说明
2. 设置突变参数

    1. 一开始先输入你需要的不同基因型个数
    2. 然后再按顺序设置基因型的百分比
    3. 接着设置二倍体的一个基因组突变情况
    4. 再设置另一个基因组突变情况
    5. 重复到最后会自动设置最后一个基因型的百分比
    6. 待设置完两个基因组后会打印设置,若没有问题,回车输出,开始运行
    7. 如果exon中有N,则会输出这个exon的序列,但不会保存到文件中查看

    下面示范参数的输入
    ```    
    E:\mywork\history4>python readout.py
    initial qphred...

    down...
    exon count:  75941
    >1. how much mutation do you want to creat(including normal type):3
    OK. then enter mutation like 'exon_num snp_num cnv_num dele_num';
            if exon_num==0 ,this type is normal type;
            dele_num need be smaller than exon_num;
            eg. 5 9 6 2 is right ; 3 2 1 4 is wrong ; 0 0 0 0 is normal type
    >2. please enter your no.1 mutation percent : 40
    >3. please enter your no.1 mutation homolog-1 : 0 0 0 0
    >4. please enter your no.1 mutation homolog-2 : 0 0 0 0
    >5. please enter your no.2 mutation percent : 30
    >6. please enter your no.2 mutation homolog-1 : 5 9 6 2
    >7. please enter your no.2 mutation homolog-2 : 7 8 15 0
    >8. please enter your no.3(last) mutation percent : 30
    >9. please enter your no.3 mutation homolog-1 : 6 4 5 2
    >10. please enter your no.3 mutation homolog-2 : 13 0 8 8
    ```
> 1. 先设置了组织样中有3个基因型.然后会打印突变参数的示范.
> 2. 第一个基因型我想设置成正常的基因型,占40%,就输入40回车,然后对两个父母双方的基因突变都设置为0 0 0 0,表示没有任何突变<br>
> 3. 第二个基因型我想设置突变的基因型,占30%,就输入30回车,如果此时输入80,显然是不合理的,程序会要求重新输入
> 4. 然后我希望总共有:<br><br>5个外显子发生了突变(如果一个外显子发生了原外显子的拷贝和原外显子突变之后又拷贝,这个算两个外显子,因为是一个外显子的两次突变)<br>9个SNP突变位点(那种突变后又拷贝的只算一次,与拷贝倍数无关)<br>6个拷贝子(只是算多出的拷贝数,snp突变不计入拷贝数,cnv突变产生2个拷贝,只记一个拷贝子,因为只多出1个)<br>2个缺失突变<br><br>就输入5 9 6 2<br>!!!需要注意缺失的突变数目要小于总外显子数,缺失的突变加上拷贝子数要大于外显子数

> 8. 反复操作到最后一个基因型的设置,会自动写剩余的百分比,不用再填写.
    
再设置完最后的基因型后回车就会随机分配突变情况(具体是哪个外显子,产生哪些突变),并打印,如下:
    
1. normal type 表示该基因组是正常的
2. mutati type 表示该基因组是发生变异的
    1. NO.XXXX 是在exonlist的文件中第XXXX个外显子发生突变
    2. cnv表示该外显子发生多拷贝突变 CNV=X 表示共X+1个该外显子
    3. origin snp+cnv 表示是原味点突变后又发生多拷贝 SNP=X 表示该外显子又X个点突变,后面的列表[]表示这些突变都分布在这个外显子百分之几的位置 . CNV=Y表示
>>>>>>> b4d11df07c2e3febb9d1bdde95e67504ebaeaa20

    XXX是突变设置文件,而后会打印突变,回车确认无误后,输入初始好的json格式碱基质量文件YYY(phred.json)和外显子序列文件ZZZ,如下
    ```
    please input mutation setting file : XXX
    **********Show Genetype**********
                ......
    **********Show Genetype**********
    press Enter to continue...
    please input qphred file : YYY
    initial qphred...
    down...
    please input exonlist file : ZZZ
    ```
    接下来等待1-2h,完成数据生成
    手动输入或编写的test.txt的突变参数参考下一步(2. 设置突变参数) 说明

2. 生成的文件

        R1_fastq+日期.fastq 和 R2_fastq+日期.fastq


## 5. 设置



3.4.可以随意更改
其余的更改可能会出问题

## 6. 例子
```
# exonlist.py
init -seq -NCBI_gh38.fna -NCBI_hg38ref.fna -NCBI
init -ano -NCBI_hg38.gff -NCBI_hg38exon.gff -NCBI
init -list -NCBI_hg38ref.fna -NCBI_hg38exon.gff -NCBI_hg38list.txt -NCBI
init -list -NCBI_hg38ref.fna -COM_hg38.bed -COM_hg38list.txt -COM
phred -source100.fastq -33
view --COM_hg38list.txt

# readout.py
3
mutations_setting.txt
phred.json
COM_hg38list.txt
```
## 7. 其他说明
1. read碱基质量:

    因为如果每次read输出的碱基质量都需要重新生成1-100位的碱基质量,开销很大,因而每次就初始化100个随机生成的碱基质量列,而后随机抽取.虽然后者可能数据不如前者,但是时间开销约是前者的1/4
2. 运行时间:

    150M左右的外显子初始化序列文件,100×的深度约一个半小时


3. 误读的碱基质量:

    减10

4. 输出的序列如果含N:

    会随机分配ATCG中的一个,但是作用在insert上,还是read上,还是exon上.有待商榷. 目前作用在insert上

5. @名称

    对于正常序列 : 外显子id:突变编号:输出的第几个read

        exonid:mutation.no:turns:

    对于突变序列: 外显子id:突变编号:母方/父方(0/1):输出的第几个read:突变信息

        exonid:mutation.no:haplot:turns:mutationinfo


  
