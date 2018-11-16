---
titl: 外显子数据生成模拟
---
# 外显子数据生成模拟
## 1. 环境准备

    python3 和对应的 numpy 扩展程序库

文件已准备了1.2.步骤生成的数据,可以跳过1.2.步,直接到第三步.

如果你要更改设置 setting.py, 1.2.步骤不可缺少
## 2. 数据初始化

这个步骤是通过注释文件提取外显子信息,在根据外显子信息从基因组序列中获得外显子组序列,供生成测序下机数据使用

需要基因组序列文件和基因组注释文件,产生NCseq.txt,ENSEMBL_exonposition.txt,ENSEMBL_exonlist.txt,

1. 基因组注释文件

    可以在 https://www.gencodegenes.org/human/ 下载.其他网站也可以,但是对基因注释文件的标注顺序有要求(因为没有写重新排序的方法),也不保证其他注释文件可以使用(因为还没测试):
    
        需要按染色体1-24+线粒体,基因从低位点到高位点,同一个基因的外显子要要是连续的,中间不能插其他外显子
    目前确定可以使用的文件是上述网站下载的'gencode.v29.annotation.gtf'.
    解压得到gtf文件

2. 基因组文件

    可以在 http://genome.ucsc.edu/ 下载,其他网站也可以,对染色体和线粒体的基因组按照同上顺序(也是没有排序的处理),也不保证其他基因组文件可以使用.

    目前确定可以使用的文件是上述网站下载的'GRCh38_latest_genomic.fna.gz'并解压得到fna文件

3. 运行exonlist.py

    在命令行输入: python3 exonlist.py genomefile(基因组文件) annotationfile(基因注释文件)

    ```
    python3 exonlist.py genomefile annotationfile
    ```
4. 生成的文件说明

    1. NCseq.txt 是完整染色体基因组序列
    2. ENSEMBL_exonposition.txt 是外显子注释文件,默认根据ENSEMBL输出,如果你的注释文件是小写的ensembl,需要更改exonlist.py的内容
    3. ENSEMBL_exonlist.txt 是外显子加上两端侧翼序列的外显子组序列

## 3. 获得碱基质量的频率分布
需要其他额外的fastq文件

产生frequencies_100.json文件

1. 碱基质量来源说明

    从 https://www.ebi.ac.uk/ 下载到的部分fastq文件是illumina hiseq 2000测到的数据,数据长度100bp,为phrea+33,因能力有限没有找到适合的分布,因而按照对>1000000条的read数据,获得100bp长度的read每个碱基质量频率分布情况,保存到'frequencie_100.json'
2. 运行qphred.py

    需要准备文件和目的read长度相同的其他fastq文件材料,比如 'SRR4238252_2.fastq'
    ```
    python3 fastqfile
    ```
3. 生成的文件
   
    frquencies_100.json

## 4. 输出fastq文件
需要ENSEMBL_exonlist.txt,frequencies_100.json

文件已准备了这两个原始数据,可以跳过前面的1.2.步

1. 运行readout.py
   
    可以在命令行执行后再手动设置突变参数
    ```
    python3 readout.py
    ```
    也可以将输入好的突变参数先写到test.txt文件中
    ```
    python3 readout.py < test.py
    ```
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

    ```
    **********Show Genetype**********
    *****percent=40*****
    homolog1:Normal type
    homolog2:Normal type

    *****percent=30*****
    homolog1:Mutation type
    NO.8400 exon : cnv
    CNV=1
    NO.24993 exon : cnv
    CNV=1
    NO.33605 exon : origin snp + cnv
    SNP=3 : [0.05678, 0.15552, 0.81389] ; CNV=2
    NO.51017 exon : origin snp + cnv
    SNP=6 : [0.0604, 0.15138, 0.2067, 0.33815, 0.70043, 0.8536] ; CNV=1
    NO.74072 exon : cnv
    CNV=1
    homolog2:Mutation type
    NO.22489 exon : new snp + cnv
    SNP=1 : [0.01707] ; CNV=1
    NO.25788 exon : new snp + cnv
    SNP=7 : [0.09404, 0.5478, 0.60091, 0.75362, 0.89251, 0.93089, 0.9969] ; CNV=2
    NO.45355 exon : cnv
    CNV=1
    NO.51624 exon : cnv
    CNV=1
    NO.52507 exon : cnv
    CNV=7
    NO.62289 exon : cnv
    CNV=2
    NO.73353 exon : cnv
    CNV=1

    *****percent=30*****
    homolog1:Mutation type
    NO.18588 exon : cnv
    CNV=1
    NO.30785 exon : cnv
    CNV=2
    NO.34853 exon : cnv
    CNV=1
    NO.36561 exon : origin snp
    SNP=4 : [0.04269, 0.30547, 0.60146, 0.7943]
    NO.47838 exon : cnv
    CNV=1
    NO.51283 exon : deletionhomolog2:Mutation type
    NO.3336 exon : cnv
    CNV=3
    NO.15283 exon : deletionNO.21948 exon : deletionNO.23101 exon : deletionNO.28918 exon : cnv
    CNV=1
    NO.29998 exon : deletionNO.41145 exon : cnv
    CNV=1
    NO.45010 exon : deletionNO.45585 exon : deletionNO.52723 exon : cnv
    CNV=1
    NO.60553 exon : deletionNO.60737 exon : cnv
    CNV=2
    NO.61738 exon : deletion
    **********Show Genetype**********
    print "enter" to continue...
    ```

3. 基因组突变参数的说明

    在你输入完百分比后回车,你需要输入4个数字:从左到右分别是

    1. 外显子突变总数:exon_num
    2. SNP突变总数:snp_sum
    3. 多拷贝新增外显子的总个数:cnc_cum (比如一个基因多了一份拷贝,总共两份拷贝,记为1)
    4. 缺失外显子总个数:dele_num

4. 生成的文件说明

    1. fastq+日期.txt 生成的fastq数据
    2. mutation_setting_log.txt 基因型设置的数据

## 5. 设置

1. INSERT=200 
2. READ =100
3. ACCURACY_RATE= 0.8 
4. DEEPTH=50
5. ROW_STEP=2
6. ROW_COLUMN=80
7. CHIP_LEN=10

3.4.可以随意更改
其余的更改可能会出问题

## 6. 例子
```
python exonlist.py GRCh38_latest_genomic.fna gencode.v29.annotation.gtf
python qphred.py SRR4238252_2.fastq
python readout.py < test.txt
```
## 7. 其他说明
1. read碱基质量:

    因为如果每次read输出的碱基质量都需要重新生成1-100位的碱基质量,开销很大,因而每次就初始化100个随机生成的碱基质量列,而后随机抽取.虽然后者可能数据不如前者,但是时间开销约是前者的1/4
2. 运行时间:

    50×的深度约50分钟,预计100×的深度约100分钟

    时间应该主要在运算上,读写时间影响不大(不确定)

3. 误读的碱基质量:

    减10

4. 输出的序列如果含N:

    会随机分配ATCG中的一个,但是作用在insert上,还是read上,还是exon上.有待商榷.

    目前作用在insert上

5. @名称

    没有设置

6. 突变的设置

    可以参考下面说明

    ```
        no. len(snp),cnv,dele,illustration
        1.  5 6 1  origin snp +cnv
        2.  0 0 1  deletion
        3.  5 0 1  origin snp
        4.  0 6 1  invalid
        ---------------------
        5.  5 6 0  new snp +cnv
        6.  0 6 0  cnv
        7.  5 0 0  invalid
        8.  0 0 0  invalid
        1-4 origin exon not exist
        5-8 origin exon exist
        one exon can have new_snp+cnv and cnv at same time
    ```
