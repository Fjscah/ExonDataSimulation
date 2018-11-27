---
titl: 外显子数据生成模拟
---
# 外显子数据生成模拟
## **介绍**

    这是一个模拟多个基因型以某种比列混杂,生成外显子测序数据的脚本.主要用于癌症比列分析软件的准确性检测
    
    癌症的变异类型多样,有的程度轻微,有的程度巨大,而且一个肿瘤组织中混杂着正常和多种变异类型的细胞.因而出现了很多分析癌变细胞比例的软件.但是这些软件分析结果差异很大,与人工分析细胞形态得出的分析结果也相差甚远.因为实验细胞样的真正比例无从事先知道,因而需从设置变异参数出发,确定变异比例,模拟生成外显子测序数据,来检测分析软件的准确性

    这个脚本的过程大致是:根据染色体碱基序列和外显子注释获得外显子序列->由原始的fastq文件获得碱基质量分布->根据外显子序列设置突变情况(各个变异型细胞比例,它们的snp,cnv突变情况)->最后根据突变文件,碱基质量分布,外显子序列生成相应fastq文件.而后用于检测分析软件的准确性

* 程序分为两部分:
  
    1. 初始化文件 initial.py :
   
        用于获得外显子序列和碱基质量分布

    2. fastq生成文件 readout.py:
        
        外显子序列和碱基质量分布

    3. 突变文件目前需要手动设置,没有相关程序
   
* 数据生成需要一些文件:

    共需要准备4个文件:

        基因组文件,基因注释文件,fastq参考文件,突变设置文件

    产生一些中间文件:
        
        染色体基因组文件,外显子基因注释文件,碱基质量分布文件,外显子初始化文件
    
    最后得到PE的fastq文件
    
## **环境准备**

    python3 

## **数据初始化**


这个步骤是通过 基因组文件 , 基因注释文件 和 fastq参考文件 , 提取外显子信息和碱基质量分布, 供生成测序下机数据使用

* 前面所提到的中间文件都是在这里生成

    1. 基因组文件

        fna格式, 可以在 UCSC ( http://genome.ucsc.edu/ ), NCBI ,下载 hg37(19)/hg38 文件 . 请确保染色体文件在同一文件中,并且按照从前往后的顺序. UCSC是每个染色体一个文件,需要将他们合并,合并完成后就是由完整染色体构成的基因组文件
        
    2. 基因组注释文件

        gff格式, 可以在 Gencode( https://www.gencodegenes.org/human/ ) , NCBI下载hg37(19)/38 的基因文件 . 注释文件的染色体顺序也是从前往后,即与基因组文件保持一致

    3. 碱基质量分布参考文件
    
        因为从现有数据分析上(还没有查资料,只是根据分布图得出),还没有找到碱基质量的匹配的分布,因而目前只能先按照碱基的频率分布获得碱基质量的分布. 这就需要有参考的fastq文件来得到碱基质量的分布情况

### **初始化步骤**
   
1. 执行 exonlist.py
    ```
    python3 exonlist.py
    ```
    会给出操作说明:
    ```
    init -seq -filex -file1 -ver    : get chromosome sequence from filex, generate file1
    init -ano -filey -file2 -ver    : get exon annotation from filey, generate file2
    init -list -file1 -file2 -file3 -ver :get exonlist from file1,file2, generate file3
    init -qph -filez -file4 -XX     : get qphred frequencies from file,generate 'phred.json'
    init -all                       :excute above all according setting,py
    view -file3                     : view exonlist's exon sequence
    ```
2. 首先将基因组文件过滤,只获得染色体的基因(如果你的文件已经是完整的染色体文件,没有参杂其他片段,可以跳过这一步).
    
    比如从NCBI下载的fna文件(NCBI_gh38.fna),既有染色体也有其他序列,要获得过滤的文件( NCBI_hg38ref.fna ), 则按如下执行,后面的'NCBI'表示是针对NCBI的格式
    ```
    >init -seq -NCBI_gh38.fna -NCBI_hg38ref.fna -NCBI
    ```
3. 接着将基因组注释文件过滤,获得外显子的注释文件(如果你的注释文件就是的每行都是你的目的序列,即外显子的标注,可以跳过这一步)

    比如从NCBI下载的gff文件(NCBI_hg38.gff),只要其中外显子的注释,生成NCBI_hg38exon.gff , 则按如下输入, 'NCBI'表示这是NCBI常用的注释格式
    ```
    >init -ano -NCBI_hg38.gff -NCBI_hg38exon.gff -NCBI
    ```
4. 根据前面生成好的两个文件,就可以获得初始化的外显子序列

    就接上面示例操作,要根据NCBI_hg38ref.fna和NCBI_hg38exon.gff获得外显子序列NCBI_hg38list.txt, 则如下输入
    ```
    >init -list -NCBI_hg38ref.fna -NCBI_hg38exon.gff -NCBI_hg38list.txt -NCBI
    ```
5. 接下来还要获得碱基的质量分布

    比如准备好了source100.fastq文件 ,这个文件的碱基质量是 ASCII+33 转换的(即碱基质量+33=它对应在fastq中的质量符号的值), 按照如下输入 , 如果文件是ASCII+64,则改为+64即可
    ```
    >phred -source100.fastq -phred.json -33
    ```
6. 对于上面四步如果你已经在 setting.py 中设置(详细参考请见:),可以跳过上面三部,直接执行这一步即可
    ```
    >init -all
    ```
    
## **突变设置**

### **突变设置的格式:**

* 格式解释
   
    首写需要给出各种突变的百分比,然后给出每种突变二倍体的基因型 :
    1. chr : 突变所在染色体
    2. exon : 在初始化的外显子文件中,外显子是染色体的第几个染色体
    3. start,end : 突变相对于外显子的终止位点
    4. ref,alt : 参考基因组序列,和突变后的序列,如果是长的缺失,ref可以为'.',alt用'.'表示
    5. cn : ?/?, 分别表示两个单倍体型, '.'表示这个单倍型没有该突变,如果为数字,则表示这样的外显子又相应数字的拷贝
    6. description : 描述突变类型,可以省略,省略请用'.'表示

* 以下面例子来说
    ```
    ===================
    >genometype:3
    mutation    percent
    1           50
    2           40
    3           10
    ------------------
    >mutation:1
    chr     exon    start   end     ref     alt     cn      description 
    1       136     .       .       .       .       ./0     hybrid,cnv
    1       136     40      44      CAGGC   .       ./5     hybrid,indel
    1       136     50      51      CT      AC      ./1     hybrid,snp
    1       136     40,50   44,51   CAGGC,CT  .,AC  ./2     hybrid,indel,snp
    1       137     .       .       .      .       ./0     hybrid,cnv  
    1       436     .       .       .       .       ./0     hybrid,CNV
    1       436     56      59      AGGC    T       2/.     hybrid,indel
    1       436     76      77      AC      AT      1/.      hybrid,snp
    1       569     156     159     CCTA    G       1/1     homo,snp
    4       376     .       .       .       .       0/0     homo,deletion 
    >mutation:2
    chr     exon    start   end     ref     alt     cn      description      
    >mutation:3
    chr     exon    start   end     ref     alt     cn      description     
    2       436     56      60      GTAGC   .       1/.     hybrid,indel    
    3       436     .       .       .       .       0/2     hybrid,deletion
    ```
    文中设置了三个突变,编号分别为1,2,3,组织块中的混合比例为50% , 40% , 10%. 
    
    用'>genometype:'做指示
    ```
    ===================
    >genometype:3
    mutation    percent
    1           50
    2           40
    3           10
    ------------------
    ```
    然后设置每种突变型前,用'>mutation:编号'做指示. 具体说明如下:
    
    1. 每个外显子最多只有一次start=end=ref=alt='.'的情况,表示这个外显子本身拷贝数是缺失(0)还是正常(.),还是多拷贝(N>0),这里的description都是cnv.比如第一行
    2. 其余的情况start都≠'.',start和end要对应,并从小到大排列,比如第二行
    3. 如果ref过长,可以用'.'表示,相应地description为indel,否则为snp.比如三到四行
    4. 如果cn两个是相同的,则为纯合子homo,比如三四行;如果不相同,则为杂合子,比如一二行
    5. description不影响最后数据生成,不写也没关系,只是为了好理解

    ```
    >mutation:1
    chr     exon    start   end     ref     alt     cn      description 
    1       136     .       .       .       .       ./0     hybrid,cnv
    1       136     40,50   44,51   CAGGC,CT  .,AC  ./2     hybrid,snp
    1       569     156     159     CCTA    G       1/1     homo,snp
    4       376     40       60       .       .       1/1     homo,indel 
    还剩下突变设置文件要准备,为辅助突变文件的书写,可用'view'查看外显子序列.

    比如查看上面获得NCBI_hg38list.txt外显子序列,如下输入,即可转到view查看
    ```
    如果没有突变,比如第二个突变型,如下即可:
    ```
    >mutation:2
    chr     exon    start   end     ref     alt     cn      description 
    ```
### **使用view 辅助写突变**

手动书写不知道外显子具体情况,view可以辅助突变的书写.

1. 在 initail.py 中 输入 view -exonlistfile (也可以直接调用 view.py ),即可进入. 具体操作如下:

    ```
    chr -n                  : show chromosome exon number
    chrx -n                 : shoe chx exon number
    chrx -pos -pos1 -pos2   : show chx from pos1 to pos2 sequence,pos is absolute
    chrx -num               : show chx no.num exon sequence
    chrx -num -pos1 -pos2   : show chx no.num from pos1 to pos2 sequence,pos is relative
    chrx -num -pis1 -pot2 -t: show chx no.num exon sequence,mark pos1 and pos2
    ```
2. 示例1
   
    ```
    >view -test.txt
    ```
    ```
    chr -n                  : show chromosome exon number
    chrx -n                 : shoe chx exon number
    chrx -pos -pos1 -pos2   : show chx from pos1 to pos2 sequence,pos is absolute
    chrx -num               : show chx no.num exon sequence
    chrx -num -pos1 -pos2   : show chx no.num from pos1 to pos2 sequence,pos is relative
    chrx -num -pis1 -pot2 -t: show chx no.num exon sequence,mark pos1 and pos2
    ```
3. 示例2
    ```
    >>chr1 -n
    ```
    ```
    chr1 = 19037
    ```
4. 示例
    ```
    >>chr1 -1
    ```
    ```
    chr1    12145   12310   CH1-1   166     200
    | TAACTTAATACCACAACCAGGCATAGGGGAAAGATTGGAGGAAAGATGAGTGAGAGCATCAACTTCTCTCACAACCTAGGCCAGTAAGTAGTGCTTGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTACCTGCCGTCTG |
    ```
5. 示例
    ```
    >>chr1 -1 -10 -20
    ```
    ```
    chr1    12145   12310   CH1-1   166     200
    | ACCACAACCAG |
    ```
6. 示例
    ```
    >>chr1 -pos -12000 -13000
    ```
    ```
    chr1    12145   12310   CH1-1   166     200
    | GTAACTTAATACCACAACCAGGCATAGGGGAAAGATTGGAGGAAAGATGAGTGAGAGCATCAACTTCTCTCACAACCTAGGCCAGTAAGTAGTGCTTGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTACCTGCCGTCTG |
    chr1    12595   12760   CH1-2   166     200
    | GCTCCTGTCTCCCCCCAGGTGTGTGGTGATGCCAGGCATGCCCTTCCCCAGCATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACTTGGATCACACTCTTGTGAGTGTCCCCAGTGTTGCAGAGGTGAGAGGAGAGTAGACAGTGAGTGGGAGTGGCGTCGCC |
    ```
7. 示例
    ```
    >>chr1 -1 -10 -20 -t
    ```
    ```
    chr1    12145   12310   CH1-1   166     200
    TAACTTAAT | ACCACAACCAG | GCATAGGGGAAAGATTGGAGGAAAGATGAGTGAGAGCATCAACTTCTCTCACAACCTAGGCCAGTAAGTAGTGCTTGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTACCTGCCGTCT
    ```

  


## 输出fastq文件

需要 外显子序列文件,碱基质量分布文件和突变设置 3个文件

1. 运行readout.py
   
    执行 readout.py ,需要外显子初始化文件,碱基质量文件和突变设置文件
    ```
    python3 readout.py exonlist phred.json mutations_setting.txt
    ```
    之后会打印突变,并开始生成数据
    ```
    **********Show Genetype**********
                ......
    **********Show Genetype**********
    press Enter to continue...
    ```
    接下来等待几小时,完成数据生成
    ```
    writing...
    ```

2. 生成的文件

        R1_fastq+日期.fastq 和 R2_fastq+日期.fastq


## 5. 设置

* 设置分为三个部分 : 全局others,初始化initial,数据生成readout
    ```
    '''others'''
    ROW_STEP=2 # 从行末端换行到下一行文本pos的移动,最好不要改动
    PHRED33=33 # ASCII+33 ,碱基质量与碱基字符的一种换算
    PHRED64=64 #  ASCII+33 ,碱基质量与碱基字符的一种换算
    PHRED=PHRED33 # 采用哪一种换算,就赋值给该变量
    CHIP_LEN=10
    ```
    ```
    '''initial.py setting'''
    MAXINSERT=200 # 外显子测序不仅会捕获外显子区域,也会捕获侧翼序列,maxinsert是允许的捕获的最大侧翼序列长度
    JOIN_GAP=200 # 如果两个外显子的侧翼区域有重叠,需要把这两个外显子合并为一个,以确保不会该区域不重复读写,join_gap默认要小于maxinsert.

    # 对于来自不同数据库的基因组和注释文件格式是有区别的,因此从'seq'获取染色体基因组,'ano'获取外显子注释文件,'list'初始化外显子,都要按照文件情况定义识别格式,下面有UCSC,NCBI,GENCODE已经测序公司给的格式,可以在CUSTOMS下自定义格式
    UCSC={
        'seq':('>chr1','>chr2','>chr3','>chr4','>chr5','>chr6','>chr7','>chr8','>chr9','>chr10',
        '>chr11','>chr12','>chr13','>chr14','>chr15','>chr16','>chr17','>chr18','>chr19',
        '>chr20','>chr21','>chr22','>chrX','>chrY','>chrM'),
        }

    NCBI={
        'seq':r'>NC',
        'ano':'exon',
        'list':r'NC_0*(\d*).*\sexon\s.*?(\d+)[\s\w]*?(\d+)'
        }
    GENCODE={
        'ano':'exon',
        'list':r'chr(\d*).*\sexon\s.*?(\d+)[\s\w]*?(\d+)'
        }
    COM={
        'list':r'chr([\w\d]*).*?(\d+)[\s\w]*?(\d+).*'
        }
    CUSTOMS={}
    # DEFAULT 是为了简化繁琐步骤设置的,您可以在这里先设置好参数,在initi.py文件中执行init -all即可全部执行
    DEFAULTS={
        'filew':'mutations_setting.txt', # 突变文件
        'filex':'NCBI_gh38.fna', # 基因组文件
        'file1':'NCBI_hg38ref.fna', # 染色体基因组文件
        'filey':'NCBI_hg38.gff', # 基因注释文件
        'file2':"NCBI_hg38exon.gff", # 外显子注释文件
        'file3':"NCBI_hg38list.txt", # 外显子初始化文件
        'filez':'qphred.fasq', # 碱基质量来源文件
        'ASCII':PHRED, # ASCII互换值
        'file4':'phred.json', # 碱基质量文件
        'ver':NCBI} # 染色体基因组合注释文件的格式类型
    ```
    ```
    '''readout.py setting'''
    ACCURACY_RATE= 0.8 # 机器准确率
    DEEPTH=100 #深度
    SUBSTITUTION={'A':('G','C','T','G'),'G':('A','C','T','A'),'C':('T','A','G','T'),'T':('C','G','A','C')} # 机器对某种碱基偏好误读情况
    INSERT_E=200 # insert的平均长度μ
    INSERT_D=15 # insert的标准差σ
    ERR_PH=3 # 误读碱基下降的质量
    ```

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

    减3

4. 输出的序列如果含N:
    
    先吐槽,我是真的不理解为什么那么多N到现在还测不出来,这些区域是有什么极其特别的部分吗,几十年了.

    如果N出现在侧翼区域,会缩小相应的两端侧翼长度以确保不读到N,检查会在initial.py中进行,修改后保存

    如果N出现在外显子区域,则会输出该外显子,不对外显子产生任何数据,检查会在readout.py中进行

    

5. insert长度
   
   insert长度是一个正态取值
   
   如果超出侧翼序列,只会取能取到的区间

   如果短于read长度,则重新获得一个insert,即输出的长度没有过短的序列

6. @名称

    对于正常序列 : 外显子id:突变编号:输出的第几个read

        exonid:mutation.no:turns:

    对于突变序列: 外显子id:突变编号:母方/父方(0/1):输出的第几个read:突变信息

        exonid:mutation.no:haplot:turns:mutationinfo


  
