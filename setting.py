
READ =100

ACCURACY_RATE= 0.8
DEEPTH=100
CHIP_LEN=10
SUBSTITUTION={'A':('G','C','T','G'),'G':('A','C','T','A'),'C':('T','A','G','T'),'T':('C','G','A','C')} 

ROW_STEP=2

PHRED=33
INSERT=200
'''exonlist_insert.py setting'''
MAXINSERT=200
JOIN_GAP=200
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
DEFAULTS={
    'filex':'NCBI_gh38.fna',
    'filey':'NCBI_hg38.gff',
    'file1':'NCBI_hg38ref.fna',
    'file2':"NCBI_hg38exon.gff",
    'file3':"NCBI_hg38list.txt",
    'fileq':'qphred.fasq',
    'ver':NCBI}