'''others'''
ROW_STEP = 1
CHIP_LEN = 10
INSERT_E = 300
MAXINSERT = 300
JOIN_GAP = 200
'''initial.py setting'''
UCSC = {
    'wgs': r'>chr\w{1,2}\s*\n',
}
NCBI = {
    'wgs': r'>NC',
    'bed': r'NC_0*(\d*).*\sexon\s.*?(\d+)[\s\w]*?(\d+)'
}
GENCODE = {
    'bed': r'chr(\d*).*\sexon\s.*?(\d+)[\s\w]*?(\d+)'
}
COM = {
    'bed': r'chr([\w\d]*).*?(\d+)[\s\w]*?(\d+).*'
}
CUSTOMS = {}

DEFAULTS = {
    'filew': 'mutations_setting.txt',
    'filex': 'NCBI_gh38.fna',
    'file1': 'NCBI_hg38ref.fna',
    'filey': 'NCBI_hg38.gff',
    'file2': "NCBI_hg38exon.gff",
    'file3': "NCBI_hg38list.txt ",  
    'filez': 'phred.fasq',
    'file4': 'phred.json',
    'ver': NCBI}
'''readout.py setting'''
ACCURACY_RATE = 1.0
ACCURACY_RATE_D = 0.0
DEPTH = 50
SUBSTITUTION = {'A': ('G', 'C', 'T', 'G'), 'G': ('A', 'C', 'T', 'A'), 'C': (
    'T', 'A', 'G', 'T'), 'T': ('C', 'G', 'A', 'C')}
INSERT_D = 15
ERR_PH = 3
