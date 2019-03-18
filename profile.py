

import time
from enum import Enum

t = time.strftime('%Y%m%d_%H_%M_%S', time.localtime(time.time()))
COMPLEMENT = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
ATCG = ('A', 'T', 'C', 'G', 'a', 't', 'c', 'g')
DEGENERATE = {'W': ['A', 'T'],
              'S': ['C', 'G'],
              'R': ['A', 'G'],
              'Y': ['C', 'T'],
              'K': ['G', 'T'],
              'M': ['A', 'C'],
              'B': ['C', 'G', 'T'],
              'D': ['A', 'G', 'T'],
              'H': ['A', 'C', 'T'],
              'V': ['A', 'C', 'G'],
              'N': ['A', 'C', 'G', 'T'],
              }
SUBSTITUTIONS = {'A': ('G', 'C', 'T', 'G'), 'G': ('A', 'C', 'T', 'A'), 'C': (
    'T', 'A', 'G', 'T'), 'T': ('C', 'G', 'A', 'C')}


FORMAT = {
    'UCSC': {
        'fna': r'>chr\w{1,2}\s*\n',
    },
    'NCBI': {
        'fna': r'>NC',
        'bed': r'NC_0*(\d*).*\sexon\s.*?(\d+)[\s\w]*?(\d+)'
    },
    'GENCODE': {
        'bed': r'chr(\d*).*\sexon\s.*?(\d+)[\s\w]*?(\d+)'
    },
    'COM': {
        'bed': r'chr([\w\d]*).*?(\d+)[\s\w]*?(\d+).*'
    },
    'CUSTOMS': {

    },
    'TEST': {
        'bed': r'([\w]*)\t(\d+)\t(\d+)'
    }
}
modes = Enum('mode', ('WES', 'WGS'))
pairs = Enum('pair', ('SE', 'PE', ))
mutation_ways = (Enum('mutation_way', ('table', 'formula', 'auto')))
qphs = Enum('qph', ('sanger', 'solexa', ))
'''input set'''
REFERENCES = [  # reference seuence , every element is one ref genome
    'GRCh38_latest_genomic.fna',
]

REFMAT = [  # regular expression for gaining data from references
    FORMAT['NCBI'],
    FORMAT['NCBI'],
]
REGIONS = [  # regions targeted by probes
    'GRCh38_latest_genomic.gff',
]
REGMAT = [  # regular expression for gaining data from regions
    FORMAT['NCBI'],
]
POLYS = [  # content
    (1, 1),
]
MUTATIONS = [  # mutation files
    'mutation.txt',
]
QUALITY = 'source100.fastq'  # row fastq file for generating quality score profile
MUTAINPUT = mutation_ways.formula
QPH = qphs.sanger
'''lenght set'''
CHIP_LEN = 10                   # probe length
SEGMENT_E = 150                 # expectation of DNA segment length
SEGMENT_D = 0                   # variance of DNA segment length
# if the distance of two neighboring regionis is smaller than join_gap,then merge them
JOIN_GAP = 10
# the lenth of irrelative region targetrd by probe (one side)
FLANK_LEN = 200

'''error set'''
ERROR_E = 0.01                  # exception of error rate
ERROR_D = 0.00                  # variance of error rate
SUBSTITUTION = 0.8              # conditional probability of substitions
DELETION = 0.1                  # conditional probability of deletions
INSERTION = 0.1                # conditional probability of deletions
# for signal mismatch , quality score reduce 'error_cut'
ERROR_CUT = 3

'''read set'''
DEPTH = 50                    # sequence depth
MODE = modes.WES
PAIR = pairs.PE


'''output set'''
CD = ''                           # path to store
SEED = 4000000
LEGAL_N = False                   # whether wipe off N
INNER_N = True                    # whether wipe off  aimed region with
COLUMN = 100
MEMORY = 10000
BED_INFO = 'bed_info'
FASTA_INFO = 'fasta_info'
ELEN = 5
R1 = CD+"R1_%s.fastq" % (t)
R2 = CD+"R2%s_%s.fastq" % (t)
R0 = CD+"R%s_%s.fastq" % (t)
