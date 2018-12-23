
VER = {
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
    'TEST':{
        'bed': r'([\w]*)\t(\d+)\t(\d+)'
    }
}


DEFAULT = {
    'filex': 'genome',
    'file1': 'genome-ini',
    'filey': 'annotation',
    'file2': 'annotation-ini',
    'filez': 'phred',
    'file3': 'phred-ini',
    'file4': 'mutation',
    'chip_len': 10,
    'insert_e': 150,
    'insert_d': 0,
    'join_gap': 200,
    'flank_len': 0,
    'accuracy_e': 1.0,
    'accuracy_d': 0.0,
    'depth': 50,
    'err': 3,
    'ver': VER['NCBI'],
    'cd': '',
    'extra': '',
    'filtrate': True,
    'arithmetic': '',
    'muta':'',
}
