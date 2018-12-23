

import sys

from read import *
from sequence import *
from setting import *

if __name__ == '__main__':
    help = '''
    -fna <filex> <file1> <ver>          : get whole genome from filex , generate file3
    -bed <filey> <file2> <ver>          : get aimed and sorted annotation from filey, generate file2
    -qph <filez> <file3>                : get qphred frequencies from file,generate file3
    -wes <file1> <file2> <file3> <file4> : get aimed regions pe fastq from file1234 (deafault flitrate 'N')
                                        : file4 is mutation file, generate fastq1+2
    -wgs <file1> <file3> <file4>        : get whole region pe fastq from file124, generate fastq1+2
    -seq <file1> <file2> <filea>        : get aimed sequence from file1,file2, generate filea
    -fli <filea>                        : flitrate 'N' sequence, generate fileb(=filea+'.fliter')
    -view <filea>                       : view  sequence
    '''
    info = sys.argv
    print(info)
    if 'help' in info[1]:
        print(help)
    elif 'fna' in info[1]:
        if len(info) >= 5:
            DEFAULT['filex'], DEFAULT['file1'], DEFAULT['ver'] = info[2:]
            Fasta.ini_fna(**DEFAULT)
    elif 'bed'in info[1]:
        if len(info) >= 5:
            DEFAULT['filey'], DEFAULT['file2'], DEFAULT['ver'] = info[2:]
            Bed.ini_bed(**DEFAULT)
    elif 'qph' in info[1]:
        if len(info) >= 4:
            filex, file3 = info[2:]
            get_fre_qph(filex, file3)
    elif 'wes' in info[1]:
        cd = extra = ''
        if len(info) >= 6:
            DEFAULT['file1'],DEFAULT['file2'],DEFAULT['file3'],DEFAULT['file4']=info[2:6]
        if len(info) >= 7:
            DEFAULT['cd']=info[6]
        if len(info)>=8:
            DEFAULT['extra']=info[7]
        readout('wes',**DEFAULT)
    elif 'wgs' in info[1]:
        cd = extra = ''
        if len(info) >= 5:
            DEFAULT['file1'],DEFAULT['file3'],DEFAULT['file4']=info[2:5]
        if len(info) >= 7:
            DEFAULT['cd']=info[6]
        if len(info)>=8:
            DEFAULT['extra']=info[7]
        readout('wgs',**DEFAULT)
    elif 'seq' in info[1]:
        if len(info) >= 5:
            file1, file2, file4, = info[2:5]
        if '-t'in info:
            DEFAULT['filtrate']=True
        if '-f'in info:
            DEFAULT['filtrate']=False
        flank_len,filtrate,insert_e=DEFAULT['flank_len'],DEFAULT['filtrate'],DEFAULT['insert_e']
        Sequen.seq_fnabed(file1, file2, file4,flank_len,filtrate,insert_e)
    elif 'view' in info[1]:
        if len(info) >= 3:
            file = info[2]
            view(file)
    elif 'fli' in info[1]:
        if len(info) >= 3:
            file4 = info[2]
            Sequen.file_filtrate(file4,DEFAULT['insert_e'])
