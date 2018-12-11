import re
import sys

import view
from exon import  get_wgs, get_bed,get_wes
from read import get_phred_fre
from setting import *

print("default:", "INSERT maxlength=", MAXINSERT, "CHIP length=", CHIP_LEN)


def get_all(defaults):
    filex, filey = defaults['filex'], defaults['filey']
    file1, file2, file3 = defaults['file1'], defaults['file2'], defaults['file3']
    filez, file4 = defaults['filez'], defaults['file4']
    ver = defaults['ver']
    get_wgs(filex, file1, ver)
    get_bed(filey, file2, ver)
    get_wes(file1, file2, file3,CHIP_LEN,
                  JOIN_GAP, ROW_STEP, MAXINSERT)
    get_phred_fre(filez, file4)


if __name__ == '__main__':
    help = '''
    -wgs filex file1 ver        : get whole genome from filex , generate file3
    -bed filey file2 ver        : get exon annotation from filey, generate file2
    -wes file1 file2 file3      : get exonlist from file1,file2, generate file3
    -qph filez file4            : get qphred frequencies from file,generate file4
    -all                        : excute above all according setting,py
    -view file3                : view exonlist's exon sequence
    '''

    info = sys.argv
    if 'help' in info:
        print(help)
    elif '-wgs' in info:
        if len(info)>=5:
            filex, file1, ver = info[2:]
            get_wgs(filex,file1,ver)
    elif '-bed'in info:
        if len(info) >= 5:
            fiely, file2, ver = info[2:]
            get_bed(fiely, file2, ver)
    elif '-wes' in info:
        if len(info) >= 5:
            file1, file2, file3,= info[2:]
            get_wes(file1, file2, file3, CHIP_LEN,
                          JOIN_GAP, ROW_STEP, MAXINSERT)
    elif '-all' in info:
        get_all(DEFAULTS)
    elif '-qph' in info:
        if len(info) >= 4:
            filex, file4 = info[2:]
            get_phred_fre(filex, file4)
    elif '-view' in info:
        if len(info) >= 2:
            file = info[1]
            view.view(file)
