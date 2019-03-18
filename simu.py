

import os
import sys
import time
from functools import reduce
from profile import *

from mutation import *
from read import *
from sequence import *

'''
import argparse
parser=argparse.ArgumentParser()
parser.parse_args()
'''
if __name__ == '__main__':
    time_start = time.time()
    help = '''
    -ref          : get whole genome from filex , generate file3
    -reg          : get aimed and sorted annotation from filey, generate file2
    -qph          : get qphred frequencies from file,generate file3
    -read         : get aimed regions pe fastq from file1234 (deafault flitrate 'N')
    -seq          : get aimed sequence from file1,file2, generate filea
    -fli          : flitrate 'N' sequence, generate fileb(=filea+'.fliter')
    -view <file>  : view  sequence
    optional:
    -wes/wgs      : set the way to generate DNA segment
    -se/pe        : set the way to generate reads
    '''
    info = sys.argv
    print(info)

    if '-wes' in info:
        MODE = modes.WES
    if '-wgs' in info:
        MODE = modes.WGS
    if '-pe' in info:
        PAIR = pairs.PE
    if '-se' in info:
        PAIR = pairs.SE
    regions = set(REGIONS)  # filtrate the same region file
    inireferences = list(map(lambda x: CD+x+'.ini', REFERENCES))
    iniregions = list(map(lambda x: CD+x+'.ini', regions))
    iniexomes = list(map(lambda x: CD+x+'.exome', REFERENCES))
    inimutations = list(map(lambda x: CD+x+'.ini', MUTATIONS))
    stanfiles = list(map(lambda x: CD+x+'.temp', REGIONS))
    iniquality = CD+QUALITY+'.ini'
    if not os.path.exists(BED_INFO):
        open(BED_INFO, 'a').close()
    if not os.path.exists(FASTA_INFO):
        open(FASTA_INFO, 'a').close()
    if '-help' in info:
        print(help)
    elif '-ref' in info:  # initialize reference genomes(single haploid)
        n = 1
        for x, y in zip(REFERENCES, inireferences):  # mark number
            Fasta.ini_ref(x, y, n)
            n += 1
    elif '-reg'in info:  # merge targed regions
        for x, y in zip(regions, iniregions):
            keys = Bed.ini_reg(x, y)
    elif '-qph' in info:
        Quality.ini_qph(QUALITY, iniquality)
    elif '-seq' in info:  # initailize exome sequence from initialized regions
        for x, y, z in zip(inireferences, iniregions, iniexomes):
            Fasta.ini_exome(x, y, z)
    elif '-mut' in info:  # initialize mutations
        for x, y, poly in zip(MUTATIONS, inimutations, POLYS):
            ini_muta(inireferences, iniregions, x, y, poly)
    elif '-read' in info:
        if pairs.PE == PAIR:
            open(R1, 'w', newline='\n').close()
            open(R2, 'w', newline='\n').close()
        elif PAIR == pairs.SE:
            open(R0, 'w', newline='\n').close()
        for x, poly in (inimutations, POLYS):
            # polyoid 's id,content,mutationseq,inireferences,mutationsbed
            # mutationsbed can get from iniregions+polys
            readout(inireferences, iniregions, iniquality, x, poly)
    elif '-view' in info:
        if len(info) >= 3:
            file = info[2]
            view(file)

    time_end = time.time()
    t = time_end-time_start
    print('totally cost: %dh : %dm : %ds' % (t//3600, (t % 3600)//60, t % 60))
