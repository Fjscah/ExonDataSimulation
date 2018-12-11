import json
import random
import re
import sys
import time
import numpy

import mutation
from exon import get_exons
from read import fastq_exon, random_qphred
from setting import (ACCURACY_RATE, ACCURACY_RATE_D, DEFAULTS, DEPTH, ERR_PH,
                     INSERT_D, INSERT_E, JOIN_GAP, MAXINSERT)

if __name__ == '__main__':
    random.seed()
    ename, qname, mname = DEFAULTS['file3'], DEFAULTS['file4'], DEFAULTS['filew']
    cd = './'
    if len(sys.argv) >= 4:
        pyname, ename, qname, mname = sys.argv[:4]
    if len(sys.argv) >= 5:
        cd = sys.argv[4]

    # set and print gene mutation type
    mutation_types = mutation.file_mutation(mname)
    mutation.show_mutation(mutation_types)

    # initial qphred frequencies
    print('initial qphred...')
    with open(qname, "r") as f:
        frequenciess = json.load(f)
        readlen = len(frequenciess)
        summ = sum(frequenciess['pos1_frequencies'].values())
    print('down...')

    # according to mutaton output reads
    time_start = time.time()
    avglen = min(2*readlen, INSERT_E)
    print('output normal reads(cnv+deletion)...')
    t = time.strftime('%Y%m%d_%H%M%S', time.localtime(time.time()))
    with open(r"%s/info_%s.txt" % (cd, t), 'w', newline='\n') as w:
        w.writelines(['ename : ', ename, '\n',
                      'qname :', qname, '\n',
                      'mname :', mname, '\n',
                      'maxinsert : ', str(MAXINSERT), '\n',
                      'joingap : ', str(JOIN_GAP), '\n',
                      'insert_e : ', str(INSERT_E), '\n',
                      'avglen : ', str(avglen), '\n',
                      'accuracy : ', str(ACCURACY_RATE), '\n',
                      'depth : ', str(DEPTH), '\n',
                      'insert_d : ', str(INSERT_D), '\n',
                      'accuracy_d : ', str(ACCURACY_RATE_D), '\n',
                      'readlen : ', str(readlen), '\n',
                      'err_ph : ', str(ERR_PH), '\n'])
    #effects = []
    with open(r"%s/R1_%s.fastq" % (cd, t), 'w', newline='\n') as w1:
        with open(r"%s/R2_%s.fastq" % (cd, t), 'w', newline='\n') as w2:
            with open(ename, 'r')as r:
                for exon in get_exons(r):
                    # creat qphread
                    n = exon.length//100
                    frequencies = []
                    for x in range(10+n):
                        frequencies.append(random_qphred(frequenciess, summ))
                    # write fastq
                    fastq_exon(w1, w2, exon, frequencies,
                                        mutation_types, readlen, avglen)
                    #effects.append(effect)
    #mean = numpy.mean(effects)
    #std = numpy.std(effects, ddof=1)
    #print('\n', 'mean = ', mean, '\n', 'std= ', std)
    print('\n***down***')
    time_end = time.time()
    t = time_end-time_start
    print('totally cost: %dh : %dm : %ds' % (t//3600, (t % 3600)//60, t % 60))
