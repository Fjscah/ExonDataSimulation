import json
import random
import re

import mutation
from exon import Exon
from filefunc import get_column_row, num_positive
from setting import (ACCURACY_RATE, ACCURACY_RATE_D, CHIP_LEN, DEPTH,
                     DEFAULTS, ERR_PH, INSERT_D, INSERT_E, JOIN_GAP, MAXINSERT,
                     SUBSTITUTION)

ATCG = ('A', 'T', 'C', 'G', 'a', 't', 'c', 'g')
COMPLEMENT = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
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
              # 'N':['A','C','G','T'],
              }


def get_phred_fre(file, write):
    print('get phred frequencies from %s...' % file)
    frequencies = {}
    readlen = get_column_row(file, r'+', 1, '==')
    for x in range(1, readlen+1):
        frequencies['pos%d_frequencies' % x] = {}
    with open(file, 'r') as f:
        i = 0
        line = f.readline()
        while(line):
            i += 1
            if i % 4 == 0:
                row_fastq = line.strip()
                x = 1
                for char in row_fastq:
                    if char in frequencies['pos%d_frequencies' % x]:
                        frequencies['pos%d_frequencies' % x][char] += 1
                    else:
                        frequencies['pos%d_frequencies' % x][char] = 1
                    x = x+1
            if i % 40000 == 0:
                print(i, end='\r')
            if i > 4000000:
                break
            line = f.readline()
    with open(write, "w") as f:
        json.dump(frequencies, f)
    print('\ndown. outfile: %s' % write)


def random_qphred(frequencies, summ):
    s = ''
    for x in range(1, len(frequencies)+1):
        char = random_weight_choice(frequencies['pos%d_frequencies' % x], summ)
        s = s+char
    return s


def random_weight_choice(dictt, summ):
    t = random.randint(1, summ)
    for char, val in dictt.items():
        t -= val
        if t <= 0:
            return char
    if t > 0:
        print('kkkkkkkkkkkkkkkkk')


def add_read_err(read, phred, accuracy=ACCURACY_RATE, accuracy_d=ACCURACY_RATE_D):
    size = len(read)
    accuracy = num_positive(random.normalvariate(accuracy, accuracy_d))
    err_num = num_positive(round((1-accuracy)*size))
    err_num_list = random.sample(range(size), err_num)
    for n in err_num_list:
        pich = random.randint(0, 2)
        read[n] = SUBSTITUTION[read[n].upper()][pich]
        phred_value = ord(phred[n])-ERR_PH
        if phred_value > 31:
            phred[n] = chr(phred_value)
        else:
            phred[n] = chr(32)


def complementation(read):
    l = []
    for char in read:
        l.append(COMPLEMENT[char.upper()])
    return l


def PEread(seq, length, frequencies, readlen, exoninsert):
    ''' get one inser's 2 reads'''
    # chip pos on exon
    cbp = random.randint(0, length-CHIP_LEN)
    while(True):
        # insertion length
        insertion = round(random.normalvariate(INSERT_E, INSERT_D))
        if insertion < readlen:
            continue
        # chip pos on inseriton
        sbp = random.randint(0, insertion-CHIP_LEN)
        # insertion pos on exon
        lbp = exoninsert+cbp-sbp
        if lbp < 0:
            insertion = lbp+insertion
            lbp = 0
        rbp = lbp+insertion-1
        if rbp > length+2*exoninsert-1:
            insertion = length+2*exoninsert-rbp-1+insertion
            rbp = length+2*exoninsert-1
        if insertion >= readlen:
            break
    l = []
    for n in range(lbp, rbp+1):
        char = seq[n]
        if char in ATCG:
            l.append(char.upper())
        elif char.upper() in DEGENERATE:
            l.append(random.choice(DEGENERATE[char.upper()]))
        else:
            raise Exception("error :can't indentify", char)
            # l.append(random.choice(ATCG))
    #effct = num_positive(min(exoninsert+length, exoninsert+readlen)-max(
    #    lbp, exoninsert)+min(rbp+1, exoninsert+length)-max(exoninsert, rbp-readlen+1))
    read1 = l[0:readlen]  # list
    read2 = complementation(l[-1:-1-readlen:-1])  # list
    phred1 = list(random.choice(frequencies))
    phred2 = list(random.choice(frequencies))
    add_read_err(read1, phred1)
    add_read_err(read2, phred2)
    return ''.join(read1), ''.join(read2), ''.join(phred1), ''.join(phred2), lbp, rbp


def write_reads(writed1, writed2, exon_id, exoninsert, seq, length, frequencies, readlen, x, y, z=2, onemutation='0'):
    read1, read2, phred1, phred2, lbp, rbp = PEread(
        seq, length, frequencies, readlen, exoninsert)
    writed1.write('@SIMU:%s:%s:%d:%d:%d:%d:%s/1\n' %
                  (exon_id, x, z, y, lbp, rbp, onemutation))
    #writed1.write('@%s/1\n'% (exon_id))
    writed1.write(read1+'\n')
    writed1.write('+\n')
    writed1.write(phred1+'\n')  # random_qphred(frequencies)
    writed2.write('@SIMU:%s:%s:%d:%d:%d:%d:%s/2\n' %
                  (exon_id, x, z, y, lbp, rbp, onemutation))
    #writed2.write('@%s/2\n'% (exon_id))
    writed2.write(read2+'\n')
    writed2.write('+\n')
    writed2.write(phred2+'\n')
    # return effect


def fastq_exon(writed1, writed2, exon, frequencies, mutation_types, readlen, avglen):
    def h(x):
        return x*(x+1)*(x+2)/3/INSERT_E/(INSERT_E-CHIP_LEN+1)

    def get_turns():
        if length+1 > INSERT_E:
            a = h(INSERT_E-CHIP_LEN)/(length-CHIP_LEN+1)
        else:
            a = (h(INSERT_E-CHIP_LEN)-h(INSERT_E-length-1))/(length-CHIP_LEN+1)
        return length*DEPTH/2/(1-a)/avglen
    exon_id = exon.exon_id
    ''' one exon's all fastq reault be writed to w'''
    print('writing chr %s - exon : %s - %s' %
          (exon.chr, exon_id, exon.begin), end='\r')
    if 'N' in exon.seq or exon.length+2*exon.insert < readlen:
        print('error: not output \n'+exon.getexon_info()+'\n'+''.join(exon.seq))
        return
    seq, length = list(exon.seq), exon.length
    turns = get_turns()
    #effect_t = 0
    for x in mutation_types:
        copy_num = mutation_types[x].haplots[0].get_normal_num(
            exon_id)+mutation_types[x].haplots[1].get_normal_num(exon_id)
        turns = round(turns*copy_num*mutation_types[x].percent/100)
        for y in range(turns):
            write_reads(writed1, writed2, exon_id, exon.insert, seq,
                        length, frequencies, readlen, x, y)
            #effect_t += effect
    for x in mutation_types:
        for z in range(2):
            if exon_id in mutation_types[x].haplots[z].mutations:
                for seq, length, copy_num, onemutation in exon.get_new_seq((mutation_types[x].haplots[z].mutations[exon_id].abnormal_mutations), CHIP_LEN):
                    turns = round(get_turns()*copy_num *
                                  mutation_types[x].percent/100)
                    seq = list(seq)
                    for y in range(turns):
                        write_reads(writed1, writed2, exon_id, exon.insert, seq, length,
                                    frequencies, readlen, x, y, z=z, onemutation=onemutation)
    #return effect_t//length
