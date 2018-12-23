import json
import random
import re
import sys
import time

import mutation
from filefunc import *
from mutation import *
from sequence import *

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
SUBSTITUTION = {'A': ('G', 'C', 'T', 'G'), 'G': ('A', 'C', 'T', 'A'), 'C': (
    'T', 'A', 'G', 'T'), 'T': ('C', 'G', 'A', 'C')}


def get_fre_qph(file, write):
    print('getting phred frequencies from %s...' % file)
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
    print('\n\tdown. outfile: %s' % write)


def get_fre_json(file):
    print('initial qphred from %s...' % file, end='')
    with open(file, "r") as filed:
        frequenciess = json.load(filed)
    print('\n\tdown.')
    return frequenciess


def get_random_qph(frequencies, num):
    qphreds = []
    summ = sum(frequencies['pos1_frequencies'].values())
    for n in range(num):
        s = ''
        for x in range(1, len(frequencies)+1):
            char = random_weight_choice(
                frequencies['pos%d_frequencies' % x], summ)
            s = s+char
        qphreds.append(s)
    return qphreds


def complementation(read):
    l = []
    for char in read:
        l.append(COMPLEMENT[char.upper()])
    return l


def add_read_err(read, phred, *, accuracy_e, accuracy_d, err, **default):
    ''' read,phred is list'''
    size = len(read)
    accuracy_e = value2unnegtive(random.normalvariate(accuracy_e, accuracy_d))
    err_num = value2unnegtive(round((1-accuracy_e)*size))
    err_num_list = random.sample(range(size), err_num)
    for n in err_num_list:
        pich = random.randint(0, 2)
        read[n] = SUBSTITUTION[read[n].upper()][pich]
        phred_value = ord(phred[n])-err
        if phred_value > 31:
            phred[n] = chr(phred_value)
        else:
            phred[n] = chr(32)


def wes_insert(sequence, frequencies, readlen, *, chip_len, insert_e, insert_d, **default):
    ''' get one inser's 2 reads'''
    # chip pos on exon
    length = sequence.length
    seq = sequence.seq
    flank_len = sequence.flank_len
    cbp = random.randint(0, length-chip_len)
    while(True):
        # insertion length
        insertion = round(random.normalvariate(insert_e, insert_d))
        if insertion < readlen:
            continue
        # chip pos on inseriton
        sbp = random.randint(0, insertion-chip_len)
        # insertion pos on exon
        lbp = flank_len+cbp-sbp
        if lbp < 0:
            insertion = lbp+insertion
            lbp = 0
        rbp = lbp+insertion-1
        if rbp > length+2*flank_len-1:
            insertion = length+2*flank_len-rbp-1+insertion
            rbp = length+2*flank_len-1
        if insertion >= readlen:
            break
    return seq[lbp:rbp], lbp, rbp


def PE_reads(writed1, writed2, insert, frequencies, readlen, sequence_id, lbp, rbp, x, **default):
    l = []
    for char in insert:
        if char in ATCG:
            l.append(char.upper())
        elif char.upper() in DEGENERATE:
            l.append(random.choice(DEGENERATE[char.upper()]))
        else:
            raise Exception("error :can't indentify", char)
            # l.append(random.choice(ATCG))
    read1 = l[0:readlen]  # list
    read2 = complementation(l[-1:-1-readlen:-1])
    phred1 = list(random.choice(frequencies))
    phred2 = list(random.choice(frequencies))
    add_read_err(read1, phred1, **default)
    add_read_err(read2, phred2, **default)
    read1, read2, phred1, phred2 = ''.join(read1), ''.join(
        read2), ''.join(phred1), ''.join(phred2)
    writed1.write('@SIMU:%s:%s:%d:%d:%s/1\n' %
                  (sequence_id, x, lbp, rbp, default['muta']))
    writed1.write(read1+'\n')
    writed1.write('+\n')
    writed1.write(phred1+'\n')  # random_qphred(frequencies)
    writed2.write('@SIMU:%s:%s:%d:%d:%s/2\n' %
                  (sequence_id, x, lbp, rbp, default['muta']))
    writed2.write(read2+'\n')
    writed2.write('+\n')
    writed2.write(phred2+'\n')


def wes_fastq(writed1, writed2, sequence, frequencies, tissue, readlen, avglen, depth, **default):
    ''' one exon's all fastq reault be writed to w'''
    insert_e, chip_len = default['insert_e'], default['chip_len']

    def h(x):
        return x*(x+1)*(x+2)/3/insert_e/insert_e

    def get_turns(length):
        if length+1 > insert_e:
            a = h(insert_e-chip_len)/(length-chip_len+1)
        else:
            a = (h(insert_e-chip_len)-h(insert_e-length-1))/(length-chip_len+1)
        return length*depth/(1-a)/avglen
    sequence_id = sequence.id
    if sequence.seq_length < readlen:
        print('\tnot output , error:seq length too short \n'+sequence.get_info())
        return
    if sequence.length < chip_len:
        print('\tnot output , errpr:aimed seq length too short \n'+sequence.get_info())
        return
    print('\twriting chr %s - exon : %s - %s' %
          (sequence.chr, sequence_id, sequence.begin), end='\r')
    border = Border(sequence.chr, sequence.begin -
                    sequence.flank_len, sequence.end+sequence.flank_len)
    if 'copy_num' in default:
        turns = round(get_turns(sequence.length)*default['copy_num'])
    else:
        copy_num = tissue.get_normal_num(border)
        length = sequence.length
        turns = round(get_turns(length)*copy_num)
    for x in range(turns):
        insert, lbp, rbp = wes_insert(
            sequence, frequencies, readlen, **default)
        PE_reads(writed1, writed2, insert, frequencies,
                 readlen, sequence_id, lbp, rbp, x, **default)


def wgs_fastq(writed1, writed2, sequence, frequencies, tissue, readlen, avglen, depth, **default):
    insert_e, insert_d = default['insert_e'], default['insert_d']
    sequence_id = sequence.id
    if sequence.seq_length < readlen:
        print('\tnot output , error:seq length too short \n'+sequence.get_info())
        return
    print('\twriting chr %s - exon : %s - %s' %
          (sequence.chr, sequence_id, sequence.begin), end='\r')
    if 'copy_num' in default:
        turns = round(depth*default['copy_num'])
    else:
        border = Border(sequence.chr, '.', '.')
        copy_num = tissue.get_normal_num(border)
        length = sequence.length
        turns = round(depth*copy_num)
    for x in range(turns):
        a = 0
        b = a+round(random.normalvariate(insert_e, insert_d))
        end = sequence.length
        seq = sequence.seq
        while(b < end):
            if b-a >= readlen:
                PE_reads(
                    writed1, writed2, seq[a:b], frequencies, readlen, sequence_id, a, b, x, **default)
            a, b = b, a+round(random.normalvariate(insert_e, insert_d))
        if end-a >= readlen:
            PE_reads(writed1, writed2,
                     seq[a:b], frequencies, readlen, sequence_id, a, end, x, **default)


def readout(arith, *, file1, file2, file3, file4, cd, extra, flank_len, join_gap, filtrate, **default):
    # print info
    insert_e = default['insert_e']
    print('illustrate :', 'cd=', cd, 'extrs=', extra)
    print('paremeter set :', default)
    if arith == 'wes':
        print('getting pe fastq from %s , %s , %s , %s ...' %
              (file1, file2, file3, file4))
    elif arith == 'wgs':
        print('getting pe fastq from %s , %s , %s ...' % (file1, file3, file4))

    # set and print gene mutation type
    tissue = Tissue.formula_mutation(file4)
    tissue.show_tissue()

    # initial qphred frequencies
    frequencies = get_fre_json(file3)
    readlen = len(frequencies)
    avglen = min(2*readlen, insert_e)
    # get mutation sequence from file4(mutation file),file2(bed file)
    fileb = file4+'.bed'
    filec = file4+'.fna'
    filee = file4+'.temp'
    chrss, column, row_step = Sequen.fasta_file_info(file1)
    info = {'chrss': chrss, 'column': column, 'step': row_step}
    # get normal aimed sequence from file1(genome file),file2(bed file)
    if arith == 'wes':
        filea = file2+'.fna'
        Sequen.seq_fnabed(file1, file2, filea, flank_len,
                          filtrate, insert_e, **info)
        # get need bed from file2
        write_wes_bed(file2, fileb, tissue, flank_len,
                      join_gap, chrss=chrss)
    elif arith == 'wgs':
        filea = file1+'.fil'
        Sequen.file_filtrate(file1, insert_e)
        write_wgs_bed(fileb, tissue, flank_len, join_gap, chrss=chrss)
    # get need fasta meterials from file1 and fileb
    Sequen.seq_fnabed(file1, fileb, filee, 0, False, **info)

    # get fastq , output reads
    # output mutation reads
    time_start = time.time()
    t = time.strftime('%Y%m%d_%H_%M_%S', time.localtime(time.time()))
    R1 = "%sR1%s_%s.fastq" % (cd, extra, t)
    R2 = "%sR2%s_%s.fastq" % (cd, extra, t)
    write_matation_sequences(filee, R1, R2, filec, arith,
                             frequencies, tissue, readlen, avglen, column, filtrate, **default)
    write_normal_sequences(filea, R1, R2, arith, frequencies,
                           tissue, readlen, avglen, **default)
    print('\nall down. outfile :')
    print(R1)
    print(R2)
    time_end = time.time()
    t = time_end-time_start
    print('totally cost: %dh : %dm : %ds' % (t//3600, (t % 3600)//60, t % 60))


def write_matation_sequences(file, write1, write2, write3, arith, frequencies, tissue, readlen, avglen, column, filtrate, **default):
    insert_e = default['insert_e']
    if arith == 'wes':
        inner = False
    elif arith == 'wgs':
        inner = True
    with open(write1, 'a', newline='\n') as w1:
        with open(write2, 'a', newline='\n')as w2:
            with open(file, 'r', newline='\n')as filed:
                print('outputting mutation reads...')
                with open(write3, 'w', newline='\n') as w3:
                    if arith == 'wes':
                        for x in sorted(tissue.polyploids.keys()):
                            info = {'no': x}
                            x = tissue.polyploids[x]
                            copy_num = x.percent/len(x.genomes)
                            for i, y in enumerate(x.genomes):
                                info['genome'] = i
                                print('\tgenome', '-', info['no'], '-', i, ':')
                                for sequence in y.seq(filed, **info):
                                    n = sequence.length//100
                                    phreads = get_random_qph(frequencies, n+10)
                                    if filtrate:
                                        sequences = sequence.filtrate(
                                            insert_e, inner)
                                    else:
                                        sequences = [sequence]
                                    for y in sequences:
                                        y.write(w3, column)
                                        wes_fastq(w1, w2, y, phreads,
                                                  tissue, readlen, avglen, **default, copy_num=copy_num)
                                        y.del_seq()
                    elif arith == 'wgs':
                        for x in sorted(tissue.polyploids.keys()):
                            info = {'no': x}
                            x = tissue.polyploids[x]
                            copy_num = x.percent/len(x.genomes)
                            for i, y in enumerate(x.genomes):
                                info['genome'] = i
                                print('\tgenome', '-', info['no'], '-', i, ':')
                                for sequence in y.seq(filed, **info):
                                    n = sequence.length//100
                                    phreads = get_random_qph(frequencies, n+10)
                                    if filtrate:
                                        sequences = sequence.filtrate(
                                            insert_e, inner)
                                    else:
                                        sequences = [sequence]
                                    for y in sequences:
                                        y.write(w3, column)
                                        wes_fastq(w1, w2, y, phreads,
                                                  tissue, readlen, avglen, **default, copy_num=copy_num)
                                        y.del_seq()
    print('\tdown. outfile : ', write3)


def write_normal_sequences(file, write1, write2, arith, frequencies, tissue, readlen, avglen, **default):
    with open(write1, 'a', newline='\n') as w1:
        with open(write2, 'a', newline='\n')as w2:
            with open(file, 'r')as filed:
                print('outputting normal reads...')
                if arith == 'wes':
                    for sequence in Sequen.sequences_iterator(filed):
                        # creat qphread
                        n = sequence.length//100
                        phreads = get_random_qph(frequencies, n+10)
                        # write fastq
                        wes_fastq(w1, w2, sequence, phreads, tissue,
                                  readlen, avglen, **default)
                elif arith == 'wgs':
                    for sequence in Sequen.sequences_iterator(filed):
                        # creat qphread
                        n = sequence.length//100
                        phreads = get_random_qph(frequencies, n+10)
                        # write fastq
                        wgs_fastq(w1, w2, sequence, phreads, tissue,
                                  readlen, avglen, **default)
    print('\tdown.                                 ')
