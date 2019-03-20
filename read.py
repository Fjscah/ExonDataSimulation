import json
import random
import re
import sys
import time

from basic import *
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
SUBSTITUTIONS = {'A': ('G', 'C', 'T', 'G'), 'G': ('A', 'C', 'T', 'A'), 'C': (
    'T', 'A', 'G', 'T'), 'T': ('C', 'G', 'A', 'C')}


def read_complement(read):
    l = []
    for char in read:
        l.append(COMPLEMENT[char.upper()])
    return ''.join(l)


def add_read_err(read, qphred, asc):
    ''' read,phred is list'''
    s = []
    for char, qua in zip(read, qphred):
        error_rate = 10**((asc-ord(qua))/10)
        if QPH==qphs.soleax.value:
            error_rate=error_rate/(1+error_rate)
        threshold = random.random()
        if threshold < error_rate:
            pich = random.randint(0, 3)
            char = SUBSTITUTIONS[char.upper()][pich]
        s.append(char)
    return ''.join(s)


def add_segment_err(segment):
    ''' segment,phred is list'''
    segment = list(segment)
    l=[]
    for char in segment:
        if char in ATCG:
            l.append(char.upper())
        elif char.upper() in DEGENERATE:
            l.append(random.choice(DEGENERATE[char.upper()]))
        else:
            print("error :can't indentify", char)
            l.append(random.choice(ATCG[:4]))
            #raise Exception("error :can't indentify", char)
    segment=l
    size = len(segment)
    error_e = ERROR_E
    error_d = ERROR_D
    error_e = value2unnegtive(random.normalvariate(error_e, error_d))
    err_num = value2unnegtive(round((error_e*size)))
    err_num_list = random.sample(range(size), err_num)  # include zero

    a = SUBSTITUTION
    b = a+DELETION
    c = a+b+INSERTION

    for n in err_num_list:
        randomnum = random.random()
        if randomnum <= b:  # substition
            pich = random.randint(0, 3)
            try:
                segment[n] = SUBSTITUTIONS[segment[n].upper()][pich]
            except Exception as e:
                print(n,pich,len(segment))
                raise e
        elif randomnum <= c:  # deletion
            segment[n]=''
        else:
            a=segment[n].upper()
            s = ATCG[:4]+[a]
            char = random.choice(s)
            segment[n]=a+s
    return ''.join(segment)


def wgsout(content, readlen, quality, reffile):
    frequencies, asc = Quality.get_qph(quality)
    if pairs.PE.value == PAIR:
        r1 = open(R1, 'w', newline='\n')
        r2 = open(R2, 'w', newline='\n')
    elif PAIR == pairs.SE.value:
        r0 = open(R0, 'w', newline='\n')
    infos = Fasta.fasta_file_info(reffile)
    column = infos['column']
    step = infos['step']
    l = len(infos)-2
    for x in range(l):
        print('write chromosome',x+1,'fasta...',end='')
        chromosome, pos = Fasta.analyse_infos(infos, x+1)
        turns = round(DEPTH*chromosome.length/readlen *content/len(REFERENCES)/PAIR)
        with open(reffile, 'r', newline='\n') as filed:
            for y in range(turns):
                if y%10000==0:
                    qphreds = Quality.get_qphred_reads(frequencies, 100)
                    print(y//10000)
                segment = wgs_segment(
                    filed, readlen, chromosome, pos, column, step)
                head = chromosome.id+'.'+str(y+1)
                if ERROR:
                    segment = add_segment_err(segment)
                if len(segment) < readlen:
                    print("segment length is too short -> ignore it")
                    continue
                if pairs.PE.value == PAIR:
                    write_fastq(segment, head, qphreds, readlen, asc, r1, r2)
                else:
                    write_fastq(segment, head, qphreds, readlen,asc, r0)
    if pairs.PE.value == PAIR:
        r1.close()
        r2.close()
    elif PAIR == pairs.SE.value:
        r0.close()
    print(' down.')


def wgs_segment(filed, readlen, chromosome, pos, column, step):
    stop = chromosome.end
    tail = max(SEGMENT_E, readlen)
    begin = random.randint(1, stop-tail)
    end = round(begin+random.normalvariate(SEGMENT_E, SEGMENT_D))
    if end > stop:
        end = stop
    length = end-begin+1
    if length < readlen:
        begin = begin-readlen+length
        length = readlen
    segment = ''.join(get_words(filed, begin, length,
                                MEMORY, column, step, pos))
    return segment


def wesout(content, readlen, quality, seqfile):

    frequencies, asc = Quality.get_qph(quality)
    if pairs.PE.value == PAIR:
        r1 = open(R1, 'w', newline='\n')
        r2 = open(R2, 'w', newline='\n')
    elif PAIR == pairs.SE.value:
        r0 = open(R0, 'w', newline='\n')
    with open(seqfile, 'r', newline='\n') as filed:
        sequence = Fasta.iterator_fasta(filed)
        for x in sequence:
            length = x.seq_length
            qphreds = Quality.get_qphred_reads(frequencies, 10+length//100)
            turns = round(DEPTH*x.mid_length/readlen *
                          content/len(REFERENCES)/PAIR)
            if length < readlen:
                print('\tseq length too short -> not ouput ')
                continue
            if x.mid_length < CHIP_LEN:
                print('\taimed seq length too short -> not output')
                continue
            print('\twriting exon : %s' % x.id, end='\r')
            for y in range(turns):
                segment = x.wes_segment(
                    readlen, CHIP_LEN, SEGMENT_E, SEGMENT_D)
                head = x.id+'.'+str(y+1)
                if ERROR:
                    segment = add_segment_err(segment)
                if pairs.PE.value == PAIR:
                    write_fastq(segment, head, qphreds, readlen, asc, r1, r2)
                else:
                    write_fastq(segment, head, qphreds, asc, r0)
    if pairs.PE.value == PAIR:
        r1.close()
        r2.close()
    elif PAIR == pairs.SE.value:
        r0.close()
    print(' down.')


def write_fastq(segment, head, qphreds, readlen, asc, *R):
    if pairs.PE.value == PAIR:
        read1 = segment[:readlen]
        read2 = segment[-1:-1-readlen:-1]
        read2 = read_complement(read2)
        qphred1 = random.choice(qphreds)
        qphred2 = random.choice(qphreds)
        if MISMATCH:
            read1 = add_read_err(read1, qphred1, asc)
            read2 = add_read_err(read2, qphred2, asc)
        R[0].write('@'+head+'/1\n+\n'+read1+'\n'+qphred1+'\n')
        R[1].write('@'+head+'/2\n+\n'+read1+'\n'+qphred1+'\n')
    else:
        read = segment[:readlen]
        qphred = random.choice(qphreds)
        if MISMATCH:
            read = add_read_err(read, qphred, asc)
        R[0].write('@'+head+'\n+\n'+read+'\n'+qphred+'\n')


def write_ref_muta(ref, inimuta, refmuta, nid):
    print("write mutation's reference from file :", ref, inimuta)
    infos = Fasta.fasta_file_info(ref)
    column = infos['column']
    step = infos['step']
    chrs = Fasta.chrs_info(infos)
    ninfos = {}
    with open(ref, 'r', newline='\n') as fd1:
        with open(inimuta, 'r', newline='\n')as fd2:
            with open(refmuta, 'w', newline='\n')as wd2:
                for chrr in chrs:
                    chromosome, pos = Fasta.analyse_infos(infos, chrr)
                    mutas = Fasta.iterator_fasta(fd2, chromosome.id)
                    print(chromosome.get_fasta_head(), end='\r')
                    chromosome.id += '.%s' % nid
                    wd2.write('%s\n' % chromosome.get_fasta_head())
                    npos = wd2.tell()
                    ninfos[chrr] = {'pos': npos,
                                    'chromosome': chromosome.get_fasta_head()}
                    begin, cur = 1, 0
                    for x in mutas:
                        seqs = get_words(fd1, begin, x.begin-begin,
                                         MEMORY, column, step, pos)
                        cur = write_big_words(wd2, seqs, COLUMN, cur)
                        cur = write_small_word(wd2, x.seq, COLUMN, cur)
                        begin = x.end+1
                    seqs = get_words(fd1, begin, chromosome.end-begin+1,
                                     MEMORY, column, step, pos=pos)
                    cur = write_big_words(wd2, seqs, COLUMN, cur)
                    wd2.write('\n')
    ninfos['column'] = COLUMN
    ninfos['step'] = 1
    with open(FASTA_INFO, 'r') as filed:
        try:
            fasta_infos = json.load(filed)
        except:
            fasta_infos = {}
    with open(FASTA_INFO, 'w') as writed:
        fasta_infos[refmuta] = ninfos
        json.dump(fasta_infos, writed)
    print('down. outfile : ', refmuta)


def readout(refs, regs, qph, inimutation, mut):
    frequencies, acs = Quality.get_qph(qph)
    readlen = len(frequencies)
    if SEGMENT_E<=readlen:
        print('segment is not much longger than readlen')
        return
    frequencies = {}
    for ref, reg in zip(refs, regs):
        regfile = reg+'.reg'+str(mut[0])
        tempref = 'mutation.ref'
        if not os.path.exists(tempref):
            write_ref_muta(ref, inimutation, tempref,mut[0])
        if MODE == modes.WES.value:
            seqfile = 'mutation.exome'
            if not os.path.exists(seqfile):
                Fasta.ini_exome(tempref, regfile, seqfile)
            os.remove(tempref)
            wesout(mut[1], readlen, qph, seqfile)
            os.remove(seqfile)
        elif MODE == modes.WGS.value:
            wgsout(mut[1], readlen, qph, tempref)
            os.remove(tempref)