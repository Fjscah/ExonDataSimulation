import json
import random
import re
import sys
import time
from profile import *

from filefunc import *
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
        error_rate = 10**(asc-ord(qua)/-10)
        threshold = random.random()
        if threshold < error_rate:
            pich = random.randint(0, 3)
            char = SUBSTITUTIONS[char.upper()][pich]
        s.append(char)
    return ''.join(s)


def wgsout(content, readlen, quality, reffile):
    frequencies, asc = Quality.get_qph(quality)
    if pairs.PE == PAIR:
        r1 = open(R1, 'w', newline='\n')
        r2 = open(R2, 'w', newline='\n')
    elif PAIR == pairs.SE:
        r0 = open(R0, 'w', newline='\n')
    infos = Fasta.fasta_file_info(reffile)
    column = infos['column']
    step = infos['step']
    l = len(infos)-2
    for x in range(l):
        chromosome, pos = Fasta.analyse_infos(infos, x+1)
        qphreds = Quality.get_qphred_reads(
            frequencies, 10+chromosome.length//100)
        turns = round(DEPTH*chromosome.length/readlen *
                      content/len(REFERENCES/PAIR.value))
        with open(reffile, 'r', newline='\n') as filed:
            for y in range(turns):
                segment = wgs_segment(
                    filed, readlen, chromosome, pos, column, step)
                head = y.id+'.'+str(y+1)
                segment = add_segment_err(segment)
                if pairs.PE == PAIR:
                    write_fastq(segment, head, qphreds, readlen, asc, r1, r2)
                else:
                    write_fastq(segment, head, qphreds, asc, r0)
    if pairs.PE == PAIR:
        r1.close()
        r2.close()
    elif PAIR == pairs.SE:
        r0.close()
    print(' down.')


def wgs_segment(filed, readlen, chromosome, pos, column, step):
    stop = chromosome.end
    tail = max(SEGMENT_E, readlen)
    begin = random.randint(1, stop-tail)
    end = begin+random.normalvariate(SEGMENT_E, SEGMENT_D)
    if end > stop:
        end = stop
    length = end-begin+1
    if length < readlen:
        begin = begin-readlen+length
        length = readlen
    segment = ''.join(get_words(filed, begin, length, MEMORY, column, step,pos))
    return segment


def wesout(content, readlen, quality, seqfile):
    sequence = Fasta.iterator_fasta(seqfile)
    frequencies, asc = Quality.get_qph(quality)
    if pairs.PE == PAIR:
        r1 = open(R1, 'w', newline='\n')
        r2 = open(R2, 'w', newline='\n')
    elif PAIR == pairs.SE:
        r0 = open(R0, 'w', newline='\n')
    for x in sequence:
        length = x.seq_length
        qphreds = Quality.get_qphred_reads(frequencies, 10+length//100)
        turns = round(DEPTH*x.mid_length/readlen *
                      content/len(REFERENCES)/PAIR.value)
        if length < readlen:
            print('\tnot output , error:seq length too short \n')
            return
        if x.mid_length < CHIP_LEN:
            print('\tnot output , errpr:aimed seq length too short \n')
            return
        print('\twriting exon : %s' % x.id, end='\r')
        for y in range(turns):
            segment = y.wes_segment(readlen, CHIP_LEN, SEGMENT_E, SEGMENT_D)
            head = y.id+'.'+str(y+1)
            segment = add_segment_err(segment)
            if pairs.PE == PAIR:
                write_fastq(segment, head, qphreds, readlen, asc, r1, r2)
            else:
                write_fastq(segment, head, qphreds, asc, r0)
    if pairs.PE == PAIR:
        r1.close()
        r2.close()
    elif PAIR == pairs.SE:
        r0.close()
    print(' down.')


def add_segment_err(segment):
    ''' segment,phred is list'''
    size = len(segment)
    segment = list(segment)
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
        if random <= b:  # substition
            pich = random.randint(0, 3)
            segment[n] = SUBSTITUTIONS[segment[n].upper()][pich]
        elif randomnum <= c:  # deletion
            segment.pop(n)
        else:
            s = ATCG[:4]+[segment[n].upper()]
            char = random.choice(s)
            segment.insert(n, char)
    return ''.join(segment)


def write_fastq(segment, head, qphreds, readlen, asc, *R):
    if pairs.PE == PAIR:
        read1 = segment[:R]
        qphred1 = random.choice(qphreds)
        add_read_err(read1, qphred1, asc)
        read2 = segment[-1:-1-readlen:-1]
        read2 = read_complement(read2)
        qphred2 = random.choice(qphreds)
        add_read_err(read2, qphred2, asc)
        R[0].write('@'+head+'/1\n+\n'+read1+'\n'+qphred1+'\n')
        R[1].write('@'+head+'/2\n+\n'+read1+'\n'+qphred1+'\n')
    else:
        read = segment[:R]
        qphred = random.choice(qphreds)
        add_read_err(read, qphred, asc)
        R[0].write('@'+head+'\n+\n'+read1+'\n'+qphred1+'\n')


def write_ref_muta(ref, inimuta, refmuta):
    infos = Fasta.fasta_file_info(ref)
    column = infos['column']
    step = infos['step']
    with open(ref, 'r', newline='\n') as fd2:
        line = get_line_text(fd2, '>', 0, 're')
        while(line):
            chromosome = Sequence.self_fasta_head(line)
            fd2.write('%s\n' % chromosome.get_fasta_head())
            mutas = Fasta.iterator_fasta(inimuta, chromosome.id)
            begin, cur = 1, 0
            for x in mutas:
                pos = infos[chromosome.chr]['pos']
                seqs = get_words(fd2, begin, x.begin-begin,
                                 MEMORY, column, step, pos)
                cur = write_big_words(fd2, seqs, COLUMN, cur)
                cur = write_small_word(fd2, x.seq, COLUMN, cur)


def readout(refs, regs, qph, inimutation, poly):
    frequencies, acs = Quality.get_qph(qph)
    readlen = len(frequencies)
    frequencies = {}
    for ref, reg in zip(refs, regs):
        regfile = reg+'.reg'+str(poly[0])
        tempref = 'mutation.ref'
        write_ref_muta(ref, inimutation, tempref)
        if MODE == modes.WES:
            seqfile = 'mutation.exome'
            Fasta.ini_exome(tempref, regfile, seqfile)
            wesout(poly[0], readlen, qph, seqfile)
        elif MODE == modes.WGS:
            wgsout(poly[0], readlen, qph, regfile)
