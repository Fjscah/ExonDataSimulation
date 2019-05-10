import json
import random
import re
import sys
import time
import os
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
modes = Enum('mode', ('WES', 'WGS'))
pairs = Enum('pair', ('SE', 'PE', ))
mutation_ways = (Enum('mutation_way', ('table', 'formula', 'auto')))
qphs = Enum('qph', ('sanger', 'soleax', ))

def read_complement(read):
    l = []
    for char in read:
        l.append(COMPLEMENT[char.upper()])
    return ''.join(l)


def add_read_err(read, qphred, asc,QPH):
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


def add_segment_err(segment,*,error_e,error_d,substitution,deletion,insertion,**kw):
    ''' segment,phred is list'''
    size = len(segment)
    error_e = value2unnegtive(random.normalvariate(error_e, error_d))
    err_num = value2unnegtive(round((error_e*size)))
    err_num_list = random.sample(range(size), err_num)  # include zero

    a = substitution
    b = a+deletion
    c = a+b+insertion

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
            segment[n]=s+char
    return segment

def check_segment(segment):
    l=[]
    segment = list(segment)
    for char in segment:
        if char in ATCG:
            l.append(char.upper())
        elif char.upper() in DEGENERATE:
            l.append(random.choice(DEGENERATE[char.upper()]))
        else:
            print("error :can't indentify", char)
            l.append(random.choice(ATCG[:4]))
            #raise Exception("error :can't indentify", char)
    return l



def wgs_segment(filed, chromosome, pos, column, step,segment_e,segment_d,MEMORY):
    stop = chromosome.end
    begin = random.randint(1, stop-segment_e)
    if segment_d:
        end = round(begin+random.normalvariate(segment_e, segment_d))
    else:
        end = round(begin+segment_e)
    if end > stop:
        end = stop
    length = end-begin+1
    segment = ''.join(get_words(filed, begin, length,
                                MEMORY, column, step, pos))       
    return segment



def motify_segment(segment,extra,effect_len,*,ERROR,**kw):
    l=len(segment)
    segment=check_segment(segment)#segment is list
    extra=check_segment(extra)
    if ERROR:
        segment = add_segment_err(segment,**kw)#list
        #print('ERROR')
    if len(segment)<l:
        segment+=extra[:l-len(segment)]
    if len(segment)<effect_len:
        print('segment is too short -> ignore it',effect_len,len(segment))
        return ''
    else:
        return ''.join(segment)
def wgsout(content,  quality, reffile,MEMORY,*,effect_len,chip_len,R0,R1,R2,PAIR,SEG,SEGMENT_E,SEGMENT_D,DEPTH,REFERENCES,**kw):
    frequencies,sums,readlens, asc = Quality.get_qph(quality)
    readlen=Quality.get_avg_readlens(readlens)
    if pairs.PE.value == PAIR:
        r1 = open(R1, 'w', newline='\n')
        r2 = open(R2, 'w', newline='\n')
    elif PAIR == pairs.SE.value:
        r0 = open(R0, 'w', newline='\n')
    if os.path.exists(SEG):
        segs,summ=segfile(SEG)[:2]
    segment_e=SEGMENT_E
    segment_d=SEGMENT_D
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
                    qphreds = Quality.get_qphred_reads(frequencies,sums,readlens,100)
                    print(y//10000)
                if os.path.exists(SEG):
                    segment_e=random_weight_choice(segs,summ)
                    segment_d=0
                segment = wgs_segment(filed, chromosome, pos, column, step,segment_e,segment_d,MEMORY)
                segment=motify_segment(segment,'',effect_len,**kw)#segment is list
                if not segment:
                    continue 
                head = chromosome.id+'.'+str(y+1)
                if pairs.PE.value == PAIR:
                    write_fastq(segment, head, qphreds, asc, PAIR,r1, r2,**kw)
                else:
                    write_fastq(segment, head, qphreds, asc,PAIR, r0,**kw)
    if pairs.PE.value == PAIR:
        r1.close()
        r2.close()
    elif PAIR == pairs.SE.value:
        r0.close()
    print(' down.')
def wesout(content, quality,seqfile,*,effect_len,chip_len,R0,R1,R2,PAIR,SEG,SEGMENT_E,SEGMENT_D,DEPTH,REFERENCES,**kw):
    ERROR=kw['ERROR']
    error_e=kw['error_e']
    frequencies,sums,readlens, asc = Quality.get_qph(quality)
    readlen=Quality.get_avg_readlens(readlens)
    if pairs.PE.value == PAIR:
        r1 = open(R1, 'w', newline='\n')
        r2 = open(R2, 'w', newline='\n')
    elif PAIR == pairs.SE.value:
        r0 = open(R0, 'w', newline='\n')
    if os.path.exists(SEG):
        print('get segment length from file :',SEG)
        segs,summ=segfile(SEG)[:2]
    sequence = Fasta.iterator_fasta(seqfile)
    segment_e=SEGMENT_E
    segment_d=SEGMENT_D
    echr=0
    for x in sequence:
        length = x.seq_length
        qphreds = Quality.get_qphred_reads(frequencies,sums,readlens, 10+length//100)
        turns = round(DEPTH*x.mid_length/readlen *
                        content/len(REFERENCES)/PAIR)
        if length < effect_len:
            print('\tseq length %s too short -> not ouput '%length)
            print(x.get_fasta_head(),'\n',x.seq)
            continue
        if echr!=x.chr:
            print('\twriting exon : %s' % x.id)
            echr=x.chr
        for y in range(turns):
            if os.path.exists(SEG):
                segment_e=random_weight_choice(segs,summ)
                segment_d=0
            segment,extra = x.wes_segment(
                    effect_len, chip_len, segment_e, segment_d,ERROR,error_e)
            segment=motify_segment(segment,extra,effect_len,**kw)
            if not segment:
                continue
            head = x.id+'.'+str(y+1)#segment is str
            if pairs.PE.value == PAIR:
                write_fastq(segment, head, qphreds, asc, PAIR,r1, r2,**kw)
            else:
                write_fastq(segment, head, qphreds, asc,PAIR,r0,**kw)
    if pairs.PE.value == PAIR:
        r1.close()
        r2.close()
    elif PAIR == pairs.SE.value:
        r0.close()
    print(' down.')


def write_fastq(segment, head, qphreds,  asc,PAIR, *R,MISMATCH,QPH,**kw):
    if pairs.PE.value == PAIR:
        qphred1 = random.choice(qphreds)
        qphred2 = random.choice(qphreds)
        l1=len(qphred1)
        l2=len(qphred2)
        read1 = segment[:l1]
        read2 = segment[-1:-1-l2:-1]
        read2 = read_complement(read2)
        s=l1-len(read1)
        if s:
            e=random.randint(0,s)
            qphred1=qphred1[e:e+len(read1)]
        s=l2-len(read2)
        if s:
            e=random.randint(0,s)
            qphred2=qphred2[e:e+len(read2)]
        if MISMATCH:
            read1 = add_read_err(read1, qphred1, asc,QPH)
            read2 = add_read_err(read2, qphred2, asc,QPH)
        if random.random()>0.5:
            R[0].write('@'+head+'\t1:Y:0:\n'+read1+'\n+\n'+qphred1+'\n')
            R[1].write('@'+head+'\t2:Y:0:\n'+read1+'\n+\n'+qphred1+'\n')
        else:
            R[0].write('@'+head+'\t1:Y:0:\n'+read2+'\n+\n'+qphred2+'\n')
            R[1].write('@'+head+'\t2:Y:0:\n'+read1+'\n+\n'+qphred1+'\n')
    else:
        qphred = random.choice(qphreds)
        l=len(qphred)
        read = segment[:l]
        s=l-len(read)
        if s:
            e=random.randint(0,s)
            qphred=qphred[e:e+len(read)]
        if MISMATCH:
            read = add_read_err(read, qphred, asc,QPH)
        R[0].write('@'+head+'\n'+read+'\n+\n'+qphred+'\n')


def write_ref_muta(ref, inimuta, refmuta, nid,COLUMN,MEMORY,FASTA_INFO):
    print("write mutation's reference from file :", ref, inimuta)
    infos = Fasta.fasta_file_info(ref,FASTA_INFO)
    column = infos['column']
    step = infos['step']
    chrs = Fasta.chrs_info(infos)
    ninfos = {}
    with open(ref, 'r', newline='\n') as fd1:
        with open(refmuta, 'w', newline='\n')as wd2:
            for chrr in chrs:
                chromosome, pos = Fasta.analyse_infos(infos, chrr)
                mutas = Fasta.iterator_fasta(inimuta, chromosome.id)
                length=0
                for x in mutas:
                    length=length-x.length+x.seq_length
                nchromosome=Sequence(chrr,chromosome.begin,chromosome.end+length,chromosome.id)
                print(nchromosome.get_fasta_head())
                chromosome.id += '.%s' % nid
                wd2.write('%s\n' % nchromosome.get_fasta_head())
                npos = wd2.tell()
                ninfos[chrr] = {'pos': npos,
                                'chromosome': chromosome.get_fasta_head()}
                begin, cur = 1, 0
                mutas = Fasta.iterator_fasta(inimuta, chromosome.id)
                for x in mutas:
                    seqs = get_words(fd1, begin, x.begin-begin,
                                        MEMORY, column, step, pos)
                    print('ref(',begin, x.begin-1,end=')-')
                    cur = write_big_words(wd2, seqs, COLUMN, cur)
                    print('mut(',x.begin,x.end,end=')-')
                    cur = write_small_word(wd2, x.seq, COLUMN, cur)
                    begin = x.end+1
                print('ref(',begin, chromosome.end,')')
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

'''
ss=[]
for x in mutas:
    ss.append((begin,x.begin))
    ss.append((x.begin+1,x.end))
    begin=x.end+1
ss.append((begin,chromosome.end))
print(ss)
'''

def readout(refs, regs, qph, inimutation, mut,label,*,MODE,COLUMN,MEMORY,FASTA_INFO,FLI_N,INNER_N,**kw):
    effect_len=kw['effect_len']
    flank_len=kw['flank_len']
    for ref, reg in zip(refs, regs):
        # haplotype + mutu - > mute reg
        regfile =os.path.splitext(reg)[0]+str(mut[0])+'.bed'
        # haplotype + mute -> mut ref
        tempref = os.path.splitext(ref)[0]+str(mut[0])+'.fna'
        tt=tempref+'.temp'
        if not os.path.exists(tempref):
            write_ref_muta(ref, inimutation, tt,mut[0],COLUMN,MEMORY,FASTA_INFO)
            os.rename(tt,tempref)
        if MODE == modes.WES.value:
            # haplotype + mute -> mute exome
            seqfile = os.path.splitext(reg)[0]+'exome'+str(mut[0])+'.fna'
            tt=seqfile+'.temp'
            if not os.path.exists(seqfile):
                Fasta.ini_exome(tempref, regfile, tt,effect_len,flank_len,FLI_N,INNER_N,COLUMN,MEMORY,FASTA_INFO)
                os.rename(tt,seqfile)
            wesout(mut[1],  qph, seqfile,**kw)
            Fasta.remove_fastainfo(FASTA_INFO,tempref)
            os.remove(tempref)
            os.remove(seqfile)
        elif MODE == modes.WGS.value:
            wgsout(mut[1], qph, tempref,MEMORY,**kw)
            Fasta.remove_fastainfo(FASTA_INFO,tempref)
            os.remove(tempref)