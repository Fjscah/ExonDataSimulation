import json
import random
import re
import shutil
import sys
import time
from random import choice

import mutation
from exon import Exon
from setting import (ACCURACY_RATE, CHIP_LEN, DEEPTH, INSERT, READ, ROW_STEP,
                     SUBSTITUTION,PHRED)
from view import chr_exon_num, show_chr_exon_num

ATCG=('A','T','C','G','a','t','c','g')
COMPLEMENT={'A':'T','G':'C','C':'G','T':'A'}


'''
def random_weight_choice(lists):
    summ=numpy.sum(lists,axis=0)[1]
    t = random.randint(1, summ)
    for i, val in lists:
        t -= val
        if t <= 0:
            return chr(i+33)
    if t>0:
        input('kkkkkkkkkkkkkkkkk')
'''


def random_qphred(rang,frequencies,plus=33):
    s=''
    for x in range(1,len(rang)+1):
        weights=frequencies['pos%d_frequencies'%x]
        num=random.choices(rang,weights=weights)
        char=chr(num[0]+plus)
        s=s+char
    return s

def add_read_err(read,phred,accuracy=ACCURACY_RATE):
    size=len(read)
    err_num=round((1-accuracy)*size)
    err_num_list=random.sample(range(size),err_num)
    for n in err_num_list:
        pich=random.randint(0,2)
        read[n]=SUBSTITUTION[read[n].upper()][pich]
        phred_value=ord(phred[n])-33
        if phred_value>10:
            phred[n]=chr(ord(phred[n])-10)
        else:
            phred[n]=chr(33)
def complementation(read):
    l=[]
    for char in read:
        l.append(COMPLEMENT[char.upper()])
    return l
def PEread(seq,length,frequencies):
    ''' get one inser's 2 reads'''
    p1=INSERT+length-CHIP_LEN-1
    cbp=random.randint(INSERT, p1)
    sbp=random.randint(-INSERT+CHIP_LEN,0)
    lbp=cbp+sbp
    l=[]
    for n in range(INSERT):
        char =seq[n+lbp]
        if char in ATCG:
            l.append(char.upper())
        else :
            l.append(choice(ATCG))
    read1=l[0:READ]#list
    read2=complementation(l[-1:-1-READ:-1])#list
    phred1=list(choice(frequencies))
    phred2=list(choice(frequencies))
    add_read_err(read1,phred1)
    add_read_err(read2,phred2)
    return ''.join(read1) , ''.join(read2),''.join(phred1),''.join(phred2)

def get_exon(file):
    ''' from every row in flie get exon info and yield exon object'''
    line=file.readline()
    while(line):
        m = re.match(r'^chr([\d]*)\t(\d*)\t(\d*)\t([\w\.\-]*)\t\d*\t(\d*)',line)
        if m:
            # if m is exon title ,then store this exon's sequence
            seq=''
            line=file.readline()
            n = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.\-]*)',line)
            while(not n):
                seq=seq+line.strip()
                line=file.readline()
                if not line:
                    break
                n = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.\-]*)',line)
            yield Exon(m.group(1),int(m.group(2)),int(m.group(3)),m.group(4),seq,insert=int(m.group(5)))
        else:
            line=file.readline()

def fastq_exon(w1,w2,exon,frequencies,mutation_types):
    ''' one exon's all fastq reault be writed to w'''
    print('writing chr %s - exon : %s  - %s' %(exon.chr,exon.exon_id,exon.begin) ,end='\r')
    if 'N' in exon.seq:
        print('\n'+''.join(exon.seq))
    exon_id=exon.exon_id
    seq,length=list(exon.seq),exon.length
    for x in mutation_types:
        copy_num=mutation_types[x].haplots[0].get_normal_num(exon_id)+mutation_types[x].haplots[1].get_normal_num(exon_id)
        turns=round(length*mutation_types[x].percent*copy_num*DEEPTH/2/INSERT/100)
        for y in range (turns):
            read1,read2,phred1,phred2=PEread(seq,length,frequencies)
            w1.write('@%s:%s:%d/1\n'% (exon_id,x,y))
            w1.write(read1+'\n')
            w1.write('+\n')
            w1.write(phred1+'\n')#random_qphred(frequencies)
            w2.write('@%s:%s:%d/2\n'% (exon_id,x,y))
            w2.write(read2+'\n')
            w2.write('+\n')
            w2.write(phred2+'\n')
    for x in mutation_types:
        for z in range(2):
            if exon_id in mutation_types[x].haplots[z].mutations:
                for seq,length,copy_num,onemutation in exon.get_new_seq((mutation_types[x].haplots[z].mutations[exon_id].snp_mutations),CHIP_LEN):
                    turns=round(length*mutation_types[x].percent*copy_num*DEEPTH/2/INSERT/100)
                    seq=list(seq)
                    for y in range(turns):
                        read1,read2,phred1,phred2=PEread(seq,length,frequencies)
                        w1.write('@%s:%s:%d:%d:%s/1\n'% (exon_id,x,z,y,onemutation))
                        w1.write(read1+'\n')
                        w1.write('+\n')
                        w1.write(phred1+'\n')#random_qphred(frequencies)
                        w2.write('@%s:%s:%d:%d:%s/2\n'% (exon_id,x,z,y,onemutation))
                        w2.write(read2+'\n')
                        w2.write('+\n')
                        w2.write(phred2+'\n')
            '''
            reads=reads+'@%s/1\n'% gene_id+read1+'\n'+'+\n'+random_qphred(frequencies)+'\n'
            reads=reads+'@%s/2\n'% gene_id+read2+'\n'+'+\n'+random_qphred(frequencies)+'\n'
            '''


if __name__ == '__main__':


    #set and print gene mutation type
    # mutafile='mutations_setting.txt'
    mutation_types=mutation.set_mutation()
    if mutation=='exit':
        exit
    input('press Enter to continue...')

    # initial qphred frequencies
    phredfile=input('please input qphred file : ')
    print('initial qphred...')
    with open(phredfile,"r") as f:
        frequencies=json.load(f)
        #frequencies['pos%d_frequencies'%x]=list(enumerate(frequencies['pos%d_frequencies'%x]))
        #frequencies['pos%d_frequencies'%x].sort(key=lambda y:y[1],reverse=True) 
    print('down...')


    # creat qphread 100
    rang=list(range(len(frequencies['pos%d_frequencies'%1])))
    phred_reads=[]
    for x in range(100):
        phred_reads.append(random_qphred(rang,frequencies,PHRED))
    frequencies=phred_reads
    phred_reads=None




    # set filein name
    filein=input('please input exonlist file : ')#'REF_exonlist(insert).txt'

    # get every chr's exon number
    #print("getting exon number on every chromosome...")
    #total_exon=countrow(filein,r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)')
    #print("exon count: ",total_exon)
    '''
    exon_nums=chr_exon_num(filein)
    show_chr_exon_num(exon_nums)
    '''
    #print('down')


    # according to mutaton output reads
    # output normal reads(not snp mutation)
    time_start=time.time()
    print('output normal reads(cnv+deletion)...')
    t=time.strftime('%Y%m%d_%H_%M_%S',time.localtime(time.time()))
    with open("R1_fastq%s.fastq" % t,'w') as w1:
        with open("R2_fastq%s.fastq" % t,'w') as w2:
            # one exon generate to mutation_types
            with open(filein,'r')as r:
                for exon in get_exon(r):
                    fastq_exon(w1,w2,exon,frequencies,mutation_types)
    print('\n***down***')
    time_end=time.time()
    t=time_end-time_start
    print('totally cost: %dh : %dm : %ds'%(t//3600,(t%3600)//60,t%60))   
