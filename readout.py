import json
import random
import re
import shutil
import sys
import time
from random import choice

import numpy

import mutationtype
from exon import Exon
from setting import (ACCURACY_RATE, CHIP_LEN, DEEPTH, INSERT, READ, ROW_COLUMN,
                     ROW_STEP)

ATCG=('A','T','C','G','a','t','c','g')
DICT={'A':('G','C','T','G'),'G':('A','C','T','A'),'C':('T','A','G','T'),'T':('C','G','A','C')} 
COMPLEMENT={'A':'T','G':'C','C':'G','T':'A'}
"""       
def bp2pos(bp):
    ''' according to bp length ,get it's position in file
    only suitable for pos at the beginning of row'''
    raw=bp//ROW_COLUMN
    return raw*(ROW_COLUMN+ROW_STEP)+bp %ROW_COLUMN-1
"""

def get_normal_exon(file):
    ''' from every row in flie get exon info and yield exon object'''
    line=file.readline()
    c=1
    while(line):
        m = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
        if m:
            # if m is exon title ,then store this exon's sequence
            seq=''
            line=file.readline()
            n = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
            while(not n):
                seq=seq+line.strip()
                line=file.readline()
                if not line:
                    break
                n = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
            yield Exon(m.group(1),int(m.group(2)),int(m.group(3)),m.group(4),seq),c
            c+=1
        else:
            line=file.readline()
def get_snp_exon(file,snp_mutations):
    count=0
    line=file.readline()
    for x in snp_mutations:
        num=x[0]
        if count==num:
            new_seq=snp2exon(seq,x[1])
            yield Exon(m.group(1),int(m.group(2)),int(m.group(3)),m.group(4),new_seq),x
            continue
        while(line):
            m = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
            if m:
                count+=1
                if count==num:
                    seq=''
                    line=file.readline()
                    n = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
                    while(not n):
                        seq=seq+line.strip()
                        line=file.readline()
                        n = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
                    new_seq=snp2exon(seq,x[1])
                    yield Exon(m.group(1),int(m.group(2)),int(m.group(3)),m.group(4),new_seq),x
                    break
                else:
                    line=file.readline()
                    n = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
                    while(not n):
                        line=file.readline()
                        n = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
            else:
                line=file.readline()
def snp2exon(seq,snps):
    seq=list(seq)
    length=len(seq)
    snp_pos=[]
    for x in snps:
        pos=round(x*length)
        while(pos in snp_pos):
            pos=pos+1
        snp_pos.append(pos)
    for x in snp_pos:
        seq[x-1]=random.choice(DICT[seq[x-1].upper()])
    return seq


def random_weight_choice(lists):
    summ=numpy.sum(lists,axis=0)[1]
    t = random.randint(1, summ)
    for i, val in lists:
        t -= val
        if t <= 0:
            return chr(i+33)
    if t>0:
        input('kkkkkkkkkkkkkkkkk')
def random_qphred(frequencies):
    s=''
    for x in range(1,READ+1):
        char=random_weight_choice(frequencies['pos%d_frequencies'%x])
        s=s+char
    return s

def add_read_err(read,phred,accuracy=ACCURACY_RATE):
    size=len(read)
    err_num=round((1-accuracy)*size)
    err_num_list=random.sample(range(size),err_num)
    for n in err_num_list:
        pich=random.randint(0,2)
        read[n]=DICT[read[n].upper()][pich]
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
def PEread(exon,frequencies):
    ''' get one inser's 2 reads'''
    p1=INSERT+exon.length()-CHIP_LEN
    cbp=random.randint(INSERT, p1)
    sbp=random.randint(-INSERT+CHIP_LEN,0)
    lbp=cbp+sbp
    l=[]
    for n in range(INSERT):
        char =exon.seq[n+lbp]
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
def fastq_exon(w1,w2,exon,frequencies,mutation_types,n,boolsnp,**snps):
    ''' one exon's all fastq reault be writed to w'''
    print('writing chr %s - exon : %s  - %s' % 
        (exon.chr,exon.gene_id,exon.begin) ,end='\r')
    if 'N' in exon.seq:
        print('\n'+''.join(exon.seq))
    gene_id=exon.gene_id
    for x in range(len(mutation_types)):
        if boolsnp==False:
            copy_num=mutation_types[x].homologs[0].get_normal_copy_num(n)+mutation_types[x].homologs[1].get_normal_copy_num(n)
        elif boolsnp==True:
            if 'muta' in snps:
                if snps['muta']!=x:
                    continue
            if 'homolog' in snps:
                copy_num=mutation_types[x].homologs[snps['homolog']].get_snp_copy_num(n)
            if copy_num==0:
                continue
        turns=round(exon.length()*mutation_types[x].percent*copy_num*DEEPTH/2/INSERT/100)
        for y in range (turns):
            read1,read2,phred1,phred2=PEread(exon,frequencies)
            w1.write('@%s:%d/1\n'% (gene_id,y))
            w1.write(read1+'\n')
            w1.write('+\n')
            w1.write(phred1+'\n')#random_qphred(frequencies)
            w2.write('@%s:%d/2\n'% (gene_id,y))
            w2.write(read2+'\n')
            w2.write('+\n')
            w2.write(phred2+'\n')
            '''
            reads=reads+'@%s/1\n'% gene_id+read1+'\n'+'+\n'+random_qphred(frequencies)+'\n'
            reads=reads+'@%s/2\n'% gene_id+read2+'\n'+'+\n'+random_qphred(frequencies)+'\n'
            '''

def countrow(file,text):
    # return certain row's number
    with open(file,'r')as r:
        count=0
        for line in r:
            m = re.match(r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)',line)
            if m:
                count+=1
        return count

# initial qphred
print('initial qphred...')
with open("frequencies_100.json","r") as f:
    frequencies=json.load(f)
for x in range(1,READ+1):
    frequencies['pos%d_frequencies'%x]=list(enumerate(frequencies['pos%d_frequencies'%x]))
    frequencies['pos%d_frequencies'%x].sort(key=lambda y:y[1],reverse=True) 
print('\ndown...')
# creat qphread 100
phred_reads=[]
for x in range(100):
    phred_reads.append(random_qphred(frequencies))
frequencies=phred_reads

# get exon total number
filein='ENSEMBL_exonlist.txt'
total_exon=countrow(filein,r'^chr([\w\d]*)\t(\d*)\t(\d*)\t([\w\.]*)')
print("exon count: ",total_exon)

#set and print gene mutation type
mutation_types=mutationtype.set_mutation(total_exon)
input('print "enter" to continue...')
'''
with open('mutation_settings.json','w') as f:
    for x in mutation_types:
        jobj=json.dumps(x, default=lambda obj: obj.__dict__)
        json.dump(jobj, f)
with open('mutation_settings.json','r') as f:
    mutation_types=json.loads(f)
    print(mutation_types)
    #mutationtype.show_mutation(mutation_types)
'''


# according to mutaton output reads
# output normal reads(not snp mutation)
time_start=time.time()

print('output normal reads(cnv+deletion)...')
t=time.strftime('%Y%m%d_%H_%M_%S',time.localtime(time.time()))
with open("R1_fastq%s.fastq" % t,'w') as w1:
    with open("R2_fastq%s.fastq" % t,'w') as w2:
        # one exon generate to mutation_types
        with open(filein,'r')as r:
            for exon,n in get_normal_exon(r):
                fastq_exon(w1,w2,exon,frequencies,mutation_types,n,False)
        print('\ndown')
        # output mutation reads(snp mutation)
        print('output snp mutation reads(cnv+snp)...')
        i=0
        for x in mutation_types:
            # one exon generate to one mutation_type
            with open(filein,'r')as r:
                for exon,mutainfo in get_snp_exon(r,x.homologs[0].snp_mutations):
                    fastq_exon(w1,w2,exon,frequencies,mutation_types,mutainfo,True,muta=i,homolog=0)
            with open(filein,'r')as r:
                for exon,mutainfo in get_snp_exon(r,x.homologs[1].snp_mutations):
                    fastq_exon(w1,w2,exon,frequencies,mutation_types,mutainfo,True,muta=i,homolog=1)
            i+=1
print('\n***down***')

time_end=time.time()
t=time_end-time_start
print('totally cost: %dh : %dm : %ds'%(t//3600,(t%3600)//60,t%60))   
