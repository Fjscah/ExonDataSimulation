import json
import random
import re
import sys
import time

import mutation
from exon import Exon
from exonlist import get_exons
from phred import random_qphred
from setting import (ACCURACY_RATE, CHIP_LEN, DEEPTH, INSERT_D,INSERT_E, PHRED, ROW_STEP,
                     SUBSTITUTION,DEFAULTS,ERR_PH)
from view import get_chr_exon_num, show_chr_exon_num

ATCG=('A','T','C','G','a','t','c','g')
COMPLEMENT={'A':'T','G':'C','C':'G','T':'A'}


def add_read_err(read,phred,accuracy=ACCURACY_RATE):
    size=len(read)
    err_num=round((1-accuracy)*size)
    err_num_list=random.sample(range(size),err_num)
    for n in err_num_list:
        pich=random.randint(0,2)
        read[n]=SUBSTITUTION[read[n].upper()][pich]
        phred_value=ord(phred[n])-33
        if phred_value>ERR_PH:
            phred[n]=chr(ord(phred[n])-ERR_PH)
        else:
            phred[n]=chr(33)

def complementation(read):
    l=[]
    for char in read:
        l.append(COMPLEMENT[char.upper()])
    return l

def PEread(seq,length,frequencies,readlen,exoninsert):
    ''' get one inser's 2 reads'''
    # chip pos on exon
    cbp=random.randint(0,length-CHIP_LEN)
    while(True):
        # insertion length
        insertion=round(random.normalvariate(INSERT_D,INSERT_E))
        # chip pos on inseriton
        sbp=random.randint(0,insertion-CHIP_LEN)
        # insertion pos on exon
        lbp=exoninsert+cbp-sbp
        if lbp<0:
            insertion=lbp+insertion
            lbp=0
        rbp=lbp+insertion-1
        if rbp+1>length+2*exoninsert:
            insertion=length+2*exoninsert-rbp-1+insertion
            rbp=length+2*exoninsert
        if insertion>=readlen:
            raise Exception("rbp+1>length+2*exoninsert")
            #break
    l=[]
    for n in range(lbp,rbp):
        char =seq[n]
        if char in ATCG:
            l.append(char.upper())
        else :
            raise Exception("error :can't indentify",char)
            #l.append(random.choice(ATCG))
    read1=l[0:readlen]#list
    read2=complementation(l[-1:-1-readlen:-1])#list
    phred1=list(random.choice(frequencies))
    phred2=list(random.choice(frequencies))
    add_read_err(read1,phred1)
    add_read_err(read2,phred2)
    return ''.join(read1) , ''.join(read2),''.join(phred1),''.join(phred2)

def fastq_exon(w1,w2,exon,frequencies,mutation_types,readlen):
    ''' one exon's all fastq reault be writed to w'''
    print('writing chr %s - exon : %s  - %s' %(exon.chr,exon.exon_id,exon.begin) ,end='\r')
    if 'N' in exon.seq:
        print('error: not output \n'+''.join(exon.seq))
        return
    exon_id=exon.exon_id
    seq,length=list(exon.seq),exon.length
    for x in mutation_types:
        copy_num=mutation_types[x].haplots[0].get_normal_num(exon_id)+mutation_types[x].haplots[1].get_normal_num(exon_id)
        turns=round(length*mutation_types[x].percent*copy_num*DEEPTH/2/INSERT_D/100)
        for y in range (turns):
            read1,read2,phred1,phred2=PEread(seq,length,frequencies,readlen,exon.insert)
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
                for seq,length,copy_num,onemutation in exon.get_new_seq((mutation_types[x].haplots[z].mutations[exon_id].abnormal_mutations),CHIP_LEN):
                    turns=round(length*mutation_types[x].percent*copy_num*DEEPTH/2/INSERT_D/100)
                    seq=list(seq)
                    for y in range(turns):
                        read1,read2,phred1,phred2=PEread(seq,length,frequencies,readlen,exon.insert)
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
    ename,qname,mname=DEFAULTS['file3'],DEFAULTS['file4'],DEFAULTS['filew']
    if len(sys.argv)==4:
        pyname,ename,qname,mname=sys.argv

    #set and print gene mutation type
    mutation_types=mutation.file_mutation(mname)
    #input('press Enter to continue...')

    # initial qphred frequencies
    print('initial qphred...')
    with open(qname,"r") as f:
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
    readlen=len(frequencies)


    # get every chr's exon number
    '''
    exon_nums=get_chr_exon_num(filein)
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
            with open(ename,'r')as r:
                for exon in get_exons(r):
                    fastq_exon(w1,w2,exon,frequencies,mutation_types,readlen)
    print('\n***down***')
    time_end=time.time()
    t=time_end-time_start
    print('totally cost: %dh : %dm : %ds'%(t//3600,(t%3600)//60,t%60))   
