import json
import os
import re
import sys
import time

import view
from exon import Exon, WholeExon
from filefunc import *
from setting import *

print("default:","INSERT maxlength=",MAXINSERT,"CHIP length=",CHIP_LEN)
def chr2num(chrr):
    if chrr in ('X','x'):
        return '23'
    elif chrr in ('Y','y'):
        return '24'
    elif chrr in ('12920','M','m'):
        return '25'
    return chrr
def char2phred(char,plus=33):
    code=ord(char)-plus
    return code
def get_exoninfo(line,word,chip_len):
    '''according ever row to genarater Exon object'''
    m = re.match(word,line)
    if m :
        if int(m.group(3))-int(m.group(2))>=chip_len:
            chrr=chr2num(m.group(1))
            return Exon(chrr,int(m.group(2)),int(m.group(3)))
        else:
            return None
    else:
        return None
def exons_generater(f,word,join_gap,chip_len,maxinsert):
    '''
    according file content to generate Exon object's itertator
    f is file stream; 
    join_gap is the shortest distance of adjaccent, ohtherwise join it
    '''
    l_exons=[]
    echr=0
    join_gap
    for line in f.readlines():
        m=get_exoninfo(line,word,chip_len)
        if not m:
            continue
        if m.chr!=echr:
            eid=1
            if len(l_exons)==1:
                yield Exon(echr,l_exons[0][0],l_exons[0][1],'CH'+str(echr)+'-'+str(eid),insert=maxinsert)
            elif len(l_exons)>1:
                lorder=sorted(l_exons)
                size =len(lorder)
                x=0
                y=1
                end=lorder[x][1]
                while(x<=size-2 and y<=size-1):
                    if end+1+join_gap>=lorder[y][0]:
                        end =max(end,lorder[y][1])
                        y+=1
                    elif end+1+join_gap<lorder[y][0]:
                        #print("yield",eid)
                        yield Exon(echr,lorder[x][0],end,'CH'+str(echr)+'-'+str(eid),insert=maxinsert)
                        eid+=1
                        end=lorder[y][1]
                        x,y=y,y+1
                yield Exon(echr,lorder[x][0],lorder[-1][1],'CH'+str(echr)+'-'+str(eid),insert=maxinsert)
            l_exons=[]
            echr=m.chr
        l_exons.append((m.begin,m.end))
    eid=1
    if len(l_exons)==1:
        yield Exon(echr,l_exons[0][0],l_exons[0][1],'CH'+str(echr)+'-'+str(eid),insert=maxinsert)
    elif len(l_exons)!=0:
        lorder=sorted(l_exons)
        size =len(lorder)
        x=0
        y=1
        end=lorder[x][1]
        while(x<=size-2 and y<=size-1):
            if end+1+join_gap>=lorder[y][0]:
                end =max(end,lorder[y][1])
                y+=1
            elif end+1+join_gap<lorder[y][0]:
                #print("yield",eid)
                yield Exon(echr,lorder[x][0],end,'CH'+str(echr)+'-'+str(eid),insert=maxinsert)
                eid+=1
                end=lorder[y][1]
                x,y=y,y+1
        yield Exon(echr,lorder[x][0],lorder[-1][1],'CH'+str(echr)+'-'+str(eid),insert=maxinsert)
            

def get_chr_seq(file,write,ver):
    word=eval(ver+"['seq']")
    # get chromosome sequence from grch38 genomic
    print("get chromosome sequence from %s..."%file)
    get_content(file,write,word,r'>.*')
    print('\ndown. outfile: %s'%write)
def get_exon_ano(file,write,ver):
    word=eval(ver+"['ano']")
    # get exons annotation from v29 annotation
    print("get exons annotation from %s..."%file)
    get_content(file,write,word,r'.*')
    print('\ndown. outfile: %s'%write)
def get_exon_list(file1,file2,write):
    word=eval(ver+"['list']")
    # get exon sequence ; write to exonlist.txt
    print("get exon initial list from %s , %s..."%(file1,file2))
    file1_column=get_row_column(file1,'>',1,'re')
    with open(file3,'w') as w:
        w.write("#chr\tbegin\tend\tgeneid\tlength\tINSERT\n")
        with open (file1,'r') as r1:
            with open (file2,'r') as r2:
                exons =exons_generater(r2,word,JOIN_GAP,CHIP_LEN,MAXINSERT)
                wholeexon = WholeExon(r1,w,file1_column,ROW_STEP)
                wholeexon.wholeexonseq(exons,file1_column)
    print('\ndown. outfile: %s'%write)
def get_all(defaults):
    filex,filey=defaults['filex'],defaults['filey']
    file1,file2,file3=defaults['file1'],defaults['file2'],defaults['file3']
    ver=defaults['ver']
    get_chr_seq(filex,file1,ver)
    get_exon_ano(filey,file2,ver)
    get_exon_list(file1,file2,file3)
def get_phred_fre(file,write,plus):
    print('get phred frequencies from %s...'%file)
    frequencies={}
    readlen=get_row_column(file,'+',1,'re')
    for x in range(1,readlen+1):
        frequencies['pos%d_frequencies'%x]=[0]*43
    with open(file,'r') as f:
        i =0
        for line in f.readlines():
            i+=1
            if i%4==0:
                row_fastq=line.strip()
                x=1
                for char in row_fastq:
                    frequencies['pos%d_frequencies'%x][char2phred(char,plus)]+=1
                    x=x+1
            if i%40000==0:
                print(i,end='\r')
    with open(write,"w") as f:
        json.dump(frequencies, f)
    print('\ndown. outfile: %s'%write)



if __name__ == '__main__':
    help='''
    init -seq -filex -file1 -ver    : get chromosome sequence from filex, generate file1
    init -ano -filey -file2 -ver    : get exon annotation from filey, generate file2
    init -list -file1 -file2 -file3 -ver :get exonlist from file1,file2, generate file3
    init -all                       :excute above all according setting,py
    phred -filez -file4 -XX         : get qphred frequencies from file,generate 'phred.json'
    view -file3                     : view exonlist's exon sequence
    '''
    print(help)
    opera=input('>')
    while(opera.strip()!='exit'):
        if opera.strip=='help':
            print(help)
        elif re.match(r'init -seq',opera):
            info=opera.split(' -')
            if len(info)>=5:
                filex,file1,ver=info[2:]
                get_chr_seq(filex,file1,ver)
        elif re.match(r'init -ano',opera):
            info=opera.split(' -')
            if len(info)>=5:
                fiely,file2,ver=info[2:]
                get_exon_ano(fiely,file2,ver)
        elif re.match(r'init -list',opera):
            info=opera.split(' -')
            if len(info)>=6:
                file1,file2,file3,ver=info[2:]
                get_exon_list(file1,file2,file3)
        elif re.match(r'init -all',opera):
            get_all(DEFAULTS)
        elif re.match(r'phred',opera):
            info=opera.split(' -')
            if len(info)>=3:
                filex,file4,plus=info[1:]
                get_phred_fre(filex,file4,PHRED)
        elif re.match(r'view',opera):
            info=opera.split(' -')
            if len(info)>=2:
                file=info[1]
                view.view(file)
        opera=input('>')
