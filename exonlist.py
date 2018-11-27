import json
import os
import re
import sys
import time


from exon import Exon,get_seq
from filefunc import write_content,write_column,get_column_row
from setting import *


def chr2num(chrr):
    if chrr in ('X','x'):
        return '23'
    elif chrr in ('Y','y'):
        return '24'
    elif chrr in ('12920','M','m'):
        return '25'
    return chrr
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
class WholeExon(object):
    # coulum is the numbeer of dNTP in a row from filein, step is usually one byte which stand for '\n'
    def __init__(self,fin,fout,column,step=2):
        self.filein=fin
        self.fileout=fout
        self.pos=0
        self.chr=0
        self.column=column
        self.step=step
    # write sequence of certain exon according to it's begin pos and it's length
    def exonseq(self,exon,pos,length,new_column):
        self.filein.seek(pos,0)
        seq=""
        sys.stdout.write(str(self.filein.tell()))
        sys.stdout.write('-->')
        line=self.filein.readline().strip()
        # sp is set for alignment, in other word, write the same number of dNTP in every row
        while(length>0):
            length=length-len(line)
            if length>=0:
                seq+=line
            else:
                seq+=line[:length]
            line=self.filein.readline().strip()
        if 'N' in seq:
            l=seq[0:exon.insert].rfind('N')
            r=seq[-1:-1-exon.insert:-1].rfind('N')
            s=max(l,r)+1#s>=1
            if s:
                seq=seq[s:-s]
                exon.set_insert(exon.insert-s)
                print(exon.getexon_info(),'\n',seq)
        self.fileout.write("%s\n" % (exon.getexon_info()))
        write_column(self.fileout,seq,new_column)
        sys.stdout.write(str(self.filein.tell())+"     ")
        sys.stdout.write('\r')
    def lentopos(self,length):
        raw=length//self.column
        return raw*(self.column+self.step)+length %self.column-1
    # according exons informatin write its sequence ; the exon must be in order from 1 to 24
    def wholeexonseq(self,exons,new_column):
        new_column=self.column
        exon=next(exons)
        while(exon):
            tempchr=exon.chr
            if not tempchr==self.chr:
                print("exon end:",self.filein.tell(),"     ")
                self.chr=tempchr
                line=self.filein.readline()
                row=0
                while(line):
                    row+=1
                    sys.stdout.write("%d"%row)
                    sys.stdout.write('\r')
                    if re.match(r'>',line):
                        print("chr",self.chr,"pos",self.filein.tell(),line.strip())
                        # change to next chromosome 
                        self.pos=self.filein.tell()
                        row=0
                        break
                    try:
                        line = self.filein.readline()
                    except StopIteration:
                        line =0
                    else:
                        pass

            # write exon title
            #print('write exon ',exon.getexon_info())
            po = self.pos+ self.lentopos(exon.begin-exon.insert)
            length=exon.length+2*exon.insert
            # write exon sequence
            self.filein.seek(self.pos,0)
            line = self.filein.readline()

            self.exonseq(exon,po,length,new_column)
            try:
                exon=next(exons)
            except StopIteration:
                exon =0
            else:
                pass            

def get_chr_seq(file,write,ver):
    word=eval(ver+"['seq']")
    if isinstance(word,str):
        match='re'
    elif isinstance(word,list,tuple):
        match='=='
    # get chromosome sequence from grch38 genomic
    print("get chromosome sequence from %s..."%file)
    write_content(file,write,word,r'>.*',0,0,match,'re')
    print('\ndown. outfile: %s'%write)
def get_exon_ano(file,write,ver):
    word=eval(ver+"['ano']")
    # get exons annotation from v29 annotation
    print("get exons annotation from %s..."%file)
    write_content(file,write,word,r'.*',0,0,'re','re')
    print('\ndown. outfile: %s'%write)
def get_exon_list(file1,file2,write,ver,chip_len,join_gap,row_step,maxinsert):
    word=eval(ver+"['list']")
    # get exon sequence ; write to exonlist.txt
    print("get exon initial list from %s , %s..."%(file1,file2))
    file1_column=get_column_row(file1,'>',1,'re')
    with open(write,'w') as w:
        w.write("#chr\tbegin\tend\tgeneid\tlength\tINSERT\n")
        with open (file1,'r') as r1:
            with open (file2,'r') as r2:
                exons =exons_generater(r2,word,join_gap,chip_len,maxinsert)
                wholeexon = WholeExon(r1,w,file1_column,row_step)
                wholeexon.wholeexonseq(exons,file1_column)
    print('\ndown. outfile: %s'%write)


def get_exons(file):
    ''' from every row in flie get exon info and yield exon object'''
    line=file.readline()
    while(line):
        m = re.match(r'^chr([\d]*)\t(\d*)\t(\d*)\t([\w\.\-]*)\t\d*\t(\d*)',line)
        if m:
            # if m is exon title ,then store this exon's sequence
            seq,line=get_seq(file)
            yield Exon(m.group(1),int(m.group(2)),int(m.group(3)),m.group(4),seq,insert=int(m.group(5)))
        else:
            line=file.readline()
