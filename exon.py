import re
import sys

from filefunc import *


class Exon(object):

    # creat exon number counter
    __count =0
    def __init__(self,chr=0,begin=0,end=0,exon_id=None,seq=None,mutations=None,insert=None):
        '''
        exon_id usually is combination of chr and it's no. on chr;
        seq is str
        insert is two side's length
        mutations is it's mutations info,type=list
        eg.
        exon_id like'CH-156'
        seq='AGTCG|TGCGATGC|GTAGT',insert=5,then ture exon seq='TGCGATGC'('|' is just mark)
        mutations=[{'cn'=3,'description'='cnv'},
        {'pos'=56,;'ref'='AGT','alt'='.','cn'=5,'description'='snp+cnv'}]
        mutations can be one dict,like {'cn'=3,'description'='cnv'}
        '''
        self.__chr=chr
        self.__begin=begin
        self.__end=end
        self.__exon_id=exon_id
        self.seq=seq
        if insert==None:
            self.__insert=0
        else:
            self.__insert=int(insert)
        self.mutations=[]
        if mutations:
            self.add_mutation(mutations)

    @property
    def chr(self):
        return self.__chr
    @property
    def begin(self): 
        return self.__begin
    @property
    def end(self):
        return self.__end
    @property
    def exon_id(self):
        return self.__exon_id
    @property
    def length(self):
        return self.__end-self.__begin
    @property
    def insert(self):
        return self.__insert
    def getexon_info(self):
        '''return chr,begin,end,exon_id,length,insert;split is tab'''
        info="chr%s\t%d\t%d\t%s\t%d\t%d" %(
                self.chr,self.begin,self.end,self.exon_id,self.length,self.insert)
        return info

    def add_mutation(self,mutation):
        '''
        'mutation like:
        {'cn'=3,'pos':52,ref='ATGG',alt='T'}
        or the list of such dic
        [{'cn'=3,'pos':52,ref='ATGG',alt='T'}]
        '''
        if isinstance(mutation, dict):
            mutation=[mutation]
        if isinstance(mutation, list):
            for x in mutation:
                if x not in self.mutations:
                    if x['pos']=='.':
                        if self.get_normal_num()==1:
                            self.mutations.append(x)
                        else:
                            print('exon has ownned deletion or cnv')
                    else:
                        self.mutations.append(x)
        
    def show(self,pos1=1,pos2=-1,all=False):
        if pos2==-1:
            pos2=self.length
        print(self.getexon_info())
        if all:
            print(self.seq[self.insert:self.insert+pos1-1],'|',self.seq[self.insert+pos1-1:self.insert+pos2],'|',self.seq[self.insert+pos2:-self.insert])
        else:
            print('|',self.seq[self.insert+pos1-1:self.insert+pos2],'|')

    def get_normal_num(self):
        for x in self.mutations:
            if x['pos']=='.':
                return int(x['cn'])
        return 1
    def show_mutations(self):
        print('pos\tref\talt\t cn\tdescription')
        for x in self.mutations:
            for y in x.values():
                print(y,end='\t')
            print('')
    @property
    def snp_mutations(self):
        '''return mutations which sequence has change'''
        snp_mutations=[]
        for x in self.mutations:
            if x['pos']!='.':
                snp_mutations.append(x)
        return snp_mutations
    def __get_one_seq(self,onemutation,chip_len):
        '''according one mutaton(dict) ,return new seq,length,copynum'''
        if onemutation['pos']=='.':
            return self.seq,self.length,1,'NORMAL'
        else:
            pos=int(onemutation['pos'])
            length=0
            mid=''
            if (onemutation['ref'])!='.':
                length=len(onemutation)
            if (onemutation['alt'])!='.':
                mid=onemutation['alt']
            newseq=self.seq[:self.insert+pos]+mid+self.seq[self.insert+pos+length-1:]
            seqlength=self.length+len(mid)-length
            if seqlength > chip_len:
                return newseq,seqlength,int(onemutation['cn']),tuple(onemutation.values())
            else:
                return '',0,0,onemutation

    def get_new_seq(self,mutation,chip_len):
        '''according mutatons(list of dicts) ,return these mutations' new seq,length,copynum'''
        if isinstance(mutation, dict):
            mutation=[mutation]
        if isinstance(mutation, list):
            for x in mutation:
                yield self.__get_one_seq(x,chip_len)



            

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
    def exonseq(self,pos,length,new_column):
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
            self.fileout.write("%s\n" % (exon.getexon_info()))
            po = self.pos+ self.lentopos(exon.begin-exon.insert)
            length=exon.length+2*exon.insert
            # write exon sequence
            self.filein.seek(self.pos,0)
            line = self.filein.readline()

            self.exonseq(po,length,new_column)
            try:
                exon=next(exons)
            except StopIteration:
                exon =0
            else:
                pass
