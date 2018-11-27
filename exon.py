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
        return self.__end-self.__begin+1
    @property
    def insert(self):
        return self.__insert
    def set_insert(self,insert):
        self.__insert=insert

    def getexon_info(self):
        '''return chr,begin,end,exon_id,length,insert;split is tab'''
        info="chr%s\t%d\t%d\t%s\t%d\t%d" %(
                self.chr,self.begin,self.end,self.exon_id,self.length,self.insert)
        return info

    def add_mutation(self,mutation):
        '''
        'mutation like:
        {'cn'='.','pos':52,ref='ATGG',alt='T'}
        {'cn'=3,'pos':52,ref='ATGG',alt='T'}
        or the list of such dic
        [{'cn'=3,'pos':52,ref='ATGG',alt='T'}]
        '''
        if isinstance(mutation, dict):
            mutation=[mutation]
        if isinstance(mutation, list):
            for x in mutation:
                if x not in self.mutations and x['cn']!='.':
                    if x['start'][0]=='.':
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
            if x['start'][0]=='.':
                return int(x['cn'])
        return 1
    
    def show_mutations(self):
        print('start\tend\tref\talt\t cn\tdescription')
        for x in self.mutations:
            for y in x.values():
                if isinstance(y,(tuple,list)):
                    print(','.join(y),end='\t')
                else:
                    print(y,end='\t')
            print('')
    
    @property
    def abnormal_mutations(self):
        '''return mutations which sequence has change'''
        abnormal_mutations=[]
        for x in self.mutations:
            if x['start'][0]!='.':
                abnormal_mutations.append(x)
        return abnormal_mutations
    
    def __get_one_seq(self,onemutation,chip_len):
        '''according one mutaton(dict) ,return new seq,length,copynum'''
        if onemutation['start'][0]=='.':
            return self.seq,self.length,1,'NORMAL'
        else:
            cur=0
            newseq=''
            for w,x,y,z in zip(onemutation['start'],onemutation['end'],onemutation['ref'],onemutation['alt']):
                try:
                    w=int(w)
                    x=int(x)
                except:
                    print(onemutation)
                y=len(y)
                if z=='.':
                    z=''
                newseq+=self.seq[cur:self.insert+w-1]+z
                cur=self.insert+x
            newseq+=self.seq[cur:]
            seqlength=len(newseq)-2*self.insert
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

   



def get_seq(filed):
    seq=''
    line=filed.readline()
    n = re.match(r'^chr',line)
    while(not n):
        seq=seq+line.strip()
        line=filed.readline()
        if not line:
            return seq
        n = re.match(r'^chr',line)
    return seq,line