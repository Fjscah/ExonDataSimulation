import random
import re
from exon import Exon

def positive_value(value):
    if value<0:
        value=0
    return value
        
class HaploType():

    
    def __init__(self,mutations=None):
        '''
        mutattions is a dict,like:
        {'CH1-563':Exon(),'CH15-53':Exon()}
        CH1-563 is exon's id,ch1 is it's chromosome,
        563 stands it's no.563 on this chromosome
        '''
        self.__mutations={}
        if mutations:
            self.add_mutation(mutations)
    @property
    def mutations(self):
        return self.__mutations
    def add_mutation(self,mutations):
        if isinstance(mutations,dict):
            mutations=[mutations]
        if isinstance(mutations[0],dict):
            for x in mutations:
                exon_id='CH'+x['chr']+'-'+x['exon']
                if exon_id not in self.__mutations:
                    self.__mutations[exon_id]=Exon(exon_id=exon_id)
                x.pop('chr')
                x.pop('exon')
                self.__mutations[exon_id].add_mutation(x)
        else:
            print("it's not avaliable type")
    
    def get_normal_num(self,exon_id):
        if exon_id in self.__mutations:
            return self.__mutations[exon_id].get_normal_num()
        return 1

    def show_mutations(self):
        if self.__mutations:
            print('mutation type')
        else:
            print('normal type')
        for x in sorted(self.__mutations.keys()):
            print(x,':')
            self.__mutations[x].show_mutations()
                


class GenomeType():
    def __init__(self,percent,haplot1,haplot2):
        self.percent=percent
        self.haplots=[haplot1,haplot2]
    def show_genetype(self):
        print('***percent=%d***'%self.percent)
        for x in range(len(self.haplots)):
            print('haplot%d : '% (x+1),end='')
            self.haplots[x].show_mutations()



def set_mutation():
    choice=input('''which mutation creating way do you want?
        1. set all typess of mutation's total number
        2. set every mutatition manually
        3. import mutation file
        >''')
    mutation_types={}
    if choice=='1':
        mutation_types=auto_mutaion()
    elif choice=='2':
        mutation_types=manual_mutation()
    elif choice=='3':
        mutation_types=file_mutation()
    elif choice=='exit':
        return 'exit'
    else:
        print("your input cannot be identified.")
        return 'exit'
    
    # show mutation setting
    show_mutation(mutation_types)
    return mutation_types


def show_mutation(mutation_types):
    print('*'*10+'Show Genetype'+'*'*10)
    for x in mutation_types.values():
        x.show_genetype()
    print('*'*10+'Show Genetype'+'*'*10)
    
def auto_mutaion():
    mutation_types={}
    return mutation_types
def manual_mutation():
    mutation_types={}
    return mutation_types

def searchline(f,text):
    f.seek(0,0)
    line =f.readline()
    while(line):
        if re.match(text,line):
            m=re.match(text,line)
            return 
        line=f.readline()

def file_mutation():
    filename=input('please input mutation setting file : ')
    mutation_types={}
    keys=['chr','exon','pos','ref','alt','cn','description']
    with open(filename,'r') as f:
        searchline(f,r'>genometype:')
        line=f.readline()
        percent=float(100)
        while(line):
            line=f.readline()
            if re.match(r'---',line) or percent<0:
                break
            mutas=line.strip().split()
            mutation_types[mutas[0]]=GenomeType(float(mutas[1]),HaploType(),HaploType())
            percent-=float(mutas[1])
        for key in mutation_types:
            searchline(f,r'>mutation:%s-1'%key)
            line=f.readline()
            while(line):
                line=f.readline()
                if re.match(r'>',line) :
                    break
                muta=line.strip().split()
                muta=dict(zip(keys,muta))
                mutation_types[key].haplots[0].add_mutation(muta)
            searchline(f,r'>mutation:%s-2'%key)
            line=f.readline()
            while(line):
                line=f.readline()
                if re.match(r'>',line) :
                    break
                muta=line.strip().split()
                if len(muta)==7:
                    muta[5]=int(muta[5])
                if muta:
                    muta=dict(zip(keys,muta))
                    mutation_types[key].haplots[1].add_mutation(muta)
    return mutation_types



        
