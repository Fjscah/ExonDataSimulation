import copy
import random
import re

from exon import Exon


class HaploType():

    def __init__(self, mutations=None):
        '''
        mutattions is a dict,like:
        {'CH1-563':Exon(),'CH15-53':Exon()}
        CH1-563 is exon's id,ch1 is it's chromosome,
        563 stands it's no.563 on this chromosome
        '''
        self.__mutations = {}
        if mutations:
            self.add_mutation(mutations)

    @property
    def mutations(self):
        return self.__mutations

    def add_mutation(self, mutations):
        if isinstance(mutations, dict):
            mutations = [mutations]
        if isinstance(mutations[0], dict):
            for x in mutations:
                exon_id = 'CH'+x['chr']+'-'+x['exon']
                if exon_id not in self.__mutations:
                    self.__mutations[exon_id] = Exon(exon_id=exon_id)
                x.pop('chr')
                x.pop('exon')
                self.__mutations[exon_id].add_mutation(x)
        else:
            print("it's not avaliable type")

    def get_normal_num(self, exon_id):
        if exon_id in self.__mutations:
            return self.__mutations[exon_id].get_normal_num()
        return 1

    def show_mutations(self):
        if self.__mutations:
            print('mutation type')
        else:
            print('normal type')
        for x in sorted(self.__mutations.keys()):
            print(x, ':')
            self.__mutations[x].show_mutations()


class GenomeType():
    def __init__(self, percent, haplot1, haplot2):
        self.percent = percent
        self.haplots = [haplot1, haplot2]

    def show_genetype(self):
        print('***percent=%d***' % self.percent)
        for x in range(len(self.haplots)):
            print('haplot%d : ' % (x+1), end='')
            self.haplots[x].show_mutations()

    def add_mutation(self, muta):
        muta1 = copy.deepcopy(muta)
        muta1['cn'] = muta['cn'][0]
        if muta1['cn'] != '.':
            self.haplots[0].add_mutation(muta1)
        muta2 = copy.deepcopy(muta)
        muta2['cn'] = muta['cn'][1]
        if muta2['cn'] != '.':
            self.haplots[1].add_mutation(muta2)


def set_mutation():
    choice = input('''which mutation creating way do you want?
        1. set all typess of mutation's total number
        2. set every mutatition manually
        3. import mutation file
        >''')
    mutation_types = {}
    if choice == '1':
        mutation_types = auto_mutaion()
    elif choice == '2':
        mutation_types = manual_mutation()
    elif choice == '3':
        mname = input('please input mutation setting file : ')
        mutation_types = file_mutation(mname)
    elif choice == 'exit':
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
    mutation_types = {}
    return mutation_types


def manual_mutation():
    mutation_types = {}
    return mutation_types


def searchline(f, text):
    f.seek(0, 0)
    line = f.readline()
    while(line):
        if re.match(text, line):
            m = re.match(text, line)
            return
        line = f.readline()


def file_mutation(mname):
    mutation_types = {}
    keys = ['chr', 'exon', 'start', 'end', 'ref', 'alt', 'cn', 'description']
    with open(mname, 'r') as f:
        searchline(f, r'>genometype:')
        line = f.readline()
        percent = float(100)
        while(line):
            line = f.readline()
            if re.match(r'---', line) or percent < 0:
                break
            mutas = line.strip().split()
            mutation_types[mutas[0]] = GenomeType(
                float(mutas[1]), HaploType(), HaploType())
            percent -= float(mutas[1])
        for key in mutation_types:
            searchline(f, r'>mutation:%s' % key)
            line = f.readline()
            while(line):
                line = f.readline()
                if re.match(r'>', line):
                    break
                muta = line.strip().split()
                if len(muta) >= 7:
                    muta = dict(zip(keys, muta))
                    muta['start'] = muta['start'].split(',')
                    muta['end'] = muta['end'].split(',')
                    muta['ref'] = muta['ref'].split(',')
                    muta['alt'] = muta['alt'].split(',')
                    muta['cn'] = muta['cn'].split('/')
                    mutation_types[key].add_mutation(muta)
    return mutation_types


if __name__ == '__main__':
    set_mutation()
