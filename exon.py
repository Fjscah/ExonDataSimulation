import json
import os
import re
import sys
import time
from collections import Iterator

from filefunc import equal_text, get_column_row, get_line_text, write_column
from setting import *


class Exon(object):

    # creat exon number counter
    __count = 0

    def __init__(self, chr=0, begin=0, end=0, exon_id=None, seq=None, mutations=None, insert=None):
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
        self.__chr = chr
        self.__begin = begin
        self.__end = end
        self.__exon_id = exon_id
        self.seq = seq
        if insert == None:
            self.__insert = 0
        else:
            self.__insert = int(insert)
        self.mutations = []
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

    def set_insert(self, insert):
        self.__insert = insert

    def getexon_info(self):
        '''return chr,begin,end,exon_id,length,insert;split is tab'''
        info = ">chr%s\t%d\t%d\t%s\t%d\t%d" % (
            self.chr, self.begin, self.end, self.exon_id, self.length, self.insert)
        return info

    def add_mutation(self, mutation):
        '''
        'mutation like:
        {'cn'='.','pos':52,ref='ATGG',alt='T'}
        {'cn'=3,'pos':52,ref='ATGG',alt='T'}
        or the list of such dic
        [{'cn'=3,'pos':52,ref='ATGG',alt='T'}]
        '''
        if isinstance(mutation, dict):
            mutation = [mutation]
        if isinstance(mutation, list):
            for x in mutation:
                if x not in self.mutations and x['cn'] != '.':
                    if x['start'][0] == '.':
                        if self.get_normal_num() == 1:
                            self.mutations.append(x)
                        else:
                            print('exon has ownned deletion or cnv')
                    else:
                        self.mutations.append(x)

    def show(self, pos1=1, pos2=-1, all=False):
        if pos2 == -1:
            pos2 = self.length
        print(self.getexon_info())
        if all:
            print(self.seq[self.insert:self.insert+pos1-1], '|', self.seq[self.insert +
                                                                          pos1-1:self.insert+pos2], '|', self.seq[self.insert+pos2:-self.insert])
        else:
            print('|', self.seq[self.insert+pos1-1:self.insert+pos2], '|')

    def get_normal_num(self):
        for x in self.mutations:
            if x['start'][0] == '.':
                return int(x['cn'])
        return 1

    def show_mutations(self):
        print('start\tend\tref\talt\t cn\tdescription')
        for x in self.mutations:
            for y in x.values():
                if isinstance(y, (tuple, list)):
                    print(','.join(y), end='\t')
                else:
                    print(y, end='\t')
            print('')

    @property
    def abnormal_mutations(self):
        '''return mutations which sequence has change'''
        abnormal_mutations = []
        for x in self.mutations:
            if x['start'][0] != '.':
                abnormal_mutations.append(x)
        return abnormal_mutations

    def __get_one_seq(self, onemutation, chip_len):
        '''according one mutaton(dict) ,return new seq,length,copynum'''
        if onemutation['start'][0] == '.':
            return self.seq, self.length, 1, 'NORMAL'
        else:
            cur = 0
            newseq = ''
            for w, x, y, z in zip(onemutation['start'], onemutation['end'], onemutation['ref'], onemutation['alt']):
                try:
                    w = int(w)
                    x = int(x)
                except:
                    print(onemutation)
                y = len(y)
                if z == '.':
                    z = ''
                newseq += self.seq[cur:self.insert+w-1]+z
                cur = self.insert+x
            newseq += self.seq[cur:]
            seqlength = len(newseq)-2*self.insert
            if seqlength > chip_len:
                return newseq, seqlength, int(onemutation['cn']), "-".join(str(i) for i in onemutation['start'])
            else:
                return '', 0, 0, onemutation

    def get_new_seq(self, mutation, chip_len):
        '''according mutatons(list of dicts) ,return these mutations' new seq,length,copynum'''
        if isinstance(mutation, dict):
            mutation = [mutation]
        if isinstance(mutation, list):
            for x in mutation:
                yield self.__get_one_seq(x, chip_len)


def get_seq(filed):
    seq = ''
    line = filed.readline()
    n = re.match(r'^>', line)
    while(not n):
        seq = seq+line.strip()
        line = filed.readline()
        if not line:
            return seq,line
        n = re.match(r'^>', line)
    return seq, line



def chr2num(chrr):
    if chrr in ('X', 'x'):
        return '23'
    elif chrr in ('Y', 'y'):
        return '24'
    elif chrr in ('12920', 'M', 'm'):
        return '25'
    return chrr


def get_exoninfo(line, chip_len):
    '''according ever row to genarater Exon object'''
    m = line.split()
    if m:
        if int(m[2])-int(m[1]) >= chip_len:
            chrr = m[0]
            return Exon(chrr, int(m[1]), int(m[2]))
        else:
            return None
    else:
        return None


def exons_generater(f, join_gap, chip_len, maxinsert):
    '''
    according file content to generate Exon object's itertator
    f is file stream; 
    join_gap is the shortest distance of adjaccent, ohtherwise join it
    '''
    l_exons = []
    echr = 0
    join_gap
    for line in f.readlines():
        m = get_exoninfo(line,chip_len)
        if not m:
            continue
        if m.chr != echr:
            eid = 1
            if len(l_exons) == 1:
                yield Exon(echr, l_exons[0][0], l_exons[0][1], 'CH'+str(echr)+'-'+str(eid), insert=maxinsert)
            elif len(l_exons) > 1:
                lorder = sorted(l_exons)
                size = len(lorder)
                x = 0
                y = 1
                end = lorder[x][1]
                while(x <= size-2 and y <= size-1):
                    if end+1+join_gap >= lorder[y][0]:
                        end = max(end, lorder[y][1])
                        y += 1
                    elif end+1+join_gap < lorder[y][0]:
                        # print("yield",eid)
                        yield Exon(echr, lorder[x][0], end, 'CH'+str(echr)+'-'+str(eid), insert=maxinsert)
                        eid += 1
                        end = lorder[y][1]
                        x, y = y, y+1
                yield Exon(echr, lorder[x][0], lorder[-1][1], 'CH'+str(echr)+'-'+str(eid), insert=maxinsert)
            l_exons = []
            echr = m.chr
        l_exons.append((m.begin, m.end))
    eid = 1
    if len(l_exons) == 1:
        yield Exon(echr, l_exons[0][0], l_exons[0][1], 'CH'+str(echr)+'-'+str(eid), insert=maxinsert)
    elif len(l_exons) != 0:
        lorder = sorted(l_exons)
        size = len(lorder)
        x = 0
        y = 1
        end = lorder[x][1]
        while(x <= size-2 and y <= size-1):
            if end+1+join_gap >= lorder[y][0]:
                end = max(end, lorder[y][1])
                y += 1
            elif end+1+join_gap < lorder[y][0]:
                # print("yield",eid)
                yield Exon(echr, lorder[x][0], end, 'CH'+str(echr)+'-'+str(eid), insert=maxinsert)
                eid += 1
                end = lorder[y][1]
                x, y = y, y+1
        yield Exon(echr, lorder[x][0], lorder[-1][1], 'CH'+str(echr)+'-'+str(eid), insert=maxinsert)


class WholeExon(object):
    # coulum is the numbeer of dNTP in a row from filein, step is usually one byte which stand for '\n'
    def __init__(self, fin, fout, column, step=1):
        self.filein = fin
        self.fileout = fout
        self.pos = 0
        self.chr = 0
        self.column = column
        self.step = step
    # write sequence of certain exon according to it's begin pos and it's length

    def exonseq(self, exon, pos, length, new_column):
        self.filein.seek(pos, 0)
        seq = ""
        # print(str(self.filein.tell()),'-->')
        # sp is set for alignment, in other word, write the same number of dNTP in every row
        line = self.filein.readline().strip()
        while(length > 0):
            length = length-len(line)
            if length >= 0:
                seq += line
            else:
                seq += line[:length]
            line = self.filein.readline().strip()
            if not line:
                raise FloatingPointError('file end..')
        if 'N' in seq:
            l = seq[0:exon.insert].rfind('N')
            r = seq[-1:-1-exon.insert:-1].rfind('N')
            s = max(l, r)+1  # s>=1
            if s:
                seq = seq[s:-s]
                exon.set_insert(exon.insert-s)
                print(exon.getexon_info())
        self.fileout.write("%s\n" % (exon.getexon_info()))
        write_column(self.fileout, seq, new_column)
        # print(str(self.filein.tell())+"\r")

    def lentopos(self, length):
        raw = length//self.column
        return raw*(self.column+self.step)+length % self.column-1
    # according exons informatin write its sequence ; the exon must be in order from 1 to 24

    def wholeexonseq(self, exons, new_column):
        new_column = self.column
        exon = next(exons)
        while(exon):
            tempchr = exon.chr
            if not tempchr == self.chr:
                print("exon end:", self.filein.tell(), "     ")
                self.chr = tempchr
                line = self.filein.readline()
                row = 0
                while(line):
                    row += 1
                    print("%d" % row, end='\r')
                    if re.match(r'>', line):
                        print("chr", self.chr, "pos",
                              self.filein.tell(), line.strip())
                        # change to next chromosome
                        self.pos = self.filein.tell()
                        row = 0
                        break
                    line = self.filein.readline()
            if not line:
                break
            # write exon title
            #print('write exon ',exon.getexon_info())
            po = self.pos + self.lentopos(exon.begin-exon.insert)
            length = exon.length+2*exon.insert
            try:
                self.exonseq(exon, po, length, new_column)
                exon = next(exons)
            except StopIteration:
                break



def get_wgs(file, write, ver):
    def write_content(file, write, stexts, etexts,smatch='==', ematch='=='):
        def text_iterator(texts):
            if isinstance(texts, str):
                while(True):
                    yield texts
            elif isinstance(texts, (tuple, list)):
                for x in texts:
                    yield x
            elif isinstance(texts, Iterator):
                return texts
        stexts = text_iterator(stexts)
        etexts = text_iterator(etexts)
        i=1
        with open(file, 'r') as filed:
            with open(write, 'w', newline='\n') as writed:
                for stext, etext in zip(stexts, etexts):
                    line = get_line_text(filed, stext, 0, smatch)
                    if not line:
                        break
                    print(line.strip()[0:40], end="\r")
                    seq=[]
                    end=0
                    while(line):
                        line = filed.readline()
                        if equal_text(line, etext, ematch):
                            filed.seek(filed.tell()-len(line)-1, 0)
                            break
                        end+=len(line.strip())
                        seq.append(line)
                    exon=Exon(i,1,end,'CH%s-1'%i)
                    i+=1
                    writed.write("%s\n" % (exon.getexon_info()))
                    writed.writelines(seq)
    word = eval(ver+"['wgs']")
    if isinstance(word, str):
        match = 're'
    elif isinstance(word, list, tuple):
        match = '=='
    # get chromosome sequence from grch38 genomic
    print("get chromosome sequence from %s..." % file)
    write_content(file, write, word, r'>.*', match, 're')
    print('\ndown. outfile: %s' % write)


def get_bed(file, write, ver):
    word = eval(ver+"['bed']")
    print("get exons annotation from %s..." % file)
    echrr=0
    with open(file, 'r') as r:
        with open(write, 'w', newline='\n') as w:
            for line in r.readlines():
                m=re.match(word,line)
                if m:
                    chrr=chr2num(m.group(1))
                    if chrr!=echrr:
                        print(line[:40],end='\r')
                        echrr=chrr
                    s=m.group(2)
                    e=m.group(3)
                    w.write('%s\t%s\t%s\n'%(chrr,s,e))
    print('\ndown. outfile: %s' % write)


def get_wes(file1, file2, write,  chip_len, join_gap, row_step, maxinsert):
    # get exon sequence ; write to exonlist.txt
    print("get exon initial list from %s , %s..." % (file1, file2))
    file1_column = get_column_row(file1, '>', 1, 're')
    with open(write, 'w', newline='\n') as w:
        w.write("#chr\tbegin\tend\tgeneid\tlength\tINSERT\n")
        with open(file1, 'r') as r1:
            with open(file2, 'r') as r2:
                exons = exons_generater(
                    r2, join_gap, chip_len, maxinsert)
                wholeexon = WholeExon(r1, w, file1_column, row_step)
                try:
                    wholeexon.wholeexonseq(exons, file1_column)
                except FloatingPointError as e:
                    print(e)
    print('\ndown. outfile: %s' % write)


def get_exons(file):
    ''' from every row in flie get exon info and yield exon object'''
    line = file.readline()
    while(line):
        m = re.match(
            r'^>chr([\d]*)\t(\d*)\t(\d*)\t([\w\.\-]*)\t\d*\t(\d*)', line)
        if m:
            # if m is exon title ,then store this exon's sequence
            seq, line = get_seq(file)
            yield Exon(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), seq, insert=int(m.group(5)))
        else:
            line = file.readline()
