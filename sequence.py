import collections
import os
import re

from filefunc import *
from setting import VER


class Edge(object):
    '''
    when begin==end==int,it a length 1 region
    when begin-1==end==int,it's a length 0 point between end and begin
    when begin==int ,end=='.',it will change to case 2
    when begin=end=='.',it is a whole region ,length -1
    when begin=='.',end==int,it will change to case 4
    when not int,change to '.'
    '''

    def __init__(self, begin, end):
        begin, end = str2int([begin, end])
        if isinstance(begin, int):
            if isinstance(end, int):
                self._begin = begin
                self._end = end
            else:
                self._end = begin
                self._begin = begin+1
        else:
            self._begin = self._end = '.'

    @property
    def begin(self):
        return self._begin

    @property
    def end(self):
        return self._end

    @property
    def length(self):
        if isinstance(self._begin, int):
            return self._end-self._begin+1
        else:
            return -1

    def check_overlap(self, edge):
        if self._begin == '.' or edge.begin == '.':
            return True
        else:
            return edge.begin <= self._end and edge.end >= self._begin


class Border(Edge):
    def __init__(self, chrr, begin, end):
        super().__init__(begin, end)
        self._chr = chrr

    @property
    def chr(self):
        return self._chr

    def get_info(self):
        info = "%s\t%s\t%s" % (self.chr, self._begin, self._end)
        return info

    @staticmethod
    def get_border(line):
        info = re.match(r'(\w+)\t(\w+)\t(\w+)', line)
        if info:
            chrr, begin, end = str2int(info.groups())
            return Border(chrr, begin, end)

    def __eq__(self, border):
        return self._begin == border.begin and self._end == self.end and self._chr == border.chr

    def __hash__(self):
        chrr = str(self._chr)
        begin = str((self._begin))
        end = str(self._end)
        return '+'.join([chrr, begin, end]).__hash__()


class Sequence(Border):
    def __init__(self, chrr=None, begin=0, end=0, seq_id=None, seq=None, flank_len=0):
        super().__init__(chrr, begin, end)
        self._id = seq_id
        self._seq = seq
        self._flank_len = flank_len

    @property
    def total_length(self):
        return 2*self._flank_len+self.length

    @property
    def length(self):
        if self._seq:
            return len(self.mid)
        else:
            return super().length

    @property
    def seq_length(self):
        return len(self._seq)

    @property
    def flank_len(self):
        return self._flank_len

    @flank_len.setter
    def flank_len(self, flank_len):
        self._flank_len = flank_len

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, seq_id):
        self._id = seq_id

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, seq):
        self._seq = seq

    def del_seq(self):
        self.seq = None

    @property
    def lflank(self):
        return self._seq[:self._flank_len]

    @property
    def rflank(self):
        return self._seq[-self._flank_len:]

    @property
    def mid(self):
        if self._flank_len:
            return self._seq[self._flank_len:-self._flank_len]
        else:
            return self._seq

    def get_part_seq(self, pos1=None, pos2=None, relative=True):
        length = self.length
        if not relative:
            pos1 = max(0, pos1-self._begin)
            pos2 = min(length, pos2-self._begin)
            if pos2 < 0:
                return ''
        if pos2 == -1:
            pos2 = self.length
        return self.mid[pos1-1:pos2]

    def show_seq(self, pos1=None, pos2=None, whole=False, relative=True):
        if relative:
            print(self.get_info(), 'absolute ', 'pos1=',
                  pos1+self._begin, 'pos2=', pos2+self._begin)
        else:
            print(self.get_info(), 'relative ', 'pos1=',
                  pos1-self._begin, 'pos2=', pos2-self._begin)
        if whole:
            print(self.get_part_seq(1, pos1, relative), '|',
                  self.get_part_seq(pos1, pos2, relative), '|', self.get_part_seq(pos2, -1, relative))
        else:
            print('|', self.get_part_seq(pos1, pos2, relative), '|')

    def get_info(self):
        '''return chr,begin,end,sequence_id,length,insert;split is tab'''
        info = ">chr%s\t%s\t%d\t%d\t%d\t%d" % (
            self._chr, self._id,  self._begin, self._end, self.length, self.flank_len)
        return info

    @staticmethod
    def get_sequence(line):
        info = re.match(
            r'>chr(\w*)\t([\w\-\.]*)\t(\d*)\t(\d*)\t(\d*)\t(\d*)', line)
        if info:
            chrr, seq_id, begin, end, length, flank_len = str2int(
                list(info.groups()))
            return Sequence(chrr, begin, end, seq_id, flank_len=flank_len)

    def search_seq(self, filed):
        if self.length:
            line = get_line_text(filed, r'>chr%s\t[\w\-\.]*\t%s\t%s' % (
                self._chr, self.begin, self.end), 0, 're', pos=0, num=1)
            if line:
                self.seq = Fasta.get_seq(filed)[0]
            else:
                self.seq = None
        else:
            self.seq = ''

    def write(self, writed, column):
        writed.write("%s\n" % (self.get_info()))
        write_small_word(writed, self.seq, column)
        writed.write('\n')

    def filtrate(self, insert_e, inner=False, char='N'):
        ranges, seqs = find_indexs(self.mid,char)
        sequences = []
        if ranges and ranges[0][1]-ranges[0][0] < self.length:
            if inner:
                print('inner filtera')
                i = 1
                for x in range(len(ranges)):
                    begin = self.begin+ranges[x][0]
                    end = self.begin+ranges[x][1]-1
                    if end-begin+1 < insert_e:
                        continue
                    sequences.append(Sequence(
                        self.chr, begin, end, self.id+'-%d' % (i), seqs[x], flank_len=0))
                    i += 1
            else:
                print('\t*ignore : inner filtera :', '\n\t',self.get_info())
                return sequences
        elif ranges:
            l = self.lflank.rfind('N')
            r = self.rflank.find('N')
            if l != -1 or r != -1:
                print('outter filtera')
                s = max(l+1, self.flank_len-l-1)
                self.flank_len -= s
                self.seq = self.seq[-s:s]
            if self.total_length > insert_e:
                sequences.append(self)
            else:
                print('\t*ignore : total length too small','\n\t',self.get_info())
        return sequences


class Fasta(object):
    '''Fasta is a layout which includes head(>name)+body(sequence), not flank_len'''
    @staticmethod
    def get_seq(filed, l=False):
        if not l:
            seqs = []
            line = filed.readline()
            while(line):
                if re.match(r'>', line):
                    break
                seqs.append(line.strip())
                line = filed.readline()
            return ''.join(seqs), line
        else:
            seqs = []
            length = 0
            line = filed.readline()
            while(line):
                if re.match(r'>', line):
                    break
                length += len(line.strip())
                seqs.append(line)
                line = filed.readline()
            return seqs, line, length

    ''' before Standardization'''

    @staticmethod
    def ini_fna(*, filex, file1, ver, **default):
        '''
        get chromosome sequence from genome, file must be in order!!
        otherwise chrx may be rename other chry
        '''
        file = filex
        write = file1
        word = VER[ver.strip()]['fna']
        print("getting chromosome sequence from %s..." % file)
        print('default: sort= true')
        chrr, line = 1, True
        with open(file, 'r') as filed:
            with open(write, 'w', newline='\n') as writed:
                line = get_line_text(filed, word, 0, 're')
                while(line):
                    print(line.strip()[0:40], end="\r")
                    seq, line, length = Fasta.get_seq(filed, True)
                    chromosome = Sequence(chrr, 1, length, 'CH%s-1' % chrr)
                    chrr += 1
                    writed.write("%s\n" % (chromosome.get_info()))
                    writed.writelines(seq)
                    if equal_text(line, word, 're'):
                        continue
                    else:
                        line = get_line_text(filed, word, 0, 're')
        print('\ndown. outfile: %s' % write)

    ''' after Standardization'''
    @staticmethod
    def sequences_iterator(filed):
        ''' from every row in flie get exon info and yield exon object'''
        line = filed.readline()
        while(line):
            sequence = Sequence.get_sequence(line)
            if sequence:
                # if sequence is exon title ,then store this exon's sequence
                seq, line = Fasta.get_seq(filed)[0:2]
                sequence.seq = seq
                yield sequence
            else:
                line = filed.readline()


class Bed(object):
    '''Bed is a layout which includes chr,begin and end'''
    '''before Standardization'''

    @staticmethod
    def ini_bed(*, filey, file2, ver, **default):
        '''
        get aimed bed from genome annotation, file must be in order!!
        otherwise chrx may be rename other chry
        '''
        file = filey
        write = file2
        join_gap = default['join_gap']
        chip_len = default['chip_len']
        word = VER[ver.strip()]['bed']
        print('dafault:', 'join_gap=', join_gap, 'chip_len=', chip_len)
        print("getting exons annotation from %s..." % file)
        ranges, echr, chrr = [], 0, 0
        with open(file, 'r') as filed:
            with open(write, 'w', newline='\n') as writed:
                for line in filed.readlines():
                    m = re.match(word, line)
                    if not m:
                        continue
                    nchr = m.group(1)
                    s = int(m.group(2))
                    e = int(m.group(3))
                    if nchr != echr:
                        print(line.strip()[:40], end='\r')
                        Bed.write_ranges(writed, ranges, chrr,
                                         join_gap, chip_len)
                        ranges = []
                        chrr += 1
                        echr = nchr
                    ranges.append((s, e))
                Bed.write_ranges(writed, ranges, chrr, join_gap, chip_len)
        print('\ndown. outfile: %s' % write)

    @staticmethod
    def write_ranges(writed, ranges, chrr, join_gap, chip_len):
        if join_gap != None:
            # join_gap is for sort
            ranges = sort_range(ranges, join_gap)
        for x in ranges:
            Bed.write_border(writed, Border(chrr, x[0], x[1]), chip_len)

    @staticmethod
    def write_border(writed, border, chip_len):
        if border.length >= chip_len:
            writed.write('%s\n' % Border.get_info(border))
        else:
            print('\t*ignore : range smaller than chip_len :','\n\t',
                  Border.get_info(border))

    @staticmethod
    def write_beds(write, beds):
        with open(write, 'w') as writed:
            for x in beds:
                writed.write('%s\n' % Border.get_info(x))

    '''after Standardization'''
    @staticmethod
    def sequences_iterator(filed, flank_len):
        beds = []
        echr = 0
        for line in filed.readlines():
            m = Border.get_border(line)
            if m:
                if m.chr != echr:
                    for i, x in enumerate(beds, start=1):
                        yield Sequence(echr, x[0], x[1], 'CH'+str(echr)+'-'+str(i), flank_len=flank_len)
                    beds = []
                    echr = m.chr
                beds.append((m.begin, m.end))
        for i, x in enumerate(beds, start=1):
            yield Sequence(echr, x[0], x[1], 'CH'+str(echr)+'-'+str(i), flank_len=flank_len)


class Sequen(object):
    '''
    Sequen is a layout extended from Fasta , which head has some information
    including flank_len.
    '''
    ''' before Standization '''
    @staticmethod
    def seq_fnabed(file1, file2, write, flank_len=0, filtrate=True, insert_e=None,**info):
        '''getting aimed sequence from chromosome sequence and bed'''
        print("getting aimed sequence list from %s , %s..." %
              (file1, file2))
        print('default:','filtrate=',filtrate,'flank_len=',flank_len,'inser_e=',insert_e)
        with open(file2, 'r') as filed2:
            sequences = Bed.sequences_iterator(filed2,  flank_len)
            try:
                Sequen.write_fasta(file1, write, sequences, filtrate, insert_e,**info)
            except StopIteration as e:
                print('\tstop.', 'seq_fnabed()', end='')
        print('\tdown. outfile: %s' % write)

    @staticmethod
    def fasta_file_info(file):
        echr = 0
        chrss = {}
        with open(file, 'r')as f:
            line = get_line_text(f, '>', 1, 're')
            column = len(line.strip())
            step = len(line)-column
            f.seek(0)
            line = f.readline()
            while(line):
                m = re.match(r'^>chr([\w\d]*)', line)
                if m:
                    if m.group(1) != echr:
                        chrss[int(m.group(1))] = {
                            'pos': f.tell(), 'sequence': Sequence.get_sequence(line)}
                        echr = m.group(1)
                line = f.readline()
        return chrss, column, step

    @staticmethod
    def write_fasta(file, write, sequences, filtrate=True, insert_e=None,**info):
        if 'chrss' in info and 'column' in info and 'step' in info:
            chrss, column, row_step=info['chrss'],info['column'],info['step']
        else:
            chrss, column, row_step = Sequen.fasta_file_info(file)
        with open(write, 'w', newline='\n') as writed:
            writed.write("#chr\texonid\tbegin\tend\tlength\tflanklen\n")
            with open(file, 'r', newline='\n') as filed:
                for x in sequences:
                    if x.chr not in chrss:
                        continue
                    if x.end > chrss[x.chr]['sequence'].end:
                        print('\t*ignore : out of range: aim end > ref end')
                        print('\taim info', x.get_info())
                        print('\tref info', chrss[x.chr]['sequence'].get_info())
                        continue
                    chr_pos = chrss[x.chr]['pos']
                    begin = x.begin-x.flank_len
                    if begin > 1:
                        pos, cur = len2pos(begin, column, row_step)
                    else:
                        pos, cur = 0, 0
                        x.flank_len = x.begin-1
                    end = x.end+x.flank_len
                    if end > chrss[x.chr]['sequence'].end:
                        x.flank_len -= (end-chrss[x.chr]['sequence'].end)
                    pos += chr_pos
                    x.seq = get_small_word(
                        filed, x.total_length, column, cur, pos=pos)[0]
                    if filtrate:
                        for y in x.filtrate(insert_e, inner=False):
                            y.write(writed, column)
                            y.del_seq()
                    else:
                        x.write(writed, column)
                    x.del_seq()

    ''' after Standardization '''
    @staticmethod
    def sequences_iterator(filed):
        ''' from every row in flie get exon info and yield exon object'''
        line = filed.readline()
        while(line):
            m = Sequence.get_sequence(line)
            if m:
                seq, line = Fasta.get_seq(filed)[:2]
                yield Sequence(m.chr, m.begin, m.end, m.id, seq, m.flank_len)
            else:
                line = filed.readline()

    @staticmethod
    def file_filtrate(file, insert_e):
        write = file+'.fil'
        column = get_column_row(file, '>', 1, 're', False)
        print('filtrate N from %s...' % file)
        print('default :','insert_e',insert_e)
        with open(file, 'r') as filed:
            with open(write, 'w') as writed:
                writed.write("#chr\texonid\tbegin\tend\tlength\tflanklen\n")
                line = filed.readline()
                while(line):
                    sequence = Sequence.get_sequence(line)
                    if sequence:
                        print(sequence.get_info(), end='\r')
                        sequence.seq, line = Fasta.get_seq(filed)[0:2]
                        sequence = sequence.filtrate(insert_e, inner=True)
                        for x in sequence:
                            x.write(writed, column)
                            x.del_seq()
                    else:
                        line = filed.readline()
        print('\ndown. outfile: %s' % write)


def get_chr_bed_num(file):
    with open(file, 'r') as f:
        l = {}
        summ = 0
        for line in f.readlines():
            m = re.match(r'^>chr([\w\d]*)\s*', line)
            if m:
                if m.group(1) not in l:
                    l[m.group(1)] = 1
                else:
                    l[m.group(1)] += 1
    summ = sum(l.values())
    l['total'] = summ
    return l


def show_chr_bed_num(bed_nums):
    for x, key in list(enumerate(bed_nums.keys())):
        if (x+1) % 5 == 0:
            print('')
        print('%s=%d' % (key, bed_nums[key]), end='\t')
    print('')


def chrs_sequence(file):
    chr_sequences = {}
    chr_poss = chrs_pos(file, 0)
    with open(file, 'r')as f:
        for x, y in chr_poss.items():
            f.seek(y, 0)
            sequence = Sequence.get_sequence(f.readline())
            chr_poss[x] = f.tell()
            chr_sequences[x] = sequence
    return chr_sequences, chr_poss


def chrs_pos(file, aim=0):
    echr = 0
    chr_poss = {}
    with open(file, 'r')as f:
        line = f.readline()
        while(line):
            m = re.match(r'^>chr([\w\d]*)', line)
            if m:
                if m.group(1) != echr:
                    if aim:
                        chr_poss[int(m.group(1))] = f.tell()
                    else:
                        chr_poss[int(m.group(1))] = f.tell()-len(line)
                    echr = m.group(1)
            line = f.readline()
    return chr_poss


def show_bed(filed, pos, chrr, info):
    print(info)
    num, pos1, pos2, whole = 1, 0, -1, False
    if len(info) >= 1:
        num = info[0]
    if len(info) >= 2:
        pos1 = info[1]
    if len(info) >= 3:
        pos2 = info[2]
    if len(info) >= 4:
        if info[3].upper() == "T":
            whole = True
    line = get_line_text(filed, r'^>', 0, 're', pos=pos, num=num)
    if line:
        bed = Sequence.get_sequence(line)
        seq = Fasta.get_seq(filed)[0]
        bed.seq = seq
        bed.show_seq(pos1, pos2, whole)


def show_pos(filed, pos, pos1, pos2, column, row_step):
    pos2 += 1
    filed.seek(pos, 0)
    print(filed.readline(), end='')
    length = pos2-pos1
    po, cur = len2pos(pos1, column, row_step)
    pos = filed.tell()+po
    print(get_small_word(filed, length, column, cur, pos=pos)[0])


def view(file):
    helps = '''
    file can be two kinds files:
    1. chromosome initial file, every length is long, use absolute pos
    2. beds initial file, every length is short, use ralative pos
    -n                              : show chromosome bed number
    -n <x>                          : shoe chromosome x bed number
    -p <x> <pos1> <pos2>            : show chx from pos1 to pos2 sequence,pos is absolute
    -b <x> <num>                    : show chx no.num bed sequence
    -b <x> <num> <pos1> <pos2>      : show chx no.num from pos1 to pos2 sequence,pos is relative
    -b <x> <num> <pis1> <pot2> -t   : show chx no.num bed sequence,mark pos1 and pos2
    '''
    column = get_column_row(file, '>', 1, 're', False)
    row_step = get_column_row(file, '>', 1, 're', True)-column
    print(helps)
    bed_nums = get_chr_bed_num(file)
    chr_pos = chrs_pos(file)
    with open(file, 'r') as f:
        line = input('>')
        while(line.strip() != 'exit'):
            line = re.findall(r'[^-\s]+', line)
            if 'help' in line:
                print(helps)
            elif'n' in line:
                if len(line) > 1:
                    if line[1] in bed_nums:
                        print('chromosome', line[1], '=', bed_nums[line[1]])
                else:
                    show_chr_bed_num(bed_nums)
            elif 'p' in line:
                if int(line[1]) in chr_pos:
                    pos = chr_pos[int(line[1])]
                    pos1 = int(line[2])
                    pos2 = int(line[3])
                    show_pos(f, pos, pos1, pos2, column, row_step)
            elif 'b' in line:
                if int(line[1]) in chr_pos:
                    pos = chr_pos[int(line[1])]
                    chrr = int(line[1])
                    info = str2int(line[2:])
                    show_bed(f, pos, chrr, info)
            line = input('>')


if __name__ == '__main__':
    file = input("please enter bedlist file you want to look for : ")
    view(file)
