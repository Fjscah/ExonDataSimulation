import collections
import json
import os
import re
from basic import *


class Line(object):
    '''
    when begin-end+1==1,it's  1 size region 
    when begin-end+1==0,it's  0 size point , point is right of end
    when begin=='.' xor end=='.',change it to case 2. end =begin-1 or begin=end+1
    when begin==end=='.', it's the whole region .assign begin=0,end=-2,begin-end+1=-1,
    eg. 1,2,3,4 . 
    (2,.)=(2,1) is ',' between 1 and 2
    (.,2)=(3,2) is ',' between 2 and 3
    (.,.)=(start,stop)=(0,-2) is whole region , but not know the length
    (2,2) is 2 not include its both ','
    (2,3) is 2,3 include one ',', 2,3

    '''

    def __init__(self, begin, end):
        self.set_line(begin, end)

    @property
    def begin(self):
        return self._begin

    @property
    def end(self):
        return self._end

    def set_line(self, begin, end):
        begin, end = str2int([begin, end])
        if isinstance(begin, int):
            self._begin = begin
            if isinstance(end, int):
                self._end = end
            else:
                self._end = begin-1
        elif isinstance(end, int):
            self._end = end
            self._begin = end+1
        else:
            self._end = -2
            self._begin = 0

    @begin.setter
    def begin(self, begin):
        try:
            begin = int(begin)
            if self.end-begin+1 >= -1 and begin >= 0:
                self._begin = begin
        except:
            self._begin = self.end+1

    @end.setter
    def end(self, end):
        try:
            end = int(end)
            if end-self.begin+1 >= -1 and end >= -2:
                self.end = end
        except:
            self._end = self._begin-1

    @property
    def length(self):
        return self._end-self._begin+1

    def check_overlap(self, edge):
        if self._end == edge.end and self._begin == edge.begin:
            return True
        else:
            return edge.begin <= self._end and edge.end >= self._begin


class Region(Line):
    def __init__(self, chrr, begin, end, idd=''):
        super().__init__(begin, end)
        try:
            self._chr = int(chrr)
        except:
            self._chr = chrr
        self.id = idd

    @property
    def chr(self):
        return self._chr

    @chr.setter
    def chr(self, chrr):
        try:
            self._chr = int(chrr)
        except:
            self._chr = chrr

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, idd):
        self._id = idd

    def get_fasta_head(self):
        '''return fasta head information'''
        info = ">%s\t%s\t%s\t%s\t%s" % (
            self.id, self.chr, self._begin, self._end, self.length)
        return info

    def get_bed(self):
        info = "%s\t%s\t%s" % (self.chr, self._begin, self._end)
        return info

    @staticmethod
    def self_bed(line):
        info = re.match(r'(\w+)\t(\w+)\t(\w+)', line)
        if info:
            chrr, begin, end = str2int(info.groups())
            return Region(chrr, begin, end)

    def __eq__(self, region):
        ''' judge both whether is the same'''
        return self._begin == region.begin and self._end == self.end and self._chr == region.chr

    def __hash__(self):
        ''' it is for set()'''
        chrr = str(self._chr)
        begin = str((self._begin))
        end = str(self._end)
        return '+'.join([chrr, begin, end]).__hash__()


class Sequence(Region):
    def __init__(self, chrr, begin, end, idd='',  flank_len=0, seq=''):
        super().__init__(chrr, begin, end, idd)
        self._seq = seq
        self._flank_len = int(flank_len)

    @property
    def all_length(self):
        return 2*self._flank_len+self.length

    @property
    def flank_len(self):
        return self._flank_len

    @property
    def mid_length(self):
        return self.seq_length-2*self._flank_len

    @property
    def seq_length(self):
        return len(self.seq)

    @property
    @flank_len.setter
    def flank_len(self, flank_len):
        flank_len = int(flank_len)
        if flank_len >= 0:
            self._flank_len = flank_len

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, seq):
        self._seq = seq

    def del_seq(self):
        self.seq = ''

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

    @staticmethod
    def self_fasta_head(line):
        info = re.match(
            r'>([\w\-\.]+)\t(\d*)\t(\d*)\t(\d*)\t(\d*)\t(\d*)', line)
        if info:
            idd, chrr,  begin, end, length, flank_len = str2int(
                list(info.groups()))
            return Sequence(chrr, begin, end, idd, flank_len)

    @staticmethod
    def self_fasta(strings):
        info = re.match(
            r'>([\w\-\.]+)\t(\d+)\t(\d+)\t(\d+)\t(\d*)\t(\d*)', strings)
        if info:
            idd, chrr,  begin, end, length, flank_len = str2int(
                list(info.groups()))
            seq = ''.join(strings[1:])
            return Sequence(chrr, begin, end, idd, flank_len, seq)


    def get_fasta_head(self):
        '''return chr,begin,end,sequence_id,length,insert;split is tab'''
        info = ">%s\t%s\t%d\t%d\t%d\t%d" % (
            self._id, self._chr,  self._begin, self._end, self.length, self._flank_len)
        return info

    def get_fasta(self, column):
        return [self.self_fasta_head]+get_words_text(self.seq, column)
 
    def wes_segment(self, readlen, chip_len=CHIP_LEN, insert_e=SEGMENT_E, insert_d=SEGMENT_D):
        ''' get one inser's 2 reads'''
        # chip pos on exon
        length = self.mid_length
        flank_len = self._flank_len
        step = length//chip_len
        cbp = random.randint(0, step)*chip_len
        while(True):
            # insertion length
            insertion = round(random.normalvariate(insert_e, insert_d))
            if insertion < readlen:
                continue
            # chip pos on inseriton
            sbp = random.randint(0, insertion-chip_len)
            # insertion pos on exon
            lbp = flank_len+cbp-sbp
            if lbp < 0:
                insertion = lbp+insertion
                lbp = 0
            rbp = lbp+insertion
            if rbp > length+2*flank_len:
                insertion = length+2*flank_len-rbp+insertion
                rbp = length+2*flank_len
            if ERROR:
                num=round(ERROR_E*insertion)
                extra=self._seq[rbp:rbp+num]
            else :
                extra=''
            if insertion >= readlen:
                break
        return self._seq[lbp:rbp],extra

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
            print(self.get_fasta_head(), 'absolute ', 'pos1=',
                  pos1+self._begin, 'pos2=', pos2+self._begin)
        else:
            print(self.get_fasta_head(), 'relative ', 'pos1=',
                  pos1-self._begin, 'pos2=', pos2-self._begin)
        if whole:
            print(self.get_part_seq(1, pos1, relative), '|',
                  self.get_part_seq(pos1, pos2, relative), '|', self.get_part_seq(pos2, -1, relative))
        else:
            print('|', self.get_part_seq(pos1, pos2, relative), '|')

    def search_seq(self, filed, stop, chr_pos, column, row_step):
        if self.end > stop:
            print('region out of reference gonome -> ignore it')
        elif self.length > -1:
            self._flank_len = min(
                self._flank_len, self._begin-1, stop-self._end)
            begin = self._begin-self._flank_len
            self._seq = ''.join(get_words(
                filed, begin, self.all_length, MEMORY, column, row_step, pos=chr_pos))

    def write_fasta(self, writed, column):
        writed.write("%s\t\n" % (self.get_fasta_head()))
        write_small_word(writed, self.seq, column)
        writed.write('\n')

    def filtrate(self, effect_len, split=False, char='N'):
        ranges, seqs = not_indexs(self.mid, char)
        sequences = []
        if ranges and ranges[0][1]-ranges[0][0]+1 < self.length:
            # if aimed sequence has N
            if not split:
                print('aimed region has N -> ignore it')
            elif split:
                print('aimed region has N -> split it')
                i = 1
                for x in range(len(ranges)):
                    begin = self.begin+ranges[x][0]
                    end = self.begin+ranges[x][1]
                    if end-begin+1 < effect_len:
                        print('splited region is short -> ingore it')
                        continue
                    idd = self.id+'.%d' % (i)
                    sequences.append(
                        Sequence(self.chr, begin, end, idd, 0, seqs[x]))
                    i += 1
            print(self.mid)
        elif ranges:
            # if aimed doesn't have N
            l = self.lflank.rfind('N')
            r = self.rflank.find('N')
            # if flank has N
            if l != -1:
                sl=l+1
            else:
                sl=0
            if r != -1:
                sr=self._flank_len-r
            else:
                sr=0
            s = max(sl,sr)
            if s:
                print('flank region has N -> short it')
                self.seq = self.seq[s:-s]
                self._flank_len -= s
            if self.all_length > effect_len:
                sequences.append(self)
            else:
                print("filtrated region's length %s is short -> ignore it" %(self.all_length-effect_len))
        return sequences


class Fasta(object):
    '''Fasta is a layout which includes head(>name)+body(sequence), not flank_len'''
    @staticmethod
    def get_before_seq(filed):
        '''get fasta sequence body, l is whether to record length'''
        seqs = []
        s, l = [], MEMORY
        length = 0
        line = filed.readline().strip()
        while(line):
            if re.match(r'>', line):
                break
            t = len(line)
            length += t
            s.append(line[:l])
            l -= t
            if l < 0:
                s = ''.join(s)
                seqs.append(s)
                s = [line[l:]]
                l = MEMORY-len(s[0])
            line = filed.readline().strip()
        if s:
            s = ''.join(s)
            seqs.append(s)
            s = []
        return seqs, line, length

    ''' before Standardization'''
    @staticmethod
    def get_small_seq(filed):
        seq = []
        line = filed.readline()
        while(line):
            if re.match(r'>', line):
                break
            seq.append(line.strip())
            line = filed.readline()
        return ''.join(seq), line

    @staticmethod
    def ini_ref(ref, inifile, n):
        '''
        get chromosomes sequence from one haploid, file must be in order
        otherwise chr x  may be rename other chr y
        '''
        reffile=ref[0]
        expression = FNA[ref[1]]
        print(re.match(expression,'>NC00'))
        print("initing reference sequence from %s..." % reffile)
        chrr = 1
        with open(reffile, 'r') as filed:
            with open(inifile, 'w', newline='\n') as writed:
                line = get_line_text(filed, expression, 0, 're')
                while(line):
                    print(line.strip()[0:40], end="\r")
                    seqs, line, length = Fasta.get_before_seq(filed)
                    chromosome = Sequence(
                        chrr, 1, length, 'CH%s.%s' % (n, chrr))
                    chrr += 1
                    writed.write("%s\n" % (chromosome.get_fasta_head()))
                    write_big_words(writed, seqs, COLUMN, 0)
                    writed.write('\n')
                    if equal_text(line, expression, 're'):
                        continue
                    else:
                        line = get_line_text(filed, expression, 0, 're')
        print('down. outfile: %s  ' % inifile)

    ''' after Standardization'''
    @staticmethod
    def iterator_fasta(file, idd=''):
        with open(file,'r',newline='\n')as filed:
            line = filed.readline()
            while(line):
                if re.match('>'+idd, line):
                    sequence = Sequence.self_fasta(line)
                    seq, line = Fasta.get_small_seq(filed)
                    sequence.seq = seq
                    yield sequence
                else:
                    line = filed.readline()

    @staticmethod
    def ini_exome(iniref, inireg, iniexome,effect_len=E_LEN):
        '''get aimed sequence from chromosome sequence and bed'''
        print("getting exome from : %s , %s" % (iniref, inireg))
        with open(inireg, 'r') as filed2:
            sequences = Bed.iterator_sequences(filed2)
            try:
                Fasta.ini_exons(iniref, iniexome, sequences,effect_len)
            except StopIteration:
                print('\tstop.', end='')
        print('down. outfile: %s ' % iniexome)

    @staticmethod
    def analyse_infos(infos, chrr):
        chrr = str(chrr)
        if chrr not in infos:
            print('chromosome', chrr, 'not in reference -> ignore it')
            return None, None
        chromosome = Sequence.self_fasta_head(infos[str(chrr)]['chromosome'])
        pos = infos[str(chrr)]['pos']
        return chromosome, pos
    @staticmethod
    def chrs_info(infos):
        chrrs=list(infos.keys())
        chrrs.remove('column')
        chrrs.remove('step')
        chrrs=str2int(chrrs)
        chrrs.sort()
        return chrrs
    @staticmethod
    def fasta_file_info(file, *info):
        log = FASTA_INFO
        print('get fasta info from :', file, end='...')
        with open(log, "r") as filed:
            try:
                fastas_info = json.load(filed)
            except:
                fastas_info = {}
        if file in fastas_info:
            infos = fastas_info[file]
        else:
            echr = 0
            infos = {}
            with open(file, 'r')as f:
                line = get_line_text(f, '>', 1, 're')
                column = len(line.strip())
                step = len(line)-column
                f.seek(0)
                line = f.readline()
                while(line):
                    if line[0]=='>':
                        m = Sequence.self_fasta_head(line)
                        if m.chr != echr:
                            infos[str(m.chr)] = {
                                'pos': f.tell(), 'chromosome': line}
                            echr = m.chr
                    line = f.readline()
            infos['column'] = column
            infos['step'] = step
            fastas_info[file] = infos
            with open(log, "w") as filed:
                json.dump(fastas_info, filed)
        sinfos = []
        for x in info:
            sinfos.append(infos[x])
        print('down.')
        return infos

    @staticmethod
    def ini_exons(file, write, sequences, effect_len=E_LEN):
        infos = Fasta.fasta_file_info(file)
        echr=0
        with open(write, 'w', newline='\n') as writed:
            writed.write("#exonid\tchr\tbegin\tend\tlength\tflanklen\n")
            with open(file, 'r', newline='\n') as filed:
                for x in sequences:
                    chromosome, pos = Fasta.analyse_infos(infos, x.chr)
                    if not pos:
                        continue
                    if echr!=chromosome.chr:
                        echr=chromosome.chr
                        print(chromosome.get_fasta_head(),end='\r')
                    column = infos['column']
                    step = infos['step']
                    x.id = chromosome.id+x.id
                    x.search_seq(filed, chromosome.end, pos, column, step)
                    if FLI_N:
                        for y in x.filtrate(effect_len, INNER_N):
                            y.write_fasta(writed, COLUMN)
                            y.del_seq()
                    else:
                        x.write_fasta(writed, column)
                    x.del_seq()


class Bed(object):
    '''Bed is a layout which includes chr,begin and end'''
    '''before Standardization'''
    @staticmethod
    def ini_reg(reg, inifile):
        '''
        get aimed regions from bed file, file needn't be in order
        '''
        regfile=reg[0]
        stanfile = inifile+'.temp'
        keys = set()
        expression=BED[reg[1]]
        echr=None
        print('getting bed from file : %s'% regfile)
        with open(regfile,'r',newline='\n')as filed:
            with open(stanfile,'w',newline='\n') as writed:
                for line in filed.readlines():
                    m = re.match(expression, line)
                    if m:
                        chrr = str2chromosome(m.group(1))
                        if echr!=chrr:
                            echr=chrr
                            keys.add(chrr)
                            print(line[:30],end='\r')
                        writed.write(str(chrr)+'\t'+m.group(2)+'\t'+m.group(3)+'\n')
        keys=list(keys)
        keys.sort()
        with open(BED_INFO,'r') as filed:
            try:
                beds_info = json.load(filed)
            except:
                beds_info={}
        with open(BED_INFO,'w') as writed:
            beds_info[inifile]=keys
            json.dump(beds_info, writed)
        Bed.sort_reg(stanfile, inifile, keys)
        os.remove(stanfile)
    @staticmethod
    def get_bed_info(file):
        log = BED_INFO
        print('get bed info from :', file)
        with open(log, "r") as filed:
            try:
                beds_info = json.load(filed)
            except:
                beds_info = {}
        if (file in beds_info) and (beds_info[file]):
            keys = beds_info[file]
        else:
            keys = set()
            echr = 0
            with open(file, 'r', newline='\n')as filed:
                for line in filed.readlines():
                    m = re.match(r'(\d+)', line)
                    if m:
                        chrr = int(m.group(1))
                        if echr != chrr:
                            echr = chrr
                            keys.add(chrr)
            keys = list(keys)
            keys.sort()
            beds_info[file] = keys
            with open(log, "w") as filed:
                json.dump(beds_info, filed)
        return keys

    @staticmethod
    def sort_reg(file, sortreg, chrrs=[]):
        '''
        sort aimed bed from filtrated regions
        '''
        if not chrrs:
            chrrs = Bed.get_bed_info(sortreg)
        #!!!!!
        with open(sortreg, 'w', newline='\n') as writed:
            pass
        for x in chrrs:
            print('sorting chromosome %s ' % x, end='\r')
            ranges = Bed.get_ranges(file, x)
            if isinstance(JOIN_GAP, int):  # if need merge
                ranges = merge_ranges(ranges, JOIN_GAP)
            with open(sortreg, 'a', newline='\n') as writed:
                for y in ranges:
                    if y[1]-y[0] > E_LEN:  # filtrate too shor region
                        s = '\t'.join([str(x), str(y[0]), str(y[1])])+'\n'
                        writed.write(s)
                    else:
                        print("region's length %s is too short -> ignore it"% (y[1]-y[0]))

    @staticmethod
    def get_ranges(file, chrr):
        ranges = []
        with open(file, 'r', newline='\n') as filed:
            for line in filed.readlines():
                m = re.match('%s\t' % chrr, line)
                if m:
                    m = line.split()
                    s = int(m[1])
                    e = int(m[2])
                    ranges.append((s, e))
        return ranges

    '''after Standardization'''
    @staticmethod
    def iterator_sequences(filed):
        i = 1  # use it for count every region on one chromosome
        echr = 0  # use it for count every chromesome
        for line in filed.readlines():
            m = Region.self_bed(line)
            if m:
                if echr != m.chr:
                    i = 1
                    echr = m.chr
                yield Sequence(echr, m.begin, m.end, '.'+str(i), flank_len=FLANK_LEN)
                i += 1


class Quality(object):
    '''before ini'''
    @staticmethod
    def ini_qph(file, inifile):
        print('getting phred frequencies')
        frequencies = {}
        readlen = get_column_row(file, r'+', 1, '==')
        for x in range(1, readlen+1):
            frequencies['pos%d_frequencies' % x] = {}
        with open(file, 'r') as f:
            i = 0
            line = f.readline()
            while(line):
                i += 1
                if i % 4 == 0:
                    row_fastq = line.strip()
                    x = 1
                    for char in row_fastq:
                        if char in frequencies['pos%d_frequencies' % x]:
                            frequencies['pos%d_frequencies' % x][char] += 1
                        else:
                            frequencies['pos%d_frequencies' % x][char] = 1
                        x = x+1
                if i % 40000 == 0:
                    print(i, end='\r')
                if i > SEED:
                    break
                line = f.readline()
        with open(inifile, "w") as f:
            json.dump(frequencies, f)
    '''after ini'''
    @staticmethod
    def get_qph(file):
        print('initial qphred from %s...' % file, end='')
        asc = 64
        with open(file, "r") as filed:
            frequencies = json.load(filed)
        for x in frequencies['pos1_frequencies'].keys():
            if x in '1234456789+-*':
                asc = 33
                break
        print('down.')
        return frequencies, asc

    @staticmethod
    def get_qphred_reads(frequencies, num):
        qphreds = []
        summ = sum(frequencies['pos1_frequencies'].values())
        l = len(frequencies)
        for n in range(num):
            s = ''
            for x in range(1, l+1):
                char = random_weight_choice(
                    frequencies['pos%d_frequencies' % x], summ)
                s = s+char
            qphreds.append(s)
        return qphreds


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



def view(file):
    helps = '''
    file need be reference file
    -p <x> <pos1> <pos2>            : show chx from pos1 to pos2 sequence
    -help
    '''
    infos=Fasta.fasta_file_info(file)
    column = infos['column']
    step = infos['step']
    print(helps)
    #bed_nums = get_chr_bed_num(file)
    with open(file, 'r') as filed:
        line = input('>')
        while(line.strip() != 'exit'):
            line = line.split()
            if '-help' in line:
                print(helps)
            elif '-p' in line and len(line)>=4:
                chrr,begin,end=str2int(line[1:4])
                length=end-begin+1
                pos=Fasta.analyse_infos(infos,chrr)[1]
                texts=''.join(get_words(filed, begin, length, MEMORY, column, step,pos=pos))
                texts=tidy_small_word(texts,100)
                for text in texts:
                    print(text,end='')
                print()
            else:
                print('cannot idendify instruction' )
            line = input('>')


if __name__ == '__main__':
    file = input("please enter bedlist file you want to look for : ")
    view(file)