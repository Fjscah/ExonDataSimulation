import collections
import json
import os
import re

import matplotlib.pyplot as plt

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

    def wes_segment(self, e_len, chip_len, insert_e, insert_d, ERROR, error_e):
        ''' get one inser's 2 reads'''
        # chip pos on exon
        length = self.mid_length
        flank_len = self._flank_len
        step = length//chip_len
        cbp = random.randint(0, step)*chip_len
        # insertion length
        if insert_d:
            insertion = round(random.normalvariate(insert_e, insert_d))
        else:
            insertion = round(insert_e)
        # chip pos on inseriton
        if insertion <= e_len:
            return '', ''
        sbp = random.randint(0, insertion-e_len)
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
            num = round(error_e*insertion)
            extra = self._seq[rbp:rbp+num]
        else:
            extra = ''
        return self._seq[lbp:rbp], extra

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

    def search_seq(self, filed, stop, chr_pos, column, row_step, MEMORY):
        if self.end > stop:
            print('region out of reference gonome -> ignore it')
            return False
        elif self.length > -1:
            self._flank_len = min(
                self._flank_len, self._begin-1, stop-self._end)
            begin = self._begin-self._flank_len
            self._seq = ''.join(get_words(
                filed, begin, self.all_length, MEMORY, column, row_step, pos=chr_pos))
            return True

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
                print('aimed region has N -> ignore it', self._chr,
                      self._begin, self._end, '            ')
            elif split:
                print('aimed region has N -> split it', self._chr,
                      self._begin, self._end, '            ')
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
            print(self.mid[:500])
        elif ranges:
            # if aimed doesn't have N
            l = self.lflank.rfind('N')
            r = self.rflank.find('N')
            # if flank has N
            if l != -1:
                sl = l+1
            else:
                sl = 0
            if r != -1:
                sr = self._flank_len-r
            else:
                sr = 0
            s = max(sl, sr)
            if s:
                print('flank region has N -> short it')
                self.seq = self.seq[s:-s]
                self._flank_len -= s
            if self.all_length > effect_len:
                sequences.append(self)
            else:
                print("filtrated region's length %s is short -> ignore it" %
                      (self.all_length-effect_len))
        return sequences


class Fasta(object):
    '''Fasta is a layout which includes head(>name)+body(sequence), not flank_len'''
    @staticmethod
    def get_before_seq(filed, MEMORY):
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
    def ini_ref(ref, inifile, n, FNA, COLUMN, MEMORY):
        '''
        get chromosomes sequence from one haploid, file must be in order
        otherwise chr x  may be rename other chr y
        '''
        reffile = ref[0]
        expression = FNA[ref[1]]
        print("initing reference sequence from %s..." % reffile)
        chrr = 1
        with open(reffile, 'r', newline='\n') as filed:
            with open(inifile, 'w', newline='\n') as writed:
                line = get_line_text(filed, expression, 0, 're')
                while(line):
                    print(line.strip()[0:40])
                    seqs, line, length = Fasta.get_before_seq(filed, MEMORY)
                    chromosome = Sequence(
                        chrr, 1, length, 'CH%s.%s' % (n, chrr))
                    chrr += 1
                    print(chromosome.get_fasta_head())
                    writed.write("%s\n" % (chromosome.get_fasta_head()))
                    write_big_words(writed, seqs, COLUMN, 0)
                    writed.write('\n')
                    if equal_text(line, expression, 're'):
                        continue
                    else:
                        #print('black line')
                        line = get_line_text(filed, expression, 0, 're')
        print('down. outfile: %s  ' % inifile)

    ''' after Standardization'''
    @staticmethod
    def iterator_fasta(file, idd=''):
        with open(file, 'r', newline='\n')as filed:
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
    def ini_exome(iniref, inireg, iniexome, effect_len, flank_len, FLI_N, INNER_N, COLUMN, MEMORY, FASTA_INFO):
        '''get aimed sequence from chromosome sequence and bed'''
        print("getting exome from : %s , %s" % (iniref, inireg))
        with open(inireg, 'r') as filed2:
            sequences = Bed.iterator_sequences(filed2, flank_len)
            try:
                Fasta.ini_exons(iniref, iniexome, sequences, effect_len,
                                FLI_N, INNER_N, COLUMN, MEMORY, FASTA_INFO)
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
        chrrs = list(infos.keys())
        chrrs.remove('column')
        chrrs.remove('step')
        chrrs = str2int(chrrs)
        chrrs.sort()
        return chrrs

    @staticmethod
    def remove_fastainfo(FASTA_INFO, file):
        log = FASTA_INFO
        if not os.path.exists(log):
            open(log, 'a').close()
        with open(log, 'r') as filed:
            try:
                fasta_infos = json.load(filed)
            except:
                fasta_infos = {}
        if file in fasta_infos:
            fasta_infos.pop(file)
            print('remove fastainfo :', file)
        else:
            print(file, 'not in fastainfo')

    @staticmethod
    def fasta_file_info(file, FASTA_INFO, *info):
        log = FASTA_INFO
        if not os.path.exists(log):
            open(log, 'a').close()
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
                    if line[0] == '>':
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
    def ini_exons(file, write, sequences, effect_len, FLI_N, INNER_N, COLUMN, MEMORY, FASTA_INFO):
        infos = Fasta.fasta_file_info(file, FASTA_INFO)
        echr = 0
        with open(write, 'w', newline='\n') as writed:
            writed.write("#exonid\tchr\tbegin\tend\tlength\tflanklen\n")
            with open(file, 'r', newline='\n') as filed:
                for x in sequences:
                    chromosome, pos = Fasta.analyse_infos(infos, x.chr)
                    if not pos:
                        continue
                    if echr != chromosome.chr:
                        echr = chromosome.chr
                        print(chromosome.get_fasta_head())
                    column = infos['column']
                    step = infos['step']
                    x.id = chromosome.id+x.id
                    if x.search_seq(filed, chromosome.end, pos, column, step, MEMORY):
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
    def ini_reg(reg, inifile, BED, BED_INFO, E_LEN, CHIP_LEN, JOIN_GAP=0):
        '''
        get aimed regions from bed file, file needn't be in order
        '''
        regfile = reg[0]
        stanfile = inifile+'.temp'
        keys = set()
        expression = BED[reg[1]]
        echr = None
        print('getting bed from file : %s' % regfile)
        with open(regfile, 'r', newline='\n')as filed:
            with open(stanfile, 'w', newline='\n') as writed:
                for line in filed.readlines():
                    m = re.match(expression, line)
                    if m:
                        chrr = str2chromosome(m.group(1))
                        if echr != chrr:
                            echr = chrr
                            keys.add(chrr)
                            print(line[:30].strip())
                        writed.write(str(chrr)+'\t'+m.group(2) +
                                     '\t'+m.group(3)+'\n')
        keys = list(keys)
        keys.sort()
        with open(BED_INFO, 'r') as filed:
            try:
                beds_info = json.load(filed)
            except:
                beds_info = {}
        with open(BED_INFO, 'w') as writed:
            beds_info[inifile] = keys
            json.dump(beds_info, writed)
        Bed.sort_reg(stanfile, inifile, BED_INFO,
                     E_LEN, CHIP_LEN, JOIN_GAP, keys)
        os.remove(stanfile)

    @staticmethod
    def get_bed_info(file, BED_INFO):
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
    def sort_reg(file, sortreg, BED_INFO, effect_len, chip_len, join_gap=0, chrrs=[]):
        '''
        sort aimed bed from filtrated regions
        '''
        print('sort reg file', file)
        if not chrrs:
            chrrs = Bed.get_bed_info(file, BED_INFO)
        #!!!!!
        with open(sortreg, 'w', newline='\n') as writed:
            pass
        for x in chrrs:
            print('sorting chromosome %s ' % x)
            ranges = Bed.get_ranges(file, x)
            ranges = merge_ranges(ranges, join_gap, effect_len)
            for i, y in enumerate(ranges):
                if y[1]-y[0] < chip_len:
                    s = (chip_len-y[1]+y[0]+1)//2
                    a = y[0]-s
                    b = y[1]+s
                    ranges[i] = (a, b)
            ranges = merge_ranges(ranges, join_gap, effect_len)
            with open(sortreg, 'a', newline='\n') as writed:
                for y in ranges:
                    s = '\t'.join([str(x), str(y[0]), str(y[1])])+'\n'
                    writed.write(s)
        print('down. outfile :', sortreg)

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
    def iterator_sequences(filed, flank_len):
        i = 1  # use it for count every region on one chromosome
        echr = 0  # use it for count every chromesome
        for line in filed.readlines():
            m = Region.self_bed(line)
            if m:
                if echr != m.chr:
                    i = 1
                    echr = m.chr
                yield Sequence(echr, m.begin, m.end, '.'+str(i), flank_len)
                i += 1


class Quality(object):
    '''before ini'''
    @staticmethod
    def ini_qph(file, inifile, SEED):
        print('getting phred frequencies')
        frequencies = {}
        readlens = {}
        with open(file, 'r') as f:
            i = 0
            line = f.readline()
            while(True and line):
                if re.match('^@', line):
                    break
                line = f.readline()
            while(line):
                i += 1
                if i % 4 == 0:
                    row_fastq = line.strip()
                    l = len(row_fastq)
                    if l in readlens:
                        readlens[l] += 1
                    else:
                        readlens[l] = 1
                    for n, char in enumerate(row_fastq, 1):
                        if ('pos%d_frequencies' % n) in frequencies:
                            if char in frequencies['pos%d_frequencies' % n]:
                                frequencies['pos%d_frequencies' % n][char] += 1
                            else:
                                frequencies['pos%d_frequencies' % n][char] = 1
                        else:
                            frequencies['pos%d_frequencies' % n] = {}
                            frequencies['pos%d_frequencies' % n][char] = 1
                if i % 40000 == 0:
                    print(i, end='\r')
                if i > SEED:
                    break
                line = f.readline()
        sums = []
        for x in range(len(frequencies)):
            sums.append(sum(frequencies['pos%d_frequencies' % (x+1)].values()))
        with open(inifile, "w") as f:
            json.dump([frequencies, sums, readlens], f)
            #json.dump(frequencies, f)
    '''after ini'''
    @staticmethod
    def get_qph(file):
        print('initial qphred from %s...' % file, end='')
        asc = 64
        with open(file, "r") as filed:
            frequencies, sums, readlens = json.load(filed)
        for x in frequencies['pos1_frequencies'].keys():
            if x in '1234456789+-*':
                asc = 33
                break
        print('down.')
        return frequencies, sums, readlens, asc

    @staticmethod
    def get_avg_readlens(readlens):
        summ = 0
        for x, y in readlens.items():
            summ += int(x)*int(y)
        return round(summ/sum(readlens.values()))

    @staticmethod
    def get_qphred_reads(frequencies, sums, readlens, num):
        qphreds = []
        summ = sum(readlens.values())
        for n in range(num):
            s = ''
            l = int(random_weight_choice(readlens, summ))
            for x in range(1, l+1):
                char = random_weight_choice(
                    frequencies['pos%d_frequencies' % x], sums[x-1])
                s = s+char
            qphreds.append(s)
        return qphreds


class Depth(object):
    @staticmethod
    def re_group(file, depth, segment, chip_len):
        f = open(file, 'r', newline='\n')
        sreg = []
        maxs = []
        x = 1
        line = f.readline()
        while(line):
            chrr, b, dep = line.strip().split()
            chrr = chrr.split('.')[-1]
            b, dep = int(b), int(dep)
            sreg.append(dep)
            if x > 10000000:
                break
            if x % 100 == 0:
                maxs.append(max(sreg))
                sreg = []
            line = f.readline()
            x += 1
        f.close()
        maxs.sort(reverse=True)
        dmax = []
        for x in maxs:
            if x-depth > -depth//2 and x-depth < depth:
                dmax.append(x)
        depth = sum(dmax)/len(dmax)
        maxs = []
        dmax = []
        group = depth*chip_len/segment
        print('depth=', depth, 'group=', group)
        return group

    @staticmethod
    def dep2bed(file, depth, segment, chip_len, effect_len, CD, join_gap, BED_INFO, slide_step):
        '''
        it use method -- chip length's total depth to identify sequence regions
        '''
        #file ='reference.depth'
        write2 = CD+file.split('/')[-1]+'.bed'
        group = Depth.re_group(file, depth, segment, chip_len)

        def search_dis(distances, dis, bsum):
            if dis >= distances[-1][0]:
                cc = distances[-1][1]
                v = len(distances)
                ahead = bsum//(v*chip_len)
                return max(round(cc-ahead), 0)
            for v, k in enumerate(distances):
                if k[0] > dis:
                    break
            cc = distances[v-1][1]
            if v > 1:
                ahead = bsum//(v*chip_len-chip_len)
            else:
                ahead = 0
            return max(round(cc-ahead), 0)  # real add
        echr = 0
        # one depth need legnth =l=a
        # total bases= l//c*r*r
        # one catch length =l//c*r=b
        # ax=by->y=ax/b
        n = 2*segment+chip_len-2*effect_len
        origins, midsum, distances = onechip(
            segment, group, chip_len, effect_len, slide_step)
        midlsum = midsum//2
        print('midsum=', midsum, 'midltsum=', midsum, 'origin', origins[0])
        origins = origins[:5]
        # print('chiplen',len(origins[s-e:e-s]))
        f = open(file, 'r', newline='\n')
        w = open(write2, 'w', newline='\n')
        ts = [0, ]*slide_step
        tsum = sum(ts)
        i = 0
        flag = False
        exonranges = []
        direction = 0
        t = 0
        rag = 0
        bb = 0
        regsums = []
        mflag = False
        mrag = 0
        nflag = False
        line = f.readline()
        mregsums = []
        bsum=0
        while(line):
            i += 1
            chrr, b, d = line.strip().split()
            chrr = chrr.split('.')[-1]
            b, d = int(b), int(d)
            if chrr != echr:
                print('repositon chromsome', chrr, '...')
                echr = chrr
                flag = False
                i = b
                ts = [0, ]*slide_step
                tsum = sum(ts)
                flag =mflag=nflag= False
                exonranges = []
                direction = 0
                t =rag=bb=mrag= 0
                regsums = []
                bsum=0
            if b-i >= slide_step and not flag:
                # mid zero
                i = b
                ts = [0, ]*slide_step
                tsum = 0
            elif b-i >= slide_step and flag:
                t = i
                i = b
                ts = [0, ]*slide_step
                tsum = 0
            elif b-i > 0:
                for count in range(b-i):
                    v = ts.pop(0)
                    ts.append(0)
                    direction-v+d
                i = b
            v = ts.pop(0)
            ts.append(d)
            tsum = tsum-v+d
            direction = direction-v+d
            if exonranges and exonranges[0][1] < b-segment+effect_len-slide_step:
                exonranges.pop(0)
            if i % slide_step == 0 or t:
                esum = 0
                for x, y in exonranges:
                    start = b-slide_step-segment+effect_len
                    x = max(round(x-start), 0)
                    y = round(y-start)
                    esum += sum(origins[x:y])
                if flag:
                    # if reads too small or zero
                    start = b-slide_step-segment+effect_len
                    tcc = search_dis(distances, b-slide_step, bsum)
                    tesum = esum + \
                        sum(origins[bb+tcc-start:b-slide_step-start])
                    if t or tsum-tesum < midsum+slide_step:
                        if t:
                            b = t
                            mm = tsum
                            tsum = 0
                        if tsum-tesum < 0:
                            ee = b-slide_step
                        else:
                            ee = b-slide_step//2
                        cc = search_dis(distances, ee-bb, bsum)
                        if ee-bb-cc > effect_len:
                            string = '\t'.join([chrr, str(bb+cc), str(ee)])
                            w.write(string+'\n')
                        flag = mflag = False
                        exonranges.append((bb+cc, ee))
                        regsums = mregsums = []
                        if t:
                            t = 0
                            #b = i
                            tsum = mm
                    elif mflag:
                        mregsums.append(tsum)
                        mrag += 1
                        # length_m=len(mregsums)
                        if tsum >= max_reg+slide_step and mrag*slide_step < chip_len:
                            mflag = False
                            regsums += mregsums
                        elif direction > slide_step and mrag*slide_step >= chip_len:

                            firstregs = first_derivate(mregsums)
                            l = firstregs.index(min(firstregs))+rag
                            #l=mrag//2+rag
                            '''
                            l=mrag//3+rag
                            r=mrag//3*2+rag
                            '''
                            ee = bb+l*slide_step
                            cc = search_dis(distances, ee-bb, bsum)
                            string = '\t'.join([chrr, str(bb+cc), str(ee)])
                            exonranges.append((bb+cc, ee))
                            w.write(string+'\n')
                            flag = mflag = nflag = False
                            mregsums = []
                            mrag = 0
                    # if reads reduce , enter overlap
                    elif direction < -slide_step and len(regsums)*slide_step >= chip_len:
                        # regsums.sort()
                        max_reg = max(regsums)
                        mflag = True
                        mregsums.append(tsum)
                        mrag = 1
                        rag = len(regsums)
                    else:
                        regsums.append(tsum)
                elif tsum >= midsum+esum+slide_step and ts[0] >= origins[0] and direction > slide_step:
                    # reads enough and left point not zero
                    if exonranges and  b-slide_step-exonranges[-1][1]<chip_len:
                        pass
                    else:
                        flag = True
                        regsums = [tsum]
                        bsum = tsum-esum
                        bb = b-slide_step
                direction = 0
            line = f.readline()
        f.close()
        w.close()
        write3 = CD+file.split('/')[-1]+'.sortbed'
        Bed.sort_reg(write2, write3, BED_INFO, effect_len, chip_len, join_gap)
        print('down. outfile', write2, write3)


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


def view_depth(depfile, reg):
    help = '''
    -d <x> <pos> <pos> :show depth,out png
    '''
    line = input('>>')
    while(line.strip() != 'exit'):
        line = line.split()
        if '-help' in line:
            print(help)
        elif '-d' in line and len(line) >= 4:
            chrr, begin, end = str2int(line[1:4])
            if '-t' in line:
                mark = 1
            else:
                mark = 0
            xs = [v for v in range(begin, end+1)]
            ys = [0, ]*(end-begin+1)
            f1 = open(depfile, 'r', newline='\n')
            for ss in f1.readlines():
                if re.match('.*.%s\t' % chrr, ss):
                    idd, pos, depth = str2int(ss.split())
                    if pos > end:
                        break
                    elif pos >= begin and pos <= end:
                        ys[pos-begin] = depth
            f1.close()
            lengtht = max((begin-end)//10000, 20)
            print('lengtht', lengtht)
            plt.figure(figsize=(lengtht, 5))
            plt.style.use('seaborn-white')
            plt.ylabel('depth')
            plt.title('the depth between %s and %s on chromosome %s' %
                      (begin, end, chrr))
            plt.xlim(begin, end)
            if ys:
                top = max(ys)
            else:
                top = 0
            zs = [0, ]*(end-begin+1)
            ws = []
            ranges = Bed.get_ranges(reg, chrr)
            for x in ranges:
                if x[0] > end:
                    break
                if x[0] >= begin:
                    zs[x[0]-begin] = top
                    ws.append(x[0]-begin)
                    if x[1] <= end:
                        for y in range(x[1]-x[0]):
                            zs[x[0]-begin+y+1] = top
                            ws.append(x[0]-begin+y+1)
                    else:
                        for y in range(end-x[0]):
                            zs[x[0]-begin+y+1] = top
                            ws.append(x[0]-begin+y+1)
                        break
                elif x[1] >= begin:
                    for y in range(x[1]-begin):
                        zs[y] = top
                        ws.append(y)
            esum = 0
            for x in ws:
                esum += ys[x]
            if esum:
                esum = esum/len(ws)
            asum = sum(ys)
            asum = asum/(end-begin+1)
            print('exon degth=', esum, 'region depth=', asum)
            plt.fill_between(xs, ys, color='pink', alpha=0.8)
            if mark:
                plt.fill_between(xs, zs, alpha=0.6, facecolor='lightblue')
            # plt.plot(xs,ys,color='grey',linewidth=0.8,alpha=1)
            # plt.plot(xs,ys,color='red',alpha=0.4)
            depfilet = depfile.split('/')[-1]
            regt = reg.split('/')[-1]
            if mark:
                plt.savefig('%s-%s-chr%s-%s-%s.png' % (begin, end,
                                                       chrr, depfilet, regt), bbox_inches='tight')
                print('outfile', '%s-%s-chr%s-%s-%s.png' %
                      (begin, end, chrr, depfilet, regt))
            else:
                plt.savefig('f%s-%s-chr%s-%s.png' %
                            (begin, end, chrr, depfilet), bbox_inches='tight')
                print('outfile', 'f%s-%s-chr%s-%s-%s.png' %
                      (begin, end, chrr, depfilet, regt))
        else:
            print('cannot idendify instruction')
        line = input('>>')


def view_seq(file, FASTA_INFO, MEMORY):
    helps = '''
    file need be reference file
    -p <x> <pos1> <pos2>            : show chx from pos1 to pos2 sequence
    -help
    '''
    infos = Fasta.fasta_file_info(file, FASTA_INFO)
    # print(infos)
    column = infos['column']
    step = infos['step']
    with open(file, 'r', newline='\n') as filed:
        line = input('>>')
        while(line.strip() != 'exit'):
            line = line.split()
            if '-help' in line:
                print(helps)
            elif '-p' in line and len(line) >= 4:
                chrr, begin, end = str2int(line[1:4])
                pos = Fasta.analyse_infos(infos, chrr)[1]
                texts = ''.join(get_words(filed, begin, end -
                                          begin+1, MEMORY, column, step, pos=pos))
                texts = tidy_small_word(texts, 100)
                for text in texts:
                    print(text, end='')
                print()
            else:
                print('cannot idendify instruction')
                break
            line = input('>>')


def view(FASTA_INFO, MEMORY):
    helps = '''
    file need be reference file
    -p reffile            : show chx from pos1 to pos2 sequence
    -d depfile regfile
    -help
    '''
    print(helps)
    line = input('>')
    while(line.strip() != 'exit'):
        line = line.split()
        if '-help' in line:
            print(helps)
        elif '-p' in line and len(line) >= 2:
            file = line[1]
            view_seq(file, FASTA_INFO, MEMORY)
        elif '-d' in line and len(line) >= 3:
            depfile, regfile = line[1:3]
            view_depth(depfile, regfile)
        else:
            print('cannot idendify instruction')
        line = input('>')
