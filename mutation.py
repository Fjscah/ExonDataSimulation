from copy import deepcopy
from functools import reduce
from basic import *
from sequence import *
import sys
COMPLEMENT = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}


def read_complement(read):
    l = []
    for char in read:
        l.append(COMPLEMENT[char.upper()])
    return ''.join(l)


class Indel(Line):
    def __init__(self, begin, end, alt):
        super().__init__(begin, end)
        self._alt = alt.strip('.')

    @property
    def alt(self):
        return self._alt

    @alt.setter
    def alt(self, alt):
        self._alt = alt

    @staticmethod
    def get_indel_formula(formula):
        alt = ''.join(re.findall(r'[A-Za-z]+', formula))
        edge = str2int(re.findall(r'[\d\.]+', formula))
        if len(edge) == 1:
            if '(' in formula:
                return Indel(edge[0], edge[0], alt)
            else:
                return Indel(edge[0]+1, edge[0], alt)
        elif len(edge) == 2:
            return Indel(edge[0], edge[1], alt)

    @property
    def theory_length(self):
        if self.length > -1:
            return len(self._alt)-self.length
        else:
            return len(self.alt)

    def get_formula(self):
        formula = ''
        if self.length == -1:
            formula = self._alt
        elif self.length > 0:
            formula = '('+str(self._begin)+'-'+str(self._end)+')'+self._alt
        elif self.length == 0:
            formula = str(self._end)+self._alt
        return formula

    def show_formula(self):
        print(self.get_formula(), end="")


class Segment(object):
    def __init__(self, sequence=None, cnv=1, indels=[]):
        self._sequence = sequence
        self._indels = indels
        self._cnv = cnv

    @property
    def theory_length(self):
        length = 0
        t = self._sequence
        if isinstance(t, Sequence):
            length = t.length
            for x in self._indels:
                if x.length > -1:
                    length += x.theory_length
                else:
                    length = x.theory_length-t.length
        if isinstance(t, list):
            for x in t:
                length += x.theory_length
        if isinstance(t, Indel):
            length = t.theory_length
        return length*abs(self._cnv)

    @property
    def cnv(self):
        return self._cnv

    @staticmethod
    def get_segment_formu(formula):
        chrr, formula, cnv = re.match(
            '(\d*)\((.*)\)\**(-*\d*)$', formula).groups()
        try:
            cnv = int(cnv)
        except:
            cnv = 1
        if formula.find(',') == -1:  # sequence
            m = re.match(r'([\d\.]+)(.*\-)([\d\.]+)$', formula.strip())
            if m:
                begin, formula, end = m.groups()
                segment = Segment(Sequence(chrr, begin, end), int(cnv))
                segment.add_var_formula(formula)
            else:
                alt = ''.join(re.findall(r'[A-Za-z]+', formula))
                segment = Segment(Indel('.', '.', alt), int(cnv))
        else:
            subsegments = []
            subformulas = split_formula(formula)
            for x in subformulas:
                subsegment = Segment.get_segment_formu(x)
                subsegments.append(subsegment)
            segment = Segment(subsegments, int(cnv))
        return segment

    def add_var_formula(self, formula):
        formula = formula.strip('-')
        subformulas = split_formula(formula, '-')
        self._indels = []
        for x in subformulas:
            self._indels.append(Indel.get_indel_formula(x))

    @property
    def seq(self):
        seq = ''
        if isinstance(self._sequence, Sequence):
            cur = 0
            seq_t = self._sequence.seq
            self._indels.sort(key=lambda x: x.begin)
            for x in self._indels:
                end = max(cur, x.begin-self._sequence.begin)
                seq += seq_t[cur:end]+x.alt
                cur = x.end-self._sequence.begin+1
            seq += seq_t[cur:]
        elif isinstance(self._sequence, list):
            for x in self._sequence:
                seq += x.seq
        elif isinstance(self._sequence, Indel):
            seq = self._sequence.alt
        if self._cnv < 0:
            seq = seq[::-1]
            seq=read_complement(seq)
        return seq*abs(self._cnv)

    def search_seq(self, filed, infos, column, row_step, MEMORY):
        if isinstance(self._sequence, Sequence):
            chromosome, pos = Fasta.analyse_infos(infos, self._sequence.chr)
            if not self._sequence.search_seq(filed, chromosome.end, pos, column, row_step, MEMORY):
                return False
        elif isinstance(self._sequence, list):
            for x in self._sequence:
                if not x.search_seq(filed, infos, column, row_step, MEMORY):
                    return False
        return True

    def del_seq(self):
        if isinstance(self._sequence, Sequence):
            self._sequence.del_seq()
        elif isinstance(self._sequence, list):
            for x in self._sequence:
                x.del_seq()
        self._sequence = None
        for x in self._indels:
            x.alt = ''

    def get_fomula(self):
        formula = '('
        if isinstance(self._sequence, Sequence):
            formula = str(self._sequence._chr)+'('+str(self._sequence.begin)
            for x in self._indels:
                formula += '-'+x.get_formula()
            formula += '-'+str(self._sequence.end)+')*'+str(self._cnv)

        elif isinstance(self._sequence, list):
            for x in self._sequence:
                formula += x.get_fomula()+','
            formula += ')*'+str(self._cnv)
        elif isinstance(self._sequence, Indel):
            formula += self._sequence.alt+')*'+str(self._cnv)
        return formula

    def get_inserts(self, length):
        t = self._sequence
        inserts = []
        if isinstance(t, Sequence):
            cur, end = t.begin, t.begin
            self._indels.sort(key=lambda x: x.begin)
            for x in self._indels:
                l = x.length
                if l > 0:
                    end = x.begin-1
                elif l > -1:
                    end = x.begin
                inserts.append((t.chr, cur, end, cur-t.begin, end-t.begin))
                cur = end+1
            if t.end > end:
                inserts.append((t.chr, cur, t.end, cur-t.begin, t.end-t.begin))
        if isinstance(t, list):
            for x in t:
                l = length
                inserts += x.get_inserts(l)
                l += x.theory_length
        l = self.theory_length//abs(self._cnv)
        if self._cnv == 1:
            return inserts
        insert = deepcopy(inserts)
        inserts = []
        if self._cnv > 0:
            for n in range(self._cnv):
                for x in insert:
                    inserts.append(
                        (x[0], x[1], x[2], x[3]+l*n+length, x[4]+l*n+length))
        elif self._cnv < 0:
            length += self.theory_length
            for n in range(1+self._cnv, 1):
                for x in insert:
                    inserts.append(
                        (x[0], x[1], x[2], -x[3]+l*n+length-1, -x[4]+l*n+length-1))
        return inserts


class Muta(Sequence):
    def __init__(self, chrr=None, begin=None, end=None,  *segments):
        super().__init__(chrr, begin, end, flank_len=0)
        self._segments = []
        self._segments += list(segments)

    @property
    def segments(self):
        return self._segments

    @ property
    def chrs(self):
        return self._chrs

    @ chrs.setter
    def chrs(self, chrs):
        self._chrs = chrs

    @property
    def theory_length(self):
        length = 0
        for x in self._segments:
            length += x.theory_length
        return length

    def add_var_formula(self, formula):
        formula = formula.strip()
        formula = re.match('\((.*)\)$', formula).group(1)
        subformulas = split_formula(formula)
        for x in subformulas:
            if x.strip():
                self._segments.append(Segment.get_segment_formu(x))

    def get_formula(self):
        formula = '('
        for x in self._segments:
            formula += x.get_fomula()+','
        formula += ')'
        return formula

    def get_self_dele(self):
        return str(self._chr)+'('+str(self.begin)+'-'+str(self.end)+')'

    def get_flank(self, flank_len):
        beds = []
        if flank_len > 0:
            if isinstance(self._end, int):
                beds.append(Sequence(self._chr, self.begin -
                                     flank_len, self._begin-1))
                beds.append(
                    Sequence(self._chr, self._end+1, self._end+flank_len))
            else:
                beds.append(Sequence(self._chr, self.begin -
                                     flank_len+1, self._begin))
                beds.append(
                    Sequence(self._chr, self._begin+1, self._end+flank_len))
        return beds

    @property
    def seq(self):
        seq = ''
        # seq=self.lflank+'.|'
        for x in self._segments:
            seq += x.seq
        # seq+='|.'+self.rflank
        return seq

    def search_seq(self, filed, infos, column, row_step, MEMORY):
        # Sequence.search_seq(self,filed,stop,chr_pos,column,row_step)
        for x in self._segments:
            if not x.search_seq(filed, infos, column, row_step, MEMORY):
                return False
        return True

    def del_seq(self):
        self._seq = ''
        for x in self._segments:
            x.del_seq()
    """
    @staticmethod
    def merge(muta1, muta2, flank_len=FLANK_LEN):
        if muta1.chr == muta2.chr:
            if muta1.check_overlap(muta2):
                print('overlap -> ignore it')
                muta = [muta1]
            elif muta1.end+flank_len >= muta2.begin:
                muta = [Muta(muta1.chr, muta1.begin, muta2.end, *muta1.segments, Segment(
                    Sequence(muta1.chr, muta1.end+1, muta2.begin-1)), *muta2.segments)]
            else:
                muta = [muta1, muta2]
        return muta
    """

    def get_inserts(self, length):
        inserts = []
        length = self._begin+length
        for x in self._segments:
            inserts += x.get_inserts(length)
            length += x.theory_length
        return inserts


class chromosome(object):
    def __init__(self, chrr):
        self._chr = chrr
        self.mutas = []

    @property
    def chr(self):
        return self._chr

    def add_var_formula(self, begin, end, formula):
        muta = Muta(self._chr, begin, end)
        muta.add_var_formula(formula)
        self.mutas.append(muta)

    def show_chromosome(self):
        for x in self.mutas:
            print(x.get_self_dele(), end='\t')
            print(x.get_formula())

    def neaten(self):
        if len(self.mutas) < 2:
            return
        self.mutas.sort(key=lambda x: (x.begin, x.end))

    def get_dele_insert(self):
        ranges = []
        for x in self.mutas:
            ranges.append((x.begin, x.end, x.theory_length))
        return ranges

    def get_inserts(self):
        length = 0
        inserts = []
        for x in self.mutas:
            inserts += x.get_inserts(length)
            length += x.theory_length
        inserts.sort(key=lambda x: (x[0], x[1]))
        return inserts


class Haploid(object):
    def __init__(self, *mutas):
        self.chromosomes = {}
        self.add_mutas(*mutas)
        self.joins = []

    def add_mutas(self, *mutas):
        pass

    def add_var_formula(self, chrr, begin, end, formula):
        if chrr in self.chromosomes:
            self.chromosomes[chrr].add_var_formula(begin, end, formula)
        else:
            self.chromosomes[chrr] = chromosome(chrr)
            self.chromosomes[chrr].add_var_formula(begin, end, formula)

    def show_haploid(self):
        if not self.chromosomes:
            print('normal haploid')
        else:
            print('mutation haploid')
        for x in sorted(self.chromosomes.values(), key=lambda y: int(y.chr)):
            x.show_chromosome()

    def neaten(self):
        for x in self.chromosomes.values():
            x.neaten()

    def write_muta_seq(self, ref, outfile, COLUMN, FASTA_INFO, MEMORY):
        infos = Fasta.fasta_file_info(ref, FASTA_INFO)
        print('write mutation sequence from ', ref)
        with open(ref, 'r', newline='\n') as filed:
            with open(outfile, 'a', newline='\n') as writed:
                writed.write("#exonid\tchr\tbegin\tend\tlength\tflanklen\n")
                for x in sorted(self.chromosomes.keys()):
                    i = 1
                    invalids = []
                    for y in self.chromosomes[x].mutas:
                        chromosome, pos = Fasta.analyse_infos(infos, y.chr)
                        y.id = chromosome.id+'.'+str(i)+'m'
                        if y.end > chromosome.end:
                            print('muta del out of range -> ignore it',
                                  y.chr, y.end, chromosome.end)
                            invalids.append(y)
                            continue
                        if not y.search_seq(filed, infos, infos['column'], infos['step'], MEMORY):
                            print('muta ins out of range -> ignore it',
                                  y.chr, y.end, chromosome.end)
                            invalids.append(y)
                            continue
                        else:
                            y.write_fasta(writed, COLUMN)
                        y.del_seq()
                    for z in invalids:
                        self.chromosomes[x].mutas.remove(z)
                        invalids = []

        print('down. outfile :', outfile)

    def write_muta_bed(self, reg, outfile, BED_INFO, effect_len, chip_len):
        print('reposition from file :', reg)
        with open(outfile, 'w', newline='\n')as f:
            pass
        chrrs = Bed.get_bed_info(reg, BED_INFO)
        for chrr in chrrs:
            print('reposition chromosome', chrr, ' regions ')
            if chrr in self.chromosomes:
                inserts = self.chromosomes[chrr].get_inserts()
                sranges = self.chromosomes[chrr].get_dele_insert()
            else:
                inserts = []
                sranges = []
            sss = []
            for schr in chrrs:
                schr_insets = [x for x in inserts if int(x[0]) == schr]
                if not schr_insets:
                    continue
                ranges = Bed.get_ranges(reg, schr)
                sss = cross_insert(schr_insets, ranges)
            ranges = Bed.get_ranges(reg, chrr)
            if sranges:
                #print('differ and insert chromosome',chrr,end='\r')
                ranges = differ_insert_ranges(ranges, sranges)
            ranges += sss
            ranges = merge_ranges(ranges, 0, effect_len)  # 0:join_gap
            '''
            if sss:
                for i,y in enumerate(ranges):
                    if y[1]-y[0] < chip_len:
                        s=(chip_len-y[1]+y[0]+1)//2
                        a=y[0]-s
                        b=y[1]+s
                        ranges[i]=(a,b)
            '''
            with open(outfile, 'a', newline='\n') as filed:
                for x in ranges:
                    filed.write('%s\t%s\t%s\n' % (chrr, x[0], x[1]))
        print('down. outfile :', outfile)


class Polyploid(object):
    def __init__(self, idd, content, *haploids):
        self._id = idd
        self._content = content
        self.haploids = haploids

    @property
    def id(self):
        return self._id

    def add_var_formula(self, genotype, chrr, begin, end, formula):
        polys = re.findall(r'[*.]+', genotype)
        for n, x in enumerate(polys):
            if x == '*':
                self.haploids[n].add_var_formula(chrr, begin, end, formula)

    def show_polyploid(self):
        print('id=', self._id, ',', 'percent=', self._content,
              ',', 'polyploid number=', len(self.haploids))
        for x in range(len(self.haploids)):
            print('*haploid %d:' % (x+1), end='')
            self.haploids[x].show_haploid()

    def neaten(self):
        for x in self.haploids:
            x.neaten()

    @staticmethod
    def self_poly_formufile(mutafile, idd, content, polynum):
        polyploid = Polyploid(idd, content, *([Haploid()]*polynum))
        with open(mutafile, 'r', newline='\n') as filed:
            for line in filed.readlines():
                if line[0] != '#' and line[0] != '\n':
                    infos = line.split('\t')
                    genotype, chrr, begin, end, formula = str2int(infos[:5])
                    polyploid.add_var_formula(
                        genotype, chrr, begin, end, formula)
        polyploid.neaten()
        return polyploid


def ini_muta(inireferences, iniregions, mut, inifile, polynum, COLUMN, MEMORY, BED_INFO, FASTA_INFO, effect_len, chip_len):
    mutation = mut[2]
    # load one polyoid mutation set
    polyploid = Polyploid.self_poly_formufile(
        mutation, mut[0], mut[1], polynum)
    # polyploid.show_polyploid()
    #!!!!!
    # output its(one polyoid) targed region on mutation reference
    for n, reg in enumerate(iniregions):
        regmutation = reg+'.reg'+str(polyploid.id)
        polyploid.haploids[n].write_muta_bed(
            reg, regmutation, BED_INFO, effect_len, chip_len)
    with open(inifile, 'w', newline='\n') as writed:
        pass
    # output its(one polyoid) all mutations sequence
    for n, ref in enumerate(inireferences):  # n is the no.n haploid
        # for every haploid , output its all mutations sequence
        polyploid.haploids[n].write_muta_seq(
            ref, inifile, COLUMN, FASTA_INFO, MEMORY)
    '''
    # output its(one polyoid) all mutation  reference
    for n, ref in enumerate(inireferences):
        mutationref = ref+str(poly[0])+'.ref'
        write_ref_muta(ref, inifile, mutationref)
    for ref, reg in zip(inireferences, iniregions):
        x = CD+ref+str(poly[0])+'.ref'
        y = CD+reg+str(poly[0])+'.reg'
        z = CD+str(poly[0])+'.seq'
        Fasta.ini_exome(x, y, z)
    '''


def vcf2formula(file):
    formulafile = file+'formula'
    with open(file, 'r', newline='\n')as filed:
        with open(formulafile, 'w', newline='\n')as writed:
            for line in filed.readlines():
                if re.match('#', line):
                    continue
                infos = line.split()
                chrr, pos, idd, ref, alt, qual, filt = infos[:7]
                chrr = chrr.strip('chr')
                end = int(pos)+len(ref.strip('.'))-1
                alt = alt.strip('.')
                if alt:
                    alt = '('+alt+')*1'
                alt = '('+alt+')'
                tumor = infos[-1]
                tumor = tumor.split(':')[0]
                tumor = tumor.replace('1', '*')
                tumor = tumor.replace('0', '.')
                if '*' in tumor:
                    s = '\t'.join([tumor, chrr, pos, str(end), alt])
                    writed.write(s+'\n')


if __name__ == '__main__':
    info = sys.argv
    vcf2formula(info[1])
