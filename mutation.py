from filefunc import *
from sequence import *


class Indel(Edge):
    def __init__(self, begin, end, alt):
        super().__init__(begin, end)
        self._alt = alt.strip('.')

    @property
    def alt(self):
        return self._alt

    @staticmethod
    def get_indel_formu(formula):
        alt = ''.join(re.findall(r'[A-Za-z]+', formula))

        edge = str2int(re.findall(r'\d+', formula))
        if len(edge) == 0:
            return Indel('.', '.', alt)
        elif len(edge) == 1:
            if '[' in formula:
                return Indel(edge[0], edge[0], alt)
            else:
                return Indel(edge[0]+1, edge[0], alt)
        elif len(edge) == 2:
            return Indel(edge[0], edge[1], alt)

    def get_formu(self):
        formula = ''
        if self.length == -1:
            formula = self._alt
        elif self.length > 0:
            formula = '['+str(self._begin)+'~'+str(self._end)+']'+self._alt
        elif self.length == 0:
            formula = str(self._end)+self._alt
        return formula

    def show_formula(self):
        print(self.get_formu(), end="")


class Segment(object):
    def __init__(self, sequence=None, cnv=1, segments=None, indels=None):
        self._sequence = sequence
        self._segments = segments
        self._indels = indels
        self._cnv = cnv

    def __check(self):
        return not(self._segments and self._sequence)

    @property
    def cnv(self):
        return self._cnv

    @staticmethod
    def __unbracket(formula):
        formula = formula.strip()
        m = re.match(r'.*[*]([\-\d]+)$', formula)
        if m:
            cnv = int(m.group(1))
            formula = formula.rstrip(' %s' % cnv)
            formula = formula.rstrip('*').strip()[1:-1]
        else:
            cnv = 1
            if formula[0] == '(':
                formula = formula[1:-1]
        formulas = re.findall(r'[^,]+', formula)
        return formulas, cnv

    @staticmethod
    def get_segment_formu(chrs, formula):
        segs, cnv = Segment.__unbracket(formula)
        if isinstance(chrs, tuple):
            csegments = []
            cur = 0
            for x in range(len(chrs)):
                if isinstance(chrs[x], (tuple)):
                    formula = ','.join(segs[cur:cur+deeplen(chrs[x])])
                    cur += deeplen(chrs[x])
                else:
                    formula = segs[cur]
                    cur += 1
                segment = Segment.get_segment_formu(chrs[x], formula)
                csegments.append(segment)
            segment = Segment(cnv=1, segments=csegments)
            return segment
        elif isinstance(chrs, int):
            formula = segs[0]
            m = re.match(r'^(\D*)(\d+)(.+?)(\d+)(\D*)$', formula)
            if m:
                formula = ''.join(m.groups()[::2])
                begin, end = str2int((m.group(2), m.group(4)))
            else:
                begin, end = '.', '.'
            segment = Segment(Sequence(chrs, begin, end), cnv, None, [])
            segment.add_formu(formula)
            return segment

    def add_formu(self, formula):
        formulas = re.findall('[^-]+', formula)
        self._indels = []
        for x in formulas:
            self._indels.append(Indel.get_indel_formu(x))

    def __sort_indels(self):
        pass

    def get_formu(self):
        formula = []
        if self._segments:
            formula = []
            for x in self._segments:
                formula.append(x.get_formu())
            formula = ','.join(formula)
            formula = '('+formula+')'
            if abs(self._cnv) != 1:
                formula += '*'+str(self._cnv)
        elif self._sequence:
            for x in self._indels:
                if x.get_formu():
                    formula.append(x.get_formu())
            if isinstance(self._sequence._begin, int):
                if self._indels and self._indels[0]._begin == '.':
                    formula.insert(1, str(self._sequence._begin))
                else:
                    formula.insert(0, str(self._sequence._begin))
                formula.append(str(self._sequence.end))
            formula = '-'.join(formula)
            if abs(self._cnv) != 1:
                formula = '('+formula+')'
                formula += '*'+str(self._cnv)
        return formula

    def insert_borders(self):
        if self._sequence:
            if self._sequence.length > 0:
                return [self._sequence]
            else:
                return []
        elif self._segments:
            beds = []
            for x in self._segments:
                beds += x.insert_borders()
            return beds
        else:
            return []

    @property
    def seq(self):
        seq = ''
        if self._sequence and self._indels:
            if self._sequence.length < 1:
                for x in self._indels:
                    seq += x.alt
            else:
                for x in self._indels:
                    a = Edge(self._sequence.begin, x.begin)
                    b = Edge(self._sequence.begin, x.end)
                    c = 0
                    # indel(10,20)
                    if a.length > 0 and b.length > 0:
                        seq += self._sequence.seq[c:a.length]+x.alt
                        c = b.length+1
                    # indel(10,0)
                    elif a.length > 0:
                        seq += self._sequence.seq[c:a.length+1]+x.alt
                        c = a.length+1
                    else:
                        seq = x.alt+seq
                seq += self._sequence.seq[c:]
        elif self._segments:
            for x in self._segments:
                seq += x.seq()
        if self._cnv < 0:
            seq = seq[::-1]
        return seq*abs(self.cnv)

    def search_seq(self, filed):
        if self._sequence:
            self._sequence.search_seq(filed)
        elif self._segments:
            for x in self._segments:
                x.search_seq(filed)

    def del_seq(self):
        self._segments = self._sequence = self._indels = None


class Muta(Border):
    def __init__(self, chrr=None, begin=None, end=None, chrs=None, cnv=1, *segments):
        super().__init__(chrr, begin, end)
        self._segments = []
        self._segments += list(segments)
        self._chrs = chrs
        self._cnv = cnv

    @property
    def cnv(self):
        return self._cnv

    @cnv.setter
    def cnv(self, cnv):
        self._cnv = cnv

    @ property
    def chrs(self):
        return self._chrs

    @ chrs.setter
    def chrs(self, chrs):
        self._chrs = chrs

    @ staticmethod
    def __unbracket(formula):
        formula = formula.strip()
        m = re.match(r'.*[*]([\-\d]+)$', formula)
        if m:
            cnv = int(m.group(1))
            formula = formula.rstrip(' %s' % cnv)
            formula = formula.rstrip('*').strip()[1:-1]
        else:
            cnv = 1
            formula = formula[1:-1]
        return formula, cnv

    @staticmethod
    def get_muta_formu(formula):
        formulas = formula.split(':')
        dele = str2int(re.findall(r'[^(,)]+', formulas[0].strip()))
        chrs = eval(formulas[1].strip())
        formulas[2], cnv = Muta.__unbracket(formulas[2])
        muta = Muta(dele[0], dele[1], dele[2], chrs, cnv)
        muta.add_formu(chrs, formulas[2])
        return muta

    def add_formu(self, chrs, formula):
        cur = 0
        segs = re.findall(r'[^,]+', formula)
        for x in chrs:
            if isinstance(x, (tuple, list)):
                formula = ','.join(segs[cur:cur+deeplen(x)])
                cur += deeplen(x)
            else:
                formula = segs[cur]
                cur += 1
            self._segments.append(Segment.get_segment_formu(x, formula))

    def get_formu(self):
        # formula1 is dele
        formula1 = str(self._chr)+'('+str(self.begin)+','+str(self.end)+')'
        # formula2 is insert chrs
        formula2 = str(self._chrs)
        # formula3 is insert content
        formula3 = []
        for x in self._segments:
            formula3.append(x.get_formu())
        formula3 = '('+','.join(formula3)+')'
        if self._cnv != 1:
            formula3 = formula3+'*'+str(self.cnv)
        formula = ':'.join([formula1, formula2, formula3])
        return formula

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
        for x in self._segments:
            seq += x.seq
        if self._cnv < 0:
            seq = seq[::-1]
        if abs(self.cnv) > 1:
            seq = seq*abs(self.cnv)
        return seq

    def insert_borders(self):
        beds = []
        for x in self._segments:
            beds += x.insert_borders()
        return beds

    def search_seq(self, filed):
        for x in self._segments:
            x.search_seq(filed)

    def dele_edge(self):
        if isinstance(self._begin, int):
            return (self._begin, self._end)

    def del_seq(self):
        for x in self._segments:
            x.del_seq()


class Genome(object):
    def __init__(self, *mutas):
        self.mutas = {}
        self.add_matas(*mutas)
        self.joins = []

    def add_matas(self, *mutas):
        for x in mutas:
            self.mutas[x.chr][x.dele_edge()] = x

    def add_formu(self, formula):
        muta = Muta.get_muta_formu(formula)
        self.mutas[muta.chr] = {muta.dele_edge(): muta}

    def show_genome(self):
        if not self.mutas:
            print('normal genome')
        else:
            print('mutation genome')
        for x in self.mutas.values():
            for y in x.values():
                print('\t'+y.get_formu())

    def insert_borders(self):
        borders = []
        for x in self.mutas.values():
            for y in x.values():
                borders += y.insert_borders()
        return borders

    def get_normal_num(self, border):
        if border.chr in self.mutas:
            for x in self.mutas[border.chr].values():
                if x.check_overlap(border):
                    return 0
        return 1

    def dele_borders(self, chromosome, ranges, flank_len, join_gap):
        '''ranges and deles[chromosome] 's ele is tuple'''
        borders = []
        if chromosome.chr not in self.mutas:
            return borders
        dele_ranges = list(self.mutas[chromosome.chr].keys())
        sort_ranges = sort_range(ranges+dele_ranges, join_gap)
        in_ranges = in_range(sort_ranges, dele_ranges)
        count_ranges = count_range(in_ranges, dele_ranges)
        for x, y in count_ranges.items():
            # y is self.mutas--muta
            joins = []
            for i in y:
                joins.append(self.mutas[chromosome.chr][i])
            split_ranges = split_range(x, y)
            for a, b in split_ranges:
                borders.append(Sequence(chromosome.chr, a, b))
            joins += borders
            joins.sort(key=lambda border: border.begin)
            flanks = []
            flank_len=min(joins[0].begin-1,flank_len,chromosome.end-joins[-1].end)
            flanks.append(Sequence(chromosome.chr, joins[0].begin-flank_len, joins[0].begin-1))
            flanks.append(Sequence(chromosome.chr, joins[-1].end+1, joins[-1].end+flank_len))
            joins.insert(0, flanks[0])
            joins.append(flanks[1])
            borders += flanks
            self.joins.append(joins)
        return borders

    def search_seq(self, filed, *, no, genome):
        seq = ''
        for i, x in enumerate(self.joins):
            try:
                for y in x:
                    y.search_seq(filed)
                    seq += y.seq
            except TypeError as e:
                print(e)
                continue
            yield Sequence(y.chr, x[0].end+1, x[-1].begin-1, 'Muta%s-%s-%s' % (no, genome, i), seq, x[0].length)

    def seq(self, filed, **info):
        f, g = self.search_seq(filed, **info), self.del_seq()
        while(True):
            try:
                yield next(f)
                next(g)
            except StopIteration as e:
                print('\tcomplete.',e)
                break

    def del_seq(self):
        for x in self.joins:
            for y in x:
                y.del_seq()
                yield


class Polyploid(object):
    def __init__(self, percent, *genomes):
        self.percent = percent
        self.genomes = genomes

    def add_formu(self, poly, formula):
        poly = re.findall(r'[*.]+', poly)
        for x in range(len(poly)):
            if poly[x] == '*':
                self.genomes[x].add_formu(formula)

    def show_polyploid(self):
        print('percent=', self.percent)
        for x in range(len(self.genomes)):
            print('\tgenome%d:' % (x+1), end='')
            self.genomes[x].show_genome()

    def get_normal_num(self, border):
        num = 0
        for x in self.genomes:
            num += x.get_normal_num(border)
        num = num/len(self.genomes)
        return num

    def dele_borders(self, chromosome, ranges, flank_len, join_gap):
        borders = []
        for x in self.genomes:
            borders += x.dele_borders(chromosome, ranges, flank_len, join_gap)
        return borders

    def insert_borders(self):
        borders = []
        for x in self.genomes:
            borders += x.insert_borders()
        borders = list(set(borders))
        return borders


class Tissue(object):
    def __init__(self, **polyploids):
        '''ployploids is dict'''
        self.polyploids = dict(polyploids)

    def show_tissue(self):
        print('\t'+'*'*10+'Show Mutations'+'*'*10)
        for x in sorted(self.polyploids.keys()):
            print('\tno.%s' % x, end=":")
            self.polyploids[x].show_polyploid()
        print('\t'+'*'*10+'Show Mutations'+'*'*10)

    def get_normal_num(self, border):
        num = 0
        for x in sorted(self.polyploids.keys()):
            num += self.polyploids[x].percent * \
                self.polyploids[x].get_normal_num(border)
        return num

    def dele_borders(self, chromosome, ranges, flank_len, join_gap):
        '''not all ,only one chromosome '''
        borders = []
        for x in sorted(self.polyploids.keys()):
            borders += self.polyploids[x].dele_borders(
                chromosome, ranges, flank_len, join_gap)
        borders = list(set(borders))
        return borders

    def insert_borders(self):
        borders = []
        for x in sorted(self.polyploids.keys()):
            borders += self.polyploids[x].insert_borders()
        borders = list(set(borders))
        return borders

    @staticmethod
    def formula_mutation(mname):
        print('getting mutation setting from %s ...' % mname)
        with open(mname, 'r') as filed:
            tissue = Tissue(**Tissue.get_head(filed))
            line = filed.readline()
            while(line):
                muta = re.findall('[^\t]+', line)
                if len(muta) == 3:
                    tissue.polyploids[muta[0]].add_formu(muta[1], muta[2])
                line = filed.readline()
        return tissue

    @staticmethod
    def get_head(filed):
        mutation_types = {}
        get_line_text(filed, r'<mutation percent>', 1, 're', pos=0)
        line = filed.readline()
        while(line):
            if re.match(r'<mutation details>', line):
                break
            mutas = line.strip().split()
            mutation_types[mutas[0]] = Polyploid(
                float(mutas[1]), Genome(), Genome())
            line = filed.readline()
        get_line_text(filed, r'<mutation details>', 1, 're', pos=0)
        return mutation_types


def set_mutation(mname):
    choice = input('''
    which mutation creating way do you want?
    1. auto set
    2. formula file set
    3. table file set
    >''')
    if choice == '1':
        tissue = Tissue()
    elif choice == '2':
        tissue = Tissue.formula_mutation(mname)
    elif choice == '3':
        mname = input('please input mutation setting file : ')
        tissue = Tissue()
    elif choice == 'exit':
        return 0
    else:
        print("your input cannot be identified.")
        return -1
    tissue.show_tissue()
    return tissue


def write_wes_bed(file, write, tissue, flank_len, join_gap,*,chrss):
    borders = []
    with open(file, 'r') as filed:
        ranges = []
        echr = 0
        for line in filed.readlines():
            m = Border.get_border(line)
            if m:
                if m.chr != echr:
                    if echr in chrss:
                        chromosome=chrss[echr]['sequence']
                        borders += tissue.dele_borders(chromosome,
                                                    ranges, flank_len, join_gap)
                    ranges = []
                    echr = m.chr
                ranges.append((m.begin, m.end))
        if echr in chrss:
            chromosome=chrss[echr]['sequence']
            borders += tissue.dele_borders(chromosome, ranges, flank_len, join_gap)
    borders += tissue.insert_borders()
    borders = list(set(borders))
    borders.sort(key=lambda border: border.chr)
    Bed.write_beds(write, borders)

def write_wgs_bed(write, tissue, flank_len, join_gap,chrss):
    borders=[]
    for x in sorted(chrss.keys()):
        chromosome=chrss[x]['sequence']
        ranges=[(chromosome.begin,chromosome.end)]
        borders += tissue.dele_borders(chromosome,ranges, flank_len, join_gap)
    borders += tissue.insert_borders()
    borders = list(set(borders))
    borders.sort(key=lambda border: border.chr)
    Bed.write_beds(write, borders)

def get_muta_seq(filed, muta, flank_len):
    muta.search(filed)
    seq = muta.get_insert_seq()
    length = len(seq)
    seq_id = 'MUTA%s-' % muta.chr+str(muta.begin)
    flanks = muta.get_flank(flank_len)
    for m in flanks:
        m.search(filed)
    sequence = Sequence(muta.chr, muta.begin, muta.begin+length,
                        seq_id, flanks[0].seq+seq+flanks[1].seq, flank_len)
    return sequence


if __name__ == '__main__':
    tissue = Tissue.formula_mutation('formula_muta.txt')
    tissue.show_tissue()
