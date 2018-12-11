import collections
import re

from exon import Exon, get_seq
from filefunc import get_line_text, get_pos_text,num_positive


def get_chr_exon_num(file):
    with open(file, 'r') as f:
        echr = 0
        l = collections.OrderedDict()
        count = 0
        summ = 0
        for line in f.readlines():
            m = re.match(r'^>(chr[\w\d]*)\s*', line)
            if m:
                if m.group(1) != echr:
                    if echr != 0:
                        l[echr] = count
                        summ += count
                    echr = m.group(1)
                    count = 1
                else:
                    count += 1
        l[echr] = count
        summ += count
    l['total'] = summ
    return l


def show_chr_exon_num(exon_nums):
    for x, key in list(enumerate(exon_nums.keys())):
        if (x+1) % 5 == 0:
            print('')
        print('%s=%d' % (key, exon_nums[key]), end='\t')
    print('')


def chr_pos(file):
    echr = 0
    chr_poss = {}
    with open(file, 'r')as f:
        line = f.readline()
        while(line):
            m = re.match(r'^>chr([\w\d]*)', line)
            if m:
                if m.group(1) != echr:
                    chr_poss[m.group(1)] = f.tell()-len(line)
                    echr = m.group(1)
            line = f.readline()
    return chr_poss


def show_exon(filed, info, pos):
    line = get_line_text(
        filed, r'^>chr([\w\d]+)\t+(\d+)\t+(\d+)\t+(%s).*' % info[1], 0, 're', pos=pos)
    m = re.match(r'^>chr([\w\d]*)\t+(\d*)\t+(\d*)\t+(%s)' % info[1], line)
    exoninfo = line.split()
    if m:
        seq, line = get_seq(filed)
        exon = Exon(m.group(1), int(m.group(2)), int(m.group(3)),
                    m.group(4), seq, insert=int(exoninfo[-1]))
    if len(info) == 2:
        exon.show(1, -1)
    elif len(info) == 3:
        exon.show(int(info[2]), -1)
    elif len(info) == 4:
        exon.show(int(info[2]), int(info[3]))
    elif len(info) == 5:
        if info[4].upper() == "T":
            info[4] = True
        else:
            info[4] = False
        exon.show(int(info[2]), int(info[3]), info[4])


def show_pos(f, info, pos):
    exons = []
    f.seek(pos, 0)
    line = f.readline()
    pos1 = int(info[2])
    pos2 = int(info[3])
    while(line):
        m = re.match(
            r'^>chr([\d]*)\t(\d*)\t(\d*)\t([\w\.\-]*)\t\d*\t(\d*)', line)
        if m:
            exoninfo = line.split()
            begin = int(m.group(2))
            end = int(m.group(3))
            if begin <= pos2 and end >= pos1:
                seq, line = get_seq(f)
                exons.append(Exon(m.group(1), int(m.group(2)), int(
                    m.group(3)), m.group(4), seq, insert=int(exoninfo[-1])))
            elif begin > pos2:
                break
            else:
                line = f.readline()
        else:
            line = f.readline()
    if exons:
        rpos1 = num_positive(pos1-exons[0].begin)
        rpos2 = exons[-1].length-num_positive(exons[-1].end-pos2)
        if len(exons) == 1:
            exons[0].show(rpos1+1, rpos2)
        elif len(exons) > 1:
            for x in exons:
                if x == exons[0]:
                    x.show(rpos1, -1)
                elif x == exons[-1]:
                    x.show(1, rpos2)
                else:
                    x.show()


def view(file):
    helps = '''
    chr -n                  : show chromosome exon number
    chrx -n                 : shoe chx exon number
    chrx -pos -pos1 -pos2   : show chx from pos1 to pos2 sequence,pos is absolute
    chrx -num               : show chx no.num exon sequence
    chrx -num -pos1 -pos2   : show chx no.num from pos1 to pos2 sequence,pos is relative
    chrx -num -pis1 -pot2 -t: show chx no.num exon sequence,mark pos1 and pos2
    '''
    print(helps)
    exon_nums = get_chr_exon_num(file)
    chr_poss = chr_pos(file)
    with open(file, 'r') as f:
        line = input('>>')
        while(line.strip() != 'exit'):
            if line.strip() == 'help':
                print(helps)
            elif line == 'chr -n':
                show_chr_exon_num(exon_nums)
            elif re.match(r'(chr[\d\w]+)\s*-n', line):
                m = re.match(r'(chr[\d\w]+)\s*-n', line)
                print(m.group(1), '=', exon_nums[m.group(1)])
            elif re.match(r'chr([\d\w]+)\s*-pos\s-(\d+)\s-(\d+)', line):
                m = re.match(r'chr([\d\w]+)\s*-pos\s-(\d+)\s-(\d+)', line)
                info = line.split(' -')
                show_pos(f, info, chr_poss[m.group(1)])
            elif re.match(r'chr([\d\w]+)\s-(\d+)', line):
                m = re.match(r'chr([\d\w]+)\s*-(\d*)', line)
                info = line.split(' -')
                info[0] = m.group(1)
                info[1] = 'CH'+(m.group(1))+'-'+info[1]
                show_exon(f, info, chr_poss[m.group(1)])

            line = input('>>')


if __name__ == '__main__':
    file = input("please enter exonlist file you want to look for : ")
    view(file)
