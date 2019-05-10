

import os
import sys
from functools import reduce

from basic import *
from mutation import *
from read import *
from sequence import *

random.seed(100)

if __name__ == '__main__':
    time_start = time.time()

    help = '''
    -ref          : get whole genome from filex , generate file3
    -reg          : get aimed and sorted annotation from filey, generate file2
    -qph          : get qphred frequencies from file,generate file3
    -read         : get aimed regions pe fastq from file1234 (deafault flitrate 'N')
    -seq          : get aimed sequence from file1,file2, generate filea
    -view         : view  sequence or depth
    -dep <file> <depth>  : repositon regions from depth file
    -clear        : clear temp files
    optional:
    -wes/wgs      : set the way to generate DNA segment
    -se/pe        : set the way to generate reads
    -cd           : set path to store
    -R            : set output fastq files extra name
    -ini          : set peofile
    if occur something wrong, try delete bed_info and fasta_info
    '''
    info = sys.argv
    print(info)

    if '-ini' in info:
        profile = info[info.index('-ini')+1]
    else:
        profile = 'profile1.ini'
    conf = configparser.ConfigParser()
    conf.read(profile)
    CHIP_LEN, E_LEN, FLANK_LEN, JOIN_GAP, SLIDE_STEP = get_value(
        conf, 'chip', int, "CHIP_LEN", "E_LEN", "FLANK_LEN", 'JOIN_GAP', 'SLIDE_STEP')
    ERROR_E, ERROR_D, SUBSTITUTION, INSERTION, DELETION = get_value(
        conf, "error", float,  "ERROR_E", "ERROR_D", "SUBSTITUTION", "INSERTION", "DELETION")
    ERROR = conf.getboolean("error", "ERROR")
    DEPTH, MODE, PAIR = get_value(conf, 'read', int, "DEPTH", 'MODE', 'PAIR')
    MISMATCH = conf.getboolean("read", 'MISMATCH')
    FLI_N, INNER_N = get_value(conf, "sequence", bool, 'FLI_N', 'INNER_N')
    QPH, QROW = get_value(conf, "quality", int, "QPH", "QROW")
    CD, QUALITY, BED_INFO, FASTA_INFO, REFERENCES, REGIONS, MUTATIONS,LABEL = get_value(
        conf, 'file', str, "CD", "QUALITY", 'BED_INFO', "FASTA_INFO", 'REFERENCES', 'REGIONS', 'MUTATIONS','LABEL')
    COLUMN, MEMORY, = get_value(conf, 'file', int, "COLUMN", "MEMORY")
    SEGMENT_E = conf.getint("segment", "SEGMENT_E")
    SEGMENT_D = conf.getfloat("segment", "SEGMENT_D")
    SEG = conf.get("segment", "SEG")
    REFERENCES = 'REFERENCES='+REFERENCES
    REGIONS = 'REGIONS='+REGIONS
    MUTATIONS = 'MUTATIONS='+MUTATIONS
    # if FLANK_LEN<SEGMENT_E-E_LEN:
    #    FLANK_LEN=SEGMENT_E
    exec(REFERENCES)
    exec(REGIONS)
    exec(MUTATIONS)
    BED_INFO = CD+BED_INFO
    FASTA_INFO = CD+FASTA_INFO
    FNA = get_dict(conf, "fna")
    BED = get_dict(conf, "bed")
    t = time.strftime('%Y%m%d_%H_%M_%S', time.localtime(time.time()))
    if LABEL:
        t=LABEL
    R1 = CD+"R1_%s.fastq" % (t)
    R2 = CD+"R2_%s.fastq" % (t)
    R0 = CD+"R_%s.fastq" % (t)
    conf = None

    if '-wes' in info:
        MODE = modes.WES.value

    if '-wgs' in info:
        MODE = modes.WGS.value

    if '-pe' in info:
        PAIR = pairs.PE.value

    if '-se' in info:
        PAIR = pairs.SE.value

    if '-cd' in info:
        CD = info[info.index('-cd')+1]

    if '-R' in info:
        string = info[info.index('-R')+1]
        R1 = CD+"R1_%s.fastq" % (string)
        R2 = CD+"R2_%s.fastq" % (string)
        R0 = CD+"R_%s.fastq" % (string)

    DEFAULT = {
        'chip_len': CHIP_LEN,
        'effect_len': E_LEN,
        "flank_len": FLANK_LEN,
        "ERROR": ERROR,
        "error_e": ERROR_E,
        "error_d": ERROR_D,
        "substitution": SUBSTITUTION,
        "insertion": INSERTION,
        "deletion": DELETION,
        "DEPTH": DEPTH,
        "MODE": MODE,
        "PAIR": PAIR,
        "MISMATCH": MISMATCH,
        "FLI_N": FLI_N,
        "INNER_N": INNER_N,
        "QPH": QPH,
        "QROW": QROW,
        "CD": CD,
        "COLUMN": COLUMN,
        "MEMORY": MEMORY,
        "REFERENCES": REFERENCES,
        "REGIONS": REGIONS,
        "QUALITY": QUALITY,
        "BED_INFO": BED_INFO,
        "FASTA_INFO": FASTA_INFO,
        'SEG': SEG,
        'SEGMENT_D': SEGMENT_D,
        'SEGMENT_E': SEGMENT_E,
        'join_gap': JOIN_GAP,
        'MUTATIONS': MUTATIONS,
        'FNA': FNA,
        "BED": BED,
        'R1': R1,
        'R2': R2,
        "R0": R0,
    }
    references = set(REFERENCES)
    regions = set(REGIONS)  # filtrate the same region file
    inireferences = list(
        map(lambda x: CD+'ini'+x[0].split('/')[-1].strip(), references))
    iniregions = list(
        map(lambda x: CD+'ini'+x[0].split('/')[-1].strip(), regions))
    iniexomes = list(
        map(lambda x: CD+'exome'+x[0].split('/')[-1].strip()+'.fna', REGIONS))
    inimutations = list(
        map(lambda x: CD+'ini'+x[2].split('/')[-1].strip(), MUTATIONS))
    stanfiles = list(
        map(lambda x: CD+x[0].split('/')[-1].strip()+'.temp', regions))
    iniquality = CD+'ini'+QUALITY.split('/')[-1].strip()

    if not os.path.exists(CD):
        os.makedirs(CD)

    if not os.path.exists(BED_INFO):
        open(BED_INFO, 'a').close()
    if not os.path.exists(FASTA_INFO):
        open(FASTA_INFO, 'a').close()

    if '-help' in info:
        print(help)

    elif '-ref' in info:  # initialize reference genomes(single haploid)
        n = 1
        for x, y in zip(REFERENCES, inireferences):  # mark number
            Fasta.ini_ref(x, y, n, FNA, COLUMN, MEMORY)
            n += 1

    elif '-reg'in info:  # merge targed regions
        for x, y in zip(regions, iniregions):
            keys = Bed.ini_reg(x, y, BED, BED_INFO, E_LEN,
                               CHIP_LEN, 0)  # 0:join_gap

    elif '-qph' in info:
        Quality.ini_qph(QUALITY, iniquality, QROW)

    elif '-seq' in info:  # initailize exome sequence from initialized regions
        for x, y, z in zip(inireferences, iniregions, iniexomes):
            Fasta.ini_exome(x, y, z, E_LEN, 0, FLI_N, INNER_N,
                            COLUMN, MEMORY, FASTA_INFO)  # 0:flank_len

    elif '-mut' in info:  # initialize mutations
        for mut, y, in zip(MUTATIONS, inimutations):
            ini_muta(inireferences, iniregions, mut, y, len(REFERENCES),
                     COLUMN, MEMORY, BED_INFO, FASTA_INFO, E_LEN, CHIP_LEN)

    elif '-read' in info:
        if pairs.PE.value == PAIR:
            open(R1, 'w', newline='\n').close()
            open(R2, 'w', newline='\n').close()
        elif PAIR == pairs.SE.value:
            open(R0, 'w', newline='\n').close()
        for x, mut in zip(inimutations, MUTATIONS):
            # polyoid 's id,content,mutationseq,inireferences,mutationsbed
            # mutationsbed can get from iniregions+polys
            readout(inireferences, iniregions, iniquality, x, mut,LABEL, **DEFAULT)
 
    elif '-dep' in info:
        depfile = info[info.index('-dep')+1]
        depth = int(info[info.index('-dep')+2])
        segment_e = SEGMENT_E
        if os.path.exists(SEG):
            print('get segment length from file :', SEG)
            segment_e = segfile(SEG)[2]
        Depth.dep2bed(depfile, depth, segment_e, CHIP_LEN,
                      E_LEN, CD, JOIN_GAP, BED_INFO, SLIDE_STEP)

    elif '-clear' in info:
        def clear(listt):
            for x in listt:
                if os.path.exists(x):
                    os.remove(x)
        clear(iniexomes)
        clear(inimutations)
        clear(iniregions)
        clear(iniexomes)
        clear(stanfiles)
        clear([iniquality, BED_INFO, FASTA_INFO])
        for x in os.listdir(CD):
            if 'ini' in x[:4]:
                os .remove(CD+x)
            elif '.temp' in x[-4:]:
                os.remove(CD+x)

    elif '-view' in info:
        view(FASTA_INFO, MEMORY)

    time_end = time.time()
    t = time_end-time_start
    print('totally cost: %dh : %dm : %ds' % (t//3600, (t % 3600)//60, t % 60))
