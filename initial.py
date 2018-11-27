import re
import view
from phred import get_phred_fre
from exonlist import get_chr_seq,get_exon_ano,get_exon_list
from setting import *

print("default:","INSERT maxlength=",MAXINSERT,"CHIP length=",CHIP_LEN)
def get_all(defaults):
    filex,filey=defaults['filex'],defaults['filey']
    file1,file2,file3=defaults['file1'],defaults['file2'],defaults['file3']
    filez,file4,plus=defaults['filez'],defaults['file4'],defaults['ASCII']
    ver=defaults['ver']
    get_chr_seq(filex,file1,ver)
    get_exon_ano(filey,file2,ver)
    get_exon_list(file1,file2,file3,ver,CHIP_LEN,JOIN_GAP,ROW_STEP,MAXINSERT)
    get_phred_fre(filez,file4,plus)




if __name__ == '__main__':
    help='''
    init -seq -filex -file1 -ver    : get chromosome sequence from filex, generate file1
    init -ano -filey -file2 -ver    : get exon annotation from filey, generate file2
    init -list -file1 -file2 -file3 -ver :get exonlist from file1,file2, generate file3
    init -qph -filez -file4 -XX     : get qphred frequencies from file,generate 'phred.json'
    init -all                       :excute above all according setting,py
    view -file3                     : view exonlist's exon sequence
    '''
    print(help)
    opera=input('>')
    while(opera.strip()!='exit'):
        if opera.strip=='help':
            print(help)
        elif re.match(r'init -seq',opera):
            info=opera.split(' -')
            if len(info)>=5:
                filex,file1,ver=info[2:]
                get_chr_seq(filex,file1,ver)
        elif re.match(r'init -ano',opera):
            info=opera.split(' -')
            if len(info)>=5:
                fiely,file2,ver=info[2:]
                get_exon_ano(fiely,file2,ver)
        elif re.match(r'init -list',opera):
            info=opera.split(' -')
            if len(info)>=6:
                file1,file2,file3,ver=info[2:]
                get_exon_list(file1,file2,file3,ver,CHIP_LEN,JOIN_GAP,ROW_STEP,MAXINSERT)
        elif re.match(r'init -all',opera):
            get_all(DEFAULTS)
        elif re.match(r'init -qph',opera):
            info=opera.split(' -')
            if len(info)>=3:
                filex,file4,plus=info[1:]
                get_phred_fre(filex,file4,plus)
        elif re.match(r'view',opera):
            info=opera.split(' -')
            if len(info)>=2:
                file=info[1]
                view.view(file)
        opera=input('>')
