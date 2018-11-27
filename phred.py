import json
import random
from filefunc import get_column_row


'''
def random_weight_choice(lists):
    summ=numpy.sum(lists,axis=0)[1]
    t = random.randint(1, summ)
    for i, val in lists:
        t -= val
        if t <= 0:
            return chr(i+33)
    if t>0:
        input('kkkkkkkkkkkkkkkkk')
'''
def char2phred(char,plus=33):
    code=ord(char)-plus
    return code
def get_phred_fre(file,write,plus):
    print('get phred frequencies from %s...'%file)
    frequencies={}
    readlen=get_column_row(file,'+',1,'re')
    for x in range(1,readlen+1):
        frequencies['pos%d_frequencies'%x]=[0]*43
    with open(file,'r') as f:
        i =0
        for line in f.readlines():
            i+=1
            if i%4==0:
                row_fastq=line.strip()
                x=1
                for char in row_fastq:
                    frequencies['pos%d_frequencies'%x][char2phred(char,plus)]+=1
                    x=x+1
            if i%40000==0:
                print(i,end='\r')
    with open(write,"w") as f:
        json.dump(frequencies, f)
    print('\ndown. outfile: %s'%write)

def random_qphred(rang,frequencies,plus=33):
    s=''
    for x in range(1,len(frequencies)+1):
        weights=frequencies['pos%d_frequencies'%x]
        num=random.choices(rang,weights=weights)
        char=chr(num[0]+plus)
        s=s+char
    return s