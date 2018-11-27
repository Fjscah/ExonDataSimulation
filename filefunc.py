import re
from collections import Iterator

def text_iterator(texts):
    if isinstance(texts,str):
        while(True):
            yield texts
    elif isinstance(texts,(tuple,list)):
        for x in texts:
            yield x
    elif isinstance(texts,Iterator):
        return texts
def equal_text(line,text,match='=='):
    if match=='==':
        return line.strip()==text
    if match=='re':
        return re.match(text,line)
def get_pos_text(filed,text,aim=0,match='==',**pos):
    '''return pos'''
    if 'pos' in pos:
        filed.seek(pos['pos'],0)
    line=filed.readline()
    while(line):
        if equal_text(line,text,match):
            if aim==0:
                return filed.tell()-len(line)-1
            elif aim==1:
                return filed.tell()
        line=filed.readline()
def get_line_text(filed,text,aim=0,match='==',**pos):
    '''return line'''
    if 'pos' in pos:
        filed.seek(pos['pos'],0)
    line=filed.readline()
    while(line):
        if equal_text(line,text,match):
            if aim==0:
                return line
            elif aim==1:
                return filed.readline()
        line=filed.readline()

def get_column_row(file,text,aim=0,match='==',**pos):
    '''
    caculate column of aimed row ,text is search content,
    aim=0 get current row ,aim=1 get next row
    match='==' is complete match, match='re' only match head ingore tail
    '''
    with open(file,'r') as f:
        line=get_line_text(f,text,aim,match,**pos)
        return len(line.strip())

def write_paragraph_content(filed,writed,stext,etext,saim=0,eaim=0,smatch='==',ematch='=='):
    '''
    stext is begin row ,text2 is and row,
    stext,etext can be list or str or iretator
    aim is whether including this text
    if saim=0 ,include begin row
    if eaim=0 ,not include end row
    '''
    line=get_line_text(filed,stext,eaim,smatch)
    print(line.strip(),end="\r")
    while(line):
        writed.write(line)
        line=writed.readline()
        if equal_text(line,etext,ematch):
            if eaim==1:
                writed.write(line)
            elif eaim==0:
                filed.seek(filed.tell()-len(line)-1,0)
            break
    if line:
        raise FloatingPointError('file end')
def write_content(file,write,stexts,etexts,saim=0,eaim=0,smatch='==',ematch='=='):
    stexts=text_iterator(stexts)
    etexts=text_iterator(etexts)
    with open(file,'r') as filed:
        with open(write,'w') as writed:
            try:
                for stext,etext in zip(stexts,etexts):
                    write_paragraph_content(filed,writed,stext,etext,saim,eaim,smatch,ematch)
            except FloatingPointError as fileend:
                print(fileend)


def write_column(writed,text,column):
    ''' write text to a opend file with coulum words every row'''
    for x in range(0,len(text),column):
        writed.write(text[x:x+column]+'\n')
def num_positive(num):
    if num<0:
        return 0
    else:
        return num