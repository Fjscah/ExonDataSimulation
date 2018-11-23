import re

def get_line_head(filed,text,match='=='):
    if match=='==':
        line=filed.readline()
        while(line):
            if line.strip()==text:
                return filed.tell()-len(line)-1
            line=filed.readline()
    elif match=='re':
        line=filed.readline()
        while(line):
            if re.match(text,line):
                return filed.tell()-len(line)-1
            line=filed.readline()
def get_line_next(filed,text,match='=='):
    if match=='==':
        line=filed.readline()
        while(line):
            if line.strip()==text:
                return filed.tell()
            line=filed.readline()
    elif match=='re':
        line=filed.readline()
        while(line):
            if re.match(text,line):
                return filed.tell()
            line=filed.readline()
def get_content(file1,w1,text1,text2):
    '''text1 is begin row ,text2 is end row'''
    state=False
    if isinstance(text1, (tuple,list)):
        texts= (x for x in text1 )
        with open(file1,'r') as r:
            with open (w1,'w') as w:
                try:
                    text=next(texts)
                    for line in r.readlines():
                        if line.strip()==text:
                            print(line.strip(),end="\r")
                            w.write(line)
                            state=True
                        elif state==True and (not re.match(text2, line)):
                            w.write(line)
                        elif re.match(text2, line) and (not line.strip()==text):
                            state=False
                            text=next(texts)
                except StopIteration:
                    pass
    elif isinstance(text1, str):
        with open(file1,'r') as r:
            with open (w1,'w') as w:
                for line in r.readlines():
                    if re.match(text1,line):
                        print(line.strip()[:40],end="\r")
                        w.write(line)
                        state=True
                    elif state==True and (not re.match(text2, line)):
                        w.write(line)
                    elif re.match(text2, line) and (not re.match(text1, line)):
                        state=False
def get_row_column(file,text,aim=0,match='=='):
    '''
    caculate column of aimed row ,text is search content,
    aim=0 get current row ,aim=1 get next row
    match='==' is complete match, match='re' only match head ingore tail
    '''
    with open(file,'r') as f:
        if aim==0:
            f.seek(get_line_head(f,text,match),0)
            line=f.readline()
            return len(line.strip())
        elif aim==1:
            f.seek(get_line_next(f,text,match),0)
            line=f.readline()
            return len(line.strip())
def write_column(writed,text,column):
    ''' write text to a opend file with coulum words every row'''
    for x in range(0,len(text),column):
        writed.write(text[x:x+column]+'\n')
def num_positive(num):
    if num<0:
        return 0
    else:
        return num