import re


def equal_text(line, text, match='=='):
    if match == '==':
        return line.strip() == text
    elif match == 're':
        return re.match(text, line)


def get_pos_text(filed, text, aim=0, match='==', **pos):
    '''return pos'''
    if 'pos' in pos:
        filed.seek(pos['pos'], 0)
    line = filed.readline()
    while(line):
        if equal_text(line, text, match):
            if aim == 0:
                return filed.tell()-len(line)-1
            elif aim == 1:
                return filed.tell()
        line = filed.readline()


def get_line_text(filed, text, aim=0, match='==', **pos):
    '''return line'''
    if 'pos' in pos:
        filed.seek(pos['pos'], 0)
    line = filed.readline()
    while(line):
        if equal_text(line, text, match):
            if aim == 0:
                return line
            elif aim == 1:
                return filed.readline()
        line = filed.readline()
    return None


def get_column_row(file, text, aim=0, match='==', **pos):
    '''
    caculate column of aimed row ,text is search content,
    aim=0 get current row ,aim=1 get next row
    match='==' is complete match, match='re' only match head ingore tail
    '''
    with open(file, 'r') as f:
        line = get_line_text(f, text, aim, match, **pos)
        if line:
            return len(line.strip())
        else:
            return 0


def write_column(writed, text, column):
    ''' write text to a opend file with coulum words every row'''
    for x in range(0, len(text), column):
        writed.write(text[x:x+column]+'\n')


def num_positive(num):
    if num < 0:
        return 0
    else:
        return num
