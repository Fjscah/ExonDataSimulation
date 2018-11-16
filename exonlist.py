import re,os,sys,time
from exon import Exon
from setting import INSERT,ROW_STEP,ROW_COLUMN,CHIP_LEN
print("default:","INSERT length=",INSERT,)
class WholeExon(object):
    # coulum is the numbeer of dNTP in a row from filein, step is usually one byte which stand for '\n'
    def __init__(self,fin,fout,column=ROW_COLUMN,step=ROW_STEP):
        self.filein=fin
        self.fileout=fout
        self.pos=0
        self.chr=0
        self.column=column
        self.step=step
    # write sequence of certain exon according to it's begin pos and it's length
    def exonseq(self,pos,length):
        self.filein.seek(pos,0)
        sys.stdout.write(str(self.filein.tell()))
        sys.stdout.write('-->')
        line=self.filein.readline().strip()
        # sp is set for alignment, in other word, write the same number of dNTP in every row
        length=length-len(line)
        if length>=0:
            self.fileout.write(line)
            sp=self.column-len(line)
            line=self.filein.readline().strip()
            while(line and length>0):
                l=len(line)
                length=length-l
                if length>=0:
                    self.fileout.write(line[0:sp])
                    self.fileout.write("\n")
                    self.fileout.write(line[sp:])
                elif length<0:
                    index=length+l
                    if index>sp:
                        self.fileout.write(line[0:sp])
                        self.fileout.write("\n")
                        self.fileout.write(line[sp:index])
                    else :
                        self.fileout.write(line[0:index])
                    break
                line=self.filein.readline().strip()
        else:
            self.fileout.write(line[0:length+len(line)])

        self.fileout.write('\n')
        sys.stdout.write(str(self.filein.tell())+"     ")
        sys.stdout.write('\r')
    def lentobyte(self,length):
        raw=length//self.column
        return raw*(self.column+self.step)+length %self.column-1
    # according exons informatin write its sequence ; the exon must be in order from 1 to 24
    def wholeexonseq(self,exons,INSERT=INSERT):
        exon=next(exons)
        #m=3
        while(exon):
            tempchr=exon.chr
            if not tempchr==self.chr:
                print("exon end:",self.filein.tell(),"     ")
                self.chr=tempchr
                line=self.filein.readline()
                row=0
                while(line):
                    row+=1
                    sys.stdout.write("%d"%row)
                    sys.stdout.write('\r')
                    if re.match(r'>NC[\s\w]*',line):
                        print("chr",self.chr,"pos",self.filein.tell(),line.strip())
                        # change to next chromosome 
                        self.pos=self.filein.tell()
                        row=0
                        break
                    try:
                        line = self.filein.readline()
                    except StopIteration:
                        line =0
                    else:
                        pass

            # write exon title
            #print('write exon ',exon.getexon_info())
            self.fileout.write("%s\t%d\n" % (exon.getexon_info(),INSERT))
            po = self.pos+ self.lentobyte(exon.begin-INSERT)
            length=exon.length()+2*INSERT
            # write exon sequence
            self.filein.seek(self.pos,0)
            try:
                line = self.filein.readline()
            except StopIteration:
                line =0
            else:
                pass
            self.exonseq(po,length)
            try:
                exon=next(exons)
            except StopIteration:
                exon =0
            else:
                pass
            
            #m=m-1
            #if m==0:
                #break
        print("dowm.")


def exonsgenerater(f):
    '''according file content to generate Exon object's itertator'''
    l_exons=[]
    eid=0
    echr=0
    for line in f.readlines():
        m=get_exoninfo(line)
        if not m:
            continue
        if m.gene_id!=eid:
            if len(l_exons)==1:
                yield Exon(echr,l_exons[0][0],l_exons[0][1],eid)
            elif len(l_exons)>1:
                lorder=sorted(l_exons)
                size =len(lorder)
                x=0
                y=1
                end=lorder[x][1]
                while(x<=size-2 and y<=size-1):
                    if end+1>=lorder[y][0]:
                        end =max(end,lorder[y][1])
                        y+=1
                    elif end+1<lorder[y][0]:
                        #print("yield",eid)
                        yield Exon(echr,lorder[x][0],end,eid)
                        end=lorder[y][1]
                        x,y=y,y+1
                yield Exon(echr,lorder[x][0],lorder[-1][1],eid)
            l_exons=[]
            echr=m.chr
            eid=m.gene_id
        l_exons.append((m.begin,m.end))
    if len(l_exons)==1:
        yield Exon(echr,l_exons[0][0],l_exons[0][1],eid)
    elif len(l_exons)!=0:
        lorder=sorted(l_exons)
        size =len(lorder)
        x=0
        y=1
        end=lorder[x][1]
        while(x<=size-2 and y<=size-1):
            if end>=lorder[y][0]:
                end =max(end,lorder[y][1])
                y+=1
            elif end<lorder[y][0]:
                #print("yield",eid)
                yield Exon(echr,lorder[x][0],end,eid)
                end=lorder[y][1]
                x,y=y,y+1
    
def exonsettle(l,echr,eid):
    lorder=sorted(l)
    size =len(lorder)
    x=0
    y=1
    end=lorder[x][1]
    while(x<=size-2 and y<=size-1):
        if end>=lorder[y][0]:
            end =max(end,lorder[y][1])
            y+=1
        elif end<lorder[y][0]:
            #print("yield",eid)
            yield Exon(echr,lorder[x][0],end,eid)
            end=lorder[y][1]
            x,y=y,y+1
    l=[]
                  
def get_exoninfo(line):
    '''according ever row to genarater Exon object'''
    m = re.match(r'[\s\w]*chr([\w]+)\s*[\s\w]*exon[\s\w]*?(\d+)[\s\w]*?(\d+)[\s\w\-\.\+]*gene_id[\s\w]*"([\w\.]+)"[\s\w]*', line)
    if m :
        if int(m.group(3))-int(m.group(2))>=CHIP_LEN:
            return Exon(m.group(1),int(m.group(2)),int(m.group(3)),m.group(4))
        else:
            return None
    else:
        return None

def getrowcontent(r1,w1,text1,text2):
    state=False
    with open(r1,'r') as r:
        with open (w1,'w') as w:
            for line in r.readlines():
                if re.match(text1, line):
                    print(line.strip()[0:60],end="\r")
                    w.write(line)
                    state=True
                elif state==True and (not re.match(text2, line)):
                    w.write(line)
                elif re.match(text2, line) and (not re.match(text1, line)):
                    state=0
    print ("down.                                                               ")


if __name__ == '__main__':
    first='GRCh38_latest_genomic.txt'
    second='gencode.v29.annotation.txt'
    if len(sys.argv)==3:
        script, first, second = sys.argv
    # get chromosome sequence from grch38 genomic
    print("get chromosome sequence from grch38 genomic...")
    getrowcontent('GRCh38_latest_genomic.txt','NCseq.txt',
        r'>NC[\w\s]*',r'>[\w\s]*')
    # get exons annotation from v29 annotation
    print("get exons annotation from v29 annotation...")
    getrowcontent('gencode.v29.annotation.gtf',"ENSEMBL_exonposition.txt",
        r'[\w\s]*ENSEMBL[\w\s]*exon[\w\s]*',r'[\w\s]*')


    # get exon sequence ; write to exonlist.txt
    print("get exon initial list...")
    with open("ENSEMBL_exonlist.txt" ,'w') as w:
        w.write("#chr\t\tbegin\tend\tgeneid\t\t\tlength\tINSERT\n")
        with open ("ENSEMBL_exonposition.txt",'r') as r1:
            with open ("NCseq.txt",'r') as r2:
                exons =exonsgenerater(r1)
                wholeexon = WholeExon(r2,w)
                wholeexon.wholeexonseq(exons)





                    
                    




                    

