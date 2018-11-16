import random

"""
class exon_mutation():
    def __init__(self,snp_num,cnv_num):
        '''
        snp_num is number of snps,cnv_num is number of this mutation exon
        if snp_num>0,cnv_num>0 ,it stands for snp+cnv ,not inclding origin exon
        if sun_num=0,cnv_num>0, it stands for only cnv,not inclding origin exon
        if snp_num>0,cnv_num=0 ,it stands for only snp, lying on origin exon
        if snp_num=0,cnv_num=-1, it stands for deletion
        '''
        self.snp_num=snp_num
        self.cnv_num=cnv_num
"""
class Homolog_type():
    def __init__(self,exon_num=0,snp_num=0,cnv_num=0,dele_num=0,total_exon=0):
        '''
        exon_num stands for mutation exon number, 
        if one exon have several mutation ,it be counted several
        snp_num stands for snp number,not including cnv resulting in several snp
        delf_num need less than exon_num
        '''
        
        assert exon_num-cnv_num-dele_num<=0 ,'need : exon <= cnv+dele'
        assert exon_num>=dele_num,'need :exon >= dele'
        self.exon_num=exon_num
        self.snp_num=snp_num
        self.cnv_num=cnv_num
        self.dele_num=dele_num
        self.exon_mutations=[]
        self.total_exon=total_exon
        self.__creat_randommutaton(total_exon)
        self.__snp_mutations=self.__get_snps()
    def __creat_randommutaton(self,total_exon):
        '''
        len(snp),cnv,dele,illustration
        1.  5 6 1  origin snp +cnv
        2.  0 0 1  deletion
        3.  5 0 1  origin snp
        4.  0 6 1  invalid
        ---------------------
        5.  5 6 0  new snp +cnv
        6.  0 6 0  cnv
        7.  5 0 0  invalid
        8.  0 0 0  invalid
        1-4 origin exon not exist
        5-8 origin exon exist
        one exon can have new_snp+cnv and cnv at same time
        '''
        self.total_exon=total_exon
        if self.exon_num==0:
            return
        # in order of dele->snp->cnv creat randommutation
        self.exon_mutations=[]
        snp=self.snp_num
        cnv=self.cnv_num
        #fist dele
        dele_no=random.sample(range(1,total_exon),self.dele_num)
        if self.exon_num==self.dele_num:
            dele_num=self.dele_num-1
        else:
            dele_num=self.dele_num
        for x in range(dele_num):
            s=random.randint(0,snp)
            c=random.randint(0,cnv-(self.exon_num-self.dele_num))
            # no is exon  identifier
            while s==0 and c>0:
                s=random.randint(0,snp)
                c=random.randint(0,cnv-(self.exon_num-self.dele_num))
            snp-=s
            cnv-=c
            l_snp=[]
            if s>0:
                for y in range(s):
                    l_snp.append(round(random.random(),5))
            l_snp.sort()
            self.exon_mutations.append([dele_no[x],l_snp,c,1])
        #second cnv
        for x in range(self.exon_num-self.dele_num-1):
            s=random.randint(0,snp)
            c=random.randint(1,cnv-(self.exon_num-self.dele_num-x-1))
            snp-=s
            cnv-=c
            l_snp=[]
            if s>0:
                for x in range(s):
                    l_snp.append(round(random.random(),5))
            l_snp.sort()
            no=random.randint(1,total_exon)
            while(no in dele_no):
                no=random.randint(1,total_exon)
            self.exon_mutations.append([no,l_snp,c,0])
        s=snp
        c=cnv
        l_snp=[]
        if s>0:
            for x in range(s):
                l_snp.append(round(random.random(),5))
        l_snp.sort()
        if self.exon_num==self.dele_num:
            self.exon_mutations.append([dele_no[-1],l_snp,c,1])
        else:
            no=random.randint(1,total_exon)
            while(no in dele_no):
                no=random.randint(1,total_exon)
            self.exon_mutations.append([no,l_snp,c,0])
        self.exon_mutations.sort()
        
    @property
    def snp_mutations(self):
        return self.__snp_mutations

    def __get_snps(self):
        '''
        len(snp),cnv,dele,illustration
        1.  5 6 1  origin snp +cnv  *
        2.  0 0 1  deletion
        3.  5 0 1  origin snp       *
        4.  0 6 1  invalid
        ---------------------
        5.  5 6 0  new snp +cnv     *
        6.  0 6 0  cnv
        7.  5 0 0  invalid
        8.  0 0 0  invalid
        '''
        snp_mutations=[]
        for x in self.exon_mutations:
            if len(x[1])>0:
                snp_mutations.append(x)
        return snp_mutations

    def show_mutation(self):
        l=''
        if len(self.exon_mutations)==0:
            l+='Normal type\n'
        else:
            l+='Mutation type\n'
        for no,x,y,z in self.exon_mutations:
            if z==1:
                if len(x)>0 and y>0:
                    l+='NO.%d exon : origin snp + cnv\n'%no+\
                        'SNP=%d : %s ; CNV=%d\n'%(len(x),str(x),y)
                elif len(x)==0 and y==0:
                    l+='NO.%d exon : deletion' %no
                elif len(x)>0 and y==0:
                    l+='NO.%d exon : origin snp\n'%no+\
                        'SNP=%d : %s\n'%(len(x),str(x))
            elif z==0:
                if len(x)==0 and y>0:
                    l+='NO.%d exon : cnv\n'%no+\
                        'CNV=%d\n'%y
                elif len(x)>0 and y>0:
                    l+='NO.%d exon : new snp + cnv\n'%no+\
                        'SNP=%d : %s ; CNV=%d\n'%(len(x),str(x),y)
        return l
                
    def set_exon_mutations(self,exon_mutations):
        pass
    def get_snp_copy_num(self,mutainfo):
        if mutainfo in self.snp_mutations:
            return mutainfo[2]
        else:
            return 0

    def get_normal_copy_num(self,n):
        '''
        len(snp),cnv,dele,illustration,return
        1.  5 6 1  origin snp +cnv  0
        2.  0 0 1  deletion         0
        3.  5 0 1  origin snp       0
        4.  0 6 1  invalid
        ---------------------
        5.  5 6 0  new snp +cnv     1
        6.  0 6 0  cnv              6
        7.  5 0 0  invalid
        8.  0 0 0  invalid
        '''
        if len(self.exon_mutations)==0:
            return 1
        elif n<self.exon_mutations[0][0]:
            return 1
        elif n>self.exon_mutations[-1][0]:
            return 1
        for x in self.exon_mutations:
            if n==x[0]:
                if x[3]==1:
                    return 0
                elif len(x[1])==0:
                    return x[2]+1
        return 1

class Genome_type():
    def __init__(self,percent,homolog1,homolog2):
        self.percent=percent
        self.homologs=[homolog1,homolog2]
    def show_genetype(self):
        l=''
        l+='*'*5+'percent=%d'%self.percent+'*'*5+'\n'
        l+='homolog1:'
        l+=self.homologs[0].show_mutation()
        l+='homolog2:'
        l+=self.homologs[1].show_mutation()
        return l
    def collect_snp(self):
        l=[]
        for x in self.homologs:
            l=l+x.get_snps()
        return l


def set_mutation(total_exon):
    mutation_types=[]
    totol_percent=100
    i=1
    state=False
    num=int(input("how much mutation do you want to creat(including normal type):"))
    print('''OK. then enter mutation like 'exon_num snp_num cnv_num dele_num';
        if exon_num==0 ,this type is normal type;
        dele_num need be smaller than exon_num;
        eg. 5 9 6 2 is right ; 3 2 1 4 is wrong ; 0 0 0 0 is normal type''')
    while(i<=num):
        if i==num:
            print("please enter your no.%d(last) mutation percent : %d"%(i,totol_percent))
            percent=totol_percent
        else:
            while(state==False):
                percent=float(input("please enter your no.%d mutation percent : "%i))
                if totol_percent > percent:
                    state=True
                    totol_percent-=percent
                else:
                    print("error: percent is too bigger!")
                    continue
            state=False
        exon_num,snp_num,cnv_num,dele_num=map(int,input(
            "please enter your no.%d mutation homolog-1 : "%i).split())
        homolog1=Homolog_type(exon_num,snp_num,cnv_num,dele_num,total_exon)
        exon_num,snp_num,cnv_num,dele_num=map(int,input(
            "please enter your no.%d mutation homolog-2 : "%i).split())
        homolog2=Homolog_type(exon_num,snp_num,cnv_num,dele_num,total_exon)  
        mutation_types.append(Genome_type(percent,homolog1,homolog2))
        i+=1
    # show mutation setting
    with open ('mutation_settings_log.txt','a') as wlog:
        wlog.write('*'*10+'Show Genetype'+'*'*10+'\n')
        for x in mutation_types:
            wlog.write(x.show_genetype())
        wlog.write('*'*10+'Show Genetype'+'*'*10+'\n')
    show_mutation(mutation_types)
    return mutation_types


def show_mutation(mutation_types):
    print('*'*10+'Show Genetype'+'*'*10)
    for x in mutation_types:
        print(x.show_genetype())
    print('*'*10+'Show Genetype'+'*'*10)
    



        
