class Exon(object):
    # stored exon info ,not include exon sequence
    geneid="None"
    num=0
    # creat exon number counter
    __count =0
    def __init__(self,chr,begin,end,gene_id,seq=None):
        self.__chr=chr
        self.__begin=begin
        self.__end=end
        self.__gene_id=gene_id
        self.seq=seq
        # print new gene ; this can be  ellipsed
        if Exon.geneid!=self.__gene_id:
            Exon.__count +=1
            Exon.geneid=self.__gene_id
            Exon.num=1
            self.num=Exon.num
            #print (Exon.__count,'chr',self.__chr,' - ',gene_id)
        else:
            Exon.num +=1
            self.num=Exon.num
    @property
    def chr(self):
        return self.__chr
    @property
    def begin(self): 
        return self.__begin
    @property
    def end(self):
        return self.__end
    @property
    def gene_id(self):
        return self.__gene_id
    def length(self):
        return self.__end-self.__begin
    def getexon_info(self):
        info="chr%s\t%d\t%d\t%s-%d\t\t%d\t" %(
                self.chr,self.begin,self.end,self.gene_id,self.num,self.length())
        return info
