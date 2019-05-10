import os
os.system("python ./simu.py -dep reference.depth 350")
os.system("python ./simu.py -reg")
os.system("python ./simu.py -mut")
os.system("python ./simu.py -read -R 111")
os.system("bwa mem ../fastq/GRCh38_latest_genomic.fna ../fastq/R1_111.fasrq ../fastq/R2_111.fastq > ../fastq/sam/444.sam")
os.system('samtool view -bS ../fastq/sam/444.sam > ../fastq/sam/444.bam')
os.system('samtool sort ../fastq/sam/444.bam > ../fastq/sam/444.sort.bam')
os.system('samtool depth ../fastq/sam/444.sort.bam > ../fastq/sam/444.depth')
