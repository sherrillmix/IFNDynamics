library(dnar)
if(FALSE){
align<-read.fa('Marvin_1-17-18-1/U508- assembled to MM33d12_BLAST 2.fasta.gz')
align$seq<-sub('N','-',align$seq)
nrow(align)
align$n<-nchar(gsub('N+','',align$seq))
}


setwd('align')
system('bwa index WEAU.fasta')
fastq1<-list.files('Marvin_1-17-18-1','_R1_001.trimmed.fastq',full.name=TRUE)
sams<-file.path('work',sub('_R1_001.trimmed.fastq','.sam',basename(fastq1)))
bams<-file.path('work',sub('_R1_001.trimmed.fastq','.bam',basename(fastq1)))
#need to cutadapt?
cmds<-sprintf('echo "bwa mem -t 40 WEAU.fasta %s %s > %s"|bash',fastq1,sub('R1_001.trimmed.fastq','R2_001.trimmed.fastq',fastq1),sams)

sapply(cmds,system)

tmp<-tempfile()
cmds<-sprintf('samtools view -bS %s>%s;samtools sort -f -m 5000000000 %s %s;samtools index %s',sams,tmp,tmp,bams,bams)
sapply(cmds,system)
