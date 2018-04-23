library(dnar)
if(FALSE){
align<-read.fa('Marvin_1-17-18-1/U508- assembled to MM33d12_BLAST 2.fasta.gz')
align$seq<-sub('N','-',align$seq)
nrow(align)
align$n<-nchar(gsub('N+','',align$seq))
}


align<-function(r1s,r2s,ref,out){
  system(sprintf('bwa index %s',ref))
  tmp<-tempfile()
  tmp2<-tempfile()
  cmd<-sprintf('echo "bwa mem -t 40 %s -k9 -r 1.1 -c 20000 <(zcat %s) <(zcat %s)> %s"|bash',ref,paste(r1s,collapse=' '),paste(r2s,collapse=' '),tmp)
  system(cmd)
  cmd<-sprintf('samtools view -bS %s>%s;samtools sort -f -m 5000000000 %s %s;samtools index %s',tmp,tmp2,tmp2,out,out)
  system(cmd)
  cmd<-sprintf('bamtoalign -s %s -e 10 %s|gzip>%s.fa.gz',ref,out,out)
  system(cmd)
}


multiAlign<-function(r1s,r2s,ref,out,nIter=5){
  message('Iterations left ',nIter)
  if(nIter>1)tmpOut<-tempfile()
  else tmpOut<-out
  align(r1s,r2s,ref,tmpOut)
  if(nIter>1){
    al<-read.fa(sprintf('%s.fa.gz',tmpOut))
    alMat<-t(apply(do.call(rbind,strsplit(al$seq,'')),1,function(xx){
      if(xx[1]=='-')xx[1:min(which(xx!='-'),length(xx))]<-'.'
      if(xx[length(xx)]=='-')xx[max(1,which(xx!='-')):length(xx)]<-'.'
      return(xx)
    }))
    newRef<-degap(paste(apply(alMat,2,function(xx)mostAbundant(xx[xx %in% c('A','C','T','G','-')])),collapse=''))
    tmpFa<-tempfile()
    write.fa('ref',newRef,tmpFa)
    newRef<-multiAlign(r1s,r2s,tmpFa,tmpOut,nIter-1)
  }
  #al<-read.fa(sprintf('%s.fa.gz',tmpOut))
  return(newRef)
}
multiAlign(fastqs[grep('_R1_',fastqs)],fastqs[grep('_R2_',fastqs)],'align/WEAU.fasta','align/test.bam')
cmd<-sprintf('samtools view -bS align/test.bam>%s;samtools sort -f -m 5000000000 %s align/test.bam;samtools index align/test.bam',tmp,tmp)
system(cmd)
system('bamtoalign -s align/WEAU.fasta -e 10 align/Marvin_sequencing_UK85.2/align.bam|gzip>align/Marvin_sequencing_UK85.2/align.fa.gz')


setwd('align')
system('bwa index WEAU.fasta')
fastq1<-list.files('Marvin_1-17-18-1','_R1_001.trimmed.fastq',full.name=TRUE)
sams<-file.path('work',sub('_R1_001.trimmed.fastq','.sam',basename(fastq1)))
bams<-file.path('work',sub('_R1_001.trimmed.fastq','.bam',basename(fastq1)))
#need to cutadapt?
cmds<-sprintf('echo "bwa mem -t 40 align/WEAU.fasta %s %s > %s"|bash',fastq1,sub('R1_001.trimmed.fastq','R2_001.trimmed.fastq',fastq1),sams)

sapply(cmds,system)

tmp<-tempfile()
cmds<-sprintf('samtools view -bS %s>%s;samtools sort -f -m 5000000000 %s %s;samtools index %s',sams,tmp,tmp,bams,bams)
sapply(cmds,system)

cmd<-sprintf('echo "bwa mem -t 40 align/WEAU.fasta -k15 -r 1.1 -c 20000 -p <(zcat align/Marvin_sequencing_UK85.2/3_documents_from_READS.fastq.gz) > align/Marvin_sequencing_UK85.2/align.sam"|bash')
system(cmd)
cmd<-sprintf('samtools view -bS align/Marvin_sequencing_UK85.2/align.sam>%s;samtools sort -f -m 5000000000 %s align/Marvin_sequencing_UK85.2/align.bam;samtools index align/Marvin_sequencing_UK85.2/align.bam',tmp,tmp)
system(cmd)
system('bamtoalign -s align/WEAU.fasta -e 10 align/Marvin_sequencing_UK85.2/align.bam|gzip>align/Marvin_sequencing_UK85.2/align.fa.gz')

system('bwa index align/acuteRef.fa')
fastqs<-list.files('align/raw','fastq.gz$',full.names=TRUE)
cmd<-sprintf('echo "bwa mem -t 40 align/acuteRef.fa -k15 -r 1.1 -c 20000 <(zcat %s) <(zcat %s)> align/raw/align.sam"|bash',paste(fastqs[grep('_R1_',fastqs)],collapse=' '),paste(fastqs[grep('_R2_',fastqs)],collapse=' '))
system(cmd)
cmd<-sprintf('samtools view -bS align/raw/align.sam>%s;samtools sort -f -m 5000000000 %s align/raw/align.bam;samtools index align/raw/align.bam',tmp,tmp)
system(cmd)
system('bamtoalign -s align/acuteRef.fa -e 10 align/raw/align.bam|gzip>align/raw/align.fa.gz')

zz<-read.fa('align/raw/align.fa.gz')
png('test3.png',width=4000,height=4000,res=250);plotDNA(zz$seq[order(regexpr('[^-]',zz$seq))]);dev.off()

xx<-read.fa('align/Marvin_sequencing_UK85.2/align.fa.gz')
png('test.png',width=4000,height=4000,res=250);plotDNA(xx$seq[order(regexpr('[^-]',xx$seq))]);dev.off()

yy<-read.fa('align/Marvin_sequencing_UK85.2/FOR SCOTT_UK85.2-P13D2.fasta.gz')
png('test2.png',width=4000,height=4000,res=250);plotDNA(yy$seq[order(regexpr('[^-]',xx$seq))]);dev.off()

if(!file.exists('align/weau.nin')){system('makeblastdb -in align/WEAU.fasta -title weau -dbtype nucl -out align/weau')}
cmd<-sprintf("zcat align/Marvin_sequencing_UK85.2/3_documents_from_READS.fastq.gz|awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}'|blastn -db align/weau -num_threads 50 -outfmt 5|gzip > align/Marvin_sequencing_UK85.2/align.blast.gz")
system(cmd)
system('echo "python align/blast2sam.py <(zcat align/Marvin_sequencing_UK85.2/align.blast.gz)>align/Marvin_sequencing_UK85.2/blast.sam"|bash')
cmd<-sprintf('samtools view -bS align/Marvin_sequencing_UK85.2/blast.sam>%s;samtools sort -f -m 5000000000 %s align/Marvin_sequencing_UK85.2/blast.bam;samtools index align/Marvin_sequencing_UK85.2/blast.bam',tmp,tmp)
system(cmd)
system('bamtoalign -s align/WEAU.fasta -e 10 align/Marvin_sequencing_UK85.2/blast.bam|gzip>align/Marvin_sequencing_UK85.2/blast.fa.gz')

xx<-read.fa('align/Marvin_sequencing_UK85.2/blast.fa.gz')
png('test.png',width=4000,height=4000,res=250);plotDNA(xx$seq[order(regexpr('[^-]',xx$seq))]);dev.off()

r1s<-list.files('align/raw/','_R1_.*.fastq.gz',full.name=TRUE)
r2s<-sub('R1_001\\.','R2_001.',r1s)
out<-'align/raw/soap'
log<-sprintf('%s.log',out)
config<-tempfile()
configText<- '
#maximal read length
max_rd_len=150
[LIB]
#average insert size
avg_ins=300
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#use only first X bps of each read
rd_len_cutoff=140
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=%s
q2=%s
q1=%s
q2=%s
q1=%s
q2=%s
'
configOut<-do.call(sprintf,c(list(configText),as.list(sort(c(r1s,r2s)))))
writeLines(configOut,config)

file.remove(list.files('align/raw/','^soap',full.name=TRUE))
cmd<-sprintf('soapdenovo2-127mer all -p 30 -s %s -o %s -K 67 -F -V -R -E 2>%s',config,out,log)
system(cmd,intern=TRUE)
zz<-readLines(log)
print(zz[grep('length',zz)])
scaffs<-read.fa('align/raw/soap.scafSeq')
print(nchar(scaffs$seq))
