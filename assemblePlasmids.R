r1s<-list.files('plasmids/Marvin_2-9-18/','R1_001.trimmed.fastq',full.name=TRUE)
r2s<-sub('R1_001\\.','R2_001.',r1s)

outs<-file.path('plasmids/soap',sub('_L001_R1_001\\..*','',basename(r1s)))
logs<-sprintf('%s.log',outs)

configs<-sapply(basename(r1s),tempfile)
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
'
mapply(function(r1,r2,config)writeLines(sprintf(configText,r1,r2),config),r1s,r2s,configs)
cmds<-sprintf('soapdenovo2-127mer all -s %s -o %s -K 67 -F -V -R -E 2>%s',configs,outs,logs)
out<-sapply(cmds,system,intern=TRUE)
zz<-sapply(logs,readLines)
sapply(zz,function(xx)print(xx[grep('length',xx)]))
