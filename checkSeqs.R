source('readSeqs.R')

#align on LANL with HMM
refName<-'B.FR.83.HXB2_LAI_IIIB_BRU_K03455'
aligned<-read.fa('work/alignHmm_2019-05-07/AlignSeq.nt.fasta.gz')
aligned$full5<-nchar(degap(substring(aligned$seq,2000,4000)))>1000
aligned$five<-nchar(degap(substring(aligned$seq,5000,6000)))>500
aligned$env<-nchar(degap(substring(aligned$seq,7000,8000)))>500
aligned$full3<-nchar(degap(substring(aligned$seq,10000,10500)))>300
aligned$seqType<-NA
aligned[aligned$full5&aligned$full3&aligned$five&aligned$env,'seqType']<-'full'
aligned[!aligned$full5&aligned$full3&aligned$five&aligned$env,'seqType']<-'threeHalf'
aligned[!aligned$full5&!aligned$full3&!aligned$five&aligned$env,'seqType']<-'env'
firstNonGap<-regexpr('[^-]',aligned$seq)
lastNonGap<-regexpr('[^-]-*$',aligned$seq)
if(any(!inCohort$newName %in% aligned$name))stop('Missing name in alignment')
if(any(!aligned$name %in% c(refName,inCohort$newName)))stop('New name found in alignment')
pdf('out/aligned_seqs.pdf',width=20,height=10)
dnaplotr::plotDNA(aligned$seq)
dnaplotr::plotDNA(aligned$seq,groups=ifelse(is.na(aligned$seqType),'bad',aligned$seqType))
withAs(aligned=aligned[is.na(aligned$seqType),],dnaplotr::plotDNA(aligned$seq,groups=apply(aligned[,c('full5','five','env','full3')],1,function(xx)paste(as.numeric(xx),collapse=''))))
#abline(v=c(1000,5500,6500,10000,10500))
dev.off()
dnaplotr::plotDNA(aligned$seq[aligned$name %in% inCohort[inCohort$name %in% corrections$name,'newName']])

select<-aligned$name %in% c(refName,'MM15.PLAS.SGS.00037.004.01');dnaplotr::plotDNA(aligned[select,'seq'],groups=aligned[select,'name'])
dnaplotr::plotDNA(substring(aligned[c(1,709),'seq'],10000),groups=aligned[c(1,709),'name'])

allProts<-readFaDir('work/alignHmm_2019-04-23','fasta.gz')
allProts<-allProts[!allProts$file %in% c('AlignSeq.nt.fasta.gz','LTR.nt.fasta.gz'),]
allProts$aa<-suppressWarnings(dna2aa(dnar::degap(allProts$seq)))
allProts$prot<-sub('\\..*','',allProts$file)
allProts$finalStop<-grepl('X',substring(allProts$aa,nchar(allProts$aa)-20))
allProts$earlyStop<-grepl('X',substring(allProts$aa,1,nchar(allProts$aa)-20))
allProts$hasM<-grepl('M',substring(allProts$aa,1,5))
hxb2Length<-dnar::withAs(xx=allProts[allProts$name==refName,],tapply(nchar(dnar::degap(xx$seq)),xx$prot,c))
hxb2Seq<-dnar::withAs(xx=allProts[allProts$name==refName,],tapply(xx$seq,xx$prot,c))
allProts$prop<-nchar(dnar::degap(allProts$seq))/hxb2Length[allProts$prot]
allProts$type<-sapply(strsplit(allProts$name,'\\.'),'[',3)
allProts$oldName<-sapply(allProts$name,function(xx)ifelse(xx==refName,'HXB2',inCohort[inCohort$newName==xx,'name']))
probs<-allProts[allProts$earlyStop & allProts$hasM&allProts$prop>.9,]
t(table(sapply(strsplit(probs$name,'\\.'),'[',3),probs$prot))
probEnv<-allProts[allProts$prot=='Env'&allProts$name %in% c(probs[probs$type=='ISO'&probs$prot=='Env','name'],refName),]
zz<-dnaplotr::seqSplit(probEnv$seq)
apply(zz,1,function(xx)cumsum(xx!='-'))
findFrameShift<-function(seq,ref,gapChars=c('-','.'),minLength=21){
  seqN<-cumsum(!strsplit(seq,'')[[1]] %in% gapChars)
  refN<-cumsum(!strsplit(ref,'')[[1]] %in% gapChars)
  offset<-(seqN-refN) %% 3
  runs<-dnar::binary2range((offset!=0))
  runs[runs$end-runs$start>minLength,]$start[1]
}
probs$shift<-mapply(findFrameShift,probs$seq,hxb2Seq[probs$prot])
probs$firstXAA<-regexpr('X',probs$aa)
probs$firstX<-mapply(dnar::noGap2Gap,probs$seq,(probs$firstXAA-1)*3+1)
probs$shiftSeq<-sprintf('%s vs %s',substring(probs$seq,probs$shift-10,probs$shift+10),substring(hxb2Seq[probs$prot],probs$shift-10,probs$shift+10))
#,substring(dnar::dna2aa(degap(substring(probs$seq,probs$firstX,probs$firstX+100))),1,1)
probs$stopSeq<-sprintf('%s vs %s',substring(probs$seq,probs$firstX-3,probs$firstX+5),substring(hxb2Seq[probs$prot],probs$firstX-3,probs$firstX+5))
write.csv(probs[,c('name','oldName','type','prot','shift','firstX','shiftSeq','stopSeq','seq','aa')],'out/problems.csv',row.names=FALSE)

findFrameShift(probs[6,'seq'],hxb2Seq[probs$prot][6])
#substring(probEnv$seq[1:2],1225,1250)


inCohort$nNt<-nchar(dnar::degap(inCohort$seq))
hist(inCohort$nNt,breaks=200)
bad<-list(
  'shorter2000'=inCohort[inCohort$nNt<2000,'newName'],
  'around4000'=inCohort[inCohort$nNt>3800&inCohort$nNt<4200,'newName'],
  '6to8000'=inCohort[inCohort$nNt>6000&inCohort$nNt<8200,'newName']
)

pdf('out/problem_seqs.pdf',width=20)
par(mar=c(4,3,1,12))
#select<-(!aligned$full&!aligned$five&!aligned$env)|aligned$name==refName;dnaplotr::plotDNA(aligned[select,'seq'],groups=aligned[select,'name'])
#select<-(aligned$full&aligned$five&!aligned$env)|aligned$name==refName;dnaplotr::plotDNA(aligned[select,'seq'],groups=aligned[select,'name'])
#select<-(!aligned$full&aligned$five&!aligned$env)|aligned$name==refName;dnaplotr::plotDNA(aligned[select,'seq'],groups=aligned[select,'name'])
#dnaplotr::plotDNA(aligned[c(1,709),'seq'],groups=aligned[c(1,709),'name'])
#dnaplotr::plotDNA(aligned$seq[aligned$name %in% unlist(bad)],groups=aligned$name[aligned$name %in% unlist(bad)])
dnaplotr::plotDNA(aligned$seq[lastNonGap>=lastNonGap[1]][c(1,2)])
dev.off()

badSelect<-(is.na(aligned$seqType)&!grepl('SGS',aligned$name))|lastNonGap>=lastNonGap[1]+10
zz<-aligned[badSelect,c('name','full5','five','env','full3','seq')]
zz$oldName<-sapply(zz$name,function(xx)inCohort[inCohort$newName==xx,'name'])
zz$reason<-ifelse(!zz$five,'Missing upstream',ifelse(!zz$env,'Missing env',ifelse(!zz$full3,'Missing downstream',ifelse(lastNonGap[badSelect]>=lastNonGap[1]+10,"Extra on 3' end",NA))))
write.csv(zz[,c('oldName','reason')],'out/weirdCoverage.csv',row.names=FALSE)

pdf('out/problem_seqs.pdf',width=22)
par(mar=c(4,3,1,15))
#select<-(!aligned$full&!aligned$five&!aligned$env)|aligned$name==refName;dnaplotr::plotDNA(aligned[select,'seq'],groups=aligned[select,'name'])
#select<-(aligned$full&aligned$five&!aligned$env)|aligned$name==refName;dnaplotr::plotDNA(aligned[select,'seq'],groups=aligned[select,'name'])
#select<-(!aligned$full&aligned$five&!aligned$env)|aligned$name==refName;dnaplotr::plotDNA(aligned[select,'seq'],groups=aligned[select,'name'])
#dnaplotr::plotDNA(aligned[c(1,709),'seq'],groups=aligned[c(1,709),'name'])
#dnaplotr::plotDNA(aligned$seq[aligned$name %in% unlist(bad)],groups=aligned$name[aligned$name %in% unlist(bad)])
dnaplotr::plotDNA(aligned$seq,groups=ifelse(is.na(aligned$seqType),'bad',aligned$seqType))
dnaplotr::plotDNA(c(zz$seq,aligned[1,'seq']),groups=c(zz$oldName,'HXB2'))
dev.off()


any(!unlist(bad) %in% zz$name & !grepl('SGS',unlist(bad)))
xx<-dnar::read.fa('tmpMM55/Consensus sequences cleaned.fasta')
#xx<-xx[nchar(gsub('[ACTGN]','',xx$seq))<100,]
mm55<-seqs[grep('MM55',seqs$name),]
dists<-levenR::leven(xx$seq,mm55$seq,nThreads=20)
write.fa(c(seqs[grep('MM55',seqs$name),'name'],xx$name),c(seqs[grep('MM55',seqs$seq),'name'],xx$seq),'tmpMM55/combined.fa')
al<-levenR::levenAlign(c(xx$seq,mm55$seq),hxb2,substring2=TRUE)
  #read.fa('tmpMM55/align20190507/AlignSeq.nt.fasta')
pdf('test.pdf');dnaplotr::plotDNA(al[[2]]);dev.off()
cuts<-substring(al[[2]],6100,8500)
splits<-dnar::seqSplit(cuts)
times<-gsub('MM55\\.([0-9]+).*','\\1',mm55$name)
xx$time<-gsub('MM55-([0-9]+)-.*','\\1',xx$name)
dists<-levenR::leven(dnar::degap(cuts[1:nrow(xx)]),dnar::degap(cuts[nrow(xx)+1:nrow(mm55)]),nThreads=50)
hams<-do.call(rbind,lapply(1:nrow(xx),function(ii){
  good<-splits[ii,] %in% c('A','C','T','G')
  hams<-sapply(nrow(xx)+1:nrow(mm55),function(jj)sum(splits[ii,good]!=splits[jj,good]))
  return(hams)
}))
out<-t(apply(hams,1,function(xx)tapply(xx,times,min)))
rownames(out)<-sub('MM55-([^-]+-[^-]+)-.*','\\1',xx$name)

hams2<-do.call(rbind,lapply(1:nrow(xx),function(ii){
  hams<-sapply(1:nrow(xx),function(jj){
    good<-splits[ii,] %in% c('A','C','T','G')&splits[jj,] %in% c('A','C','T','G')
    sum(splits[ii,good]!=splits[jj,good])
  })
  return(hams)
}))

data.frame('id'=sub('MM55-([^-]+-[^-]+)-.*','\\1',xx$name),'ambig'=nchar(gsub('[ACGTN-]+','',xx$seq)))


xx<-read.fa('~/Downloads/MM39 MM40 isolate documents To Scott.fasta')
xx$sample<-sub('\\.([0-9])$','.0\\1',sub('UK85','MM39',sub('_.*','',sub('^_','',sub('_con.*','',sub('U[0-9]+_','',sub('-.*','',xx$name)))))))
xx[xx$sample=='MM39.10','seq'] %in% sga[sga$pat=='MM39'&sga$timeString=='10','seq']
xx$name %in% inCohort[,'name']
withAs(xxxx=inCohort[inCohort$mm %in% c('MM39','MM40')&!is.na(inCohort$isolate),'name'], xxxx[!xxxx %in% xx$name])


lapply(sga[sga$pat=='MM39'&sga$timeString=='10','seq'],function(xxx)grep(xxx,xx[xx$sample=='MM39.10','seq'],perl=TRUE))

