library(xlsx)
library(dnar)
#iyer<-read.xlsx("data/Final all data master 101116 copy 2.xlsx",sheetIndex=1,stringsAsFactors=FALSE)
if(!exists('hiv'))source('../hivPair/readData.R',chdir=TRUE)
if(!exists('dat'))source('readNewData.R')
seqs<-data.frame(
  'id'=c(hiv$Renamed,dat$id),
  'study'=rep(c('iyer','gondim'),c(nrow(hiv),nrow(dat))),
  'alphaIc50'=c(hiv$IFNa2.Pooled.Donor.cells.IC50..pg..ml,dat$ic50),
  'alphaVres'=c(hiv$Residual.Pooled.Donor.cells..1500U..UT,dat$vres),
  'betaIc50'=c(hiv$IFNbeta.Pooled.Donor.cells.IC50..pg.ml,dat$beta),
  'betaVres'=c(hiv$vres,dat$betaVres),
  'patID'=c(hiv$baseName,dat$pat),
  'pair'=c(hiv$Pair.ID,dat$pat),
  'donor'=c(hiv$donor,rep(NA,nrow(dat))),
  'gender'=c(hiv$Gender,dat$Gender),
  'select'=c(hiv$select,rep('UT',nrow(dat))),
  'fluid'=c(hiv$fluid,rep('PL',nrow(dat))),
  'subtype'=c(hiv$Subtype,rep('B',nrow(dat))),
  'time'=c(rep(NA,nrow(hiv)),dat$time),
  'bulk'=c(rep(FALSE,nrow(hiv)),dat$bulk),
  'seq'=c(hiv$Sequence,dat$Sequence..Full.Length.or.3.half.),
  'replicativeCapacity'=c(hiv$Replicative.capacity.Pooled.Donor.cells.p24.d7,dat$replication),
  'envPerRT'=c(hiv$Env.RT,rep(NA,nrow(dat))),
  'infectivity'=c(hiv$Infectivity.RLU.pg.RT...T1249,rep(NA,nrow(dat))),
  'p24Release'=c(hiv$p24.release.No.IFN,rep(NA,nrow(dat))),
  stringsAsFactors=FALSE
)

seqs$isNa<-is.na(seqs$seq)|seqs$seq==''|seqs$seq=='NA'
dists<-withAs(xx=seqs[!seqs$isNa,'seq'],cacheOperation('work/combineDists.Rdat',cbind,levenR::leven(degap(xx),degap(xx[1]),substring2=TRUE,nThreads=3),levenR::leven(revComp(degap(xx)),degap(xx[1]),substring2=TRUE,nThreads=3)))
seqs$revComp[!seqs$isNa]<-dists[,1]>dists[,2]&dists[,2]<2000
seqs$raw<-degap(ifelse(seqs$revComp,revComp(seqs$seq),seqs$seq),c('*','.','-','\n'))
withAs(xx=seqs[!seqs$isNa,],write.fa(xx$id,seqs$raw,'combined/combined.fa'))
rownames(seqs)<-seqs$id


hmmer<-read.fa('combined/hmmer/AlignSeq.nt.fasta')
hxb2StartEnd<-noGap2Gap(hmmer[1,'seq'],c(1,nchar(degap(hmmer[1,'seq']))))
hmmer$trim<-substring(hmmer$seq,hxb2StartEnd[1],hxb2StartEnd[2])
coverRange<-range(which(apply(seqSplit(hmmer$trim)=='-',2,mean)<.25))
substring(hmmer$trim[-1],1,coverRange[1]-1)<-gsub('[^-]','-',substring(hmmer$trim[-1],1,coverRange[1]-1))
substring(hmmer$trim[-1],coverRange[2]+1)<-gsub('[^-]','-',substring(hmmer$trim[-1],coverRange[2]+1))
allAligns<-readFaDir('combined/hmmer/')
allAligns<-allAligns[!allAligns$file %in% c('AlignSeq.nt.fasta','LTR.nt.fasta'),]
allAligns$prot<-sub('.nt.fasta','',allAligns$file)

mafft<-read.fa('combined/mafft/AlignSeq.nt.fasta')
png('combined/hmmer.png',width=6000,height=3000)
dnaplotr::plotDNA(hmmer$seq)
dev.off()
png('combined/hmmerTrim.png',width=6000,height=3000)
dnaplotr::plotDNA(hmmer$trim)
dev.off()

png('combined/mafft.png',width=6000,height=3000)
dnaplotr::plotDNA(mafft$seq)
#axis(2,1:nrow(xx),mgp=c(-2,-3,-4))
dev.off()

rownames(hmmer)<-hmmer$name
seqs$align<-hmmer[seqs$id,'trim']
if(any(is.na(seqs$align)&!seqs$isNa))stop('Lost sequence during alignment')

for(ii in unique(allAligns$prot)){
  tmp<-allAligns[allAligns$prot==ii,]
  rownames(tmp)<-tmp$name
  seqs[,ii]<-tmp[seqs$id,'seq']
  if(any(is.na(seqs[,ii])&!seqs$isNa))stop('Lost genes sequence during alignment')
}


out<-seqs[!seqs$isNa,!colnames(seqs) %in% c('isNa','seq','revComp')]
tmp<-out[1,]
tmp[,]<-NA
tmp$id<-hmmer[1,'name']
tmp$raw<-degap(hmmer$seq[1])
tmp$align<-hmmer$trim[1]
for(ii in unique(allAligns$prot))tmp[,ii]<-allAligns[allAligns$name==tmp$id&allAligns$prot==ii,'seq']
out<-rbind(out,tmp)

for(ii in unique(allAligns$prot)){
  message(ii)
  seqMat<-seqSplit(out[,ii])
  ref<-seqMat[out$id==tmp$id,]
  aas<-rep('',nrow(seqMat))
  frame<-1
  refPos<-which(ref!='-')
  dummy<-paste(rep('-',1000),collapse='')
  for(jj in seq(3,length(refPos),3)){
    cat('.')
    subset<-(ifelse(jj-3==0,0,refPos[jj-3])+1):refPos[jj]
    subSeqs<-seqMat[,subset]
    subAA<-dna2aa(degap(apply(subSeqs,1,paste,collapse='')),warn=FALSE)
    n<-nchar(subAA)
    subAA<-sprintf('%s%s',subAA,substring(dummy,1,max(n)-n))
    if(any(nchar(subAA)!=nchar(subAA[1])))browser()
    aas<-sprintf('%s%s',aas,subAA)
  }
  out[,sprintf('%sAA',ii)]<-aas
}

write.csv(out,gzfile('combined/combinedIC50.csv.gz'),row.names=FALSE)
withAs(zz=seqs[!seqs$isNa,],write.fa(zz$id,zz$align,'combined/trimmed.fa'))

xx<-read.csv('combined/IFNbeta - DR.csv',header=FALSE,stringsAsFactors=FALSE)[2:12,2:3]
yy<-read.csv('combined/IFNbeta - DR.csv',header=FALSE,stringsAsFactors=FALSE)[15:23,2:3]
pdf('combined/ifnCompare.pdf')
plot(c(xx[,1],yy[,1]),c(xx[,2],yy[,2]),col=rep(c('red','blue'),c(nrow(xx),nrow(yy))),log='x',xaxt='n',las=1,ylab='p24',xlab='IFNb')
dnar::logAxis(1)
abline(h=50,lty=2)
x50<-approx(xx[,2],xx[,1],50)$y
y50<-approx(yy[,2],yy[,1],50)$y
abline(v=x50,col='red')
abline(v=y50,col='blue')
dev.off()

probs<-lapply(unique(allAligns$prot),function(prot){
  isNa<-is.na(out[,prot])|out[,prot]==''|!grepl('[ACTG]',out[,prot])
  aa<-dna2aa(degap(out[!isNa,prot]),warn=FALSE)
  probs<-substring(aa,nchar(aa))!='X'|grepl('X',substring(aa,1,nchar(aa)-1))
  out[!isNa,][probs,c('id',prot)]
})
names(probs)<-unique(allAligns$prot)
