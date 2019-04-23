if(!exists('mmLookup'))source('readNewData.R')

hxb2<-read.fa('work/hxb2.fa')$seq
cohortRegex<-paste(c(sub('EJ','UK',ejLookup),ejLookup,mmLookup),collapse='|')
seqs<-dnar::readFaDir('seqs')
corrections<-dnar::readFaDir('seqs/correct')
corrections$name<-sub('_fixed','',corrections$name)
rownames(corrections)<-corrections$name
if(any(!corrections$name %in% seqs$name))stop('Unmatched correction name')
corrections[seqs$name[seqs$name %in% corrections$name],'rev']<-levenR::leven(seqs[seqs$name %in% corrections$name,'seq'],dnar::revComp(corrections[seqs$name[seqs$name %in% corrections$name],'seq']),oneToOne=TRUE,substring2=TRUE)<levenR::leven(seqs[seqs$name %in% corrections$name,'seq'],corrections[seqs$name[seqs$name %in% corrections$name],'seq'],oneToOne=TRUE,substring2=TRUE)
corrections$seq[corrections$rev]<-dnar::revComp(corrections$seq[corrections$rev])
seqs[seqs$name %in% corrections$name,'seq']<-corrections[seqs$name[seqs$name %in% corrections$name],'seq'] 
sampleRegex<-'^.*(UK|MM)([0-9]+)[.-]([0-9]+).*$'
if(any(!grepl(sampleRegex,seqs$name)))stop('UK/MM not found in all seqs')
seqs$inCohort<-grepl(cohortRegex,seqs$name)
seqs$sample<-sub(sampleRegex,'\\1\\2.\\3',seqs$name)
seqs$pat<-sub('\\..*','',seqs$sample)
seqs$visit<-sub('.*\\.','',seqs$sample)
seqs$mm<-ifelse(grepl('MM',seqs$pat),seqs$pat,c(mmLookup,structure(mmLookup,.Names=sub('EJ','UK',names(mmLookup))))[seqs$pat])
idRegex<-'^.*(UK|MM)[0-9]+[.-][0-9]+[._-]([0-9]*[A-Z][A-Z0-9]+).*$'
seqs$hasId<-grepl(idRegex,seqs$name)
seqs$id[seqs$hasId]<-sub(idRegex,'\\2',seqs$name[seqs$hasId])
seqs$sga<-grepl('\\.SGS|sgs|sga|SGA|env\\.',seqs$name)
#MM33.XX.XX.XX are sga
seqs$sga<-seqs$sga|grepl('MM33\\.[0-9]+\\.[0-9]+\\.[0-9]+',seqs$name)
seqs$pbmc<-grepl('\\.pbmc\\.',seqs$name)
seqs$voa<-grepl('(VOA|voa)[._-]',seqs$name)
inCohort<-seqs[seqs$inCohort,]
inCohort$time<-apply(inCohort[,c('mm','visit')],1,function(xx)compiledMeta[compiledMeta$mm==xx[1]&suppressWarnings(as.numeric(compiledMeta$visit))==as.numeric(xx[2])&!is.na(suppressWarnings(as.numeric(compiledMeta$visit))),'time'])
inCohort$internal[inCohort$hasId]<-dnar::withAs(xx=inCohort[inCohort$hasId,],paste(xx$mm,sprintf('%02d',as.numeric(xx$visit)),sub('^P','',xx$id),sep='.'))
inCohort$isolate[inCohort$hasId]<-sapply(inCohort$internal[inCohort$hasId],function(xx){
  hits<-grepl(xx,dat$id)&!grepl('bulk|BULK',dat$id)
  if(sum(hits)>1){
    warning(sprintf('Ambiguous id match: %s',xx))
    return(NA)
  }
  if(sum(hits)==1)return(dat$id[hits])
  return(NA)
})
if(any(!grepl('VOA',inCohort[inCohort$voa,'isolate'])))stop('VOA not matched with VOA')
#have ambiguities in sequence so excluded for non-limiting
pruned<-c('MM40.07.11C4','MM33.17.1D5') 
inCohort<-inCohort[!inCohort$internal %in% pruned,]
if(any(is.na(inCohort$isolate)&inCohort$hasId&!inCohort$internal %in% pruned))stop('Unmatched ID')
if(any(dat$Sequence!='' & !dat$id %in% inCohort$isolate))stop('Problem finding known sequence in all')
if(any(!apply(cbind(inCohort[,c('voa','sga','pbmc')],!is.na(inCohort$isolate)),1,any)))stop('Unassigned sequence')
if(any(inCohort$pbmc & inCohort$sga))stop('Overlapping pbmc and sga')
if(any(!is.na(inCohort$isolate) & inCohort$sga))stop('Overlapping isolate and sga')
if(any(!is.na(inCohort$isolate) & inCohort$pbmc))stop('Overlapping isolate and pbmc')
inCohort$baseName<-sprintf('%s.%s.%s.%05d',inCohort$mm,ifelse(inCohort$voa|inCohort$pbmc,'PBMC','PLAS'),ifelse(inCohort$voa,'VOA',ifelse(!is.na(inCohort$isolate),'ISO','SGS')),inCohort$time)
inCohort$id<-ave(1:nrow(inCohort),inCohort$baseName,FUN=function(xx)1:length(xx))
inCohort$newName<-sprintf('%s.%03d.01',inCohort$baseName,inCohort$id)
if(any(table(inCohort$newName)>1))stop('Multiple sequences with same name')
write.csv(inCohort[,c('name','isolate','newName','seq')],'out/newNames.csv',row.names=FALSE)
write.fa(inCohort$newName,inCohort$seq,'out/mmLongitudinalSequences.fa')

#align on LANL with HMM
refName<-'B.FR.83.HXB2_LAI_IIIB_BRU_K03455'
aligned<-read.fa('work/alignHmm_2019-04-23/AlignSeq.nt.fasta.gz')
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
