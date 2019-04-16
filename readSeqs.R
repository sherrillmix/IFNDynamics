if(!exists('mmLookup'))source('readNewData.R')

cohortRegex<-paste(c(sub('EJ','UK',ejLookup),ejLookup,mmLookup),collapse='|')
seqs<-dnar::readFaDir('seqs')
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
write.fa(inCohort$newName,inCohort$seq,'out/newNames.fa')
#align on LANL with HMM
refName<-'B.FR.83.HXB2_LAI_IIIB_BRU_K03455'
aligned<-read.fa('work/alignHmm_2019-04-16/AlignSeq.nt.fasta.gz')
if(any(!inCohort$newName %in% aligned$name))stop('Missing name in alignment')
if(any(!aligned$name %in% c(refName,inCohort$newName)))stop('New name found in alignment')
pdf('out/aligned_seqs.pdf',width=20)
dnaplotr::plotDNA(aligned$seq)
dev.off()
allProts<-readFaDir('work/alignHmm_2019-04-16','fasta.gz')
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



