if(!exists('mmLookup'))source('readNewData.R')

hxb2<-dnar::read.fa('work/hxb2.fa')$seq
cohortRegex<-paste(c(sub('EJ','UK',ejLookup),ejLookup,mmLookup),collapse='|')
seqs<-dnar::readFaDir('seqs')
corrections<-dnar::readFaDir('seqs/correct')
corrections$name<-sub('_fixed','',corrections$name)
rownames(corrections)<-corrections$name
if(any(!corrections$name %in% seqs$name))stop('Unmatched correction name')
corrections[seqs$name[seqs$name %in% corrections$name],'rev']<-levenR::leven(seqs[seqs$name %in% corrections$name,'seq'],dnar::revComp(corrections[seqs$name[seqs$name %in% corrections$name],'seq']),oneToOne=TRUE,substring2=TRUE)<levenR::leven(seqs[seqs$name %in% corrections$name,'seq'],corrections[seqs$name[seqs$name %in% corrections$name],'seq'],oneToOne=TRUE,substring2=TRUE)
corrections$seq[corrections$rev]<-dnar::revComp(corrections$seq[corrections$rev])
seqs[seqs$name %in% corrections$name,'seq']<-corrections[seqs$name[seqs$name %in% corrections$name],'seq'] 
deletions<-read.csv('seqs/correct/delete.csv',stringsAsFactors=FALSE)$delete
if(any(!deletions %in% seqs$name))stop('Deletion not found')
seqs<-seqs[!seqs$name %in% deletions,]
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
#if(any(dat$Sequence!='' & !dat$id %in% inCohort$isolate))stop('Problem finding known sequence in all')
if(any(!apply(cbind(inCohort[,c('voa','sga','pbmc')],!is.na(inCohort$isolate)),1,any)))stop('Unassigned sequence')
if(any(inCohort$pbmc & inCohort$sga))stop('Overlapping pbmc and sga')
if(any(!is.na(inCohort$isolate) & inCohort$sga))stop('Overlapping isolate and sga')
if(any(!is.na(inCohort$isolate) & inCohort$pbmc))stop('Overlapping isolate and pbmc')
inCohort$baseName<-sprintf('%s.%s.%s.%05d',inCohort$mm,ifelse(inCohort$voa|inCohort$pbmc,'PBMC','PLAS'),ifelse(inCohort$voa,'VOA',ifelse(!is.na(inCohort$isolate),'ISO','SGS')),inCohort$time)
inCohort$id<-ave(1:nrow(inCohort),inCohort$baseName,FUN=function(xx)1:length(xx))
inCohort$newName<-sprintf('%s.%03d.01',inCohort$baseName,inCohort$id)
if(any(table(inCohort$newName)>1))stop('Multiple sequences with same name')

sga<-dnar::read.fa('seqs/new/SGA and PBMC SGA for London Cohort.fasta')
sga$clean<-sub('B\\.US\\.199.\\.WEAUd?','WEAU.',sga$name)
sga$pat<-sub('\\..*','',sga$clean)
sga$timeString<-sapply(strsplit(sga$clean,'[_.]'),'[[',2)
timeConvert<-read.csv('seqs/new/London cohort date decoder.170509.csv')
timeConvert<-timeConvert[timeConvert$Time.point!=''&timeConvert$days.post.infection!='',]
timeConvert$ID[grep('EJ',timeConvert$ID)]<-''
timeConvert$ID<-dnar::fillDown(sub('\\*','',timeConvert$ID))
timeConvert$visit<-sub('(MM|EJ|EJ |)[0-9]+\\.','',timeConvert$Time.point)
rownames(timeConvert)<-paste(timeConvert$ID,timeConvert$days.post.infection)
sga[sga$timeString=='12MW'&sga$pat=='MM39','timeString']<-'13'
isD<-grep('^d[0-9]+$',sga$timeString)
if(any(!paste(sga[isD,'pat'],sub('^d0*','',sga[isD,'timeString'])) %in% rownames(timeConvert)))stop('Problem matching dxxx and consersion')
sga[isD,'timeString']<-timeConvert[paste(sga[isD,'pat'],sub('^d0*','',sga[isD,'timeString'])),'visit']
#symptoms started 20 days after infection
sga[sga$pat=='WEAU','time']<-as.numeric(sga[sga$pat=='WEAU','timeString'])-20
sga[sga$pat!='WEAU','time']<-apply(sga[sga$pat!='WEAU',c('pat','timeString')],1,function(xx){
  as.numeric(compiledMeta[compiledMeta$mm==xx[1]&sub('^0','',compiledMeta$visit)==sub('^0','',xx[2])&!is.na(compiledMeta$visit),'time'])
})
sga$pbmc<-grepl('pbmc',sga$name)
sga$mm<-sga$pat
sga$voa<-FALSE
sga$isolate<-FALSE
 
combine<-inCohort[!inCohort$sga,c('name','seq','mm','time','pbmc','voa','isolate')]
combine<-rbind(combine,sga[,colnames(combine)])
colnames(combine)

combine$baseName<-sprintf('%s.%s.%s.%05d',combine$mm,ifelse(combine$voa|combine$pbmc,'PBMC','PLAS'),ifelse(combine$voa,'VOA',ifelse(!is.na(combine$isolate),'ISO','SGS')),combine$time)
combine$id<-ave(1:nrow(combine),combine$baseName,FUN=function(xx)1:length(xx))
combine$newName<-sprintf('%s.%03d.01',combine$baseName,combine$id)
if(any(table(combine$newName)>1))stop('Multiple sequences with same name')

write.csv(combine[,c('name','isolate','newName','seq')],'out/newNames.csv',row.names=FALSE)
write.fa(inCohort$newName,inCohort$seq,'out/mmLongitudinalSequences.fa')
