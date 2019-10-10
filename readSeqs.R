if(!exists('mmLookup'))source('readNewData.R')
bulk<-c('MM14','MM15','MM55','MM62','WEAU')

hxb2<-dnar::read.fa('work/hxb2.fa')$seq
cohortRegex<-paste(c('WEAU',sub('EJ','UK',ejLookup),ejLookup,mmLookup),collapse='|')
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
patRegex<-'WEAU|EJ[0-9]+|UK[0-9]+|MM[0-9]+'
sampleRegex<-sprintf('^.*(%s)[.-]([0-9]+).*$',patRegex)
seqs$noAss<-sub('ass?emble.*','',seqs$name)
if(any(!grepl(sampleRegex,seqs$name)))stop('UK/MM not found in all seqs')
seqs$inCohort<-grepl(cohortRegex,seqs$noAss)
seqs$sample<-sub(sampleRegex,'\\1.\\2',seqs$noAss)
seqs$pat<-sub('\\..*','',seqs$sample)
seqs$visit<-sub('.*\\.','',seqs$sample)
seqs$mm<-ifelse(grepl('MM|WEAU',seqs$pat),seqs$pat,c(mmLookup,structure(mmLookup,.Names=sub('EJ','UK',names(mmLookup))))[seqs$pat])
idRegex<-sprintf('^.*(%s)[.-][0-9]+[._-]([0-9]*[A-Z][A-Z0-9]+).*$',patRegex)
seqs$hasId<-grepl(idRegex,seqs$noAss)
seqs$id[seqs$hasId]<-sub(idRegex,'\\2',seqs$noAss[seqs$hasId])
seqs$sga<-grepl('\\.SGS|sgs|sga|SGA|env\\.',seqs$noAss)
#MM33.XX.XX.XX are sga
seqs$sga<-seqs$sga|grepl('MM33\\.[0-9]+\\.[0-9]+\\.[0-9]+',seqs$noAss)
seqs$pbmc<-grepl('\\.pbmc\\.',seqs$noAss)
seqs$voa<-grepl('(VOA|voa)[._-]',seqs$noAss)

inCohort<-seqs[seqs$inCohort,]
inCohort$time<-apply(inCohort[,c('mm','visit')],1,function(xx)compiledMeta[compiledMeta$mm==xx[1]&suppressWarnings(as.numeric(compiledMeta$visit))==as.numeric(xx[2])&!is.na(suppressWarnings(as.numeric(compiledMeta$visit))),'time'])
inCohort$internal[inCohort$hasId]<-dnar::withAs(xx=inCohort[inCohort$hasId,],paste(xx$mm,sprintf('%02d',as.numeric(xx$visit)),sub('^P','',xx$id),sep='.'))
inCohort$isolate[inCohort$hasId]<-sapply(inCohort$internal[inCohort$hasId],function(xx){
  hits<-grepl(xx,dat$id)&(!grepl('bulk|BULK',dat$id)|dat$pat %in% bulk)
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
if(any(table(inCohort$internal)>1))stop('Duplicate isolate')
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


timeConvert<-read.csv('seqs/new/London cohort date decoder.170509.csv')
timeConvert<-timeConvert[timeConvert$Time.point!=''&timeConvert$days.post.infection!='',]
timeConvert$ID[grep('EJ',timeConvert$ID)]<-''
timeConvert$ID<-dnar::fillDown(sub('\\*','',timeConvert$ID))
timeConvert$visit<-sub('(MM|EJ|EJ |)[0-9]+\\.','',timeConvert$Time.point)
rownames(timeConvert)<-paste(timeConvert$ID,timeConvert$days.post.infection)

sga<-dnar::read.fa('seqs/new/SGA and PBMC SGA for London Cohort.fasta')
sga$clean<-sub('B\\.US\\.199.\\.WEAUd?','WEAU.',sga$name)
sga$pat<-sub('\\..*','',sga$clean)
sga$pbmc<-grepl('pbmc',sga$name)
bakSga<-sga[sga$pat=='MM33'|sga$pat=='MM23',]
sga<-sga[sga$pbmc|(sga$pat!='MM33'&sga$pat!='MM23'),]
#MM23.18.pbmc.72: weird pol sequence in middle of env
sgaExclude<-c('MM23.18.pbmc.72')
sga<-sga[!sga$name %in% sgaExclude,]
mm23mm33<-rbind(dnar::read.fa('seqs/new/MM33 SGA ORIGINAL.fasta'),dnar::read.fa('seqs/new/MM23 SGA original sequences.fasta'))
mm23mm33$clean<-sub('UK61\\.','MM23.',mm23mm33$name)
mm23mm33$pat<-sub('\\..*','',mm23mm33$clean)
mm23mm33$pbmc<-FALSE
sga<-rbind(sga,mm23mm33)
sga$timeString<-sapply(strsplit(sga$clean,'[_.]'),'[[',2)
sga[sga$timeString=='12MW'&sga$pat=='MM39','timeString']<-'13'
isD<-grep('^d[0-9]+$',sga$timeString)
if(any(!paste(sga[isD,'pat'],sub('^d0*','',sga[isD,'timeString'])) %in% rownames(timeConvert)))stop('Problem matching dxxx and conversion')
sga[isD,'timeString']<-timeConvert[paste(sga[isD,'pat'],sub('^d0*','',sga[isD,'timeString'])),'visit']
#symptoms started 20 days after infection
sga[sga$pat=='WEAU','time']<-as.numeric(sga[sga$pat=='WEAU','timeString'])-20
sga[sga$pat!='WEAU','time']<-apply(sga[sga$pat!='WEAU',c('pat','timeString')],1,function(xx){
  as.numeric(compiledMeta[compiledMeta$mm==xx[1]&sub('^0','',compiledMeta$visit)==sub('^0','',xx[2])&!is.na(compiledMeta$visit),'time'])
})
sga$pbmc<-grepl('pbmc',sga$name)
sga$mm<-sga$pat
sga$voa<-FALSE
sga$isolate<-NA
sga$visit<-sga$timeString
 
combine<-inCohort[!is.na(inCohort$isolate),c('name','seq','mm','time','pbmc','voa','isolate','visit')]
combine<-rbind(combine,sga[,colnames(combine)])

combine$desc<-sprintf('%s.%s',ifelse(combine$voa|combine$pbmc,'PBMC','PLAS'),ifelse(combine$voa,'VOA',ifelse(!is.na(combine$isolate),'ISO','SGS')))
combine$baseName<-sprintf('%s.%s.%05d',combine$mm,combine$desc,combine$time)
combine$id<-ave(1:nrow(combine),combine$baseName,FUN=function(xx)1:length(xx))
combine$newName<-sprintf('%s.%03d.01',combine$baseName,combine$id)
if(any(table(combine$newName)>1))stop('Multiple sequences with same name')
table(combine$mm,combine$desc)

write.csv(table(sprintf('%s.%05d (%s)',combine$mm,combine$time,combine$visit),combine$desc),'out/countByTime.csv')

write.csv(combine[,c('name','isolate','newName','seq')],'out/newNames.csv',row.names=FALSE)
dnar::withAs(combine=combine[!combine$mm %in% bulk,],write.fa(combine$newName,combine$seq,'out/mmLongitudinalSequences_noBulk.fa'))
dnar::withAs(combine=combine[combine$mm %in% bulk,],write.fa(combine$newName,combine$seq,'out/mmLongitudinalSequences_bulk.fa'))
dnar::withAs(combine=combine,write.fa(combine$newName,combine$seq,'out/mmLongitudinalSequences.fa.gz'))

for(ii in unique(combine$mm)){
  dnar::withAs(combine=combine[combine$mm==ii,],write.fa(combine$newName,combine$seq,sprintf('out/seqSplit/%s.fa.gz',ii)))
}
#align<-read.fa('')
#pdf('out/allSeqs.pdf')
#dev.off()

xx<-read.fa('out/mmLongitudinalSequences_20190702.fa.gz')
xx$mm<-sapply(strsplit(xx$name,'\\.'),'[',1)
xx$source<-sapply(strsplit(xx$name,'\\.'),'[',2)
xx$type<-sapply(strsplit(xx$name,'\\.'),'[',3)
xx$time<-(sapply(strsplit(xx$name,'\\.'),'[',4))
write.csv('out/mmDayCount.csv',table(paste(xx$mm,xx$time),paste(xx$source,xx$type)))
