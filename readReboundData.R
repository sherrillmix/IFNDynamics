acuteCut<-30

dat<-read.csv('out/allLongitudinal.csv',stringsAsFactors=FALSE)
acuteChronicMM<-dat[!dat$qvoa&(dat$time<30|dat$time>=min(dat$time[dat$isYear])),]
acuteChronicMM$study<-'MM'
acuteChronicMM$virus<-acuteChronicMM$id
acuteChronicMM$ic50_IFNa2<-acuteChronicMM$ic50
acuteChronicMM$ic50_IFNb<-acuteChronicMM$beta
acuteChronicMM$class<-ifelse(acuteChronicMM$time<=acuteCut,'Acute','Chronic')
acuteChronicMM$displayClass<-ifelse(acuteChronicMM$isNadir,'Nadir',ifelse(acuteChronicMM$isYear,'1 Year',ifelse(acuteChronicMM$isLast,'Last',ifelse(acuteChronicMM$time<=acuteCut,'Acute',NA))))
acuteChronicMM$virus<-acuteChronicMM$id
#distinct "patient" for each time point (doesn't work because acute)
#acuteChronicMM$pat<-acuteChronicMM$sample

#voaMM<-dat[dat$qvoa,]
#voaMM$class<-'QVOA'
#voa$study<-'MM'
#voa$virus<-voa$id

rebound<-read.csv('data/Data Master 2020_ReboundandQVOA_20200504.csv',stringsAsFactors=FALSE)
rebound<-rebound[!grepl('_BE$',rebound$ID),]
table(rebound[!is.na(rebound$ic50_IFNa2),'Study'])
rebound$study<-sub(' / ','/',rebound$Study)
rebound$displayClass<-rebound$class<-ifelse(rebound$Type=='Rebound','Rebound',ifelse(grepl('Post',rebound$Type),'Post-ATI',ifelse(grepl('Pre',rebound$Type),'Pre-ATI','Outgrowth')))
rebound$pat<-sub(' \\(R-?[0-9]+\\)','',rebound$Patient)
rebound<-rebound[!is.na(rebound$Type),]
rebound[rebound$study=='OUTGROWTH','study']<-'MM'
rebound$virus<-rebound$ID
#table(comboA$study)

pair<-read.csv('rebound/donorRecipient.csv',stringsAsFactors=FALSE)
pair$class<-ifelse(pair$donor,'Donor','Recipient')
pair$type<-'CHAVI cohort'
#pair$label<-sprintf('CHAVI %s',pair$class)
pair$displayClass<-sprintf('%s %s',ifelse(pair$class=='Donor','Chronic','Acute'),pair$class)
pair$ic50_IFNa2<-pair$IFNa2.Pooled.Donor.cells.IC50..pg..ml
pair$ic50_IFNb<-pair$IFNbeta.Pooled.Donor.cells.IC50..pg.ml
#pair[,'class']<-c('Recipient'='Acute Recipient','Donor'='Chronic Donor')[pair$class]
pair$pat<-pair$sample
pair$repCap<-pair$Replicative.capacity.Pooled.Donor.cells.p24.d7
pair[,colnames(dat)[!colnames(dat) %in% colnames(pair)]]<-NA
pair$virus<-pair$Renamed
pair$source<-'shilpa'
pair$study<-'Transmission'

#pat,simpleClass,study,ic50,class
targetCols<-c('pat','class','study','ic50_IFNa2','ic50_IFNb','displayClass','virus')
combined<-rbind(acuteChronicMM[,targetCols],rebound[,targetCols],pair[,targetCols])
combined$simpleClass<-ifelse(combined$class %in% c('Chronic Donor','Donor','1 Year','Nadir','Last'),'Chronic',combined$class)
combined$simpleClass[combined$simpleClass=='Recipient']<-'Acute'
combined$simpleClass[combined$simpleClass=='Pre-ATI']<-'Outgrowth' #may want own class
#combined$simpleClass[combined$simpleClass=='Pre-ATI']<-'QVOA' #may want own class
combined$simplePat<-ifelse(combined$simpleClass %in% c('Acute','Chronic','Donor','Recipient'),'',combined$pat)
combined$label<-ifelse(combined$simplePat=='',combined$displayClass,sprintf('%s%s%s',combined$displayClass,ifelse(combined$simplePat!='',' ',''),combined$simplePat))
combined<-combined[!is.na(combined$ic50_IFNb)|!is.na(combined$ic50_IFNa2),]

fastRegex<-'MM15|WEAU'
slowRegex<-'MM55|MM62'
standardRegex<-'MM14|MM23|MM33|MM34|MM39|MM40'
combined$speed<-ifelse(grepl(fastRegex,combined$virus),'Fast',ifelse(grepl(slowRegex,combined$virus),'Slow',ifelse(grepl(standardRegex,combined$virus),'Standard','Other')))


if(FALSE){
zz<-table(comboA$pat)
zz2<-table(combined$pat[!is.na(combined$ic50_IFNa2)])
zz<-table(comboB$pat)
zz2<-table(combined$pat[!is.na(combined$ic50_IFNb)])
all(names(zz2) == names(zz))
cbind(zz,zz2)
all(zz2>=zz)
}
