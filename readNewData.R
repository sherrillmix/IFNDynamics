library(dnar)
library(lubridate)

if(!exists('compiledMeta'))source('readMeta.R')


#dat<-read.csv('data/Data Master MG, Sept_2.csv')
#dat<-read.csv('data/MM cohort cata master 11.29.2017.csv')
#allDats<-lapply(c('data/for Scott_2017_12_13.csv','data/For Scott - All bulk isol. alpha and beta.csv'),read.csv,stringsAsFactors=FALSE)
#allDats<-lapply(c('data/for Scott_Data Marster.csv'),read.csv,stringsAsFactors=FALSE)
#allDats<-lapply(c('data/For Scott Jan.2019.csv'),read.csv,stringsAsFactors=FALSE)
#allDats<-lapply(c('data/Data Master Marvin_2019-04-10.csv'),read.csv,stringsAsFactors=FALSE)
allDats<-lapply(c('data/Data Master 20200213.csv'),read.csv,stringsAsFactors=FALSE)
idCol<-'ID.for.statistical.analysis.....Scott.March..2018.'
allDats<-lapply(allDats,function(dat)dat[!is.na(dat[,idCol])&dat[,idCol]!='',])
allCols<-unique(unlist(lapply(allDats,colnames)))
allDats<-lapply(allDats,function(dat){dat[,allCols[!allCols %in% colnames(dat)]]<-NA;dat[,allCols]})
dat<-do.call(rbind,allDats)

#dat<-read.csv('data/for Scott_2017_12_13.csv')
#dat2<-read.csv('data/For Scott - All bulk isol. alpha and beta.csv')

#Standardize BULK naming and catch _BULK with no .
dat$qvoa<-grepl('VOA',dat[,idCol])
dat$id<-sub('[ _.][Bb]ulk-([IVX]+)','.\\1.bulk',sub('VOA[_ ]','',sub(' +[Bb]ulk|_Bulk|_BULK','.bulk',dat[,idCol])))
dat$id<-sprintf('%s%s',dat$id,ifelse(dat$qvoa&!grepl('\\.VOA$',dat$id),'.VOA',''))
dat$id<-sub('^Mm','MM',dat$id)
dat$id<-sub('^(MM[0-9]+)\\.([0-9])\\.','\\1.0\\2.',dat$id)
dat$bulk<-grepl('[Bb]ulk',dat$id)
probs<-grepl('^MM[0-9]+\\.[0-9]+\\.bulk',dat$id)
dat[probs,'id']<-withAs(xx=dat[probs,],ave(sub('\\.bulk','',xx$id),sub('\\.bulk','',xx$id),FUN=function(yy)sprintf('%s.%d.bulk',yy,1:length(yy))))
splits<-strsplit(dat$id,'\\.')
probs<-!(sapply(splits,length)==3 & !dat$qvoa & !dat$bulk) & !(sapply(splits,length)==4 & (dat$qvoa|dat$bulk))
if(any(probs))stop('Problem id found')
tmp<-dat[,c(idCol,'id')]
tmp$diff<-tmp[,idCol]!=tmp$id
write.csv(tmp,'newIds.csv',row.names=FALSE)
dat$sample<-sapply(strsplit(dat$id,'\\.'),function(xx)paste(xx[1:2],collapse='.'))
dat$visit<-sapply(strsplit(dat$sample,'\\.'),'[',2)
dat$virusId<-sapply(strsplit(dat$id,'\\.'),'[',3)
dat$pat<-sub('\\.[^.]+$','',dat$sample)
dat$time<-as.numeric(compiledMeta[dat$sample,'time'])
dat$timeBeforeArt<-compiledMeta[dat$sample,'daysBeforeArt']
dat$timeBefore350<-compiledMeta[dat$sample,'daysBefore350']
dat$vl<-compiledMeta[dat$sample,'vl']
dat$CD4<-compiledMeta[dat$sample,'cd4']
#dat<-dat[!is.na(dat$ic50)|!is.na(dat$vres)|!is.na(dat$beta)|!is.na(dat$betaVres),]
if(any(is.na(dat$time)))stop('Missing time')
#if(any(is.na(dat$vl)))stop('Missing vl')
#if(any(is.na(dat$CD4)))stop('Missing CD4')
dat$time2<-dat$time^2
dat$logTime<-log(dat$time)
dat$logTime2<-log(dat$time)^2
dat$logTime3<-log(dat$time)^3
dat$logTime4<-log(dat$time)^4
dat$logVl<-log(dat$vl)
dat$rt<-as.numeric(sub('N/A','',dat$RT.activity..ηg.μl.))

dat$fillVl<-dat$vl
dat$fillVl[is.na(dat$vl)]<-sapply(which(is.na(dat$vl)),function(xx){
  thisTime<-dat[xx,'time']
  thisPat<-dat[xx,'pat']
  thisDat<-compiledMeta[compiledMeta$mm==thisPat&!is.na(compiledMeta$vl),c('time','vl')]
  #print(thisDat)
  #print(thisTime)
  if(thisTime>max(thisDat$time))return(thisDat$vl[thisDat$time==max(thisDat$time)])
  else return(approx(thisDat$time,thisDat$vl,thisTime)$y)
})


dat$fillCD4<-dat$CD4
dat$fillCD4[is.na(dat$CD4)]<-sapply(which(is.na(dat$CD4)),function(xx){
  thisTime<-dat[xx,'time']
  thisPat<-dat[xx,'pat']
  thisDat<-compiledMeta[compiledMeta$mm==thisPat&!is.na(compiledMeta$cd4),c('time','cd4')]
  #if(thisTime>max(thisDat$time))return(thisDat$vl[thisDat$time==max(thisDat$time)])
  approx(thisDat$time,thisDat$cd4,thisTime)$y
})


cols<-rainbow.lab(length(unique(dat$pat)))
names(cols)<-unique(dat$pat)

ifna2_ic50<-colnames(dat)[grep('IFNa2.*IC50',colnames(dat))]
ifna2_vres<-colnames(dat)[grep('IFNa2.*Vres',colnames(dat))]
ifnb_ic50<-colnames(dat)[grep('IFNb.*IC50',colnames(dat))]
ifnb_vres<-colnames(dat)[grep('IFNb.*Vres',colnames(dat))]
dat$ic50<-apply(dat[,ifna2_ic50],1,function(xx)exp(mean(log(xx),na.rm=TRUE)))
dat$vres<-apply(dat[,ifna2_vres],1,mean,na.rm=TRUE)
dat$beta<-apply(dat[,ifnb_ic50],1,function(xx)exp(mean(log(xx),na.rm=TRUE)))
#dat$IFNbeta..Pooled.Donor.cells.IC50..pg.ml.
dat$betaVres<-apply(dat[,ifnb_vres],1,mean,na.rm=TRUE)
if(any(dat$betaVres==0)){
  warning('Beta Vres equal to zero. Setting to lowest value/10.')
  dat[dat$betaVres==0&!is.na(dat$betaVres),'betaVres']<-min(dat[dat$betaVres!=0&!is.na(dat$betaVres),'betaVres'])/10
}
dat$replication<-dat$Replicative.capacity.Pooled.Donor.cells.p24.d7..from.June.2017.repeat.
ifnVars<-c('Interferon alpha 2 IC50'='ic50','Interferon beta IC50'='beta','Interferon alpha 2 Vres'='vres','Interferon beta Vres'='betaVres','Replication capacity'='replication')

#dat$time<-sapply(dat$sample,function(xx)comboMeta[paste(compiledMeta$mm,compiledMeta$visit,sep='.')==xx,'time'])
#if(any(is.na(dat$time)))stop('Missing time metadata')

rownames(dat)<-dat$id

weau<-dat[dat$pat=='WEAU',]
#dat<-dat[dat$pat!='WEAU',]

infect<-read.csv('out/ius.csv',stringsAsFactors=FALSE)
infect$id[infect$id=='VOA_MM23.18.1A1']<-'MM23.18.1A1.VOA'
infect$id[infect$id=='MM33.13 Bulk-I']<-'MM33.13.1.bulk'
infect$id[infect$id=='MM33.14_Bulk-I']<-'MM33.14.1A1.bulk'
infect$id[infect$id=='MM33.17_Bulk-I']<-'MM33.17.1A1.bulk'
infect$id<-sub(' bulk','.bulk',infect$id)
if(any(!infect$id %in% dat$id))stop('Problem associating infectivity with isolates')
dat$infectivityMedia<-dat$infectivityDextran<-NA
dat[infect$id,'infectivityMedia']<-infect$media/dat[infect$id,'rt']/1000
dat[infect$id,'infectivityDextran']<-infect$dextran/dat[infect$id,'rt']/1000

replicative<-read.csv('out/weau3_ic50.csv',row.names=1)
replicative<-replicative[grep('Alpha|alpha',rownames(replicative)),]
replicative$id<-sub('(UK[0-9]+|EJ[0-9]+|WEAU)\\.([0-9])[._-]','\\1.0\\2.',sub(' .*$','',rownames(replicative)))
for(ii in names(mmLookup))replicative$id<-sub(sprintf('%s|%s',ii,sub('EJ','UK',ii)),mmLookup[ii],replicative$id)
replicative$id<-sub('[._-]P','.',replicative$id)
replicative$id<-sprintf('%s%s',replicative$id,ifelse(grepl('WEAU|MM55|MM62|MM15',replicative$id),'.bulk',''))
inDat<-replicative$id %in% dat$id
dat[replicative$id[inDat],'replication']<-replicative$replication[inDat]

ice<-read.csv('out/iceHalf.csv',row.names=1)
ice$id<-sub('VOA_(.*)','\\1.VOA',sub(' bulk','.bulk',rownames(ice)))
ice$id[ice$id=='MM33.13 Bulk-I']<-'MM33.13.1.bulk'
inDat<-ice$id %in% dat$id
dat[ice$id[inDat],'iceHalf']<-ice$half[inDat]


#withAs(xx=dat[dat$time>35*7,],plot(ave(xx$vl,xx$pat,FUN=function(xx)(xx-min(xx,na.rm=TRUE))/max(xx-min(xx,na.rm=TRUE),na.rm=TRUE)),xx$ic50,bg=patCols[xx$pat],log='y',pch=21,cex=2))

# regenerate old plot
# show old plot + new data
## new data plots with fits 
## viral load plot (split by 35 weeks)
## cd4 plot (split by 35 weeks)
# split out bulks "_Bulk"
## plot Replicative over time
# plot everything for beta
# plot vl and cd4 vs beta
# plot residuals vs vl and cd4
# plot raw data

# bayesian model 

nameLookup<-read.csv('out/newNames.csv')
nameLookup<-nameLookup[!is.na(nameLookup$isolate),]
rownames(nameLookup)<-nameLookup$isolate
dat$finalId<-sub('\\.01$','',nameLookup[dat$id,'newName'])
dat$count<-ave(sapply(strsplit(dat$finalId,'\\.'),'[',5),paste(dat$pat,dat$time),FUN=function(xx){
  if(any(is.na(xx)))xx[is.na(xx)]<-sprintf('%03d',max(c(0,as.numeric(xx)),na.rm=TRUE)+1:sum(is.na(xx)))
  xx
})
dat$finalId[is.na(dat$finalId)]<-withAs(xx=dat[is.na(dat$finalId),],sprintf('%s.%s.%s.%05d.%s',xx$pat,ifelse(xx$qvoa,'PBMC','PLAS'),ifelse(xx$qvoa,'VOA','ISO'),xx$time,xx$count))
if(any(table(dat$finalId)>1))stop('Duplicate name created')
if(any(!sub('\\.01$','',nameLookup$newName) %in% dat$finalId))stop('Missing isolate data')
if(any(!tapply(dat$count,paste(dat$pat,dat$time),function(xx)all(sort(as.numeric(xx))==1:length(xx)))))stop('Missing ID')



message('Time after infection')
print(mean(withAs(dat=dat[!dat$qvoa,],tapply(dat$time,dat$pat,FUN=max))/7))
zz<-withAs(dat=dat[!dat$qvoa,],tapply(dat$ic50,list(dat$pat,dat$time),mean,na.rm=TRUE))
message('Fold range of ic50')
print(mean(apply(zz,1,function(xx)max(xx,na.rm=TRUE)/min(xx,na.rm=TRUE))))
if(any(as.numeric(colnames(zz))!=sort(as.numeric(colnames(zz)))))stop('need sorted cols')
message('Fold increase from nadir increase of ic50')
print(mean(apply(zz,1,function(xx)tail(xx[!is.na(xx)],1)/min(xx,na.rm=TRUE))))
print(mean(apply(zz[!rownames(zz) %in% c('MM55','MM62','MM15','WEAU'),],1,function(xx)tail(xx[!is.na(xx)],1)/min(xx,na.rm=TRUE))))
message('Weeks after nadir')
print(mean(sapply(rownames(zz),function(xx){tmp<-zz[xx,!is.na(zz[xx,])];diff(as.numeric(names(tmp))[c(which.min(tmp),length(tmp))])}))/7)
message('Median CD4 at last time point')
print(median(withAs(xx=dat[!dat$qvoa&!is.na(dat$CD4),][dat$time[!dat$qvoa&!is.na(dat$CD4)]==ave(dat$time[!dat$qvoa&!is.na(dat$CD4)],dat$pat[!dat$qvoa&!is.na(dat$CD4)],FUN=max),c('pat','CD4')],tapply(xx$CD4,xx$pat,unique))))

dat$meanIc50<-ave(dat$ic50,paste(dat$pat,dat$time),FUN=function(xx)mean(xx,na.rm=TRUE))
dat$meanBeta<-ave(dat$beta,paste(dat$pat,dat$time),FUN=function(xx)mean(xx,na.rm=TRUE))
dat$isNadir<-dat$meanIc50==ave(dat$meanIc50,dat$pat,FUN=function(xx)min(xx,na.rm=TRUE))
dat$isBetaNadir<-dat$meanBeta==ave(dat$meanBeta,dat$pat,FUN=function(xx)min(xx,na.rm=TRUE))
dat$isFirst<-dat$time==ave(dat$time,dat$pat,FUN=function(xx)min(xx,na.rm=TRUE))
#dat$isSix<-dat$time==ave(dat$time,dat$pat,FUN=function(xx){diff<-abs(180-xx);out<-xx[which.min(diff)][1];if(diff[which.min(diff)][1]>30)return(FALSE);out})
dat$isSix<-abs(dat$time-180)<30
dat$isYear<-abs(dat$time-365)<60
dat$isLast<-dat$time==ave(dat$time*ifelse(dat$qvoa,0,1),dat$pat,FUN=max)
write.csv(dat[dat$isFirst|dat$isNadir|dat$qvoa|dat$isSix|dat$isLast|dat$isBetaNadir,c('pat','time','ic50','isFirst','isNadir','isBetaNadir','isSix','isLast','qvoa','beta','replication')],'out/firstNadir.csv')
write.csv(dat,'out/allLongitudinal.csv')

timeMeans<-lapply(unique(dat$pat),function(xx)withAs(zz=dat[dat$pat==xx,],tapply(zz$beta,zz$time,function(yy)exp(mean(log(yy),na.rm=TRUE)))))
cbind(unique(dat$pat),sapply(timeMeans,which.min),names(sapply(timeMeans,which.min)))

dat$patTime<-paste(dat$pat,dat$time)
bulkTimes<-unique(dat[grepl('bulk',dat$id)&(!is.na(dat$ic50)|!is.na(dat$beta)),'patTime'])
bulkTimes<-bulkTimes[sapply(bulkTimes,function(xx)any(!dat$bulk&dat$patTime==xx))]
outside<-sapply(bulkTimes,function(xx,var='beta'){
  thisDat<-dat[dat$patTime==xx,]
  r<-range(thisDat[!thisDat$bulk,var],na.rm=TRUE)
  s<-sd(log(thisDat[!thisDat$bulk,var]),na.rm=TRUE)
  m<-mean(log(thisDat[!thisDat$bulk,var]),na.rm=TRUE)
  thisDat[thisDat$bulk,var]<r[1]|thisDat[thisDat$bulk,var]>r[2]
  #table(c(TRUE,FALSE,thisDat[thisDat$bulk,var]<r[1]|thisDat[thisDat$bulk,var]>r[2]))-1
  #max(abs(thisDat[thisDat$bulk,var]-m)/s)
  #abs(log(thisDat[thisDat$bulk,var])-m)/s
  #pt(abs(log(thisDat[thisDat$bulk,var])-m)/s,sum(!thisDat$bulk)-1)
  #t.test(log(thisDat[!thisDat$bulk,var]),log(thisDat[thisDat$bulk,var]),var.equal=TRUE)$p.value
  2^abs(log2(thisDat[thisDat$bulk,var]/exp(m)))
})

summary(lm(I(log(ic50))~patTime+bulk,dat[dat$patTime %in% bulkTimes,]))
summary(lm(I(log(beta))~patTime+bulk,dat[dat$patTime %in% bulkTimes,]))
zz<-withAs(xx=dat[dat$patTime %in% bulkTimes,],tapply(log(xx$beta),list(xx$patTime,xx$bulk),mean,na.rm=TRUE))
cor.test(zz[,1],zz[,2])

