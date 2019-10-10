library(dnar)
library(lubridate)

if(!exists('compiledMeta'))source('readMeta.R')


patCols<-c('MM23'='#e41a1c','MM33'='#4daf4a','MM34'='#984ea3','MM39'='#377eb8','MM40'='#FF7f00','MM14'='#FFD700','MM15'='#f781bf','MM55'='#a65628','MM62'='#00CED1')
patCols<-c(patCols,'WEAU'='#708090')
patCols2<-sprintf('%s33',patCols)
patCols3<-sprintf('%s11',patCols)
names(patCols2)<-names(patCols3)<-names(patCols)
lay<-matrix(0,nrow=7,ncol=4)
lay[2:6,2:3]<-matrix(1:10,nrow=5,byrow=TRUE)

lay2<-matrix(0,nrow=7+2,ncol=4)
lay2[c(2,3,4,6,8),2:3]<-matrix(1:10,nrow=5,byrow=TRUE)
lowerP24Limit<-60
patOrder<-c("MM14","MM23","MM33","MM34","MM39","MM40","MM55","MM62","MM15","WEAU")


#dat<-read.csv('data/Data Master MG, Sept_2.csv')
#dat<-read.csv('data/MM cohort cata master 11.29.2017.csv')
#allDats<-lapply(c('data/for Scott_2017_12_13.csv','data/For Scott - All bulk isol. alpha and beta.csv'),read.csv,stringsAsFactors=FALSE)
#allDats<-lapply(c('data/for Scott_Data Marster.csv'),read.csv,stringsAsFactors=FALSE)
#allDats<-lapply(c('data/For Scott Jan.2019.csv'),read.csv,stringsAsFactors=FALSE)
allDats<-lapply(c('data/Data Master Marvin_2019-04-10.csv'),read.csv,stringsAsFactors=FALSE)
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
dat$time<-as.numeric(compiledMeta[dat$sample,'DFOSx'])
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

cols<-rainbow.lab(length(unique(dat$pat)))
names(cols)<-unique(dat$pat)

ifna2_ic50<-colnames(dat)[grep('IFNa2.*IC50',colnames(dat))]
ifna2_vres<-colnames(dat)[grep('IFNa2.*Vres',colnames(dat))]
ifnb_ic50<-colnames(dat)[grep('IFNb.*IC50',colnames(dat))]
ifnb_vres<-colnames(dat)[grep('IFNb.*Vres',colnames(dat))]
dat$ic50<-apply(dat[,ifna2_ic50],1,mean,na.rm=TRUE)
dat$vres<-apply(dat[,ifna2_vres],1,mean,na.rm=TRUE)
dat$beta<-apply(dat[,ifnb_ic50],1,mean,na.rm=TRUE)
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
dat[infect$id,'infectivityMedia']<-infect$media/dat[infect$id,'rt']
dat[infect$id,'infectivityDextran']<-infect$dextran/dat[infect$id,'rt']

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


message('Time after infection')
print(mean(withAs(dat=dat[!dat$qvoa,],tapply(dat$time,dat$pat,FUN=max))/7))
zz<-withAs(dat=dat[!dat$qvoa,],tapply(dat$beta,list(dat$pat,dat$time),mean,na.rm=TRUE))
message('Fold range of ic50')
print(mean(apply(zz,1,function(xx)max(xx,na.rm=TRUE)/min(xx,na.rm=TRUE))))
if(any(as.numeric(colnames(zz))!=sort(as.numeric(colnames(zz)))))stop('need sorted cols')
message('Fold increase from nadir increase of ic50')
print(mean(apply(zz,1,function(xx)tail(xx[!is.na(xx)],1)/min(xx,na.rm=TRUE))))
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
dat$isLast<-dat$time==ave(dat$time*ifelse(dat$qvoa,0,1),dat$pat,FUN=max)
write.csv(dat[dat$isFirst|dat$isNadir|dat$qvoa|dat$isSix|dat$isLast|dat$isBetaNadir,c('pat','time','ic50','isFirst','isNadir','isBetaNadir','isSix','isLast','qvoa','beta','replication')],'out/firstNadir.csv')
write.csv(dat,'out/allLongitudinal.csv')

