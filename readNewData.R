library(dnar)
library(lubridate)

if(!exists('meta'))source('readMeta.R')


patCols<-c('MM23'='#e41a1c','MM33'='#4daf4a','MM34'='#984ea3','MM39'='#377eb8','MM40'='#FF7f00','MM14'='#FFD700','MM15'='#f781bf','MM55'='#a65628','MM62'='#00CED1')
patCols<-c(patCols,'WEAU'='#708090')
patCols2<-sprintf('%s33',patCols)
patCols3<-sprintf('%s11',patCols)
names(patCols2)<-names(patCols3)<-names(patCols)


#dat<-read.csv('data/Data Master MG, Sept_2.csv')
#dat<-read.csv('data/MM cohort cata master 11.29.2017.csv')
#allDats<-lapply(c('data/for Scott_2017_12_13.csv','data/For Scott - All bulk isol. alpha and beta.csv'),read.csv,stringsAsFactors=FALSE)
#allDats<-lapply(c('data/for Scott_Data Marster.csv'),read.csv,stringsAsFactors=FALSE)
allDats<-lapply(c('data/For Scott November-Data Master.csv'),read.csv,stringsAsFactors=FALSE)
allDats<-lapply(allDats,function(dat)dat[!is.na(dat$ID.for.Publications)&dat$ID.for.Publications!='',])
allCols<-unique(unlist(lapply(allDats,colnames)))
allDats<-lapply(allDats,function(dat){dat[,allCols[!allCols %in% colnames(dat)]]<-NA;dat[,allCols]})
dat<-do.call(rbind,allDats)

#dat<-read.csv('data/for Scott_2017_12_13.csv')
#dat2<-read.csv('data/For Scott - All bulk isol. alpha and beta.csv')

#Standardize BULK naming and catch _BULK with no .
dat$qvoa<-grepl('VOA',dat$ID.for.Publications)
dat$id<-sub('[ _.][Bb]ulk-([IVX]+)','.\\1.bulk',sub('VOA[_ ]','',sub(' +[Bb]ulk|_Bulk|_BULK','.bulk',dat$ID.for.Publications)))
dat$id<-sprintf('%s%s',dat$id,ifelse(dat$qvoa&!grepl('\\.VOA$',dat$id),'.VOA',''))
dat$id<-sub('^Mm','MM',dat$id)
dat$id<-sub('^(MM[0-9]+)\\.([0-9])\\.','\\1.0\\2.',dat$id)
dat$bulk<-grepl('[Bb]ulk',dat$id)
probs<-grepl('^MM[0-9]+\\.[0-9]+\\.bulk',dat$id)
dat[probs,'id']<-withAs(xx=dat[probs,],ave(sub('\\.bulk','',xx$id),sub('\\.bulk','',xx$id),FUN=function(yy)sprintf('%s.%d.bulk',yy,1:length(yy))))
splits<-strsplit(dat$id,'\\.')
probs<-!(sapply(splits,length)==3 & !dat$qvoa & !dat$bulk) & !(sapply(splits,length)==4 & (dat$qvoa|dat$bulk))
if(any(probs))stop('Problem id found')
tmp<-dat[,c('ID.for.Publications','id')]
tmp$diff<-tmp$ID.for.Publications!=tmp$id
write.csv(tmp,'newIds.csv',row.names=FALSE)
dat$sample<-sapply(strsplit(dat$id,'\\.'),function(xx)paste(xx[1:2],collapse='.'))
dat$visit<-sapply(strsplit(dat$sample,'\\.'),'[',2)
dat$virusId<-sapply(strsplit(dat$id,'\\.'),'[',3)
dat$pat<-sub('\\.[^.]+$','',dat$sample)
dat$time<-as.numeric(meta[dat$sample,'DFOSx'])
dat$timeBeforArt<-meta[dat$sample,'daysBeforeArt']
dat$vl<-meta[dat$sample,'vl']
dat$CD4<-meta[dat$sample,'cd4']
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
dat$isLast<-dat$time==ave(dat$time,dat$pat,FUN=max)
write.csv(dat[dat$isFirst|dat$isNadir|dat$qvoa|dat$isSix|dat$isLast|dat$isBetaNadir,c('pat','time','ic50','isFirst','isNadir','isBetaNadir','isSix','isLast','qvoa','beta')],'out/firstNadir.csv')

