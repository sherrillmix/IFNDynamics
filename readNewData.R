library(dnar)
library(lubridate)

if(!exists('meta'))source('readMeta.R')

lay<-matrix(0,nrow=5,ncol=5)
lay[2:4,2:4]<-matrix(1:9,nrow=3,byrow=TRUE)

patCols<-c('MM23'='#e41a1c','MM33'='#4daf4a','MM34'='#984ea3','MM39'='#377eb8','MM40'='#FF7f00','MM14'='#FFD700','MM15'='#f781bf','MM55'='#a65628','MM62'='#708090')
patCols2<-sprintf('%s33',patCols)
patCols3<-sprintf('%s11',patCols)
names(patCols2)<-names(patCols3)<-names(patCols)

#dat<-read.csv('data/Data Master MG, Sept_2.csv')
#dat<-read.csv('data/MM cohort cata master 11.29.2017.csv')
#allDats<-lapply(c('data/for Scott_2017_12_13.csv','data/For Scott - All bulk isol. alpha and beta.csv'),read.csv,stringsAsFactors=FALSE)
allDats<-lapply(c('data/for scott feb.15.2018-4.csv'),read.csv,stringsAsFactors=FALSE)
allDats<-lapply(allDats,function(dat)dat[!is.na(dat$ID.for.Publications)&dat$ID.for.Publications!='',])
allCols<-unique(unlist(lapply(allDats,colnames)))
allDats<-lapply(allDats,function(dat){dat[,allCols[!allCols %in% colnames(dat)]]<-NA;dat[,allCols]})
dat<-do.call(rbind,allDats)

#dat<-read.csv('data/for Scott_2017_12_13.csv')
#dat2<-read.csv('data/For Scott - All bulk isol. alpha and beta.csv')

#EJ79/MM33 not in big spreadsheet
#EJ85/MM39 not listed
#EJ86/MM40 not listed
artDates<-c('MM14'='','MM15'='','MM23'='29/06/06','MM33'='','MM34'='14/09/09','MM39'='','MM40'='','MM55'='','MM62'='')
artDays<-sapply(names(artDates)[artDates!=''],function(ii)withAs(zz=meta[grep(sprintf('%s\\*',ii),meta$id)[1],],dmy(artDates[ii])-mdy(zz$Date)+as.numeric(zz$DFOSx)))


#Standardize BULK naming and catch _BULK with no .
dat$qvoa<-grepl('VOA',dat$ID.for.Publications)
dat$id<-sub('VOA ','',sub(' bulk|_Bulk|_BULK','-BULK',sub('  +bulk','.ZZZ-BULK',sub('\\.([0-9]+)_Bulk','.\\1.BULK',dat$ID.for.Publications))))
dat$sample<-sub('[._][^._]*$','',dat$id)
dat$sample<-sub('\\.([0-9])$','.0\\1',dat$sample)
dat$sample<-sub('\\.([0-9])\\.','.0\\1.',dat$sample)
dat$sample<-sub('Mm','MM',dat$sample)
dat$pat<-sub('\\.[^.]+$','',dat$sample)
dat$time<-as.numeric(meta[dat$sample,'DFOSx'])
dat$vl<-as.numeric(sub('<','',gsub(',','',meta[dat$sample,'VL'])))
dat$CD4<-as.numeric(gsub(',','',meta[dat$sample,'CD4']))
dat$bulk<-grepl('BULK',dat$id)
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
