library(dnar)
dat<-read.csv('data/ForScottSep.2017-MG.csv')

meta<-read.csv('data/For Scot, Complete Master table AUG.2017_meta.csv')[,-1:-2]
meta$id<-fillDown(meta$ID)
meta<-meta[meta$Time.Points!='Total number of sequences',]
rownames(meta)<-meta$Time.Points

#Standardize BULK naming and catch _BULK with no .
dat$id<-sub('_Bulk|_BULK','-BULK',sub('\\.([0-9]+)_Bulk','.\\1.BULK',dat$ID.for.Publications))
dat$sample<-sub('[._][^._]*$','',dat$id)
dat$sample<-sub('\\.([0-9])$','.0\\1',dat$sample)
dat$sample<-sub('\\.([0-9])\\.','.0\\1.',dat$sample)
dat$pat<-sub('\\.[^.]+$','',dat$sample)
dat$time<-meta[dat$sample,'DFOSx']
dat$vl<-as.numeric(gsub(',','',meta[dat$sample,'VL']))
dat$CD4<-as.numeric(gsub(',','',meta[dat$sample,'CD4']))
dat$bulk<-grepl('BULK',dat$id)
if(any(is.na(dat$vl)))stop('Missing vl')
if(any(is.na(dat$time)))stop('Missing time')
if(any(is.na(dat$CD4)))stop('Missing CD4')


dat$ic50<-apply(dat[,c('IFNa2..Oct.2016..Pooled.Donor.cells.IC50..pg.ml.','REPEAT....June.2017.....IFNa2.Pooled.Donor.cells.IC50..pg.ml.')],1,mean,na.rm=TRUE)
dat$vres<-apply(dat[,c('IFNa2.....Oct..2016..Pooled.Donor.Vres.at.5.5pg.ml...UT.','REPEAT....June.2017...IFNa2.Pooled.Donor.Vres.at.5.5pg.ml...UT.')],1,mean,na.rm=TRUE)
dat$beta<-dat$IFNbeta..Pooled.Donor.cells.IC50..pg.ml.
dat$betaVres<-dat$IFNbeta..Pooled.Donor.Vres.at.4.400.pg.ml...UT.
if(any(dat$betaVres==0)){
  warning('Beta Vres equal to zero. Setting to lowest value/10.')
  dat[dat$betaVres==0&!is.na(dat$betaVres),'betaVres']<-min(dat[dat$betaVres!=0&!is.na(dat$betaVres),'betaVres'])/10
}
dat$replication<-dat$Replicative.capacity.Pooled.Donor.cells.p24.d7..from.June.2017.repeat.

patCols<-c('MM23'='#FF0000','MM33'='#00FF00','MM34'='#8000FF','MM39'='#0000FF','MM40'='#FF8000')
patCols2<-c('MM23'='#FF000033','MM33'='#00FF0033','MM34'='#8000FF33','MM39'='#0000FF33','MM40'='#FF800033')
cols<-rainbow.lab(length(unique(dat$pat)))
names(cols)<-unique(dat$pat)
pdf('out/CD4_vs_IC50.pdf',width=5,height=5)
  par(mar=c(4,4,.1,.1))
  plot(dat$CD4,dat$ic50,log='y',pch=21,bg=patCols[dat$pat],xlab='CD4 count',ylab='Interferon alpha IC50',las=1,yaxt='n')
  legend('bottomright',names(patCols),pch=21,pt.bg=patCols,inset=.01)
  logAxis(las=1)
  par(mfrow=c(3,2),mar=c(3.6,4.2,1.1,.1))
  for(ii in unique(dat$pat)){
    plot(dat$CD4[dat$pat==ii],dat$ic50[dat$pat==ii],log='y',bg=patCols[ii],pch=21,las=1,main=ii,xlab='',ylab='INFa IC50',yaxt='n',xlim=range(dat$CD4),ylim=range(dat$ic50,na.rm=TRUE))
    title(xlab='CD4 count',mgp=c(2,1,0))
    logAxis(las=1)
  }
  par(mar=c(4,4,1.1,.1),mfrow=c(1,1))
  withAs(xx=dat[dat$time>35*7,],plot(xx$CD4,xx$ic50,log='y',pch=21,bg=patCols[xx$pat],xlab='CD4 count',ylab='Interferon alpha IC50',las=1,yaxt='n',main="Only >35 weeks",xlim=range(dat$CD4),ylim=range(dat$ic50,na.rm=TRUE)))
  legend('topright',names(patCols),pch=21,pt.bg=patCols,inset=.01)
  logAxis(las=1)
dev.off()

pdf('out/time_vs_CD4.pdf',width=5,height=5)
  par(mar=c(4,4,.1,.1))
  par(mfrow=c(1,1))
  plot(dat$time/7,dat$CD4,type='n',bg=patCols[dat$pat],pch=21,las=1,xlab='Time (weeks)',ylab='CD4 count')
  sapply(unique(dat$pat),function(xx)withAs(d=dat[dat$pat==xx,],lines(d$time/7,d$CD4,col=patCols[xx],lwd=2)))
  legend('topright',names(patCols),col=patCols,lty=1,lwd=2,inset=.01)
dev.off()

pdf('out/time_vs_VL.pdf',width=5,height=5)
  par(mar=c(4,4,.1,.1))
  plot(dat$time/7,dat$vl,type='n',bg=patCols[dat$pat],pch=21,las=1,log='y',yaxt='n',xlab='Time (weeks)',ylab='Viral load')
  sapply(unique(dat$pat),function(xx)withAs(d=dat[dat$pat==xx,],lines(d$time/7,d$vl,col=patCols[xx],lwd=2)))
  legend('topright',names(patCols),col=patCols,lty=1,lwd=2)
  logAxis(las=1)
  par(mfrow=c(3,2),mar=c(3.6,4.2,1.1,.1))
dev.off()

plot3vars<-function(var,lab,dat,logX=FALSE){
  logAdd<-ifelse(logX,'x','')
  logAddY<-sprintf('%sy',logAdd)
  par(mfrow=c(3,2),mar=c(3.6,4.2,1.1,11))
  for(ii in unique(dat$pat)){
    withAs(xx=dat[dat$pat==ii,],plot(xx$time/7,xx[,var],bg=patCols[xx$pat],pch=21,las=1,log=logAddY,yaxt='n',xlab='',ylab=lab,main=ii,xlim=range(dat$time/7),ylim=range(dat[,var],na.rm=TRUE)))
    title(xlab='Time (weeks)',mgp=c(2,1,0))
    logAxis(las=1)
    thisDat<-dat[dat$pat==ii,]
    par(new=TRUE)
    withAs(xx=unique(dat[dat$pat==ii,c('time','vl')]),plot(xx$time/7,xx$vl,type='l',log=logAddY,yaxt='n',xlab='',ylab='',xlim=range(dat$time/7),ylim=range(dat$vl,na.rm=TRUE),xaxt='n'))
    logAxis(4,las=1)
    text(convertLineToUser(3,4),10^mean(par('usr')[3:4]),'Viral load',srt=-90,xpd=NA)
    par(new=TRUE)
    withAs(xx=unique(dat[dat$pat==ii,c('time','CD4')]),plot(xx$time/7,xx$CD4,type='l',log=logAdd,yaxt='n',xlab='',ylab='',xlim=range(dat$time/7),ylim=range(dat$CD4,na.rm=TRUE),xaxt='n',col='red'))
    axis(4,pretty(dat$CD4),mgp=c(6,5,4),las=1,col='red',col.axis='red')
    text(convertLineToUser(3,4),10^mean(par('usr')[3:4]),'Viral load',srt=-90,xpd=NA)
    text(convertLineToUser(8,4),mean(par('usr')[3:4]),'CD4 count',srt=-90,xpd=NA,col='red')
  }
  #changing scales
  par(mfrow=c(3,2))
  for(ii in unique(dat$pat)){
    withAs(xx=dat[dat$pat==ii,],plot(xx$time/7,xx[,var],bg=patCols[xx$pat],pch=21,las=1,log=logAddY,yaxt='n',xlab='',ylab=lab,main=ii))
    title(xlab='Time (weeks)',mgp=c(2,1,0))
    logAxis(las=1)
    thisDat<-dat[dat$pat==ii,]
    par(new=TRUE)
    withAs(xx=unique(dat[dat$pat==ii,c('time','vl')]),plot(xx$time/7,xx$vl,type='l',log=logAddY,yaxt='n',xlab='',ylab='',xaxt='n'))
    logAxis(4,las=1)
    text(convertLineToUser(3,4),10^mean(par('usr')[3:4]),'Viral load',srt=-90,xpd=NA)
    par(new=TRUE)
    withAs(xx=unique(dat[dat$pat==ii,c('time','CD4')]),plot(xx$time/7,xx$CD4,type='l',log=logAdd,yaxt='n',xlab='',ylab='',xaxt='n',col='red'))
    axis(4,pretty(dat$CD4[dat$pat==ii]),mgp=c(6,5,4),las=1,col='red',col.axis='red')
    text(convertLineToUser(3,4),10^mean(par('usr')[3:4]),'Viral load',srt=-90,xpd=NA)
    text(convertLineToUser(8,4),mean(par('usr')[3:4]),'CD4 count',srt=-90,xpd=NA,col='red')
  }
  #sapply(unique(dat$pat),function(xx)withAs(d=dat[dat$pat==xx,],lines(d$time/7,d$vl,col=patCols[xx],lwd=2)))
}

pdf('out/time_vs_all.pdf',width=8,height=8)
  plot3vars('ic50','IFNa IC50',dat)
  plot3vars('ic50','IFNa IC50',dat,logX=TRUE)
dev.off()
pdf('out/time_vs_all_beta.pdf',width=8,height=8)
  plot3vars('beta','IFNb IC50',dat)
  plot3vars('beta','IFNb IC50',dat,logX=TRUE)
dev.off()
pdf('out/time_vs_all_rep.pdf',width=8,height=8)
  plot3vars('replication','Replicative capacity',dat)
  plot3vars('replication','Replicative capacity',dat,logX=TRUE)
dev.off()
pdf('out/time_vs_all_vres.pdf',width=8,height=8)
  plot3vars('vres','IFNa Vres',dat)
  plot3vars('vres','IFNa Vres',dat,logX=TRUE)
dev.off()
pdf('out/time_vs_all_betaVres.pdf',width=8,height=8)
  plot3vars('betaVres','IFNb Vres',dat)
  plot3vars('betaVres','IFNb Vres',dat,logX=TRUE)
dev.off()



pdf('out/VL_vs_IC50.pdf',width=5,height=5)
  par(mar=c(4,4,.1,.1))
  plot(dat$vl,dat$ic50,log='xy',bg=patCols[dat$pat],pch=21,las=1,xaxt='n',yaxt='n',xlab='Viral load',ylab='Interferon alpha IC50')
  logAxis(1)
  logAxis(las=1)
  legend('bottomright',names(patCols),pch=21,pt.bg=patCols,inset=.01)
  par(mfrow=c(3,2),mar=c(3.6,4.2,1.1,.1))
  for(ii in unique(dat$pat)){
    plot(dat$vl[dat$pat==ii],dat$ic50[dat$pat==ii],log='xy',bg=patCols[ii],pch=21,las=1,main=ii,xlab='',ylab='INFa IC50',xaxt='n',yaxt='n',xlim=range(dat$vl),ylim=range(dat$ic50,na.rm=TRUE))
    title(xlab='Viral load',mgp=c(2,1,0))
    logAxis(1)
    logAxis(las=1)
    times<-withAs(xx=dat[dat$pat==ii,],cbind('vl'=tapply(xx$vl,xx$time,mean),'ic50'=tapply(xx$ic50,xx$time,mean,na.rm=TRUE)))
    print(times)
    arrows(times[-nrow(times),'vl'],times[-nrow(times),'ic50'],times[-1,'vl'],times[-1,'ic50'],col='#00000033',length=.075)
  }
  par(mar=c(4,4,1.1,.1),mfrow=c(1,1))
  withAs(xx=dat[dat$time>35*7,],plot(xx$vl,xx$ic50,log='xy',pch=21,bg=patCols[xx$pat],xlab='Viral load',ylab='Interferon alpha IC50',las=1,yaxt='n',xaxt='n',main="Only >35 weeks",xlim=range(dat$vl),ylim=range(dat$ic50,na.rm=TRUE)))
  legend('topright',names(patCols),pch=21,pt.bg=patCols,inset=.01)
  logAxis(1)
  logAxis(las=1)
  par(mfrow=c(3,2),mar=c(3.6,4.2,1.1,.1))
  for(ii in unique(dat$pat)){
    withAs(xx=dat[dat$time>35*7&dat$pat==ii,],plot(xx$vl,xx$ic50,log='xy',pch=21,bg=patCols[xx$pat],xlab='',ylab='Interferon alpha IC50',las=1,yaxt='n',xaxt='n',main=sprintf("%s Only >35 weeks",ii),xlim=range(dat$vl),ylim=range(dat$ic50,na.rm=TRUE)))
    title(xlab='Viral load',mgp=c(2,1,0))
    logAxis(1)
    logAxis(las=1)
  }
dev.off()

dat$time2<-dat$time^2
dat$logTime<-log(dat$time)
dat$logTime2<-log(dat$time)^2
dat$logTime3<-log(dat$time)^3
dat$logTime4<-log(dat$time)^4
dat$logVl<-log(dat$vl)
fits<-lapply(unique(dat$pat),function(xx)lm(I(log(ic50))~logTime+logTime2+logVl+CD4,dat=dat[dat$pat==xx&!is.na(dat$ic50),]))
fits<-lapply(unique(dat$pat),function(xx)lm(I(log(ic50))~time+time2+time*logVl+time*CD4,dat=dat[dat$pat==xx&!is.na(dat$ic50),]))
simpleFits<-lapply(unique(dat$pat),function(xx)lm(I(log(ic50))~time+time2,dat=dat[dat$pat==xx&!is.na(dat$ic50),]))
names(simpleFits)<-names(fits)<-unique(dat$pat)
dat$pred<-NA
for(ii in names(fits)){
  preds<-predict(simpleFits[[ii]])
  dat[dat$pat==ii&!is.na(dat$ic50),'pred']<-preds
}


pdf('out/predictions.pdf',width=5,height=5)
  plot(dat$time/7,dat$ic50,yaxt='n',log='y',bg=patCols[dat$pat],pch=21,type='n',xlab='Time following onset of symptoms (weeks)',ylab='Interferon alpha IC50')
  #logAxis(1,addExtra=TRUE,exponent=FALSE)
  logAxis(2,las=2)
  patTimeMeans<-tapply(dat$ic50,list(dat$pat,dat$time),mean,na.rm=TRUE)
  patTimeSd<-tapply(dat$ic50,list(dat$pat,dat$time),sd,na.rm=TRUE)
  patTimeSe<-patTimeSd/tapply(dat$ic50,list(dat$pat,dat$time),function(xx)sum(!is.na(xx)))
  for(xx in rownames(patTimeMeans)){
    thisFit<-simpleFits[[xx]]
    fakeDays<-(min(dat[dat$pat==xx,'time'])):(max(dat[dat$pat==xx,'time'])+50)
    predIc50<-predict(thisFit,data.frame('time'=fakeDays,'time2'=fakeDays^2,'logTime'=log(fakeDays),'logTime2'=log(fakeDays)^2,'logTime3'=log(fakeDays)^3,'logTime4'=log(fakeDays)^4),interval='confidence')
    #predIc50<-predict(thisFit,data.frame('time'=fakeDays,'time2'=fakeDays^2),interval='confidence')
    lines(fakeDays/7,exp(predIc50[,'fit']),col=patCols[xx])
    polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols2[xx],border=NA)
  }
  for(xx in rownames(patTimeMeans)){
    points(as.numeric(colnames(patTimeMeans))/7,patTimeMeans[xx,],pch=21,bg=patCols[xx],cex=1)
    #segments(as.numeric(colnames(patTimeMeans))/7/1.02,patTimeMeans[xx,],as.numeric(colnames(patTimeMeans))/7*1.02,patTimeMeans[xx,],pch='-',col=patCols[xx],lwd=3)
    #segments(as.numeric(colnames(patTimeMeans))/7,patTimeMeans[xx,]+patTimeSe[xx,]*1.96,as.numeric(colnames(patTimeMeans))/7,patTimeMeans[xx,]+patTimeSe[xx,]*-1.96,col=patCols[xx],lwd=2)
  }
  legend('top',names(patCols),pch=21,pt.bg=patCols,inset=.01,ncol=2)
dev.off()


#withAs(xx=dat[dat$time>35*7,],plot(ave(xx$vl,xx$pat,FUN=function(xx)(xx-min(xx,na.rm=TRUE))/max(xx-min(xx,na.rm=TRUE),na.rm=TRUE)),xx$ic50,bg=patCols[xx$pat],log='y',pch=21,cex=2))

# regenerate old plot
# show old plot + new data
# new data plots with fits 
## viral load plot (split by 35 weeks)
## cd4 plot (split by 35 weeks)
# split out bulks "_Bulk"
# plot Replicative over time
# plot everything for beta

# bayesian model 
