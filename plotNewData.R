library(vipor)
if(!exists('dat'))source('readNewData.R')
lay<-matrix(0,nrow=5,ncol=5)
lay[2:4,2:4]<-matrix(1:9,nrow=3,byrow=TRUE)
lowerP24Limit<-60

condenseArrows<-function(dat,ic50Name,clinicalName,xlab,ylab,xlog=FALSE){
  par(mar=c(0,0,0,0))
  layout(lay,width=c(.34,rep(1,3),.01),height=c(.01,rep(1,3),.4))
  counter<-1
  for(ii in sort(unique(dat$pat))){
    ylim<-range(dat[,ic50Name],na.rm=TRUE)*10^(c(-.02,.2)*diff(log10(range(dat[,ic50Name],na.rm=TRUE))))
    plot(dat[dat$pat==ii,clinicalName],dat[dat$pat==ii,ic50Name],log=sprintf('%sy',ifelse(xlog,'x','')),bg=patCols[ii],pch=21,las=1,xlab='',ylab=ylab,yaxt='n',xlim=range(dat[,clinicalName],na.rm=TRUE)*c(1,1.1),ylim=ylim,xaxt='n',cex=1.3)
    title(ii,line=-1)
    if(counter>6&!xlog)axis(1,pretty(dat[,clinicalName],n=3),cex.axis=1.2,mgp=c(2.75,.7,0))
    if(counter>6&xlog)logAxis(1,las=1,cex.axis=1.1,mgp=c(3,.7,0))
    if(counter%%3==1)logAxis(2,las=1,cex.axis=1.1,mgp=c(3,.7,0))
    if(counter==4)text(ifelse(xlog,function(xx)10^xx,function(xx)xx)(par('usr')[1]-.27*diff(par('usr')[1:2])),10^mean(par('usr')[3:4]),ylab,srt=90,xpd=NA,cex=2)
    if(counter==8)text(ifelse(xlog,function(xx)10^xx,function(xx)xx)(mean(par('usr')[1:2])),10^(par('usr')[3]-.3*diff(par('usr')[3:4])),xlab,xpd=NA,cex=2)
    times<-withAs(xx=dat[dat$pat==ii,],cbind('clinic'=tapply(xx[,clinicalName],xx$time,mean),'ic50'=10^tapply(log10(xx[,ic50Name]),xx$time,mean,na.rm=TRUE)))
    withAs(times=times[!is.na(times[,'clinic'])&!is.na(times[,'ic50']),],arrows(times[-nrow(times),'clinic'],times[-nrow(times),'ic50'],times[-1,'clinic'],times[-1,'ic50'],col='#00000033',length=.09))
    #if(counter==9)legend(par('usr')[2]-diff(par('usr')[1:2])*.05,10^(par('usr')[3]-diff(par('usr')[3:4])*.26),c('Quadratic regression','95% confidence interval','95% prediction interval'),col=c(patCols[1],NA,NA),pt.bg=c(NA,patCols2[1],patCols3[1]),lty=c(1,NA,NA),pch=c(NA,22,22),border=NA,pt.cex=3.2,cex=1.2,xjust=1,yjust=1,xpd=NA)
    counter<-counter+1
  }
}
pdf('out/CD4_vs_IC50_condense.pdf',width=7,height=4)
  condenseArrows(dat[!dat$qvoa,],'ic50','CD4','CD4 count','Interferon alpha 2 IC50 (pg/ml)') 
  condenseArrows(dat[!dat$qvoa,],'beta','CD4','CD4 count','Interferon beta IC50 (pg/ml)') 
dev.off()
pdf('out/VL_vs_IC50_condense.pdf',width=7,height=4)
  condenseArrows(dat[!dat$qvoa,],'ic50','vl','Viral load','Interferon alpha 2 IC50 (pg/ml)',xlog=TRUE) 
  condenseArrows(dat[!dat$qvoa,],'beta','vl','Viral load','Interferon beta IC50 (pg/ml)',xlog=TRUE) 
dev.off()



pdf('out/CD4_vs_IC50.pdf',width=6,height=4)
  par(mar=c(4,4,.1,.1))
  plot(dat$CD4,dat$ic50,log='y',pch=21,bg=patCols[dat$pat],xlab='CD4 count',ylab='Interferon alpha IC50',las=1,yaxt='n')
  legend('bottomright',names(patCols),pch=21,pt.bg=patCols,inset=.01)
  logAxis(las=1)
  #layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0),nrow=2,byrow=TRUE))
  par(mar=c(3.6,4.2,1.1,.1),mfrow=c(3,3))
  for(ii in unique(dat$pat)){
    plot(dat$CD4[dat$pat==ii],dat$ic50[dat$pat==ii],log='y',bg=patCols[ii],pch=21,las=1,main=ii,xlab='',ylab='INFa IC50',yaxt='n',xlim=range(dat$CD4,na.rm=TRUE),ylim=range(dat$ic50,na.rm=TRUE))
    title(xlab='CD4 count',mgp=c(2,1,0))
    logAxis(las=1)
    times<-withAs(xx=dat[dat$pat==ii,],cbind('CD4'=tapply(xx$CD4,xx$time,mean),'ic50'=tapply(xx$ic50,xx$time,mean,na.rm=TRUE)))
    arrows(times[-nrow(times),'CD4'],times[-nrow(times),'ic50'],times[-1,'CD4'],times[-1,'ic50'],col='#00000033',length=.075)
    withAs(times=times[!is.na(times[,'CD4']),],arrows(times[-nrow(times),'CD4'],times[-nrow(times),'ic50'],times[-1,'CD4'],times[-1,'ic50'],col='#00000033',length=.075))
  }
  #layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0),nrow=2,byrow=TRUE))
  par(mar=c(3.6,4.2,1.1,.1),mfrow=c(3,3))
  for(ii in unique(dat$pat)){
    plot(dat$CD4[dat$pat==ii],dat$beta[dat$pat==ii],log='y',bg=patCols[ii],pch=21,las=1,main=ii,xlab='',ylab='INFb IC50',yaxt='n',xlim=range(dat$CD4,na.rm=TRUE),ylim=range(dat$beta,na.rm=TRUE))
    title(xlab='CD4 count',mgp=c(2,1,0))
    logAxis(las=1)
    times<-withAs(xx=dat[dat$pat==ii,],cbind('CD4'=tapply(xx$CD4,xx$time,mean),'beta'=tapply(xx$beta,xx$time,mean,na.rm=TRUE)))
    withAs(times=times[!is.na(times[,'CD4']),],arrows(times[-nrow(times),'CD4'],times[-nrow(times),'beta'],times[-1,'CD4'],times[-1,'beta'],col='#00000033',length=.075))
  }
  par(mar=c(4,4,1.1,.1),mfrow=c(1,1))
  withAs(xx=dat[dat$time>35*7,],plot(xx$CD4,xx$ic50,log='y',pch=21,bg=patCols[xx$pat],xlab='CD4 count',ylab='Interferon alpha IC50',las=1,yaxt='n',main="Only >35 weeks",xlim=range(dat$CD4,na.rm=TRUE),ylim=range(dat$ic50,na.rm=TRUE)))
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
  par(mar=c(0,0,0,0))
  layout(lay,width=c(.225,rep(1,3),.45),height=c(.01,rep(1,3),.25))
  counter<-1
  for(ii in sort(unique(dat$pat))){
    #,bg=sprintf('%s33',patCols[xx$pat])
    withAs(xx=dat[dat$pat==ii,],plot(xx$time/7,xx[,var],bg='#0000FF33',pch=ifelse(xx$qvoa,23,ifelse(xx$bulk,22,21)),las=1,log=logAddY,yaxt='n',xlab='',ylab='',xlim=range(dat$time/7),ylim=range(dat[,var],na.rm=TRUE),col='#00000033',xaxt='n',yaxt='n'))
    title(ii,line=-1)
    timeMeans<-10^withAs(xx=dat[dat$pat==ii,],tapply(log10(xx[,var]),xx$time,mean,na.rm=TRUE))
    timeMeans<-timeMeans[order(as.numeric(names(timeMeans)))]
    #,col=patCols[ii]
    lines(as.numeric(names(timeMeans[!is.na(timeMeans)]))/7,timeMeans[!is.na(timeMeans)],col='#0000FF99',lwd=2)
    thisDat<-dat[dat$pat==ii,]
    if(counter>6)axis(1,pretty(range(dat$time/7)),cex.axis=1.2,mgp=c(2.75,.7,0))
    if(counter%%3==1)logAxis(2,las=1,cex.axis=1.1,mgp=c(3,.7,0),col.ticks='blue',col.axis='blue')
    if(counter==4)text(par('usr')[1]-.18*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),lab,srt=90,xpd=NA,cex=2,col='Blue')
    if(counter==8)text(mean(par('usr')[1:2]),10^(par('usr')[3]-.2*diff(par('usr')[3:4])),'Weeks after onset of symptoms',xpd=NA,cex=2)
    par(new=TRUE)
    withAs(xx=unique(dat[dat$pat==ii&!is.na(dat$vl),c('time','vl')]),plot(xx$time/7,xx$vl,type='l',log=logAddY,yaxt='n',xlab='',ylab='',xlim=range(dat$time/7),ylim=range(dat$vl,na.rm=TRUE),xaxt='n',col='#00000077'))
    withAs(xx=unique(dat[dat$pat==ii&!is.na(dat$vl),c('time','vl')]),points(xx$time/7,xx$vl,col='black',pch='.',cex=5))
    if(counter%%3==0)logAxis(4,las=1)
    if(counter==6)text(par('usr')[2]+.14*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'Viral load (copies/ml)',srt=-90,xpd=NA,cex=2)
    if(counter==6)text(par('usr')[2]+.42*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'CD4 count (cells/mm3)',srt=-90,xpd=NA,cex=2,col='red')
    par(new=TRUE)
    withAs(xx=unique(dat[dat$pat==ii&!is.na(dat$CD4),c('time','CD4')]),plot(xx$time/7,xx$CD4,type='l',log=logAdd,yaxt='n',xlab='',ylab='',xlim=range(dat$time/7),ylim=range(dat$CD4,na.rm=TRUE),xaxt='n',col='red'))
    withAs(xx=unique(dat[dat$pat==ii&!is.na(dat$CD4),c('time','CD4')]),points(xx$time/7,xx$CD4,col='red',pch='.',cex=5))
    if(counter%%3==0)axis(4,pretty(dat$CD4),mgp=c(0,7,6),las=1,col='red',col.axis='red')
    counter<-counter+1
  }
}
pdf('out/time_vs_alpha.pdf',width=12,height=8)
  plot3vars('ic50','IFNa IC50',dat[!dat$qvoa,])
dev.off()
pdf('out/time_vs_all_beta.pdf',width=12,height=5)
  plot3vars('beta','IFNb IC50',dat[!dat$qvoa,])
dev.off()
pdf('out/time_vs_all_rep.pdf',width=12,height=5)
  plot3vars('replication','Replicative capacity',dat[!dat$qvoa,])
dev.off()
pdf('out/time_vs_all_alphaVres.pdf',width=12,height=5)
  plot3vars('vres','IFNa Vres',dat[!dat$qvoa,])
dev.off()
pdf('out/time_vs_all_betaVres.pdf',width=12,height=5)
  plot3vars('betaVres','IFNb Vres',dat[!dat$qvoa,])
dev.off()


pdf('out/VL_vs_IC50.pdf',width=6,height=4)
  par(mar=c(4,4,.1,.1))
  plot(dat$vl,dat$ic50,log='xy',bg=patCols[dat$pat],pch=21,las=1,xaxt='n',yaxt='n',xlab='Viral load',ylab='Interferon alpha IC50')
  logAxis(1)
  logAxis(las=1)
  legend('bottomright',names(patCols),pch=21,pt.bg=patCols,inset=.01)
  #layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0),nrow=2,byrow=TRUE))
  par(mar=c(3.6,4.2,1.1,.1),mfrow=c(3,3))
  for(ii in unique(dat$pat)){
    plot(dat$vl[dat$pat==ii],dat$ic50[dat$pat==ii],log='xy',bg=patCols[ii],pch=21,las=1,main=ii,xlab='',ylab='INFa IC50',xaxt='n',yaxt='n',xlim=range(dat$vl,na.rm=TRUE),ylim=range(dat$ic50,na.rm=TRUE))
    title(xlab='Viral load',mgp=c(2,1,0))
    logAxis(1)
    logAxis(las=1)
    times<-withAs(xx=dat[dat$pat==ii,],cbind('vl'=tapply(xx$vl,xx$time,mean),'ic50'=tapply(xx$ic50,xx$time,mean,na.rm=TRUE)))
    withAs(times=times[!is.na(times[,'vl']),],arrows(times[-nrow(times),'vl'],times[-nrow(times),'ic50'],times[-1,'vl'],times[-1,'ic50'],col='#00000033',length=.075))
  }
  #layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0),nrow=2,byrow=TRUE))
  par(mar=c(3.6,4.2,1.1,.1),mfrow=c(3,3))
  for(ii in unique(dat$pat)){
    plot(dat$vl[dat$pat==ii],dat$beta[dat$pat==ii],log='xy',bg=patCols[ii],pch=21,las=1,main=ii,xlab='',ylab='INFb IC50',xaxt='n',yaxt='n',xlim=range(dat$vl,na.rm=TRUE),ylim=range(dat$beta,na.rm=TRUE))
    title(xlab='Viral load',mgp=c(2,1,0))
    logAxis(1)
    logAxis(las=1)
    times<-withAs(xx=dat[dat$pat==ii,],cbind('vl'=tapply(xx$vl,xx$time,mean),'beta'=tapply(xx$beta,xx$time,mean,na.rm=TRUE)))
    withAs(times=times[!is.na(times[,'vl']),],arrows(times[-nrow(times),'vl'],times[-nrow(times),'beta'],times[-1,'vl'],times[-1,'beta'],col='#00000033',length=.075))
  }
  par(mar=c(4,4,1.1,.1),mfrow=c(1,1))
  withAs(xx=dat[dat$time>35*7,],plot(xx$vl,xx$ic50,log='xy',pch=21,bg=patCols[xx$pat],xlab='Viral load',ylab='Interferon alpha IC50',las=1,yaxt='n',xaxt='n',main="Only >35 weeks",xlim=range(dat$vl,na.rm=TRUE),ylim=range(dat$ic50,na.rm=TRUE)))
  legend('topright',names(patCols),pch=21,pt.bg=patCols,inset=.01)
  logAxis(1)
  logAxis(las=1)
  par(mfrow=c(3,3),mar=c(3.6,4.2,1.1,.1))
  for(ii in unique(dat$pat)){
    withAs(xx=dat[dat$time>35*7&dat$pat==ii,],plot(xx$vl,xx$ic50,log='xy',pch=21,bg=patCols[xx$pat],xlab='',ylab='Interferon alpha IC50',las=1,yaxt='n',xaxt='n',main=sprintf("%s Only >35 weeks",ii),xlim=range(dat$vl,na.rm=TRUE),ylim=range(dat$ic50,na.rm=TRUE)))
    title(xlab='Viral load',mgp=c(2,1,0))
    logAxis(1)
    logAxis(las=1)
  }
dev.off()

simpleFits<-withAs(dat=dat[!dat$qvoa,],lapply(unique(dat$pat),function(xx)lm(I(log(ic50))~time+time2,dat=dat[dat$pat==xx&!is.na(dat$ic50),])))
simpleFitsBeta<-withAs(dat=dat[!dat$qvoa,],lapply(unique(dat$pat),function(xx)lm(I(log(beta))~time+time2,dat=dat[dat$pat==xx&!is.na(dat$beta),])))
names(simpleFits)<-names(simpleFitsBeta)<-unique(dat$pat)

pdf('out/predictions.pdf',width=4.5,height=4.5)
for(pats in list(unique(dat$pat),unique(dat$pat[!dat$pat %in% meta2$id]),unique(dat$pat[dat$pat %in% meta2$id]))){
  par(mar=c(3.5,3.7,.1,.1))
  for(ii in names(ifnVars)){
    isVres<-grepl('Vres',ii)
    plot(dat$time/7,dat[,ifnVars[ii]],las=1,yaxt='n',log=ifelse(isVres,'','y'),type='n',xlab='',ylab=ii,yaxt=ifelse(isVres,'s','n'),mgp=c(2.5,1,0))
    title(xlab='Time following onset of symptoms (weeks)',mgp=c(2.2,1,0))
    if(!isVres)logAxis(2,las=1)
    #logAxis(1,addExtra=TRUE,exponent=FALSE)
    patTimeMeans<-tapply(dat[,ifnVars[ii]],list(dat$pat,dat$time),mean,na.rm=TRUE)
    patTimeSd<-tapply(dat[,ifnVars[ii]],list(dat$pat,dat$time),sd,na.rm=TRUE)
    patTimeSe<-patTimeSd/tapply(dat[,ifnVars[ii]],list(dat$pat,dat$time),function(xx)sum(!is.na(xx)))
    for(xx in pats){
      thisDat<-dat[dat$pat==xx,]
      thisDat$target<-thisDat[,ifnVars[ii]]
      if(!isVres)thisDat$target<-log10(thisDat$target)
      else thisDat$target<-log(thisDat$target/100/(1-thisDat$target/100))
      thisFit<-lm(target~time+time2,data=thisDat)
      fakeDays<-(min(dat[dat$pat==xx,'time'])):(max(dat[dat$pat==xx,'time'])+50)
      predIc50<-predict(thisFit,data.frame('time'=fakeDays,'time2'=fakeDays^2,'logTime'=log(fakeDays),'logTime2'=log(fakeDays)^2,'logTime3'=log(fakeDays)^3,'logTime4'=log(fakeDays)^4),interval='confidence')
      if(!isVres)predIc50<-10^predIc50
      else predIc50<-1/(1+exp(-predIc50))*100
      #predIc50<-predict(thisFit,data.frame('time'=fakeDays,'time2'=fakeDays^2),interval='confidence')
      lines(fakeDays/7,predIc50[,'fit'],col=patCols[xx])
      polygon(c(fakeDays/7,rev(fakeDays)/7),c(predIc50[,'lwr'],rev(predIc50[,'upr'])),col=patCols2[xx],border=NA)
    }
    for(xx in pats){
      #segments(as.numeric(colnames(patTimeMeans))/7,patTimeMeans[xx,]+1.96*patTimeSd[xx,],as.numeric(colnames(patTimeMeans))/7,patTimeMeans[xx,]-1.96*patTimeSd[xx,],pch='-',col=patCols[xx],lwd=3)
      segments(as.numeric(colnames(patTimeMeans))/7,patTimeMeans[xx,]+patTimeSe[xx,]*1.96,as.numeric(colnames(patTimeMeans))/7,patTimeMeans[xx,]+patTimeSe[xx,]*-1.96,col=patCols[xx],lwd=2)
      points(as.numeric(colnames(patTimeMeans))/7,patTimeMeans[xx,],pch=21,bg=patCols[xx],cex=1)
    }
    legend('top',names(patCols),col=patCols,inset=.01,ncol=3,lty=1,lwd=2,bty='n')
  }
}
dev.off()

plotPredictions<-function(pat,var,dat,thisFit,artDate,customColors,ylab,bulkHighlight=FALSE){
    xlim<-range(dat$time[dat$pat==pat]/7)
    ylim<-range(dat[,var],na.rm=TRUE)
    plot(1,1,yaxt='n',log='y',type='n',xlab='',ylab=ylab,mgp=c(2.75,1,0),main=pat,xlim=xlim,ylim=ylim)
    if(!is.na(artDate)){
      rect(srtDate/7,10^par('usr')[3],par('usr')[2],10^par('usr')[4],col='#00000022',border=NA)
      text(mean(c(artDate/7,par('usr')[2])),10^(par('usr')[4]-diff(par('usr')[3:4])*.05),'ART treatment',xpd=NA)
    }
    title(xlab='Time following onset of symptoms (weeks)',mgp=c(2.4,1,0))
    logAxis(2,las=2)
    fakeDays<-(min(dat[dat$pat==pat,'time'])):(max(dat[dat$pat==pat&!dat$qvoa,'time'])+50)
    fakeDf<-data.frame('time'=fakeDays,'time2'=fakeDays^2,'logTime'=log(fakeDays),'logTime2'=log(fakeDays)^2,'logTime3'=log(fakeDays)^3,'logTime4'=log(fakeDays)^4)
    predIc50<-predict(thisFit,fakeDf,interval='confidence')
    #predIc50<-predict(thisFit,data.frame('time'=fakeDays,'time2'=fakeDays^2),interval='confidence')
    lines(fakeDays/7,exp(predIc50[,'fit']),col=patCols[pat])
    polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols2[pat],border=NA)
    thisDat<-dat[dat$pat==pat,]
    thisDat<-thisDat[order(thisDat$bulk),]
    predIc50<-predict(thisFit,fakeDf,interval='prediction')
    polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols3[pat],border=NA)
    xPos<-thisDat$time/7+5*offsetX(log(thisDat[,var]),thisDat$time/7)
    points(xPos,thisDat[,var],pch=21+thisDat$bulk,bg=if(customColor)customCols[thisDat$sample,'color'] else sprintf('%sDD',patCols[pat]),cex=1.5+ifelse(bulkHighlight&thisDat$bulk,.5,0),lwd=ifelse(bulkHighlight&thisDat$bulk,4,1),col=ifelse(bulkHighlight&thisDat$bulk,'yellow','black'))
}
pdf('out/indivPredictions.pdf',width=5,height=4)
  par(mar=c(3.5,3.8,1.5,.6))
  for(customColor in c(FALSE,TRUE)){
    for(xx in rownames(patTimeMeans)){
      for(bulkHighlight in 0:2){
        selector<-!dat$qvoa&!is.na(dat$ic50)
        if(bulkHighlight==0)selector<-selector&!dat$bulk
        if(any(selector&dat$pat==xx))plotPredictions(xx,'ic50',dat[selector,],simpleFits[[xx]],NA,customColor,'Interferon alpha 2 IC50 (pg/ml)',bulkHighlight=bulkHighlight==1)
        #plotPredictions(xx,'ic50',dat,simpleFits[[xx]],artDfosx[xx],customColor,'Interferon alpha 2 IC50')
    }
  }
  }
dev.off()
pdf('out/indivPredictions_beta.pdf',width=5,height=4)
  par(mar=c(3.5,3.8,1.5,.6))
  for(customColor in c(FALSE,TRUE)){
    for(xx in rownames(patTimeMeans)){
      for(bulkHighlight in 0:2){
        selector<-!dat$qvoa&!is.na(dat$beta)
        if(bulkHighlight==0)selector<-selector&!dat$bulk
        if(any(selector&dat$pat==xx))plotPredictions(xx,'beta',dat[selector,],simpleFitsBeta[[xx]],NA,customColor,'Interferon beta IC50 (pg/ml)',bulkHighlight=bulkHighlight==1)
      }
    }
  }
dev.off()

plotCondenseIfn<-function(dat,ic50,ylab,simpleFits,showLegend=TRUE){
  par(mar=c(0,0,0,0))
  layout(lay,width=c(.25,rep(1,3),.01),height=c(.01,rep(1,3),1.04))
  counter<-1
  for(ii in sort(unique(dat$pat))){
    plot(dat$time/7,ic50,yaxt='n',log='y',bg=patCols[dat$pat],pch=21,type='n',xlab='',ylab=ylab,xaxt='n',cex=1.4)
    title(ii,line=-1)
    if(counter>6)axis(1,pretty(dat$time/7),cex.axis=1.2,mgp=c(2.75,.7,0))
    if(counter%%3==1)logAxis(2,las=1,cex.axis=1.1,mgp=c(3,.7,0))
    title(xlab='Time (weeks)',mgp=c(2,1,0))
    thisDat<-dat[dat$pat==ii&!dat$qvoa,]
    thisIc50<-ic50[dat$pat==ii&!dat$qvoa]
    thisFit<-lm(I(log(thisIc50))~time+time2,dat=thisDat)
    fakeDays<-(min(thisDat$time)):(max(thisDat$time)+50)
    fakeDf<-data.frame('time'=fakeDays,'time2'=fakeDays^2,'logTime'=log(fakeDays),'logTime2'=log(fakeDays)^2,'logTime3'=log(fakeDays)^3,'logTime4'=log(fakeDays)^4)
    predIc50<-predict(thisFit,fakeDf,interval='confidence')
    lines(fakeDays/7,exp(predIc50[,'fit']),col=patCols[ii])
    polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols2[ii],border=NA)
    predIc50<-predict(thisFit,fakeDf,interval='prediction')
    polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols3[ii],border=NA)
    points(thisDat$time/7,thisIc50,pch=21+thisDat$bulk,bg=patCols[ii])
    if(counter==4)text(par('usr')[1]-.19*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),ylab,srt=90,xpd=NA,cex=2)
    if(counter==8)text(mean(par('usr')[1:2]),10^(par('usr')[3]-.4*diff(par('usr')[3:4])),'Weeks after onset of symptoms',xpd=NA,cex=2)
    if(counter==9&showLegend)legend(par('usr')[2]-diff(par('usr')[1:2])*.05,10^(par('usr')[3]-diff(par('usr')[3:4])*.26),c('Quadratic regression','95% confidence interval','95% prediction interval','Limiting dilution isolate','Bulk isolate'),col=c(patCols[1],NA,NA,'black','black'),pt.bg=c(NA,patCols2[1],patCols3[1],patCols[1],patCols[1]),lty=c(1,NA,NA,NA,NA),pch=c(NA,22,22,21,22),border=NA,pt.cex=c(3.2,3.2,3.2,1.4,1.4),cex=1.2,xjust=1,yjust=1,xpd=NA)
    counter<-counter+1
  }
}
pdf('out/indivPredict_alpha_condense.pdf',width=9,height=5,useDingbats=FALSE)
plotCondenseIfn(dat[!dat$qvoa,],dat$ic50[!dat$qvoa],ylab='Interferon alpha 2 IC50 (pg/ml)',simpleFits)
dev.off()
pdf('out/indivPredict_beta_condense.pdf',width=9,height=5,useDingbats=FALSE)
  plotCondenseIfn(dat[!dat$qvoa,],dat$beta[!dat$qvoa],ylab='Interferon beta IC50 (pg/ml)',simpleFitsBeta)
  #plotCondenseIfn(dat[!dat$qvoa,],dat$beta[!dat$qvoa],ylab='Interferon beta IC50 (pg/ml)',simpleFitsBeta,showLegend=FALSE)
dev.off()


pdf('out/subjects.pdf',width=11,height=8)
  par(mar=c(4.3,4,1.5,4.5),mfrow=c(3,3))
  #layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0),nrow=2,byrow=TRUE))
  for(ii in unique(dat$pat)){
    withAs(xx=dat[dat$pat==ii,],plot(xx$time/7,xx$CD4,pty='l',las=1,log='',xlab='',ylab='',main=ii,xlim=range(dat$time/7),ylim=range(dat$CD4,na.rm=TRUE),col='blue',type='l',lwd=2))
    text(convertLineToUser(3.3,2),mean(par('usr')[3:4]),'CD4 count',srt=90,xpd=NA,col='blue')
    title(xlab='Time (weeks)',mgp=c(2,1,0))
    thisDat<-dat[dat$pat==ii,]
    par(new=TRUE)
    withAs(xx=unique(dat[dat$pat==ii,c('time','vl')]),plot(xx$time/7,xx$vl,type='l',log='y',yaxt='n',xlab='',ylab='',xlim=range(dat$time/7),ylim=range(dat$vl,na.rm=TRUE),xaxt='n',col='red',lwd=2))
    logAxis(4,mgp=c(3,1,0),las=1)
    text(convertLineToUser(3,4),10^mean(par('usr')[3:4]),'Viral load',srt=-90,xpd=NA,col='red')
  }
dev.off()

pdf('out/subjects_condense_new.pdf',width=9,height=4)
  par(mar=c(0,0,0,0))
  layout(lay,width=c(.25,rep(1,3),.27),height=c(.01,rep(1,3),.42))
  counter<-1
  for(ii in sort(unique(dat$pat))){
    xlim<-range(c(dat$time/7,lastDfosx/7))
    withAs(xx=compiledMeta[compiledMeta$mm==ii&!is.na(compiledMeta$cd4)&compiledMeta$time<=lastDfosx[ii],],plot(xx$time/7,xx$cd4,pty='l',las=1,log='',xlab='',ylab='',xlim=xlim,ylim=range(dat$CD4,na.rm=TRUE)+c(-30,90),col='blue',type='l',lwd=2,xaxt='n',yaxt='n'))
    #title(sprintf('%s %s',ii,ifelse(ii %in% rownames(founders),sprintf(' (%s)',founders[ii,'tf']),'')),line=-1)
    title(ii,line=-1)
    if(counter>6)(axis(1,pretty(compiledMeta$time/7),cex.axis=1.2))
    title(xlab='Time (weeks)',mgp=c(2,1,0))
    if(counter%%3==0)(axis(4,pretty(compiledMeta$cd4,n=5),las=1,col.axis='blue',cex.axis=1.1))
    par(new=TRUE)
    thisDat<-unique(compiledMeta[compiledMeta$mm==ii&!is.na(compiledMeta$vl)&compiledMeta$time<=lastDfosx[ii],c('time','vl')])
    plot(thisDat$time/7,thisDat$vl,type='n',log='y',yaxt='n',xlab='',ylab='',xlim=xlim,ylim=range(dat$vl,na.rm=TRUE),xaxt='n',col='red',lwd=2)
    reduceDat<-thisDat[c(TRUE,!sapply(2:(nrow(thisDat)-1),function(zz)all(thisDat[zz+-1:1,'vl']<=50)),TRUE),]
    #connects two <50 or big gap to <50
    isDashed<-(reduceDat$vl[-nrow(reduceDat)]<=lowerP24Limit&reduceDat$vl[-1]<=lowerP24Limit)|(reduceDat$vl[-1]<=lowerP24Limit&reduceDat$time[-1]-reduceDat$time[-nrow(reduceDat)]>120)
    segments(reduceDat$time[-nrow(reduceDat)]/7,reduceDat$vl[-nrow(reduceDat)],reduceDat$time[-1]/7,reduceDat$vl[-1],col='red',lwd=2,lty=ifelse(isDashed,2,1))
    if(counter==6)text(par('usr')[2]+.22*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'CD4 count (cells/mm3)',srt=-90,xpd=NA,col='blue',cex=2)
    if(counter%%3==1)logAxis(2,mgp=c(3,1,0),las=1,col.axis='red',cex.axis=1.3)
    if(counter==4)text(par('usr')[1]-.2*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'Viral load (copies/ml)',srt=90,xpd=NA,col='red',cex=2)
    if(counter==8)text(mean(par('usr')[1:2]),10^(par('usr')[3]-.32*diff(par('usr')[3:4])),'Weeks after onset of symptoms',xpd=NA,cex=2)
    counter<-counter+1
    if(ii %in% names(artDfosx)&&!is.na(artDfosx[ii])){
      rect(artDfosx[ii]/7,10^par('usr')[3],lastDfosx[ii]/7,10^par('usr')[4],col='#00000022',border=NA)
      #text(mean(c(artDfosx[ii]/7,par('usr')[2])),10^(par('usr')[4]-diff(par('usr')[3:4])*.3),'ART treatment',xpd=NA)
    }
    abline(v=lastDfosx[ii]/7,lty=2)
  }
dev.off()


pdf('out/data.pdf',width=6,height=4.5)
  par(mar=c(3.3,4,0.1,0.1))
  for(ii in names(ifnVars)){
    isVres<-grepl('Vres',ii)
    for(jj in 1:3){
      plot(dat$time/7,dat[,ifnVars[ii]],las=1,log=ifelse(isVres,'','y'),xlab='',ylab=ii,pch=21+dat$bulk,bg=patCols[dat$pat],yaxt=ifelse(isVres,'s','n'))
      title(xlab='Time (weeks)',mgp=c(2.2,1,0))
      if(!isVres)logAxis(2,las=1)
      legend('top',names(patCols),pch=21,pt.bg=patCols,inset=.01,ncol=3)
      thisDat<-dat
      thisDat$target<-thisDat[,ifnVars[ii]]
      if(!isVres)thisDat$target<-log10(thisDat$target)
      thisDat$time3<-thisDat$time^3
      thisDat$time4<-thisDat$time^4
      #thisFit<-lm(target~time+time2+time3+time4,data=thisDat)
      thisFit<-lm(target~time+time2,data=thisDat)
      print(summary(step(thisFit)))
      fakeDays<-1:(max(dat$time)+50)
      fakeDat<-data.frame('time'=fakeDays,'time2'=fakeDays^2,'time3'=fakeDays^3,'time4'=fakeDays^4,'logTime'=log(fakeDays),'logTime2'=log(fakeDays)^2,'logTime3'=log(fakeDays)^3,'logTime4'=log(fakeDays)^4)
      predIc50<-predict(thisFit,fakeDat,interval='confidence')
      if(!isVres)predIc50<-10^predIc50
      if(jj>=2){
        lines(fakeDays/7,(predIc50[,'fit']))
        polygon(c(fakeDays/7,rev(fakeDays)/7),c((predIc50[,'lwr']),rev((predIc50[,'upr']))),col='#00000033',border=NA)
        predIc50<-predict(thisFit,fakeDat,interval='prediction')
        if(!isVres)predIc50<-10^predIc50
        if(jj==3)polygon(c(fakeDays/7,rev(fakeDays)/7),c((predIc50[,'lwr']),rev((predIc50[,'upr']))),col='#00000011',border=NA)
      }
    }
  }
dev.off()


pdf('test.pdf')
selector<-!apply(is.na(dat[,ifnVars]),1,any)
tmp<-dat[selector,ifnVars]
pca<-prcomp(tmp,scale.=TRUE)
biplot(pca)
dev.off()



plotVoa<-function(dat,ic50,ylab,addLegend=TRUE){
  par(mar=c(2.9,3.8,.1,.1))
  pats<-sort(unique(dat$pat[dat$qvoa]))
  for(ii in pats){
    plot(dat$time/7,ic50,yaxt='n',log='y',bg=patCols[dat$pat],pch=21,type='n',xlab='',ylab=ylab,mgp=c(2.8,.6,0),cex=1.2)
    logAxis(las=1)
    title(ii,line=-1)
    title(xlab='Time (weeks)',mgp=c(1.8,1,0))
    thisDat<-dat[dat$pat==ii&!dat$qvoa,]
    thisVoa<-dat[dat$pat==ii&dat$qvoa,]
    thisIc50<-ic50[dat$pat==ii&!dat$qvoa]
    thisVoaIc50<-ic50[dat$pat==ii&dat$qvoa]
    thisFit<-lm(I(log(thisIc50))~time+time2,dat=dat[dat$pat==ii&!dat$qvoa,])
    fakeDays<-(min(thisDat[,'time'])):(max(thisDat[,'time'])+50)
    fakeDf<-data.frame('time'=fakeDays,'time2'=fakeDays^2,'logTime'=log(fakeDays),'logTime2'=log(fakeDays)^2,'logTime3'=log(fakeDays)^3,'logTime4'=log(fakeDays)^4)
    predIc50<-predict(thisFit,fakeDf,interval='confidence')
    #predIc50<-predict(thisFit,data.frame('time'=fakeDays,'time2'=fakeDays^2),interval='confidence')
    lines(fakeDays/7,exp(predIc50[,'fit']),col=patCols[ii])
    polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols2[ii],border=NA)
    predIc50<-predict(thisFit,fakeDf,interval='prediction')
    polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols3[ii],border=NA)
    points(thisDat$time/7,thisIc50,pch=21+thisDat$bulk,bg=patCols[ii])
    rect(artDfosx[ii]/7,10^par('usr')[3],par('usr')[2],10^par('usr')[4],col='#00000022',border=NA)
    text(mean(c(artDfosx[ii]/7,par('usr')[2])),10^(par('usr')[4]-diff(par('usr')[3:4])*.15),'ART treatment',xpd=NA)
    abline(h=thisVoaIc50,lty=2,col='#00000055')
    points(thisVoa$time/7,thisVoaIc50,cex=1.3,bg=patCols[ii],pch=23,lwd=3)
    if(ii==tail(pats,1)&addLegend)legend('bottomright',c('Limiting dilution','Viral outgrowth'),pch=c(21,23),pt.cex=c(1.2,1.3),pt.lwd=c(1,2),pt.bg=patCols[ii],inset=c(.01,.02),bg='white',cex=.9)
  }
}
pdf('out/voa.pdf',width=6,height=2.4)
  plotVoa(dat,dat$ic50,ylab='Interferon alpha 2 IC50',addLegend=FALSE)
  plotVoa(dat,dat$beta,ylab='Interferon beta IC50')
dev.off()

nadirs<-sapply(simpleFits,function(xx)-xx$coefficients['time']/xx$coefficients['time2']/2)
nadirsBeta<-sapply(simpleFitsBeta,function(xx)-xx$coefficients['time']/xx$coefficients['time2']/2)
comboNadir<-apply(cbind(nadirs,nadirsBeta),1,min)
names(comboNadir)<-names(nadirs)<-names(nadirsBeta)<-names(simpleFits)

condenseArrows<-function(dat,ic50Name,clinicalName,xlab,ylab,xlog=FALSE){
  par(mar=c(0,0,0,0))
  layout(lay,width=c(.34,rep(1,3),.01),height=c(.01,rep(1,3),.4))
  counter<-1
  for(ii in sort(unique(dat$pat))){
    ylim<-range(dat[,ic50Name],na.rm=TRUE)*10^(c(-.02,.2)*diff(log10(range(dat[,ic50Name],na.rm=TRUE))))
    plot(dat[dat$pat==ii,clinicalName],dat[dat$pat==ii,ic50Name],log=sprintf('%sy',ifelse(xlog,'x','')),bg=patCols[ii],pch=21,las=1,xlab='',ylab=ylab,yaxt='n',xlim=range(dat[,clinicalName],na.rm=TRUE)*c(1,1.1),ylim=ylim,xaxt='n',cex=1.3)
    title(ii,line=-1)
    if(counter>6&!xlog)axis(1,pretty(dat[,clinicalName],n=3),cex.axis=1.2,mgp=c(2.75,.7,0))
    if(counter>6&xlog)logAxis(1,las=1,cex.axis=1.1,mgp=c(3,.7,0))
    if(counter%%3==1)logAxis(2,las=1,cex.axis=1.1,mgp=c(3,.7,0))
    if(counter==4)text(ifelse(xlog,function(xx)10^xx,function(xx)xx)(par('usr')[1]-.27*diff(par('usr')[1:2])),10^mean(par('usr')[3:4]),ylab,srt=90,xpd=NA,cex=2)
    if(counter==8)text(ifelse(xlog,function(xx)10^xx,function(xx)xx)(mean(par('usr')[1:2])),10^(par('usr')[3]-.3*diff(par('usr')[3:4])),xlab,xpd=NA,cex=2)
    times<-withAs(xx=dat[dat$pat==ii,],cbind('clinic'=tapply(xx[,clinicalName],xx$time,mean),'ic50'=10^tapply(log10(xx[,ic50Name]),xx$time,mean,na.rm=TRUE)))
    if(sum(!is.na(times[,'clinic'])&!is.na(times[,'ic50']))>1)withAs(times=times[!is.na(times[,'clinic'])&!is.na(times[,'ic50']),],arrows(times[-nrow(times),'clinic'],times[-nrow(times),'ic50'],times[-1,'clinic'],times[-1,'ic50'],col='#00000033',length=.09))
    #if(counter==9)legend(par('usr')[2]-diff(par('usr')[1:2])*.05,10^(par('usr')[3]-diff(par('usr')[3:4])*.26),c('Quadratic regression','95% confidence interval','95% prediction interval'),col=c(patCols[1],NA,NA),pt.bg=c(NA,patCols2[1],patCols3[1]),lty=c(1,NA,NA),pch=c(NA,22,22),border=NA,pt.cex=3.2,cex=1.2,xjust=1,yjust=1,xpd=NA)
    counter<-counter+1
  }
}
pdf('out/CD4_vs_IC50_condense.pdf',width=7,height=4)
  condenseArrows(dat[!dat$qvoa,],'ic50','CD4','CD4 count','Interferon alpha 2 IC50 (pg/ml)') 
  condenseArrows(dat[!dat$qvoa,],'beta','CD4','CD4 count','Interferon beta IC50 (pg/ml)') 
  condenseArrows(dat[!dat$qvoa&dat$time>nadirs[dat$pat]-100,],'ic50','CD4','CD4 count','Interferon alpha 2 IC50 (pg/ml)') 
  condenseArrows(dat[!dat$qvoa&dat$time>nadirsBeta[dat$pat]-100,],'beta','CD4','CD4 count','Interferon beta IC50 (pg/ml)') 
dev.off()
pdf('out/VL_vs_IC50_condense.pdf',width=7,height=4)
  condenseArrows(dat[!dat$qvoa,],'ic50','vl','Viral load','Interferon alpha 2 IC50 (pg/ml)',xlog=TRUE) 
  condenseArrows(dat[!dat$qvoa,],'beta','vl','Viral load','Interferon beta IC50 (pg/ml)',xlog=TRUE) 
  condenseArrows(dat[!dat$qvoa&dat$time>nadirs[dat$pat]-100,],'ic50','vl','Viral load','Interferon alpha 2 IC50 (pg/ml)',xlog=TRUE) 
  condenseArrows(dat[!dat$qvoa&dat$time>nadirsBeta[dat$pat]-100,],'beta','vl','Viral load','Interferon beta IC50 (pg/ml)',xlog=TRUE) 
dev.off()

pdf('out/ifna2_vs_ifnb.pdf',width=5,height=5)
  par(mar=c(3,3.6,.1,.1))
  plot(dat$ic50,dat$beta,log='xy',bg=sprintf('%sDD',patCols[dat$pat]),xaxt='n',yaxt='n',pch=21,ylab='IFN beta IC50 (pg/ml)',xlab='',mgp=c(2.7,1,0),cex=1.2)
  title(xlab='IFN alpha 2 IC50 (pg/ml)',mgp=c(2,1,0))
  logAxis(1,mgp=c(2.5,.8,0))
  logAxis(2,las=1,mgp=c(2.5,.8,0))
  legend('bottomright',names(patCols),pch=21,pt.bg=patCols,inset=.01,ncol=3,x.intersp=.5,pt.cex=1.2)
dev.off()

tmp<-dat[order(dat$pat,dat$time,dat$id,decreasing=TRUE),]
nAlpha<-length(ifna2_ic50)
nBeta<-length(ifnb_ic50)
alphaCols<-rainbow.lab(nAlpha)
pdf('out/ifn_reproducibility.pdf',width=15,height=30)
  par(mar=c(3,5,.1,.1))
  plot(1,1,type='n',ylim=c(1,nrow(tmp))+c(-1,1),xlim=range(tmp$ic50,na.rm=TRUE),log='x',yaxt='n',xlab='IFNa2 IC50 (pg/ml)',xaxt='n',yaxs='i',ylab='')
  axis(2,1:nrow(tmp),tmp$id,las=2,cex.axis=.5)
  abline(h=1:1000,col='#00000011')
  logAxis(1,las=1)
  for(ii in 1:nrow(tmp)){
    ic50s<-as.vector(tmp[ii,ifna2_ic50])
    if(any(!is.na(ic50s)))segments(min(ic50s,na.rm=TRUE),ii,max(ic50s,na.rm=TRUE),ii)
    points(ic50s,rep(ii,nAlpha),pch=21,bg=alphaCols)
  }
dev.off()

