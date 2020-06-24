library(vipor)
if(!exists('dat'))source('readNewData.R')
source('functions.R')

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

simpleFits<-withAs(dat=dat[!dat$qvoa,],lapply(unique(dat$pat),function(xx)lm(I(log(ic50))~time+time2,dat=dat[dat$pat==xx&!is.na(dat$ic50),])))
simpleFitsBeta<-withAs(dat=dat[!dat$qvoa,],lapply(unique(dat$pat),function(xx)lm(I(log(beta))~time+time2,dat=dat[dat$pat==xx&!is.na(dat$beta),])))
simpleFitsRep<-withAs(dat=dat[!dat$qvoa,],lapply(unique(dat$pat),function(xx)lm(I(log(replication))~time+time2,dat=dat[dat$pat==xx&!is.na(dat$replication),])))
names(simpleFits)<-names(simpleFitsBeta)<-unique(dat$pat)



pdf('out/indivPredict_alpha_condense.pdf',width=4,height=8,useDingbats=FALSE)
  plotCondenseIfn(dat[!dat$qvoa,],dat$ic50[!dat$qvoa],ylab='Interferon alpha 2 IC50 (pg/ml)')
dev.off()
pdf('out/indivPredict_beta_condense.pdf',width=4,height=8,useDingbats=FALSE)
  plotCondenseIfn(dat[!dat$qvoa,],dat$beta[!dat$qvoa],ylab='Interferon beta IC50 (pg/ml)')
  #plotCondenseIfn(dat[!dat$qvoa,],dat$beta[!dat$qvoa],ylab='Interferon beta IC50 (pg/ml)',simpleFitsBeta,showLegend=FALSE)
dev.off()
pdf('out/indivPredict_infectivityMedia.pdf',width=4,height=8,useDingbats=FALSE)
  plotCondenseIfn(dat[!dat$qvoa,],dat$infectivityMedia[!dat$qvoa],ylab='Infectivity (IU/ul)')
dev.off()
pdf('out/indivPredict_infectivityDextran.pdf',width=4,height=8,useDingbats=FALSE)
  plotCondenseIfn(dat[!dat$qvoa,],dat$infectivityDextran[!dat$qvoa],ylab='Infectivity with dextran (IU/ul)')
dev.off()
pdf('out/indivPredict_mediaDextran.pdf',width=4,height=8,useDingbats=FALSE)
  plotCondenseIfn(dat[!dat$qvoa,],dat$infectivityDextran[!dat$qvoa]/dat$infectivityMedia[!dat$qvoa],ylab='Enhancement with dextran (fold change)')
dev.off()
pdf('out/indivPredict_replication.pdf',width=4,height=8,useDingbats=FALSE)
  plotCondenseIfn(dat[!dat$qvoa,],dat$replication[!dat$qvoa],ylab='Replication capacity')
dev.off()

pdf('out/indivPredict_ice.pdf',width=4,height=8,useDingbats=FALSE)
  plotCondenseIfn(dat[!dat$qvoa,],dat$iceHalf[!dat$qvoa],ylab='Ice half life',logY=FALSE)
dev.off()

plotVlCd4<-function(thisMeta,main,xlim,cd4Lim,vlLim,xAxis=TRUE,vlAxis=TRUE,cd4Axis=TRUE){
  withAs(xx=thisMeta[!is.na(thisMeta$cd4),],plot(xx$time/7,xx$cd4,pty='l',las=1,log='',xlab='',ylab='',xlim=xlim,ylim=cd4Lim+c(-80,90),col='blue',type='l',lwd=2,xaxt='n',yaxt='n'))
  title(main,line=-1)
  #abline(h=c(350),col='orange',lty=3)
  #abline(v=less350Time[thisMeta$mm[1]]/7,col='orange')
  thisDat<-unique(thisMeta[!is.na(thisMeta$vl),c('time','vl')])
  if(cd4Axis)axis(4,pretty(compiledMeta$cd4,n=5),las=1,col.axis='blue',cex.axis=1.1,mgp=c(3,.6,0),tcl=-.3)
  par(new=TRUE)
  plot(thisDat$time/7,thisDat$vl,type='n',log='y',yaxt='n',xlab='',ylab='',xlim=xlim,ylim=vlLim,xaxt='n',col='red',lwd=2)
  reduceDat<-thisDat[c(TRUE,!sapply(2:(nrow(thisDat)-1),function(zz)all(thisDat[zz+-1:1,'vl']<=50)),TRUE),]
  #connects two <50 or big gap to <50
  isDashed<-(reduceDat$vl[-nrow(reduceDat)]<=lowerP24Limit&reduceDat$vl[-1]<=lowerP24Limit)|(reduceDat$vl[-1]<=lowerP24Limit&reduceDat$time[-1]-reduceDat$time[-nrow(reduceDat)]>120)
  segments(reduceDat$time[-nrow(reduceDat)]/7,reduceDat$vl[-nrow(reduceDat)],reduceDat$time[-1]/7,reduceDat$vl[-1],col='red',lwd=ifelse(isDashed,1.5,2),lty=ifelse(isDashed,3,1))
  if(xAxis)axis(1,pretty(xlim),cex.axis=1.2,tcl=-.3,mgp=c(3,.5,0))
  if(vlAxis)logAxis(2,mgp=c(3,.8,0),las=1,col.axis='red',cex.axis=1.3)
}
for(showLegend in c(TRUE,FALSE)){
pdf(sprintf('out/subjects_condense_new%s.pdf',ifelse(showLegend,'_legend','')),width=3.54,height=8)
  par(mar=c(0,0,0,0))
  layout(lay2,width=c(.47,rep(1,2),.53),height=c(.13,c(1,1,1,.2,1,.2,1),ifelse(showLegend,.8,.28)))
  #layout(lay,width=c(.5,rep(1,2),.5),height=c(.01,rep(1,5),.35))
  counter<-1
  xlim<-range(c(dat$time/7,lastDfosx/7))
  for(ii in patOrder){
    plotVlCd4(compiledMeta[compiledMeta$mm==ii,],ii,xlim,range(compiledMeta$cd4,na.rm=TRUE),range(compiledMeta$vl,na.rm=TRUE),counter>4,counter%%2==1,counter%%2==0)
    #if(counter>=5&counter<9)axis(1,pretty(xlim),rep('',length(pretty(xlim))))
    #title(sprintf('%s %s',ii,ifelse(ii %in% rownames(founders),sprintf(' (%s)',founders[ii,'tf']),'')),line=-1)
    labCex<-1.7
    if(counter==6)text(par('usr')[2]+.41*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),expression(paste('CD4 count (cells/mm'^3,')')),srt=-90,xpd=NA,col='blue',cex=labCex)
    if(counter==5)text(par('usr')[1]-.38*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'Viral load (copies/ml)',srt=90,xpd=NA,col='red',cex=labCex)
    par(lheight=.7)
    if(counter==9)text(max(par('usr')[1:2]),10^(par('usr')[3]-ifelse(showLegend,.33,.22)*diff(par('usr')[3:4])),ifelse(showLegend,'Weeks from onset\nof symptoms','Weeks from onset of symptoms'),xpd=NA,cex=labCex)
    if(counter==10&showLegend){
      xRight<-grconvertX(.785,'ndc','user')
      annotWidth<-grconvertX(.99,'ndc','user')-grconvertX(.96,'ndc','user')
      annotYPos<-rev(grconvertY(seq(.012,.082,length.out=6),'ndc','user'))
      yStep<-abs(diff(log10(annotYPos[2:3])))
      segments(xRight-annotWidth,annotYPos[1:2],xRight,annotYPos[1:2],xpd=NA,col=c('red','blue'))
      #VOA
      segments(rep(xRight-annotWidth/2,3),rep(annotYPos[3],3)/10^(yStep*.4),rep(xRight-annotWidth/2,3)+c(-6,0,6),rep(annotYPos[3],3)*c(1.7,10^(yStep*.8),1.7)/10^(yStep*.4),lwd=1.2,col='purple',xpd=NA)
      #Super
      segments(rep(xRight-annotWidth/2,3)+c(-9,0,9),rep(annotYPos[4],3)*10^c(yStep*.7,0,yStep*.7)/10^(yStep*.35),rep(xRight-annotWidth/2,3)+c(0,9,-9),rep(annotYPos[4],3)*10^c(0,yStep*.7,yStep*.7)/10^(yStep*.35),lwd=1.2,xpd=NA)
      #ARTs
      rect(xRight-annotWidth,annotYPos[c(5:6,6)]/10^(yStep*.3),xRight,c(annotYPos[5:6],annotYPos[6]/10^(yStep*.3))*10^(yStep*.3),xpd=NA,border=NA,col=c('#00000022','#00000011'))
      text(xRight+annotWidth/5,annotYPos,c('Viral load','CD4 count','VOA isolation','Superinfection','cART treatment','AZT treatment'),xpd=NA,adj=0)
      rect(xRight-annotWidth-24,min(annotYPos)/10^(yStep*.7),grconvertX(.995,'ndc','user'),max(annotYPos)*10^(yStep*.7),xpd=NA,border='#00000099')
    }
    if(counter==1)text(grconvertX(-.31,from='npc'),grconvertY(1.1,from='npc'),'A',xpd=NA,adj=c(0,1),cex=2.5)
    if(counter==7)text(grconvertX(-.31,from='npc'),grconvertY(1.12,from='npc'),'B',xpd=NA,adj=c(0,1),cex=2.5)
    if(counter==9)text(grconvertX(-.31,from='npc'),grconvertY(1.12,from='npc'),'C',xpd=NA,adj=c(0,1),cex=2.5)
    counter<-counter+1
    #thisLast<-lastDfosx[ii]
    thisLast<-max(compiledMeta[compiledMeta$mm==ii&(!is.na(compiledMeta$cd4)|!is.na(compiledMeta$vl)),'time'])
    if(ii %in% names(artDfosx)&&!is.na(artDfosx[ii])){
      rect(artDfosx[ii]/7,10^par('usr')[3],thisLast/7,10^par('usr')[4],col='#00000011',border=NA)
      #text(mean(c(artDfosx[ii]/7,par('usr')[2])),10^(par('usr')[4]-diff(par('usr')[3:4])*.3),'ART treatment',xpd=NA)
    }
    if(!is.na(thisSuper<-founders[ii,'superTime'])){
      superTime<-compiledMeta[compiledMeta$mm==ii&compiledMeta$time==thisSuper,]
      baseY<-20
      segments(rep(thisSuper/7,3)+c(-9,0,9),rep(baseY,3)*c(6,2,6),rep(thisSuper/7,3)+c(0,9,-9),rep(baseY,3)*c(2,6,6),lwd=1.2)
    }
    if(nrow(thisVoa<-dat[dat$pat==ii&dat$qvoa,])>0){
      baseY<-6e5 #superTime$vl*3
      for(jj in 1:nrow(thisVoa)){
        segments(rep(thisVoa$time[jj]/7,3),rep(baseY,3),rep(thisVoa$time[jj]/7,3)+c(-7,0,7),rep(baseY,3)*c(2,6,2),lwd=1.2,col='purple')
      }
    }
    if(ii %in% names(aztDfosx)){
      thisAzt<-aztDfosx[[ii]]
      #abline(v=thisAzt[1]/7,lty=2)
      nRects<-31
      vertBreaks<-seq(par('usr')[3],par('usr')[4],length.out=nRects+1)
      rect(thisAzt[1]/7,10^vertBreaks[seq(1,nRects,2)],thisAzt[2]/7,10^vertBreaks[seq(2,nRects+1,2)],col='#00000008',border=NA)
      rect(thisAzt[1]/7,10^vertBreaks[1],thisAzt[2]/7,10^vertBreaks[nRects+1],col='#00000004',border=NA)
      #abline(v=thisAzt[1]/7,lty='F1',col='#999999')
      #points(rep(thisAzt[1]/7,30),10^seq(par('usr')[3],par('usr')[4],length.out=30),pch=21,col=NA,bg='#00000055',cex=.5)
    }
    abline(v=thisLast/7,lty=2,col='#00000099')
  }
dev.off()
}
file.copy('out/subjects_condense_new.pdf','out/Fig._1.pdf',overwrite=TRUE)

pdf('out/weau_cd4Vl_alpha_beta.pdf',width=3*1.1,height=5,useDingbats=FALSE)
  par(mar=c(0,0,0,0))
  layout(matrix(c(0,0,0,0,1,0,0,2,0,0,3,0,0,0,0),nrow=5,byrow=TRUE),width=c(.32,1,.3),height=c(.01,1,1,1,.25))
  #CD4 VL
  plotVlCd4(compiledMeta[compiledMeta$mm=='WEAU',],'WEAU',range(c(dat$time/7)),range(compiledMeta$cd4,na.rm=TRUE),range(compiledMeta$vl,na.rm=TRUE),xAxis=FALSE)
  text(par('usr')[2]+.25*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'CD4 count (cells/mm3)',srt=-90,xpd=NA,col='blue',cex=1.2)
  text(par('usr')[1]-.22*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'Viral load (copies/ml)',srt=90,xpd=NA,col='red',cex=1.2)
  #ic50s
  par(lheight=.75)
  plotPointsLine(rbind(weau,dat[,colnames(weau)]),c(weau$ic50,dat$ic50),'WEAU',ylab='',addTitle=FALSE)
  mtext('Interferon alpha 2\nIC50 (pg/ml)',2,2.7,cex=.9)
  logAxis(2,las=1,cex.axis=1.1,mgp=c(3,.7,0))
  plotPointsLine(rbind(weau,dat[,colnames(weau)]),c(weau$beta,dat$beta),'WEAU',ylab='',addTitle=FALSE)
  axis(1,pretty(dat$time/7),cex.axis=1.2,mgp=c(2.75,.7,0))
  mtext('Interferon beta\nIC50 (pg/ml)',2,2.7,cex=.9)
  logAxis(2,las=1,cex.axis=1.1,mgp=c(3,.7,0))
  mtext('Weeks after onset of symptoms',1,1.9,cex=.75)
dev.off()

pdf('out/subjects_condense_new_noArt.pdf',width=4,height=8.5)
  par(mar=c(0,0,0,0))
  layout(lay2,width=c(.5,rep(1,2),.5),height=c(.01,c(1,1,1,.15,1,.15,1),.35))
  counter<-1
  xlim<-c(min(compiledMeta$time),max(artDfosx)/7+100)
  for(ii in patOrder){
    selector<-compiledMeta$mm==ii&compiledMeta$time<ifelse(is.na(artDfosx[ii]),Inf,artDfosx[ii])
    plotVlCd4(compiledMeta[selector,],ii,xlim,range(compiledMeta$cd4,na.rm=TRUE),range(compiledMeta$vl,na.rm=TRUE),counter>8,counter%%2==1,counter%%2==0)
    #title(sprintf('%s %s',ii,ifelse(ii %in% rownames(founders),sprintf(' (%s)',founders[ii,'tf']),'')),line=-1)
    if(counter==6)text(par('usr')[2]+.4*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'CD4 count (cells/mm3)',srt=-90,xpd=NA,col='blue',cex=2)
    if(counter==5)text(par('usr')[1]-.38*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),'Viral load (copies/ml)',srt=90,xpd=NA,col='red',cex=2)
    if(counter==9)text(max(par('usr')[1:2]),10^(par('usr')[3]-.28*diff(par('usr')[3:4])),'Weeks after onset of symptoms',xpd=NA,cex=2)
    counter<-counter+1
    #thisLast<-lastDfosx[ii]
    thisLast<-max(compiledMeta[compiledMeta$mm==ii&(!is.na(compiledMeta$cd4)|!is.na(compiledMeta$vl)),'time'])
    if(ii %in% names(artDfosx)&&!is.na(artDfosx[ii])){
      rect(artDfosx[ii]/7,10^par('usr')[3],thisLast/7,10^par('usr')[4],col='#00000022',border=NA)
      #text(mean(c(artDfosx[ii]/7,par('usr')[2])),10^(par('usr')[4]-diff(par('usr')[3:4])*.3),'ART treatment',xpd=NA)
    }
    if(!is.na(thisSuper<-founders[ii,'superTime'])){
      superTime<-compiledMeta[compiledMeta$mm==ii&compiledMeta$time==thisSuper,]
      baseY<-20
      #segments(rep(thisSuper/7,3)+c(-9,0,9),rep(baseY,3)*c(6,2,6),rep(thisSuper/7,3)+c(0,9,-9),rep(baseY,3)*c(2,6,6),lwd=1.2)
    }
    abline(v=thisLast/7,lty=2)
  }
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
    #points(thisDat$time/7,thisIc50,pch=21+thisDat$bulk,bg=patCols[ii])
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

targetCols<-c('ic50','vres','beta','betaVres','replication','infectivityDextran','infectivityMedia')
selector<-!apply(is.na(dat[,targetCols]),1,any)
pca<-prcomp(dat[selector,targetCols],scale.=TRUE)
pcaPoints<-t(t(pca$x)/pca$sdev/sqrt(nrow(pca$x))) #figure out the point positions based on scores scaled by standard deviations
importance<-summary(pca)$importance[2,]
cols<-rainbow.lab(max(dat$time),alpha=.6)
pdf('out/pca.pdf',width=10,height=10)
  for(select in list(1:2)){#,2:3,3:4,4:5)){
    xlim <- range(pcaPoints[,select[1]])
    ylim <- range(pcaPoints[,select[2]])
    pcaArrows<-t(t(pca$rotation[,select])*pca$sdev[select]*sqrt(nrow(pca$x)))	#figure out the arrow positions based on loadings scaled by sdev
    xlimArrow<-range(pcaArrows[,1])
    ylimArrow<-range(pcaArrows[,2])
    ratio <- max(xlimArrow/xlim, ylimArrow/ylim)
    par(mar=c(3.5,3.5,.1,.1))
    plot(
      1,1,type='n',las=1,
      #xlim=range(-pcaPoints[,select[1]]),ylim=range(-pcaPoints[,select[2]]),
      xlim=xlim*1.1,ylim=ylim*1.1,mgp=c(2.4,.7,0),
      xlab=sprintf('Principal component %d (%d%% of variance)',select[1],round(importance[select[1]]*100)),
      ylab=sprintf('Principal component %d (%d%% of variance)',select[2],round(importance[select[2]]*100))
    )
    arrows(0,0,pcaArrows[,1]/ratio,pcaArrows[,2]/ratio,length=.1) #draw arrows
    arrowText<-dimnames(pcaArrows)[[1]] #get rownames of loadings for labels
    par(lheight=.7)
    text(pcaArrows/ratio,rownames(pcaArrows),cex=.8) #label the arrows
    points(pcaPoints[,select[1]],pcaPoints[,select[2]],bg=cols[ifelse(dat[selector,'time']<1,1,dat[selector,'time'])],pch=21,cex=1.5)
  }
dev.off()


xx<-withAs(dat=dat[!is.na(dat$replication)&!dat$qvoa,],as.data.frame(do.call(rbind,tapply(dat$replication,paste(dat$pat,dat$time),function(xx)c(length(xx),mean(xx,na.rm=TRUE),sd(xx,na.rm=TRUE),tryCatch(t.test(xx)$conf.int,error=function(e)return(c(NA,NA)))))),stringsAsFactors=FALSE))
colnames(xx)<-c('n','mean','sd','lower','upper')
xx$mm<-sub(' .*$','',rownames(xx))
xx$time<-as.numeric(sub('^[^ ]* ','',rownames(xx)))
xx<-xx[order(xx$mm,xx$time),]
xx<-xx[,c('mm','time','n','mean','sd','lower','upper')]
write.csv(xx,'out/repCap.csv',row.names=FALSE)


flInf<-lapply(patOrder,function(xx){
  thisDat<-dat[dat$pat==xx&!is.na(dat$infectivityDextran),]
  list(first=thisDat[thisDat$time==min(thisDat$time),'infectivityDextran'],last=thisDat[thisDat$time==max(thisDat$time[!thisDat$qvoa]),'infectivityDextran'])
})
firsts<-lapply(flInf,'[[',1)
lasts<-lapply(flInf,'[[',2)
pos<-((1:length(flInf))-1)*3
fPos<-rep(pos,sapply(firsts,length))+1
lPos<-rep(pos,sapply(lasts,length))+2
pdf('out/infectivity_firstLast.pdf',width=7,height=5)
  par(mar=c(1.5,3.5,.1,.1))
  fOff<-vipor::offsetX(unlist(firsts),fPos)
  lOff<-vipor::offsetX(unlist(lasts),lPos)
  plot(c(fPos+fOff,lPos+lOff),c(unlist(firsts),unlist(lasts)),xlab='',ylab='Infectivity (IU/ul with dextran)',log='y',pch=21,xaxt='n',yaxt='n',mgp=c(2.5,.5,0),bg=rep(c('#FF000033','#0000FF33'),c(length(fPos),length(lPos))),cex=1.5)
  dnar::logAxis(las=1)
  for(ii in 1:length(patOrder))axis(1,pos[ii]+1.5,patOrder[ii],mgp=c(3,.5,0))
  abline(v=pos+3,col="#00000033")
dev.off()


