if(!exists('dat'))source('readNewData.R')

metaCols<-list(
  'CD4 count'="cd4",
  'Viral load'="vl",
  'IFNa2'=c("ifna1", "ifna2"),
  'IFNb'=c("ifnb1", "ifnb2"),
  'IFNo'=c("ifno1", "ifno2"),
  'IFNg'=c("ifng1", "ifng2"), 
  'sCD14'=c("sCD14.well.1", "sCD14.well.2"),
  'LBP'=c("LBP.well.1", "LBP.well.2"), 
  'CRP'=c("CRP.well.1", "CRP.well.2"),
  'BST'="BST2.MFI", 
  'HLA'="HLA.DR", 
  'CD38'="CD38", 
  'IFITM1'=c("IFITM1.RQ__1", "IFITM1.RQ__2", "IFITM1.RQ__3"),
  'IFITM2'=c("IFITM2.RQ__1", "IFITM2.RQ__2", "IFITM2.RQ__3"), 
  'IFITM3'=c("IFITM3.RQ__1", "IFITM3.RQ__2", "IFITM3.RQ__3"),
  'MX2'=c("MX2.RQ__1", "MX2.RQ__2", "MX2.RQ__3"),
  'PKR'=c("PKR.RQ__1", "PKR.RQ__2", "PKR.RQ__3"),
  'TRIM5a'=c("TRIM5a.RQ__1", "TRIM5a.RQ__2", "TRIM5a.RQ__3")
)

virusCols<-c(
  'IFNa2 IC50'='ic50',
  'IFNa2 Vres'='vres',
  'IFNb IC50'='beta',
  'IFNb Vres'='betaVres'
)


nVar<-length(metaCols)
nVirusVar<-length(virusCols)
pdf('out/allVars.pdf',height=20,width=8)
for(ii in sort(unique(dat$pat))){
  par(mar=c(0,4,0,.1))
  layout(matrix(1:(nVirusVar+nVar+1),ncol=1),height=rep(c(1,.6),c(nVar+nVirusVar,1)))
  thisDat<-dat[dat$pat==ii,]
  thisMeta<-comboMeta[comboMeta$mm==ii,]
  qvoaTime<-min(thisDat[thisDat$qvoa,'time'],Inf)
  thisMeta<-thisMeta[thisMeta$time<qvoaTime,]
  thisDat<-thisDat[thisDat$time<qvoaTime,]
  xlabs<-pretty(comboMeta$time[!comboMeta$qvoa]/7)
  for(jj in names(virusCols)){
    plot(thisDat$time,thisDat[,virusCols[jj]],log='y',xaxt='n',ylim=range(dat[,virusCols[jj]],na.rm=TRUE),las=1,cex.axis=.5,ylab=jj,cex.lab=.9,xlab='',yaxt='n',xlim=range(comboMeta$time[!comboMeta$qvoa]))
    means<-10^tapply(log10(thisDat[,virusCols[jj]]),thisDat$time,mean,na.rm=TRUE)
    hasMean<-!is.nan(means)
    lines(as.numeric(names(means)[hasMean]),means[hasMean])
    logAxis(2,las=1,cex.axis=.5)
    text(mean(par('usr')[1:2]),10^(par('usr')[4]-.1*diff(par('usr')[3:4])),jj,col='#00000033',xpd=NA,adj=c(.5,1))
    if(jj==names(virusCols)[1])text(sum(par('usr')[1:2]*c(.02,.98)),10^(par('usr')[4]-.1*diff(par('usr')[3:4])),ii,col='#00000033',xpd=NA,adj=c(1,1),cex=3)
    axis(1,xlabs*7,rep('',length(xlabs)))
  }
  for(jj in names(metaCols)){
    plot(rep(thisMeta$time,length(metaCols[[jj]])),unlist(thisMeta[,metaCols[[jj]],drop=FALSE]),ylim=range(comboMeta[,metaCols[[jj]]],na.rm=TRUE),las=1,cex.axis=.5,ylab=jj,cex.lab=.9,xlab='',xaxt='n',xlim=range(comboMeta$time[!comboMeta$qvoa]))
    means<-apply(thisMeta[,metaCols[[jj]],drop=FALSE],1,mean,na.rm=TRUE)
    hasMean<-!is.nan(means)
    if(any(hasMean))lines(thisMeta$time[hasMean],means[hasMean])
    text(mean(par('usr')[1:2]),(par('usr')[4]-.1*diff(par('usr')[3:4])),jj,col='#00000033',xpd=NA,adj=c(.5,1))
    axis(1,xlabs*7,rep('',length(xlabs)))
  }
  axis(1,xlabs*7,xlabs)
  text(mean(par('usr')[1:2]),par('usr')[3]-diff(par('usr')[3:4])*.35,'Time after onset of symptoms (weeks)',cex=2,adj=c(.5,1),xpd=NA)
}
dev.off()


for(dayCut in c(0,180)){
pdf(sprintf('out/virus_vs_meta%s.pdf',ifelse(dayCut==0,'',sprintf('_after%ddays',dayCut))),height=10,width=10)
  for(jj in names(metaCols)){
    for(ii in names(virusCols)){
      xlim<-range(10^apply(log10(comboMeta[,metaCols[[jj]],drop=FALSE]),1,mean,na.rm=TRUE),na.rm=TRUE)
      ylim<-range(10^tapply(log10(dat[,virusCols[[ii]]]),list(dat$pat,dat$time),mean,na.rm=TRUE),na.rm=TRUE)
      #par(mar=c(0,4,0,.1))
      #layout(matrix(1:,ncol=1),width=rep(c(1,.6),c(3,1)),height=rep(c(1,.6),c(3,1)))
      plot(1,1,type='n',xlim=xlim,ylim=ylim,xlab=jj,ylab=ii,log='yx',yaxt='n',xaxt='n')
      logAxis(2,las=1)
      logAxis(1,addExtra=TRUE)
      allMeans<-data.frame('mMean'=-99,'vMean'=-99)[0,]
      for(pat in sort(unique(dat$pat))){
        thisDat<-dat[dat$pat==pat&!dat$qvoa&dat$time>dayCut,]
        thisMeta<-comboMeta[comboMeta$mm==pat&!comboMeta$qvoa&comboMeta$time>dayCut,]
        vMeans<-10^tapply(log10(thisDat[!thisDat$qvoa,virusCols[ii]]),thisDat[!thisDat$qvoa,'time'],mean,na.rm=TRUE)
        mMeans<-10^apply(log10(thisMeta[,metaCols[[jj]],drop=FALSE]),1,mean,na.rm=TRUE)
        vMeans<-vMeans[as.character(thisMeta$time)]
        isGood<-!is.na(mMeans)&!is.na(vMeans)
        #,col='#00000011'
        arrows(mMeans[isGood][-sum(isGood)],vMeans[isGood][-sum(isGood)],mMeans[isGood][-1],vMeans[isGood][-1],length=.1,col=sprintf('%s33',patCols[pat]))
        points(mMeans[isGood],vMeans[isGood],pch=21,bg=patCols[pat])
        legend(convertLineToUser(2,4),convertLineToUser(2,1),names(patCols),pt.bg=patCols,pch=21,inset=c(.01,-.01),ncol=5,xpd=NA,xjust=1,bty='n')
        allMeans<-rbind(allMeans,data.frame('mMean'=mMeans[isGood],'vMean'=vMeans[isGood]))
      }
      title(main=sprintf('Correlation r2=%0.3f, p=%0.3f',cor(allMeans$mMean,allMeans$vMean),suppressWarnings(cor.test(allMeans$mMean,allMeans$vMean,method='spearman'))$p.value))
    }
  }
dev.off()
}

