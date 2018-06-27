library(flowCore)
library(dnar)

which.maxN<-function(x,n=1){
  which(-x==sort(-x,partial=n)[n])
}

findMostExtremePoint<-function(red,green,nOutlier=0,slopeNOutlier=nOutlier){
  redBase<-quantile(red,.25)
  #greenBase<-quantile(green[red<redBase],.95)
  greenBase<-sort(green[red<redBase],decreasing=TRUE)[nOutlier+1]
  id<-which.maxN(atan2(green-greenBase,red),slopeNOutlier+1)
  return(c('slope'=(green[id]-greenBase)/red[id],'yint'=greenBase,'id'=id))
}

findMostExtremePoint2<-function(red,green,nOutlier=0){
  redBase<-quantile(red,.25)
  #greenBase<-quantile(green[red<redBase],.95)
  greenBase<-sort(green[red<redBase],decreasing=TRUE)[nOutlier+1]
  mod<-lm(green~red)
  id<-0
  return(c('slope'=mod$coef[['red']],'yint'=greenBase,'id'=id))
}
wellToRowCol<-function(wells){
  cols<-(wells %%12)
  cols[cols==0]<-12
  rows<-ceiling(wells/12)
  data.frame('col'=cols,'row'=rows)
}



readFcsDir<-function(fcsDir){
  fcsFiles<-list.files(fcsDir,'[0-9]+.FCS$')
  fcs<-lapply(fcsFiles,function(xx){
    fcs<-tryCatch(flowCore::read.FCS(file.path(fcsDir,xx)),error=function(e)return(NULL))
    if(is.null(fcs))return(NULL)
    out<-exprs(fcs)
    timeReq<-as.numeric(diff(as.difftime(unlist(fcs@description[c("$BTIM","$ETIM")])))*60*60)
    nCells<-nrow(out)
    out<-cbind(out,'density'=nCells/timeReq)
    return(out)
  })
  wellNums<-as.numeric(sub('.*_([0-9]+).FCS','\\1',fcsFiles))
  samples<-sub('(.*)_[0-9]+.FCS','\\1',fcsFiles)
  names(fcs)<-names(wellNums)<-names(samples)<-fcsFiles
  return(list('fcs'=fcs,'well'=wellNums,'sample'=samples))
}
readFcsDf<-function(fcsDir){
  fcs<-readFcsDir(fcsDir)
  out<-as.data.frame(do.call(rbind,fcs[['fcs']]))
  out$well<-rep(fcs[['well']],sapply(fcs[['fcs']],function(xx)ifelse(is.null(xx),0,nrow(xx))))
  out$sample<-rep(fcs[['sample']],sapply(fcs[['fcs']],function(xx)ifelse(is.null(xx),0,nrow(xx))))
  return(out)
}


fcs<-readFcsDir('ice/2018-04-23_jltr_spin/')
wellNums<-fcs[['well']]
samples<-fcs[['sample']]
fcs<-fcs[['fcs']]
wellNums[samples=='cd4_spin_dilute2']<-ceiling(wellNums[samples=='cd4_spin_dilute2']/3)

treats<-rep(rep(c("DEAE wash", "DEAE no wash","Poly wash","Poly no wash","CD44 beads","Neat"),2),8)
virus<-list('jltr'=rep(c(rbind(
    c("NL43+GFP (WEAU) 003x", "NL43+GFP (WEAU)150x", "YU2 002x", "YU2 100x", "SG3 004x", "SG3 200x", "CH236 020x", "No virus"),
    c("NL43+GFP (TRO11) 003x", "NL43+GFP (TRO11) 150x", "YU2 004x", "YU2 200x", "SG3 008x", "SG3 400x", "CH236 200x", "No virus")
  )),each=6),
  'cd4'=rep(c(rbind(
      c("NL43+GFP (WEAU) 003x","NL43+GFP (WEAU) 006x", "NL43+GFP (WEAU) 024x", "NL43+GFP (WEAU) 048x", "NL43+GFP (WEAU) 192x","NL43+GFP (WEAU) 384x", "No virus", "No virus"),
    c("NL43+GFP (TRO11) 003x",  "NL43+GFP (TRO11) 006x", "NL43+GFP (TRO11) 024x", "NL43+GFP (TRO11) 048x", "NL43+GFP (TRO11) 192x", "NL43+GFP (TRO11) 384x", "No virus", "No virus")
  )),each=6)
)

desiredPlates<-c('cd4_culture_dilute','cd4_spin_dilute2','jltr_spin_dilute','jltr_culture_dilute')

allPlates<-lapply(desiredPlates,function(plate){
  thisFcs<-fcs[samples==plate]
  allDat<-as.data.frame(do.call(rbind,thisFcs))
  allDat$study<-sub('_.*','',plate)
  allDat$well<-rep(wellNums[samples==plate],sapply(thisFcs,function(xx)ifelse(is.null(nrow(xx)),0,nrow(xx))))
  allDat$virus<-virus[[allDat$study[1]]][allDat$well]
  allDat$treat<-treats[allDat$well]
  #mod<-lm(`GRN-HLog`~`RED-HLog`*treat,data=allDat[allDat[,'RED-HLog']<250&allDat$virus=='No virus',])
  nOutlier<-round(.0001*sum(allDat[,'RED-HLog']<250&allDat$virus=='No virus'))
  slopeInt<-withAs(xx=allDat[allDat[,'RED-HLog']<250&allDat$virus=='No virus',],findMostExtremePoint(xx[,'RED-HLog'],xx[,'GRN-HLog'],nOutlier))
  #allDat$residual<-allDat[,'GRN-HLog']-predict(mod,allDat) #allDat$outlier<-allDat[,'GRN-HLog']-predict(mod,allDat,interval='prediction',level=.99999)[,'upr']
  allDat$outlier<-allDat[,'GRN-HLog']-slopeInt['yint']-slopeInt['slope']*allDat[,'RED-HLog']
  return(allDat)
})
names(allPlates)<-desiredPlates

lapply(allPlates,function(xx)round(tapply(xx$outlier>0,list(xx$virus,xx$treat),mean)*100,2))




for(plate in desiredPlates){
  message(plate)
  allDat<-allPlates[[plate]]
  mod<-lm(`GRN-HLog`~`RED-HLog`*treat,data=allDat[allDat[,'RED-HLog']<250&allDat$virus=='No virus',])
  nOutlier<-round(.0001*sum(allDat[,'RED-HLog']<250&allDat$virus=='No virus'))
  slopeInt<-withAs(xx=allDat[allDat[,'RED-HLog']<250&allDat$virus=='No virus',],findMostExtremePoint(xx[,'RED-HLog'],xx[,'GRN-HLog'],nOutlier))
  selectPoint<-withAs(xx=allDat[allDat[,'RED-HLog']<250&allDat$virus=='No virus',],xx[slopeInt['id'],])
  png(sprintf('out/flow/redGreen_%s.png',plate),height=8*400,width=12*500)
    par(mfrow=c(8,12),mar=c(0,0,0,0))
    for(ii in 1:96){
      thisDat<-allDat[allDat$well==ii,]
      if(nrow(thisDat)==0){
        plot(1,1,type='n',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
        message('Missing ',ii)
      }else{
        #xlab='Green-HLog',ylab='Red-HLog',
        plot(thisDat[,'RED-HLog'],thisDat[,'GRN-HLog'],ylim=range(allDat[,'GRN-HLog']),xaxt='n',yaxt='n',xlim=range(allDat[allDat[,'RED-HLog']<250,'RED-HLog']),cex=.5,col='#0000FF66')
        fakeDat<-data.frame('RED-HLog'=seq(1:10000),'treat'=thisDat$treat[1],check.names=FALSE)
        #predGrn<-predict(mod,fakeDat,interval='prediction',level=.999999)[,'upr']
        predGrn<-slopeInt['yint']+fakeDat$`RED-HLog`*slopeInt['slope']
        lines(fakeDat$`RED-HLog`,predGrn,lty=2)
        if(selectPoint$well==ii)points(selectPoint$`RED-HLog`,selectPoint$`GRN-HLog`,col='red',cex=3)
        if(thisDat$virus[1]=='No virus'&&any(thisDat$outlier>0))withAs(xx=thisDat[thisDat$outlier>0,],points(xx$`RED-HLog`,xx$`GRN-HLog`,col='red',cex=1))
        #abline(mod,lty=2)
      }
      title(main=sprintf('%s\n%s\n%0.2f%%\n%0.1f cells/s',thisDat$virus[1],thisDat$treat[1],mean(thisDat$outlier>0)*100,thisDat$density[1]),line=-9,cex.main=3)
    }
  dev.off()
}

if(FALSE){
  png(sprintf('out/flow/hist_%s.png',plate),height=8*400,width=12*500)
    par(mfrow=c(8,12),mar=c(0,0,0,0))
    for(ii in 1:96){
      thisDat<-allDat[allDat$well==ii,]
      if(nrow(thisDat)==0){
        plot(1,1,type='n',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
      }else{
        thisDat<-as.data.frame(thisDat[thisDat[,'RED-HLog']<250,])
        thisDat$residual<-thisDat[,'GRN-HLog']-predict(mod,list(allDat=thisDat))
        hist(thisDat$residual,xlim=range(allDat$residual),breaks=100,main='',xaxt='n',yaxt='n')
      }
    }
  dev.off()
}





plotPlate<-function(allDat,outFile){
  mod<-lm(`GRN-HLog`~`RED-HLog`*treat,data=allDat[allDat[,'RED-HLog']<250&allDat$virus=='No virus',])
  nOutlier<-round(.001*sum(allDat[,'RED-HLog']<250&allDat$virus=='No virus'))
  slopeInt<-withAs(xx=allDat[allDat[,'RED-HLog']<250&allDat$virus=='No virus',],findMostExtremePoint(xx[,'RED-HLog'],xx[,'GRN-HLog'],nOutlier,round(nOutlier/5)))
  selectPoint<-withAs(xx=allDat[allDat[,'RED-HLog']<250&allDat$virus=='No virus',],xx[slopeInt['id'],])
  allDat$outlier<-allDat[,'GRN-HLog']-slopeInt['yint']-slopeInt['slope']*allDat[,'RED-HLog']
  png(outFile,height=8*400,width=12*500)
    par(mfrow=c(8,12),mar=c(0,0,0,0))
    for(ii in 1:96){
      thisDat<-allDat[allDat$well==ii,]
      if(nrow(thisDat)==0){
        plot(1,1,type='n',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
        message('Missing ',ii)
      }else{
        #max is a bit crazy in CD4
        ylim<-quantile(allDat[,'GRN-HLog'],c(0,.999))
        xlim<-quantile(allDat[allDat[,'RED-HLog']<250,'RED-HLog'],c(0,1))
        plot(thisDat[,'RED-HLog'],thisDat[,'GRN-HLog'],ylim=ylim,xaxt='n',yaxt='n',xlim=xlim,cex=.5,col='#0000FF66')
        fakeDat<-data.frame('RED-HLog'=seq(1:10000),'treat'=thisDat$treat[1],check.names=FALSE)
        predGrn<-slopeInt['yint']+fakeDat$`RED-HLog`*slopeInt['slope']
        lines(fakeDat$`RED-HLog`,predGrn,lty=2)
        if(nrow(selectPoint)>0&&selectPoint$well==ii)points(selectPoint$`RED-HLog`,selectPoint$`GRN-HLog`,col='red',cex=3)
        if(thisDat$virus[1]=='No virus'&&any(thisDat$outlier>0))withAs(xx=thisDat[thisDat$outlier>0,],points(xx$`RED-HLog`,xx$`GRN-HLog`,col='red',cex=1))
      }
      title(main=sprintf('%s\n%s\n%0.2f%%\n%0.1f cells/sec',thisDat$virus[1],thisDat$treat[1],mean(thisDat$outlier>0)*100,thisDat$density[1]),line=-11,cex.main=3)
    }
  dev.off()
}

fcs2<-readFcsDf('ice/2018-05-14_infectionAids')
#t293<-fcs2[fcs2$sample=='293T',]
fcs2<-fcs2[fcs2$sample!='293T',]
fcs2$sample[fcs2$sample=='cd4_noSpin2']<-'cd4_noSpin'
treats<-rep(sprintf('%s %s',rep(c('JLTR','CEM'),each=6),rep(c("DEAE wash", "DEAE no wash","Poly wash","Poly no wash","CD44 beads","Neat"),2)),8)
treats2<-rep(rep(c("DEAE wash", "DEAE no wash","Poly wash","Poly no wash","CD44 beads","Neat"),2),8)
virus<-list('jltr'=rep(c(rbind(
    c("SG3 4x", "SG3 10x", "SG3 100x", "896 sup 3x", "896 sup 8x", "896 sup 100x", "NL43+GFP (WEAU) 3x", "No virus"),
    c("SG3 4x", "SG3 10x", "SG3 100x", "896 sup 3x", "896 sup 8x", "896 sup 100x", "NL43+GFP (WEAU) 3x", "No virus")
  )),each=6),
  'cd4'=rep(c(rbind(
      c("NL43+GFP (WEAU) 003x","NL43+GFP (WEAU) 006x", "NL43+GFP (WEAU) 012x", "NL43+GFP (WEAU) 024x", "NL43+GFP (WEAU) 048x","NL43+GFP (WEAU) 96x","NL43+GFP (WEAU) 192x", "No virus"),
      c("NL43+GFP (TRO11) 003x","NL43+GFP (TRO11) 006x", "NL43+GFP (TRO11) 012x", "NL43+GFP (TRO11) 024x", "NL43+GFP (TRO11) 048x","NL43+GFP (TRO11) 96x","NL43+GFP (TRO11) 192x", "No virus")
  )),each=6)
)
fcs2$virus<-ifelse(grepl('cd4',fcs2$sample),virus[['cd4']][fcs2$well],virus[['jltr']][fcs2$well])
fcs2$treat<-ifelse(grepl('cd4',fcs2$sample),treats2[fcs2$well],treats[fcs2$well])

tapply(fcs2[,'GRN-HLog'],list(fcs2$treat,fcs2$virus,fcs2$sample),mean)
for(ii in unique(fcs2$sample)){
  message(ii)
  plotPlate(fcs2[fcs2$sample==ii,],sprintf('out/flow/infectionAid_%s.png',ii))
}


plotGuava<-function(allDat,outFile,slope=1,int=0){
  allDat$outlier<-allDat[,'GRN-HLog']-int-slope*allDat[,'RED-HLog']
  png(outFile,height=8*400,width=12*500)
    par(mfrow=c(8,12),mar=c(0,0,0,0))
    for(ii in 1:96){
      thisDat<-allDat[allDat$well==ii,]
      if(nrow(thisDat)==0){
        plot(1,1,type='n',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
        message('Missing ',ii)
      }else{
        #max is a bit crazy in CD4
        ylim<-quantile(allDat[,'GRN-HLog'],c(0,1))
        xlim<-quantile(allDat[allDat[,'RED-HLog']<250,'RED-HLog'],c(0,1))
        plot(thisDat[,'RED-HLog'],thisDat[,'GRN-HLog'],ylim=ylim,xaxt='n',yaxt='n',xlim=xlim,cex=.5,col='#0000FF66')
        fakeDat<-data.frame('RED-HLog'=seq(1:10000),check.names=FALSE)
        predGrn<-int+fakeDat$`RED-HLog`*slope
        lines(fakeDat$`RED-HLog`,predGrn,lty=2)
        if(thisDat$virus[1]=='No virus'&&any(thisDat$outlier>0))withAs(xx=thisDat[thisDat$outlier>0,],points(xx$`RED-HLog`,xx$`GRN-HLog`,col='red',cex=1))
      }
      title(main=sprintf('%s\n%s\n%0.2f%%\n%0.1f cells/sec',thisDat$virus[1],thisDat$treat[1],mean(thisDat$outlier>0)*100,thisDat$density[1]),line=-11,cex.main=3)
    }
  dev.off()
}


virus<-rep(c('SG3','89.6','WEAU','YU2','MM15.14.2C4','MM23.07.2B2'),2)
treat<-rep(c('No drug','AMD','Tak','AMD+Tak'),2)
trop<-readFcsDf('ice/2018-06-18_tropism/')
trop<-cbind(trop,wellToRowCol(trop$well))
trop$virus<-virus[trop$col]
trop$treat<-treat[trop$row]
plotGuava(trop[trop$sample=='2018-06-18-jltr',],'out/flow/2018-06-18-tropism.png',1,200)
plotGuava(trop[trop$sample=='2018-06-18-jltr_dilute',],'out/flow/2018-06-18-tropism2.png',1,200)
plotGuava(trop[trop$sample=='2018-06-19_jltr',],'out/flow/2018-06-18-tropism3.png',1,200)





#xx<-read.FCS('ice/2018-04-23_jltr_spin/cd4_culture_dilute_0001.FCS')
