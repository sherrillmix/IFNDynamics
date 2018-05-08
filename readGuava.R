library(flowCore)

which.maxN<-function(x,n=1){
  which(-x==sort(-x,partial=n)[n])
}

findMostExtremePoint<-function(red,green,nOutlier=0){
  redBase<-quantile(red,.25)
  #greenBase<-quantile(green[red<redBase],.95)
  greenBase<-sort(green[red<redBase],decreasing=TRUE)[nOutlier+1]
  id<-which.maxN(atan2(green-greenBase,red),nOutlier+1)
  return(c('slope'=(green[id]-greenBase)/red[id],'yint'=greenBase,'id'=id))
}

fcsDir<-'ice/2018-04-23_jltr_spin/'
fcsFiles<-list.files(fcsDir,'[0-9]+.FCS$')
fcs<-lapply(fcsFiles,function(xx){
  fcs<-tryCatch(read.FCS(file.path(fcsDir,xx)),error=function(e)return(NULL))
  if(is.null(fcs))return(NULL)
  out<-exprs(fcs)
  return(out)
})
wellNums<-as.numeric(sub('.*_([0-9]+).FCS','\\1',fcsFiles))
samples<-sub('(.*)_[0-9]+.FCS','\\1',fcsFiles)
wellNums[samples=='cd4_spin_dilute2']<-ceiling(wellNums[samples=='cd4_spin_dilute2']/3)
names(fcs)<-names(wellNums)<-names(samples)<-fcsFiles

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
  #allDat$residual<-allDat[,'GRN-HLog']-predict(mod,allDat)
  #allDat$outlier<-allDat[,'GRN-HLog']-predict(mod,allDat,interval='prediction',level=.99999)[,'upr']
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
      title(main=sprintf('%s\n%s\n%0.2f%%',thisDat$virus[1],thisDat$treat[1],mean(thisDat$outlier>0)*100),line=-9,cex.main=3)
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


#xx<-read.FCS('ice/2018-04-23_jltr_spin/cd4_culture_dilute_0001.FCS')
