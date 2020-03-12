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
    out<-cbind(out,'density'=if(nCells>0)nCells/timeReq else NULL)
    return(list('dat'=out,'wellId'=description(fcs)$`GTI$WELL`))
  })
  wellNums<-as.numeric(sub('.*_([0-9]+).FCS','\\1',fcsFiles))
  samples<-sub('(.*)_[0-9]+.FCS','\\1',fcsFiles)
  names(fcs)<-names(wellNums)<-names(samples)<-fcsFiles
  wellIds<-sapply(fcs,'[[','wellId')
  return(list('fcs'=lapply(fcs,'[[','dat'),'well'=wellNums,'sample'=samples,'wellId'=wellIds))
}
readFcsDf<-function(fcsDir){
  fcs<-readFcsDir(fcsDir)
  out<-as.data.frame(do.call(rbind,fcs[['fcs']]))
  out$well<-rep(fcs[['well']],sapply(fcs[['fcs']],function(xx)ifelse(is.null(xx),0,nrow(xx))))
  out$sample<-rep(fcs[['sample']],sapply(fcs[['fcs']],function(xx)ifelse(is.null(xx),0,nrow(xx))))
  out$wellId<-rep(fcs[['wellId']],sapply(fcs[['fcs']],function(xx)ifelse(is.null(xx),0,nrow(xx))))
  return(out)
}
plotGuava<-function(allDat,outFile,slope=1,int=0,wells=1:96,col1='GRN-HLog',col2='RED-HLog',nrow=8,ncol=12){
  allDat$outlier<-allDat[,col1]-int-slope*allDat[,col2]
  png(outFile,height=nrow*400,width=ncol*500)
    par(mfrow=c(nrow,ncol),mar=c(0,0,0,0))
    for(ii in wells){
      thisDat<-allDat[allDat$well==ii,]
      if(nrow(thisDat)==0){
        plot(1,1,type='n',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
        message('Missing ',ii)
      }else{
        #max is a bit crazy in CD4
        ylim<-quantile(allDat[,col1],c(0,1))
        xlim<-quantile(allDat[allDat[,col2]<250,col2],c(0,1))
        plot(thisDat[,col2],thisDat[,col1],ylim=ylim,xaxt='n',yaxt='n',xlim=xlim,cex=.5,col='#0000FF66')
        fakeDat<-data.frame('tmp'=seq(1:10000),check.names=FALSE)
        colnames(fakeDat)<-col2
        predGrn<-int+fakeDat[,col2]*slope
        lines(fakeDat[,col2],predGrn,lty=2)
        if(!is.null(thisDat$virus)&&thisDat$virus[1]=='No virus'&&any(thisDat$outlier>0))withAs(xx=thisDat[thisDat$outlier>0,],points(xx[,col2],xx[,col1],col='red',cex=1))
      }
      if(!is.null(thisDat$virus)&&!is.null(thisDat$treat))desc<-sprintf('%s\n%s\n',thisDat$virus[1],thisDat$treat[1])
      else desc<-''
      title(main=sprintf('%s%0.2f%%\n%0.1f cells/sec',desc,mean(thisDat$outlier>0)*100,thisDat$density[1]),line=-11,cex.main=3)
    }
  dev.off()
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



virus<-rep(c('SG3','89.6','WEAU','YU2','MM15.14.2C4','MM23.07.2B2'),2)
treat<-rep(c('No drug','AMD','Tak','AMD+Tak'),2)
trop<-readFcsDf('ice/2018-06-18_tropism/')
trop<-cbind(trop,wellToRowCol(trop$well))
trop$virus<-virus[trop$col]
trop$treat<-treat[trop$row]
plotGuava(trop[trop$sample=='2018-06-18-jltr',],'out/flow/2018-06-18-tropism.png',1,200)
plotGuava(trop[trop$sample=='2018-06-18-jltr_dilute',],'out/flow/2018-06-18-tropism2.png',1,200)
plotGuava(trop[trop$sample=='2018-06-19_jltr',],'out/flow/2018-06-18-tropism3.png',1,200)




#int<-200
#slope<-1
#virus<-c('CH40_TF','CH40_6mo','NL4.3','89.6','MM23.1','MM23.8','MM23.13','89.6+CH40_TF')

plotEvo<-function(fcsEvo,main='',maxChange=max(abs(scaledGreen),na.rm=TRUE),int=200,slope=1,virus=c('CH40_TF','CH40_6mo','NL4.3','89.6','MM23.1','MM23.8','MM23.13','89.6+CH40_TF')){
  fcsEvo$outlier<-fcsEvo[,'GRN-HLog']-int-slope*fcsEvo[,'RED-HLog']
  propGreen<-matrix(tapply(fcsEvo$outlier>0,fcsEvo$well,mean)[as.character(1:96)],nrow=8,ncol=12,byrow=TRUE)
  cols<-colorRampPalette(c('blue','white','red'))(100)
  scaledGreen<-log2(propGreen/as.vector(propGreen[,4:6]))
  breaks<-seq(-maxChange-1e-6,maxChange+1e-6,length.out=101)
  par(mar=c(5,8,1.1,.1),lheight=.7)
  image(1:12,1:8,t(scaledGreen),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='',main=main)
  axis(1,1:12,c('UT','A2','BE',rep(c('UT','A2','BE'),each=3)))
  axis(2,1:8,virus,las=1)
  scaleTicks<--floor(maxChange):floor(maxChange)
  dnar::insetScale(breaks,cols,at=scaleTicks,labels=sapply(scaleTicks,function(xx)as.expression(bquote(2^.(xx)))),main='Relative amount of\nGFP-bright JLTR')
  text(rep(1:12,each=8),rep(1:8,12),as.vector(round(2^scaledGreen,2)),cex=.75)
  abline(v=.5+c(3,6,9))
  box()
  return(maxChange)
}

fcsEvo<-readFcsDf('ice/2018-07-09_evoJLTR/')
fcsEvo2<-readFcsDf('ice/2018-08-04_evoJLTR/')
fcsEvo3<-readFcsDf('ice/2018-08-17_evoJLTR/')
fcsEvo4<-readFcsDf('ice/2018-08-24_evoJLTR/')
plotGuava(fcsEvo,'out/evoSort_0709.png',1,200)
plotGuava(fcsEvo2,'out/evoSort_0804.png',1,200)
plotGuava(fcsEvo3,'out/evoSort_0817.png',1,200)
plotGuava(fcsEvo3,'out/evoSort_0824.png',1,200)

#lazy way to figure out maxchange
pdf('out/evoAnalyze.pdf')
  maxChange<-max(
    plotEvo(fcsEvo,'July 9, 2018'),
    plotEvo(fcsEvo2,'August 4, 2018'),
    plotEvo(fcsEvo3,'August 17, 2018'),
    plotEvo(fcsEvo4,'August 24, 2018')
  )
dev.off()
pdf('out/evoAnalyze.pdf')
  plotEvo(fcsEvo,'July 9, 2018',maxChange=maxChange)
  plotEvo(fcsEvo2,'August 4, 2018',maxChange=maxChange)
  plotEvo(fcsEvo3,'August 17, 2018',maxChange=maxChange)
  plotEvo(fcsEvo4,'August 24, 2018',maxChange=maxChange)
dev.off()

xx<-rbind(
  readFcsDf('ice/2019-12-16_p24Test/'),
  readFcsDf('ice/2019-12-13_p24TestAfterClog/'),
  readFcsDf('ice/2019-12-13_p24Test/')
)
xx$well<-xx$wellId
xx$virus<-xx$well
xx$treat<-''
xx$GYHLog<-xx$`YEL-HLog`+xx$`GRN-HLog`
plotGuava(xx,'out/20191216_p24.png',1,200,wells=sort(unique(xx$well[!grepl('E',xx$well)])))
plotGuava(xx,'out/20191216_p24_NIR.png',1,200,wells=sort(unique(xx$well[!grepl('E',xx$well)])),col2='NIR2-HLog')
plotGuava(xx,'out/20191216_p24_NIR_red.png',1,200,wells=sort(unique(xx$well[!grepl('E',xx$well)])),col2='NIR2-HLog',col1='RED-HLog')
plotGuava(xx,'out/20191216_p24_NIR_yellow.png',1,200,wells=sort(unique(xx$well[!grepl('E',xx$well)])),col2='NIR2-HLog',col1='YEL-HLog')
plotGuava(xx,'out/20191216_p24_NIR_greenYellow.png',1,200,wells=sort(unique(xx$well[!grepl('E',xx$well)])),col2='NIR2-HLog',col1='GYHLog')

xx<-readFcsDf('ice/2020-01-10_ifnIntSites/')
xx$well<-xx$wellId
xx$virus<-xx$well
xx$treat<-''
xx$GYHLog<-xx$`YEL-HLog`+xx$`GRN-HLog`
xx$isUninfected<-xx$well %in% c('A07','B07')|grepl('A0[1-6]$',xx$well)
xx$col<-as.numeric(sub('[A-H]','',xx$well))
xx$row<-sub('[0-9]+','',xx$well)
censor<-xx
for(ii in colnames(censor)[grep('GRN|YEL|RED|NIR|RED2|NIR2',colnames(censor))]){
  print(ii)
  cut<-quantile(censor[,ii],.999,na.rm=TRUE)
  print(cut)
  print(summary(censor[,ii]))
  censor[censor[,ii]>cut,ii]<-cut
  print(summary(censor[,ii]))
}

#GYslopeInt<-withAs(yy=censor[censor$isUninfected,],findMostExtremePoint(yy[,'NIR2-HLog'],yy[,'GYHLog'],.01*nrow(yy)))
plotGuava(censor,'out/20200110_p24.png',1,200,wells=sort(unique(censor$well)),ncol=7,nrow=8)
plotGuava(censor,'out/20200110_p24_NIR_greenYellow.png',1,200,wells=sort(unique(censor$well)),col1='NIR2-HLog',col2='GYHLog',ncol=7,nrow=8)
plotGuava(censor,'out/20200110_p24_NIR_greenYellow.png',0,10,wells=sort(unique(censor$well)),col2='FSC-HLin',col1='GRN-HLog',ncol=7,nrow=8)
withAs(censor=censor[sample(which(censor$well %in% c('H07','B04','B01')),20000),],plot(censor$`RED-HLog`,censor$`GRN-HLog`+censor$`YEL-HLog`,col=censor$isUninfected+1))
withAs(censor=censor[sample(which(censor$isUninfected),20000),],points(censor$`RED-HLog`,censor$`GRN-HLog`+censor$`YEL-HLog`,col=censor$isUninfected+1))
withAs(censor=censor[sample(which(censor$well %in% c('H07','B04','B01')),20000),],plot(censor$`FSC-HLog`,censor$`GRN-HLog`,col=censor$isUninfected+1))
withAs(censor=censor[sample(which(censor$isUninfected),20000),],points(censor$`FSC-HLog`,censor$`GRN-HLog`,col=censor$isUninfected+1))
