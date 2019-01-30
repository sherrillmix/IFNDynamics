library(dnar)
out<-read.csv('combined/combinedIC50.csv.gz',stringsAsFactors=FALSE)
isCompleteEnv<-regexpr('[^-]-*$',out$Env)>regexpr('[^-]-*$',out$Env[out$id=='B.FR.83.HXB2_LAI_IIIB_BRU_K03455'])-50
envCheck<-out[isCompleteEnv&!is.na(out$alphaIc50),]
envCheck$cg<-sapply(gregexpr('CG',degap(envCheck$Env)),length)
envCheck$lm<-substring(envCheck$Env,noGap2Gap(out$Env[out$id=='B.FR.83.HXB2_LAI_IIIB_BRU_K03455'],87),noGap2Gap(out$Env[out$id=='B.FR.83.HXB2_LAI_IIIB_BRU_K03455'],1086))
envCheck$lmCg<-sapply(gregexpr('CG',degap(envCheck$lm)),length)
uniqPair<-sort(unique(envCheck$pair))
pdf('env_cg_vs_ic50.pdf',width=10,height=10)
  thisDat<-envCheck
  spreadX<-vipor::offsetX(log(thisDat$alphaIc50),thisDat$lmCg,width=.2)
  xlim<-range(envCheck$lmCg)
  ylim<-range(envCheck$alphaIc50)
  withAs(xx=thisDat,plot(xx$lmCg+spreadX,xx$alphaIc50,log='y',main='All viral isolates',ylab='IFNa2 IC50',xlab='Number of CG in LM region of Env',xlim=xlim,ylim=ylim,yaxt='n',pch=21,bg='#FF000055',col='#00000033'))
  logAxis(las=1)
  thisMod<-lm(log(alphaIc50)~lmCg,thisDat)
  fakeX<-seq(min(thisDat$lmCg),max(thisDat$lmCg),.01)
  predY<-exp(predict(thisMod,data.frame('lmCg'=fakeX),int='conf'))
  lines(fakeX,predY[,'fit'])
  polygon(c(fakeX,rev(fakeX)),c(predY[,'lwr'],rev(predY[,'upr'])),col='#00000011',border=NA)
  par(mfrow=c(3,ceiling(length(uniqPair)/3)))
  for(ii in uniqPair){
    thisDat<-envCheck[envCheck$pair==ii,]
    spreadX<-vipor::offsetX(log(thisDat$alphaIc50),thisDat$lmCg,width=.2)
    withAs(xx=thisDat,plot(xx$lmCg+spreadX,xx$alphaIc50,log='y',main=ifelse(grepl('^[0-9]+$',ii),paste(unique(thisDat[order(!thisDat$donor),'patID']),collapse='-'),ii),ylab='IFNa2 IC50',xlab='Number of CG in LM region of Env',xlim=xlim,ylim=ylim,yaxt='n',pch=21,bg='#FF000055',col='#00000033'))
    logAxis(las=1)
    thisMod<-lm(log(alphaIc50)~lmCg,thisDat)
    fakeX<-seq(min(thisDat$lmCg),max(thisDat$lmCg),.01)
    predY<-exp(predict(thisMod,data.frame('lmCg'=fakeX),int='conf'))
    lines(fakeX,predY[,'fit'])
    polygon(c(fakeX,rev(fakeX)),c(predY[,'lwr'],rev(predY[,'upr'])),col='#00000011',border=NA)
  }
dev.off()

pdf('env_cg_vs_time.pdf',width=10,height=10)
  xlim<-range(envCheck$time,na.rm=TRUE)
  ylim<-range(envCheck$lmCg)
  uniqPair<-sort(unique(envCheck$pair[!is.na(envCheck$time)]))
  par(mfrow=c(2,ceiling(length(uniqPair)/2)))
  for(ii in uniqPair){
    thisDat<-envCheck[envCheck$pair==ii,]
    spreadX<-vipor::offsetX(thisDat$time,thisDat$lmCg,width=.2)
    withAs(xx=thisDat,plot(xx$time,xx$lmCg+spreadX,main=ii,xlab='Days after onset of symptoms',ylab='Number of CG in LM region of Env',xlim=xlim,ylim=ylim,pch=21,bg='#FF000055',col='#00000033',las=1))
  }
dev.off()


library(dnar)
arupPos<-noGap2Gap(out$EnvAA[grep('HXB2',out$id)],c(149,177))
lapply(arupPos,function(xx)table(substring(out$EnvAA,xx,xx)))

