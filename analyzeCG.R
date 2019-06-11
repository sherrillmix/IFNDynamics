library(dnar)
out<-read.csv('combined/combinedIC50.csv.gz',stringsAsFactors=FALSE)
out$isRef<-out$id=='B.FR.83.HXB2_LAI_IIIB_BRU_K03455'
isCompleteEnv<-regexpr('[^-]-*$',out$Env)>regexpr('[^-]-*$',out$Env[out$isRef])-50
envCheck<-out[isCompleteEnv&!is.na(out$alphaIc50),]
envCheck$cg<-sapply(gregexpr('CG',degap(envCheck$Env)),length)
envCheck$lm<-substring(envCheck$Env,noGap2Gap(out$Env[out$isRef],87),noGap2Gap(out$Env[out$isRef],1086))
envCheck$lmCg<-sapply(gregexpr('CG',degap(envCheck$lm)),length)
uniqPair<-sort(unique(envCheck$pair))

pdf('out/env_cg_vs_ic50.pdf',width=10,height=10)
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

pdf('out/env_cg_vs_time.pdf',width=10,height=10)
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


pdf('out/env_cg_vs_vl.pdf',width=10,height=10)
  xlim<-range(envCheck$viralLoad,na.rm=TRUE)
  ylim<-range(envCheck$lmCg)
  uniqPair<-sort(unique(envCheck$pair[!is.na(envCheck$viralLoad)]))
  par(mfrow=c(2,ceiling(length(uniqPair)/2)))
  for(ii in uniqPair){
    thisDat<-envCheck[envCheck$pair==ii,]
    spreadX<-vipor::offsetX(thisDat$viralLoad,thisDat$lmCg,width=.2)
    withAs(xx=thisDat,plot(xx$viralLoad,xx$lmCg+spreadX,main=ii,xlab='Viral load',ylab='Number of CG in LM region of Env',xlim=xlim,ylim=ylim,pch=21,bg='#FF000055',col='#00000033',las=1,log='x',xaxt='n'))
    dnar::logAxis(1)
  }
dev.off()



library(dnar)
gappedPos<-c(150,177,178,447)+c(1,1,1,1)
ungappedPos<-gap2NoGap(out$EnvAA[out$isRef],gappedPos)
lapply(ungappedPos,function(xx)table(substring(out$EnvAA,xx,xx)))
plotOut<-out
plotOut$niceLab<-sprintf('%s%s',ifelse(plotOut$isRef,'HXB2',plotOut$patID),ifelse(grepl('^[0-9]+$',plotOut$pair)&!is.na(plotOut$pair),sprintf(' (%s %s)',ifelse(plotOut$donor,'D','R'),plotOut$pair),''))
patOrder<-unique(plotOut$niceLab[order(plotOut$isRef,plotOut$pair,plotOut$donor,plotOut$niceLab)])
plotOut$groupLabel<-factor(plotOut$niceLab,levels=c(patOrder))
pdf('out/envAA.pdf',height=10,width=10)
par(mar=c(3.5,3,2,6))
dnaplotr::plotAA(rep(plotOut$EnvAA,ifelse(plotOut$isRef,10,1)),refSeq=plotOut$EnvAA[plotOut$isRef],groups=rep(plotOut$groupLabel,ifelse(plotOut$isRef,10,1)),xlab='HXB2 Env position')
axis(3,151,150,cex.axis=.8)
axis(3,178,177,cex.axis=.8)
axis(3,448,447,cex.axis=.8)
#dnaplotr::plotAA(rep(substring(plotOut$EnvAA,1,200),ifelse(plotOut$isRef,10,1)),refSeq=plotOut$EnvAA[plotOut$isRef],groups=rep(plotOut$groupLabel,ifelse(plotOut$isRef,10,1)))
#axis(3,c(150,177))
xStart<-gap2NoGap(out$EnvAA[out$isRef],125)
dnaplotr::plotAA(rep(substring(plotOut$EnvAA,125,202),ifelse(plotOut$isRef,10,1)),refSeq=substring(plotOut$EnvAA[plotOut$isRef],125,202),groups=rep(plotOut$groupLabel,ifelse(plotOut$isRef,10,1)),xStart=xStart,xlab='HXB2 Env position')
axis(3,c(150,177)+1-125+xStart,as.character(c(150,177)))
xStart<-gap2NoGap(out$EnvAA[out$isRef],400)
dnaplotr::plotAA(rep(substring(plotOut$EnvAA,400,500),ifelse(plotOut$isRef,10,1)),refSeq=substring(plotOut$EnvAA[plotOut$isRef],400,500),groups=rep(plotOut$groupLabel,ifelse(plotOut$isRef,10,1)),xStart=xStart,xlab='HXB2 Env position')
axis(3,c(447+1)-400+xStart,'447')
dev.off()

pdf('out/position_vs_ic50.pdf',width=10)
for(ii in gappedPos){
  xx<-out[!is.na(out$alphaIc50),]
  aas<-substring(xx$EnvAA,ii,ii)
  aaPos<-structure(1:length(unique(aas)),.Names=sort(unique(aas)))
  xOffset<-vipor::offsetX(log(xx$alphaIc50),aas,width=.35)
  pairCols<-structure(rainbow.lab(length(unique(xx$pair)),alpha=.5),.Names=sort(unique(xx$pair)))
  plot(aaPos[aas]+xOffset,xx$alphaIc50,log='y',xaxt='n',main=sprintf('Position %d (HXB2 %d)',ii-1,ungappedPos[gappedPos==ii]),xlab='',ylab='IFNa2 IC50',yaxt='n',xlim=c(.5,max(aaPos)+.5),xaxs='i',pch=ifelse(out$study=='gondim',21,22),bg=pairCols[xx$pair])
  axis(1,aaPos,names(aaPos))
  logAxis(las=1)
  abline(v=.5+0:100,col='#00000033')
  legend(grconvertX(0.01,'device','user'),grconvertY(0.01,'device','user'),names(pairCols),pch=ifelse(sapply(names(pairCols),function(yy)xx[xx$pair==yy,'study'][1]=='gondim'),21,22),pt.bg=pairCols,xpd=NA,xjust=0,yjust=0,ncol=6,bty='n')
}
dev.off()

marvinTargets<-c('NGKDSAIRLT','LTGRANE','RDNSTAG');shilpaTargets<-c('GKAFLSNW','TSDEGNVPLC','ARMDNTPCKV')
test<-lapply(1:3,function(ii){
  sapply(1:600,function(jj){
    thisAA<-substring(out$EnvAA,jj,jj)
    marvin<-all(sapply(strsplit(marvinTargets[ii],'')[[1]],function(xx)any(grepl(xx,thisAA[out$study=='gondim']))))
    shilpa<-all(sapply(strsplit(shilpaTargets[ii],'')[[1]],function(xx)any(grepl(xx,thisAA[out$study=='iyer']))))
    return(c(marvin,shilpa))
  })
})
sapply(test,which,arr.ind=TRUE)
gap2NoGap(out$EnvAA[out$isRef],c(151,178,179))


table(substring(out[,'GagAA'],146+50-1,146+50-1),out$subtype)
table(substring(out[,'GagAA'],146+120-1,146+120-1))
lanl<-dnar::read.fa('lanl/HIV1_ALL_2017_gag_PRO.fasta')
lanl$subtype<-sub('\\..*','',lanl$name)
lanl2<-dnar::read.fa('lanl/HIV2_ALL_2017_gag_PRO.fasta')
lanl3<-dnar::read.fa('lanl/SIV_ALL_2017_gag_PRO.fasta')
regexpr('PIVQN',lanl$seq)
regexpr('PIVQN',lanl3$seq)
regexpr('PVQQ',lanl3$seq)
regexpr('PVQQ',lanl2$seq)
lanl$cap<-substring(lanl$seq,171)
lanl2$cap<-substring(lanl2$seq,136)
noGap2Gap(lanl2$cap[1],50)
lanl3$cap<-substring(lanl3$seq,163)
noGap2Gap(lanl3$cap[1],50)
noGap2Gap(lanl$cap[grep('HXB2',lanl$name)],50)
table(substring(lanl$cap,51,51))
noGap2Gap(lanl$cap[grep('HXB2',lanl$name)],120)
t(t(sort(table(substring(lanl$cap,137,137)),decreasing=TRUE)))

t(t(sort(table(substring(lanl2$cap,50,50)),decreasing=TRUE)))

substring(lanl3$cap,50-2,50+2)
