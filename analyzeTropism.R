library(dnar)
source('iuStan.R')
source('functions.R')
 
if(!exists('dat'))source('readNewData.R')

counts<-readCounts('ice/out/coreceptor_iceAssay_18-03-15/counts.csv')
receptor<-counts[grepl('coreceptor',counts$plate),]
counts<-counts[grepl('^iceAssay_[0-9]+h$',counts$plate),]
counts$time<-sapply(strsplit(basename(counts$dir),'_'),tail,1)
counts$time<-sprintf('%s%s',ifelse(grepl('^[0-9]h$',counts$time),'0',''),counts$time)
plateIds<-readPlateViruses('ice/ice_2018-03-09.csv','ice/ice_2018-03-09_virus.csv')
counts$virus<-plateIds[counts$well,'virus']
if(any(is.na(counts$virus)))stop('Unidentified well')

round(tapply(counts$n,list(counts$virus,counts$time),mean),1)

counts<-readCounts('ice/out/03.06.2018_recount/counts.csv')
counts$time<-basename(counts$dir)
control<-counts[counts$time=='control',]
counts<-counts[grepl('h$',counts$time),]
counts$time<-sprintf('%s%s',ifelse(grepl('^[0-9]h$',counts$time),'0',''),counts$time)
plateIds<-readPlateViruses('ice/iceAssay1.csv','ice/iceAssay1_viruses.csv')
plateIds$virus<-sub('WildTipe','WildType',plateIds$virus)
plateIds$virus<-sub('(MM[0-9]+)\\.([0-9])\\.','\\1.0\\2.',plateIds$virus)
counts$virus<-plateIds[counts$well,'virus']
iceMeans<-tapply(counts$n,list(counts$virus,counts$time),mean)
iceSd<-tapply(counts$n,list(counts$virus,counts$time),sd)
counts$h<-as.numeric(sub('h','',counts$time))
iceProps<-t(apply(iceMeans,1,function(xx)xx/xx[1]))
calculateHalf<-function(counts,times){
  #adding 1 pseudocount
  tmp<-data.frame('x'=times,'y'=log(counts+1))
  mod<-lm(y~x,tmp)
  halfLife<-log(.5)/coefficients(mod)[2]  
  return(halfLife)
}
iceHalfs<-sapply(unique(counts$virus),function(xx)calculateHalf(counts[counts$virus==xx,'n'],counts[counts$virus==xx,'h']))
iceHalfs2<-sapply(unique(counts$virus),function(xx)calculateHalf(counts[counts$virus==xx&counts$h<25,'n'],counts[counts$virus==xx&counts$h<25,'h']))
names(iceHalfs)<-names(iceHalfs2)<-unique(counts$virus)

icePats<-unique(na.omit(dat[rownames(iceMeans),'pat']))
times<-as.numeric(sub('h','',colnames(iceMeans)))
pdf('out/iceMeans.pdf')
  plot(1,1,type='n',xlab='Hours on ice',ylab='Proportion of 0h TZM-BL infectivity',xlim=range(times),ylim=range(iceProps),log='y')
  for(ii in 1:nrow(iceProps))lines(times,iceProps[ii,],col=1+grepl('JRFL',rownames(iceProps)[ii]))
  par(mar=c(3.5,4,.1,4.5))
  plot(dat[rownames(iceMeans),'ic50'],iceProps[,'08h'],ylab='Proportion TZM-BL infectivity after 8h on ice',xlab='IFNa2 IC50',log='x',xaxt='n',pch=21,bg=patCols[dat[rownames(iceMeans),'pat']],cex=1.4,las=1,mgp=c(2.5,.7,0))
  logAxis(1,las=1)
  axis(4,iceProps[!rownames(iceMeans) %in% rownames(dat),'08h'],sub('CoJRFL_','',rownames(iceProps)[!rownames(iceMeans) %in% rownames(dat)]),las=1)
  par(mar=c(3.5,4,.1,4.5))
  legend('topleft',icePats,pt.bg=patCols[icePats],pch=21,inset=.01)
  plot(dat[rownames(iceMeans),'time'],iceProps[,'08h'],ylab='Proportion TZM-BL infectivity after 8h on ice',xlab='Days after onset of symptoms',pch=21,bg=patCols[dat[rownames(iceMeans),'pat']],las=1,cex=1.4,mgp=c(2.5,.7,0))
  axis(4,iceProps[!rownames(iceMeans) %in% rownames(dat),'08h'],sub('CoJRFL_','',rownames(iceProps)[!rownames(iceMeans) %in% rownames(dat)]),las=1)
  legend('top',icePats,pt.bg=patCols[icePats],pch=21,inset=.01)
  plot(dat[rownames(iceMeans),'ic50'],iceProps[,'48h'],ylab='Proportion TZM-BL infectivity after 48h on ice',xlab='IFNa2 IC50',log='x',xaxt='n',pch=21,bg=patCols[dat[rownames(iceMeans),'pat']],cex=1.4,las=1,mgp=c(2.5,.7,0))
  logAxis(1,las=1)
  axis(4,iceProps[!rownames(iceMeans) %in% rownames(dat),'48h'],sub('CoJRFL_','',rownames(iceProps)[!rownames(iceMeans) %in% rownames(dat)]),las=1)
  par(mar=c(3.5,4,.1,4.5))
  legend('topleft',icePats,pt.bg=patCols[icePats],pch=21,inset=.01)
  plot(dat[rownames(iceMeans),'time'],iceProps[,'48h'],ylab='Proportion TZM-BL infectivity after 48h on ice',xlab='Days after onset of symptoms',pch=21,bg=patCols[dat[rownames(iceMeans),'pat']],las=1,cex=1.4,mgp=c(2.5,.7,0))
  axis(4,iceProps[!rownames(iceMeans) %in% rownames(dat),'48h'],sub('CoJRFL_','',rownames(iceProps)[!rownames(iceMeans) %in% rownames(dat)]),las=1)
  legend('top',icePats,pt.bg=patCols[icePats],pch=21,inset=.01)
  #half life
  plot(dat[names(iceHalfs),'time'],iceHalfs,ylab='Half life on ice (h)',xlab='Days after onset of symptoms',pch=21,bg=patCols[dat[names(iceHalfs),'pat']],las=1,cex=1.4,mgp=c(2.5,.7,0))
  axis(4,iceHalfs[!names(iceHalfs) %in% rownames(dat)],sub('CoJRFL_','',names(iceHalfs)[!names(iceHalfs) %in% rownames(dat)]),las=1)
  legend('top',icePats,pt.bg=patCols[icePats],pch=21,inset=.01)
dev.off()


controlVirus<-rep(c('UK61.1-P13A4','Finzi WT'),each=4)
names(controlVirus)<-LETTERS[1:8]
control<-control[order(control$row,control$col),]
control$virus<-controlVirus[control$row]
virusCol<-c('red','blue')
names(virusCol)<-unique(controlVirus)
virusInput<-c(3^(0:10),3^11) #c(50*.33^(0:10),0)
pdf('out/iceControl.pdf')
  #plot(virusInput[control$col],control$n,xlab='Virus input',ylab='TZM-BL spot count',las=1)
  #for(ii in unique(control$row))lines(virusInput[control[control$row==ii,'col']],control[control$row==ii,'n'])
  #plot(virusInput[control$col],control$area,xlab='Virus input',ylab='TZM-BL spot area',las=1)
  #for(ii in unique(control$row))lines(virusInput[control[control$row==ii,'col']],control[control$row==ii,'area'])
  plot(virusInput[control$col],control$n+1,xlab='Fold virus dilution',ylab='TZM-BL spot count',las=1,log='xy',xaxt='n',yaxt='n',bg=virusCol[control$virus],pch=21)
  logAxis(1,axisMax=3^10)
  logAxis(2,las=1,offset=1,axisMax=3^10)
  axis(1,tail(virusInput,1),0)
  for(ii in unique(control$row))lines(virusInput[control[control$row==ii,'col']],control[control$row==ii,'n']+1,col=virusCol[controlVirus[ii]])
  legend('topright',names(virusCol),lty=1,col=virusCol,pt.bg=virusCol,inset=.01)
  plot(virusInput[control$col],control$area+1,xlab='Fold virus dilution',ylab='TZM-BL spot area',las=1,log='xy',xaxt='n',yaxt='n',bg=virusCol[control$virus],pch=21)
  logAxis(1,axisMax=3^10)
  axis(1,tail(virusInput,1),0)
  logAxis(2,las=1,offset=1)
  for(ii in unique(control$row))lines(virusInput[control[control$row==ii,'col']],control[control$row==ii,'area']+1,col=virusCol[controlVirus[ii]])
  legend('topright',names(virusCol),lty=1,col=virusCol,pt.bg=virusCol,inset=.01)
dev.off()



drugViruses<-c("SG3 (X4-trop. ctrl)","YU2 (R5-trop. ctrl)","CH77 (Dual-trop. ctrl)","MM15.02.2B5 bulk","MM15.03.2B6 bulk","MM15.03.2A3 bulk","MM15.06.2A1 bulk","MM15.07.2A3 bulk","MM15.08.2B1 bulk","MM15.11.2B4 bulk","MM15.13.2C1 bulk","MM15.14.2C3 bulk")
drugs<-c('No-drug', 'AMD (X4-inhib.)', 'TAK (R5-inhib.)', 'Mara(R5-inhib.)','TAK + AMD','Mara + AMD')
names(drugs)<-LETTERS[1:length(drugs)]
receptor$virus<-drugViruses[receptor$col]
receptor$drug<-drugs[receptor$row]
receptor<-receptor[!is.na(receptor$drug),]
drugMeans<-tapply(receptor$n,list(receptor$virus,receptor$drug),mean)
drugMeans[drugViruses,drugs]

iceTrop<-readCounts('ice/out/2018-04-02-iceTropism2/counts.csv')
iceTropMats<-tapply(iceTrop$n,list(iceTrop$row,iceTrop$col,iceTrop$plate),c)

iceLuc<-read.csv('ice/2018-04-02-iceTropism/2018-04-02_ice.csv',check.names=FALSE,row.names=1)[,-13]
#ran out of virus
iceLuc[c('B','D','F','H'),c(3:4,7:8)]<-NA
tropLuc<-read.csv('ice/2018-04-02-iceTropism/2018-04-02_tropism.csv',check.names=FALSE,row.names=1)[,-13]

viruses<-rep(c("WT","L193A"),4)
times<-rep(c(0,8,24,48),each=2)
virCol<-c('WT'='#0000FF','L193A'='#FF0000')
thawViruses<-rep(rep(c("WT","L193A"),4),2)
freezeThaws<-rep(1:4,each=2)
#abusing globals
plotIce<-function(xx,main='',ylab='TZM-BL B-Gal spots',normalize=FALSE){
  ylim<-c(0,max(xx,na.rm=TRUE))
  if(normalize){
    zeros<-apply(xx[1:2,],1,mean,na.rm=TRUE)
    xx<-xx/zeros
    ylim<-range(xx,na.rm=TRUE)
  }
  plot(rep(times,ncol(xx)),as.vector(unlist(xx)),main=main,xlab='Time on ice (h)',ylab=ylab,las=1,bg=virCol[viruses],pch=21,cex=1.5,ylim=ylim,log=ifelse(normalize,'y',''),yaxt=ifelse(normalize,'n','s'))
  if(normalize) logAxis(2,las=1,exponent=FALSE)
  means<-apply(xx,1,median,na.rm=TRUE)
  for(ii in 1:2){
    if(normalize){
      tmp<-data.frame('x'=rep(times[seq(ii,nrow(xx),2)],ncol(xx)),'y'=log(as.vector(unlist(xx[seq(ii,nrow(xx),2),]))))
      mod<-lm(y~x,tmp)
      fakeTime<-seq(0,60,.01) 
      preds<-predict(mod,data.frame(x=fakeTime),interval='confidence')
      lines(fakeTime,exp(preds[,1]),col=virCol[ii])
      polygon(c(fakeTime,rev(fakeTime)),exp(c(preds[,2],rev(preds[,3]))),border=NA,col=sprintf('%s11',virCol[ii]))
      halfLife<-log(.5)/coefficients(mod)[2]
      text(par('usr')[1]+.01*diff(par('usr')[1:2]),10^(par('usr')[3]+(.01+(ii-2)*-0.06)*diff(par('usr')[3:4])),sprintf('%s half-life=%0.1fhr',names(virCol)[ii],halfLife),adj=c(0,0))
    }else{
      lines(times[seq(ii,nrow(xx),2)],means[seq(ii,nrow(xx),2)],col=virCol[ii])
    }
  }
}
plotThaw<-function(xx,ylab='TZM-BL B-Gal spots'){
  plot(rep(freezeThaws,ncol(xx)),as.vector(unlist(xx)),xlab='Number of freeze thaws',ylab=ylab,las=1,bg=virCol[thawViruses],pch=21,cex=1.5,xaxt='n',log='y')
  axis(1,1:4)
  for(ii in 1:2){
    thisFreezes<-rep(freezeThaws[seq(ii,nrow(xx),2)],ncol(xx))
    vals<-as.vector(unlist(xx[seq(ii,nrow(xx),2),]))
    mod<-lm(log(vals)~thisFreezes) 
    fakeThaws<-seq(0,6,.01)
    preds<-predict(mod,data.frame(thisFreezes=fakeThaws),interval='confidence')
    lines(fakeThaws,exp(preds[,1]),col=virCol[ii])
    polygon(c(fakeThaws,rev(fakeThaws)),exp(c(preds[,2],rev(preds[,3]))),border=NA,col=sprintf('%s11',virCol[ii]))
  }
}
pdf('out/iceWtMut.pdf',width=10,height=10)
  par(mfrow=c(3,2))
  for(normalize in c(FALSE,TRUE)){
    if(normalize)ylabs<-sprintf('Proportion of 0-hour infectivity %s',c('(B-Gal)','(Luciferase)'))
    else ylabs<-c('TZM-BL B-Gal spots','TZM-BL luciferase luminescence')
    withAs(xx=iceTropMats[,1:4,'ice'],plotIce(xx,main='50ul in 96 well plate',ylab=ylabs[1],normalize=normalize))
    withAs(xx=iceLuc[,1:4],plotIce(xx,main='50ul in 96 well plate',ylab=ylabs[2],normalize=normalize))
    legend('topright',names(virCol),pt.bg=virCol,pch=21,inset=.01)
    withAs(xx=iceTropMats[,5:8,'ice'],plotIce(xx,main='50ul in a 1.5ml tube',ylab=ylabs[1],normalize=normalize))
    withAs(xx=iceLuc[,5:8],plotIce(xx,main='50ul in a 1.5ml tube',ylab=ylabs[2],normalize=normalize))
    withAs(xx=iceTropMats[,9:12,'ice'],plotIce(xx,main='500ul in a 1.5ml tube (no freeze-thaw)',ylab=ylabs[1],normalize=normalize))
    withAs(xx=iceLuc[,9:12],plotIce(xx,main='500ul in a 1.5ml tube (no freeze-thaw)',ylab=ylabs[2],normalize=normalize))
  }
  par(mfrow=c(2,2))
  withAs(xx=iceTropMats[,11:12,'tropism'],plotThaw(xx,ylab='TZM-BL B-Gal spots'))
  withAs(xx=tropLuc[,11:12],plotThaw(xx,ylab='TZM-BL luciferase luminescence'))
  legend('topright',names(virCol),pt.bg=virCol,pch=21,inset=.01)
dev.off()

calcIc50<-function(n,time,pseudocount=0){
  if(length(n)==0)return(NA)
  n<-log10(n+pseudocount)
  mod<-lm(n~time)
  log10(.5)/coefficients(mod)[2]
}

iceSpin<-readCounts('ice/out/2018-04-27_iceSpin/counts.csv')
iceSpinMats<-tapply(iceSpin$n,list(iceSpin$row,iceSpin$col,iceSpin$plate),c)
iceTimes<-matrix(rep(rep(c(48,24,8,0),each=4),6),nrow=8,dimnames=list(LETTERS[1:8],1:12))
iceVirusesCulture<-structure(rep(c('SG3+WT','SG3+Mut','NL43+WT','NL43+Mut','NL43+Mut 4x','236'),each=2),.Names=1:12)
iceVirusesSpin<-structure(rep(c('SG3+WT','SG3+Mut','NL43+WT','NL43+Mut','NL43+WT 4x','NL43+Mut 4x'),each=2),.Names=1:12)
iceSpin$time<-sapply(iceSpin$well,function(xx)iceTimes[sub('[0-9]+','',xx),sub('[A-Z]','',xx)])
iceSpin$virus<-ifelse(iceSpin$plate=='ice_spin',iceVirusesSpin[iceSpin$col],iceVirusesCulture[iceSpin$col])
zeroMeans<-withAs(xx=iceSpin[iceSpin$time==0,],tapply(xx$n,paste(xx$virus,xx$plate),mean))
iceSpin$prop<-iceSpin$n/zeroMeans[paste(iceSpin$virus,iceSpin$plate)]
iceSpinIc50<-outer(unique(iceSpin$plate),unique(iceSpin$virus),function(xx,yy)mapply(function(xxx,yyy){
    withAs(zzz=iceSpin[iceSpin$virus==yyy&iceSpin$plate==xxx,],calcIc50(zzz$n,zzz$time))
},xx,yy))
dimnames(iceSpinIc50)<-list(unique(iceSpin$plate),unique(iceSpin$virus))
round(iceSpinIc50[,order(colnames(iceSpinIc50))],1)


pdf('out/iceSpin.pdf')
for(ii in unique(iceSpin$virus)){
  for(jj in unique(iceSpin$plate)){
    thisDat<-iceSpin[iceSpin$virus==ii&iceSpin$plate==jj,]
    if(nrow(thisDat)==0)next()
    plot(1,1,type='n',xlim=range(iceSpin$time),ylim=range(iceSpin$prop),log='y',main=paste(ii,jj),yaxt='n',ylab='TZM-Bl infectivity (proportion of 0 hour)',xlab='Time on ice')
    logAxis(las=1)
    points(thisDat$time,thisDat$prop,pch=21,bg='#0000FFBB',cex=1.5)
    tmp<-data.frame('x'=thisDat$time,'y'=log(thisDat$prop))
    mod<-lm(y~x,tmp)
    fakeTime<-seq(0,60,.01) 
    preds<-predict(mod,data.frame(x=fakeTime),interval='confidence')
    lines(fakeTime,exp(preds[,1]),col='black')
    polygon(c(fakeTime,rev(fakeTime)),exp(c(preds[,2],rev(preds[,3]))),border=NA,col='#00000033')
  }
}
dev.off()


magenta<-readCounts('ice/out/2018-05-04_mmImcIceTropism/counts.csv')
magentaViruses<-c('YU2','896','NL43','SG3')
magentaMats<-tapply(magenta$n,list(magenta$row,magenta$col,magenta$plate),c)
magCompare<-cbind('24hr'=magentaMats[,5:8,'24hr'],'48hr'=magentaMats[,8:11,'magentaGal'],'magenta_48hr'=magentaMats[,2:5,'magentaGal'])
colnames(magCompare)<-paste(rep(magentaViruses,3),rep(c('24hr','48hr','magenta'),each=4))
rownames(magCompare)<-sprintf('%dx',3^(1:8))
magCompare[,order(colnames(magCompare))]

iceMM<-readCounts('ice/out/2018-05-07_iceTropism/counts.csv')
iceMMMats<-tapply(iceMM$n,list(iceMM$row,iceMM$col,iceMM$plate),c)
iceMM<-iceMM[iceMM$plate!='tropism',]
iceMMVirus<-list('ice1'=rep(c('896','NL43','SG3','YU2','470 TF','470 6mo'),each=2),'ice2'=rep(c('40 TF','40 6mo','MM23.1','MM23.7','MM23.11','MM23.13'),each=2),'ice3'=rep(c('NL43+WT','NL43+Mut','SG3+WT','SG3+Mut','SG3+WEAU','No virus'),each=2))
iceMMTime<-matrix(rep(c(50,26,10,0),each=4),ncol=12,nrow=8,dimnames=list(LETTERS[1:8],1:12))
iceMM$time<-indexMatrix(iceMM$row,iceMM$col,iceMMTime)
iceMM$virus<-mapply(function(xx,yy)xx[yy],iceMMVirus[iceMM$plate],iceMM$col)
halfs<-sapply(unique(iceMM$virus),function(xx){withAs(zzz=iceMM[iceMM$virus==xx,],calcIc50(zzz$n+1,zzz$time))})
names(halfs)<-unique(iceMM$virus)
t(t(round(halfs[orderIn(names(halfs),unique(unlist(iceMMVirus)))],1)))


zeroMeans<-withAs(xx=iceMM[iceMM$time==0,],tapply(xx$n+1,xx$virus,mean))
iceMM$prop<-(iceMM$n+1)/zeroMeans[iceMM$virus]
pdf('out/iceMM.pdf')
for(ii in unique(unlist(iceMMVirus))){
  thisDat<-iceMM[iceMM$virus==ii,]
  if(nrow(thisDat)==0)next()
  plot(1,1,type='n',xlim=range(iceMM$time),ylim=range(iceMM$prop),log='y',main=sprintf('%s\nhalf=%0.1f',ii,halfs[ii]),yaxt='n',ylab='TZM-Bl infectivity (proportion of 0 hour)',xlab='Time on ice')
  logAxis(las=1)
  points(thisDat$time,thisDat$prop,pch=21,bg='#0000FFBB',cex=1.5)
  tmp<-data.frame('x'=thisDat$time,'y'=log(thisDat$prop))
  mod<-lm(y~x,tmp)
  fakeTime<-seq(0,60,.01) 
  preds<-predict(mod,data.frame(x=fakeTime),interval='confidence')
  lines(fakeTime,exp(preds[,1]),col='black')
  polygon(c(fakeTime,rev(fakeTime)),exp(c(preds[,2],rev(preds[,3]))),border=NA,col='#00000033')
  abline(h=1,lty=2)
}
dev.off()

library(abind)
virus<-list(
  'A'=c('MM23.01.2A3','UK61.1-P2C2','MM23.07.2B2','UK61.7-P2A4','MM23.11.1D6','UK61.11-P2B3','UK61.13-P1B4','UK61.13-P21A3','MM15.02.2B5 bulk','MM15.02.2A2 bulk','MM15.03.2B6 bulk','MM15.08.2B2 bulk'),
  'E'=c('MM15.11.2B4 bulk','MM15.14.2C4 bulk','WEAU.2.1A2','WEAU.3.1A6','WEAU.11.1B6','WEAU.15.1C6','WEAU.19.1A1','89.6','SG3','YU2','SG3+WEAU','No virus')
)
treats<-factor(structure(rep(c('No drug','AMD','Tak','AMD+Tak'),2),.Names=LETTERS[1:8]))
trop<-readCounts('ice/out/2018-06-09_tropismTitration/counts.csv')
trop<-trop[trop$plate!='titration',]
trop$conc<-sub('^([0-9]x)$','0\\1',trop$plate)
trop$virus<-mapply(function(xx,yy)xx[yy],virus[sub('[FGH]','E',sub('[BCD]','A',trop$row))],trop$col)
trop$treat<-treats[trop$row]
tropCounts<-tapply(trop$n,list(trop$virus,trop$treat,trop$conc),c)
props<-do.call(abind,c(lapply(dimnames(tropCounts)[[3]],function(xx){out<-tropCounts[,,xx]/ifelse(tropCounts[,'No drug',xx]>500,tropCounts[,'No drug',xx],NA)}),list(along=3)))[unlist(virus),treats[1:4],]
propMean<-apply(props,c(1,2),mean,na.rm=TRUE)
pdf('tropism.pdf',height=12)
  par(mar=c(3,10,.1,.1))
  cols=rev(heat.colors(1000))
  breaks=seq(0,max(propMean,na.rm=TRUE),length.out=1001)
  image(1:ncol(propMean),1:nrow(propMean),t(propMean),xaxt='n',yaxt='n',xlab='',ylab='',col=cols,breaks=breaks)
  axis(2,1:nrow(propMean),rownames(propMean),las=1)
  axis(1,1:ncol(propMean),colnames(propMean))
  insetScale(breaks,cols,insetPos = c(0.015, 0.015, 0.025, 0.25),main='Proportion of no drug')
  box()
  abline(h=1:nrow(propMean)-.5,col='#00000033')
  abline(v=1:ncol(propMean)-.5,col='#00000033')
dev.off()
#tak = ccr5 inhibitor, amd = cxcr4 inhibitor


tit<-readCounts('ice/out/2018-06-23_titration/counts.csv')
virus<-read.csv('ice/2018-06-23_titration.csv',stringsAsFactors=FALSE)
virus$final<-virus$Dilution*3
virus$well<-paste(rep(LETTERS[1:4],each=12),rep(1:12,4),sep='')
rownames(virus)<-virus$well
pdf('out/neutralRed.pdf');withAs(xx=tit[order(tit$well),],plot(xx[xx$plate=='1x','n'],xx[xx$plate=='1x_NR2','n'],las=1,xlab='Count without NR',ylab='Count with NR'));dev.off()
tit<-tit[grep('^[0-9]+x$',tit$plate),]
tit$dilution<-as.numeric(sub('x','',tit$plate))*ifelse(tit$row %in% c('E','F','G','H'),3,1)
tit<-tit[order(tit$col,tit$dilution),]
tit$virusNum<-ifelse(tit$row %in% LETTERS[5:8],paste(structure(LETTERS[1:4],.Names=LETTERS[5:8])[tit$row],tit$col,sep=''),tit$well)
tit$baseDilute<-virus[tit$virusNum,'final']
tit$virus<-virus[tit$virusNum,'Virus']
tit$totalDilute<-tit$baseDilute*tit$dilution
tit$filter<-ave(tit$n,tit$virus,FUN=function(xx)ifelse(1:length(xx)>=which.max(xx)&xx>25,xx,NA))
virus$iuPerUl<-NA
pdf('out/titration.pdf')
plot(tit$totalDilute,tit$n+1,xlab='Dilution',ylab='TZMBL count',las=1,type='n',log='yx')
for(ii in unique(tit$virus))withAs(zz=tit[tit$virus==ii,],lines(zz$totalDilute,zz$n))
abline(100,1,col='#FF0000')
for(ii in virus$Virus){
  thisDat<-tit[tit$virus==ii,]
  thisDat$logDil<--log(thisDat$totalDilute)
  withAs(zz=thisDat,plot(zz$totalDilute,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=ii,ylim=range(tit$n+1),xlim=range(tit$totalDilute),pch=21,bg=is.na(thisDat$filter)+1,cex=2))
  if(sum(!is.na(thisDat$filter))>0){
    mod<-glm(n~offset(logDil),thisDat[!is.na(thisDat$filter),],family='poisson')
    fakeDat<-data.frame(logDil=-log(1:50000))
    preds<-predict(mod,fakeDat)
    lines(exp(-fakeDat$logDil),exp(preds)+1)
    iuPerUl<-exp(coef(mod)['(Intercept)'])/150
    virus$iuPerUl[virus$Virus==ii]<-iuPerUl
    mtext(sprintf('%0.1f IU/ul (polybrene+spin)',iuPerUl),3)
  }
}
dev.off()
write.csv(virus,'ice/2018-06-23_titration_calc.csv')

viruses<-read.csv('ice/2018-06-17_bigTropism.csv',stringsAsFactors=FALSE)
viruses$virusNum<-paste(viruses$Plate,substring(viruses$Row,1,1),viruses$Col,sep='_')
rownames(viruses)<-viruses$virusNum
uniqDrugs<-c('No drug','AMD','Tak','AMD+Tak')
drugs<-structure(rep(uniqDrugs,2),.Names=LETTERS[1:8])

viruses3<-read.csv('ice/2018-07-03_tropism.csv',stringsAsFactors=FALSE)
viruses3$virusNum<-paste(viruses3$Plate,viruses3$Col,sep='_')
rownames(viruses3)<-viruses3$virusNum
moreDrugs<-structure(c(uniqDrugs,'Mar','HiMar','Mar+AMD','HiMar+AMD'),.Names=LETTERS[1:8])
trop3<-readCounts('ice/out/2018-07-07_tropism/counts.csv')
nr3<-trop3[grepl("NR",trop3$plate),]
trop3<-trop3[!grepl("NR",trop3$plate),]
#fix a flipped plate
trop3[trop3$plate=='P1_27x','col']<-rev(1:12)[trop3[trop3$plate=='P1_27x','col']]
trop3[trop3$plate=='P1_27x','row']<-structure(LETTERS[8:1],.Names=LETTERS[1:8])[trop3[trop3$plate=='P1_27x','row']]
trop3$drug<-moreDrugs[trop3$row]
trop3$virusNum<-paste(sub('^P','',sub('_.*','',trop3$plate)),trop3$col,sep='_')
trop3$virus<-viruses3[trop3$virusNum,'Virus']
trop3$date<-'7-3'
trop3$dilution<-viruses3[trop3$virusNum,'Dilution']*as.numeric(sub('x$','',sub('P[12]_','',trop3$plate)))

viruses2<-read.csv('ice/2018-06-21_tropismRedo.csv',stringsAsFactors=FALSE)
trop2<-readCounts('ice/out/2018-06-21-drugTitration/counts.csv')
trop2<-trop2[grep('redo',trop2$plate),]
trop2$virus<-viruses2[trop2$col,'Virus']
trop2$dilution<-ifelse(trop2$plate=='redo_1',9,1)*ifelse(trop2$row %in% LETTERS[5:8],3,1)*viruses2[trop2$col,'Dilution']
trop2$drug<-drugs[trop2$row]
trop2$date<-'6-21'

trop<-readCounts('ice/out/2018-06-17_tropism/counts.csv')
trop$rowGroup<-ifelse(trop$row %in% LETTERS[5:8],'E','A')
trop<-trop[trop$plate!='BGAL',]
trop$pId<-sub('^P([0-9]).*','\\1',trop$plate)
trop$virusNum<-paste(trop$pId,trop$rowGroup,trop$col,sep='_')
trop$dilution<-as.numeric(sub('^P[0-9]_([0-9]+)x$','\\1',trop$plate))*viruses[trop$virusNum,'Dilution']
trop$virus<-viruses[trop$virusNum,'Virus']
trop$drug<-drugs[trop$row]
trop$date<-'6-17'

desiredCols<-c('virus','dilution','n','drug','date')
combo<-rbind(trop[,desiredCols],trop2[,desiredCols],trop3[trop3$drug %in% uniqDrugs,desiredCols])
combo<-combo[order(combo$virus,combo$dilution),]
combo$dateF<-as.factor(combo$date)
drugColors<-structure(rainbow.lab(length(uniqDrugs)),.Names=uniqDrugs)
viruses$bestDilution<-sapply(viruses$Virus,function(xx){
  thisDat<-combo[combo$drug=='No drug'&combo$virus==xx,]
  tail(thisDat$dilution[thisDat$n==max(thisDat$n)],1)
})
datePch<-c('6-17'=21,'6-21'=22,'7-3'=23)
#best<-combo[combo$dilution==viruses[combo$virusNum,'bestDilution'],]
mods<-list()
pdf('out/bigTropism.pdf')
for(ii in unique(c(viruses$Virus,viruses3$Virus))){
  thisDat<-combo[combo$virus==ii,]
  thisDrugs<-uniqDrugs
  #if(all(moreDrugs %in% thisDat$drug))thisDrugs<-unname(moreDrugs)
  #else thisDrugs<-uniqDrugs
  firstDilute<-do.call(rbind,lapply(unique(combo$date),function(day){sapply(thisDrugs,function(dd){
    thisTime<-thisDat[thisDat$date==day,]
    if(all(thisTime[thisTime$drug==dd,'n']<50))return(-Inf)
    else return(tail(thisTime[(thisTime$n>.8*max(thisTime[thisTime$drug==dd,'n'])) & thisTime$drug==dd,'dilution'],1))
  })}))
  rownames(firstDilute)<-unique(combo$date)
  thisDat$filter<-thisDat$n
  thisDat$drug<-factor(thisDat$drug,levels=thisDrugs)
  for(dd in thisDrugs)thisDat[thisDat$drug==dd&thisDat$dilution<firstDilute[,dd][thisDat$date],'filter']<-NA
  thisDat$logDil<--log(thisDat$dilution)
  withAs(zz=thisDat,plot(zz$dilution,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=ii,ylim=range(combo$n+1)*c(1,2),xlim=range(combo$dilution),type='n'))
  #abline(v=viruses[viruses$Virus==ii,'bestDilution'],lty=2)
  for(dd in thisDrugs){
    withAs(zz=thisDat[thisDat$drug==dd,],points(zz$dilution,zz$n+1,pch=datePch[zz$date],cex=2,bg=drugColors[dd],col=c('black','#00000000')[is.na(zz$filter)+1],lwd=3))
  }
  legend('topright',names(drugColors),col=drugColors,lty=1,inset=.01)
  if(length(unique(thisDat$date))>1)mod<-glm(filter~offset(logDil)+drug+dateF,thisDat,family='poisson')
  else mod<-glm(filter~offset(logDil)+drug,thisDat,family='poisson')
  for(dd in thisDrugs){
    fakeDat<-data.frame(logDil=-log(2^seq(0,16,.1)),'drug'=dd,'dateF'='6-17')
    preds<-predict(mod,fakeDat)
    lines(exp(-fakeDat$logDil),exp(preds)+1,col=drugColors[dd])
  }
  amd<-1-exp(coef(mod)['drugAMD'])
  tak<-1-exp(coef(mod)['drugTak'])
  amdTak<-1-exp(coef(mod)['drugAMD+Tak'])
  mtext(sprintf('Inhibition:\nAMD: %0.2f%%\nTAK: %0.2f%%\nAMD+TAK: %0.2f%%',amd*100,tak*100,amdTak*100),3,at=exp(log(max(combo$dilution))*.8))
  mods[[ii]]<-mod
}
dev.off()

estimates<-do.call(rbind,lapply(mods,function(mod){
  est<-coef(summary(mod))[sprintf('drug%s',uniqDrugs[-1]),'Estimate']
  se<-coef(summary(mod))[sprintf('drug%s',uniqDrugs[-1]),'Std. Error']
  c('est'=exp(est),'lower'=exp(est-se*2),'upper'=exp(est+se*2))
}))
estimates<-estimates[rownames(estimates)!='No virus',]
mms<-sub('(MM[0-9]+|WEAU\\.)?.*','\\1',rownames(estimates))
times<-as.numeric(ifelse(mms=='',NA,sub('(MM[0-9]+|WEAU)\\.([0-9]+).*$','\\2',rownames(estimates))))
estimates<-estimates[order(mms!='',mms,times),]
mms<-sub('(MM[0-9]+|WEAU\\.)?.*','\\1',rownames(estimates))
mmSwitch<-c(FALSE,mms[-1]!=mms[-length(mms)])
virusPos<-cumsum(rep(1,nrow(estimates))+mmSwitch)
rectPos<-seq(-.4,.4,length.out=ncol(estimates)/3+2)
rectCols<-structure(rainbow.lab(4),.Names=uniqDrugs)
pdf('out/bigTropismSummary.pdf',width=20)
  for(isLog in c(FALSE,TRUE)){
  par(mar=c(7,4,.2,4))
  #plot(1,1,type='n',ylim=range(estimates),xlab='',las=1,ylab='TZMBL BGal spots (proportion relative to no drug)',xlim=c(1,nrow(estimates))+c(-.5,.5),yaxt='n',xaxt='n',xaxs='i')
  plot(1,1,type='n',ylim=range(estimates[rownames(estimates)!='No virus',1:3]),xlab='',las=1,ylab='TZMBL BGal spots (proportion relative to no drug)',xlim=range(virusPos)+c(-1,1),xaxt='n',xaxs='i',log=ifelse(isLog,'y',''),yaxt=ifelse(isLog,'n','s'))
  slantAxis(1,virusPos,rownames(estimates),las=2)
  if(!isLog){
    rect(rectPos[1]+virusPos,1,rectPos[2]+virusPos,0,col=rectCols['No drug'])
  } else {
    abline(h=1,lty=2)
    logAxis(las=1)
  }
  for(ii in 1:3)rect(rectPos[1+ii]+virusPos,estimates[,ii],rectPos[2+ii]+virusPos,ifelse(isLog,1,0),col=rectCols[sub('est.drug','',colnames(estimates)[ii])])
  for(ii in 1:3)segments(mean(rectPos[1:2+ii])+virusPos,estimates[,ii+3],mean(rectPos[1:2+ii])+virusPos,estimates[,ii+6])
  legend(ifelse(isLog,'bottomright','topleft'),names(rectCols),fill=rectCols,bty='n')
  }
dev.off()

drugColors<-structure(rainbow.lab(length(moreDrugs)),.Names=moreDrugs)
mods<-list()
allDrugs<-list()
pdf('out/tropismMaraviroc.pdf')
for(ii in c(viruses3$Virus[viruses3$Virus %in% c(viruses2$Virus,viruses$Virus)],'89.6 (293T)')){
  thisDat<-trop3[trop3$virus==ii,]
  thisDrugs<-unname(moreDrugs)
  firstDilute<-do.call(rbind,lapply(unique(combo$date),function(day){sapply(thisDrugs,function(dd){
    thisTime<-thisDat[thisDat$date==day,]
    if(all(thisTime[thisTime$drug==dd,'n']<50))return(-Inf)
    else return(tail(thisTime[(thisTime$n>.8*max(thisTime[thisTime$drug==dd,'n'])) & thisTime$drug==dd,'dilution'],1))
  })}))
  rownames(firstDilute)<-unique(combo$date)
  thisDat$filter<-thisDat$n
  thisDat$drug<-factor(thisDat$drug,levels=thisDrugs)
  for(dd in thisDrugs)thisDat[thisDat$drug==dd&thisDat$dilution<firstDilute[,dd][thisDat$date],'filter']<-NA
  thisDat$logDil<--log(thisDat$dilution)
  withAs(zz=thisDat,plot(zz$dilution,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=ii,ylim=range(combo$n+1)*c(1,2),xlim=range(combo$dilution),type='n'))
  for(dd in thisDrugs){
    withAs(zz=thisDat[thisDat$drug==dd,],points(zz$dilution,zz$n+1,pch=21,cex=2,bg=drugColors[dd],col=c('black','#00000000')[is.na(zz$filter)+1],lwd=3))
  }
  legend('topright',names(drugColors),col=drugColors,lty=1,inset=.01)
  if(length(unique(thisDat$date))>1)mod<-glm(filter~offset(logDil)+drug+dateF,thisDat,family='poisson')
  else mod<-glm(filter~offset(logDil)+drug,thisDat,family='poisson')
  for(dd in thisDrugs){
    fakeDat<-data.frame(logDil=-log(2^seq(0,16,.1)),'drug'=dd,'dateF'='6-17')
    preds<-predict(mod,fakeDat)
    lines(exp(-fakeDat$logDil),exp(preds)+1,col=drugColors[dd])
  }
  drugs<-1-exp(coef(mod))[-1]
  drugText<-sprintf('%s: %0.2f',sub('^drug','',names(drugs)),drugs*100)
  mtext(sprintf('Inhibition\n%s',paste(drugText[1:3],collapse='\n')),3,at=exp(log(max(combo$dilution))*.2))
  mtext(sprintf('Inhibition\n%s',paste(drugText[4:7],collapse='\n')),3,at=exp(log(max(combo$dilution))*.8))
  mods[[ii]]<-mod
  allDrugs[[ii]]<-drugs
}
dev.off()

out<-round(do.call(rbind,allDrugs)*100,2)
colnames(out)<-sub('^drug','',colnames(out))

nrCompare<-cbind(trop3[order(trop3$plate,trop3$col,trop3$row),c('plate','col','row','n')],nr3[order(nr3$plate,nr3$col,nr3$row),c('plate','col','row','n')])
colnames(nrCompare)<-paste(rep(c('no','NR'),each=4),c('plate','col','row','n'),sep='_')
if(any(nrCompare$no_plate!=sub('_NR$','',nrCompare$NR_plate)))stop('Mismatched plates')
pdf('out/neutralRed2.pdf')
  plot(nrCompare$no_n+1,nrCompare$NR_n+1,las=1,xlab='Count without NR',ylab='Count with NR',log='xy')
  abline(0,1,lty=2)
dev.off()


tit<-readCounts('ice/out/2018-07-21_titration/counts.csv')
virus<-read.csv('ice/2018-07-21_titration.csv',stringsAsFactors=FALSE)
virus$final<-1
virus$well<-paste(rep(LETTERS[1:8],4),rep(1:4,each=8),sep='')
rownames(virus)<-virus$well
wellLookup<-structure(rep(paste(rep(LETTERS[1:8],4),rep(1:4,each=8),sep=''),3),.Names=paste(rep(LETTERS[1:8],12),rep(1:12,each=8),sep=''))
tit$plateBase<-as.numeric(sapply(strsplit(tit$plate,'_'),'[[',2))
tit$totalDilute<-virus[wellLookup[tit$well],'final'] * tit$plateBase * 4^floor((tit$col-1)/4)
tit$virus<-virus[wellLookup[tit$well],'Virus']
tit$treat<-sapply(strsplit(tit$plate,'_'),'[[',1)
tit$rep<-sapply(strsplit(tit$plate,'_'),'[[',3)
nr<-tit[grepl('NR$',tit$plate),]
tit<-tit[!grepl('NR$',tit$plate),]
tit<-tit[order(tit$col,tit$totalDilute),]
tit$filter<-ave(tit$n,tit$virus,tit$treat,tit$rep,FUN=function(xx)ifelse(1:length(xx)>=which.max(xx)&xx>20,xx,NA))
tit$treat<-sapply(strsplit(tit$plate,'_'),'[[',1)
tit$treat[tit$treat=='Spin']<-'Spin+Poly'
pdf('out/2018-07-21_titration.pdf',width=10,height=8)
par(mfrow=c(1,2))
for(ii in virus$Virus){
  for(treat in unique(tit$treat)){
    thisDat<-tit[tit$virus==ii&tit$treat==treat,]
    thisDat$logDil<--log(thisDat$totalDilute)
    withAs(zz=thisDat,plot(zz$totalDilute,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s',ii,treat),ylim=range(tit$n+1),xlim=range(tit$totalDilute),pch=21,bg=is.na(thisDat$filter)+1,cex=2))
    if(sum(!is.na(thisDat$filter))>0){
      mod<-glm(n~offset(logDil),thisDat[!is.na(thisDat$filter),],family='poisson')
      fakeDat<-data.frame(logDil=-log(1:50000))
      preds<-predict(mod,fakeDat)
      lines(exp(-fakeDat$logDil),exp(preds)+1)
      iuPerUl<-exp(coef(mod)['(Intercept)'])/150
      virus[virus$Virus==ii,sprintf('iuPerUl%s',treat)]<-iuPerUl
      mtext(sprintf('%0.1f IU/ul (%s)',iuPerUl,treat),3)
    }
  }
}
dev.off()

write.csv(virus[,!colnames(virus) %in% c('well','final')],'ice/2018-07-21_titration_calc.csv')

tit<-readCounts('ice/out/2018-11-15_infectionAids/counts.csv')
tit<-tit[tit$plate!='MM39_spin_red_flipped',]
tit$virusPlate<-sapply(strsplit(tit$plate,'_'),'[[',1)
tit$moat<-grepl('moat',tit$plate)
tit$red<-grepl('red',tit$plate)
tit$spin<-grepl('spin',tit$plate)
virus<-list('MM33'=c('MM33.01.13C1','MM33.13.2A4','MM33.19.3B4','CH040 TF'),'MM39'=c('MM39.02.2B5','MM39.10.11B4','MM39.16.2C3','CH040 6mo'))
if(any(!tit$virusPlate %in% names(virus)))stop('Unknown plate')
tit$virus<-mapply(function(xx,yy)xx[ceiling(yy/3)],ifelse(tit$virusPlate=='MM33',virus['MM33'],virus['MM39']),tit$col)
tit$treat<-rep(c('poly','dextran','neat'),4)[tit$col]
tit$treatSpin<-paste(ifelse(tit$spin,'spin','noSpin'),tit$treat)
tit$totalDilute<-1/c(150/550*(110/550)^(0:6),Inf)
tit$filter<-ave(tit$n,tit$virus,tit$treatSpin,tit$moat,tit$red,FUN=function(xx)ifelse(1:length(xx)>=which.max(xx)&xx>20,xx,NA))
tit[tit$virus=='MM39.02.2B5'&tit$treatSpin=='spin poly'&!tit$red&!is.na(tit$filter)&tit$totalDilute>5000,]
pdf('out/2018-11-15-infectionAid.pdf',width=10,height=8)
par(mfrow=c(2,3))
for(ii in unique(tit$virus)){
  for(treat in unique(tit$treatSpin)){
    thisDat<-tit[tit$virus==ii&tit$treatSpin==treat&tit$totalDilute>0&!tit$red,]
    thisDat$logDil<--log(thisDat$totalDilute)
    withAs(zz=thisDat,plot(zz$totalDilute,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s',ii,treat),ylim=range(tit$n+1),xlim=range(thisDat$totalDilute),pch=21,bg=is.na(thisDat$filter)+1,cex=2))
    if(sum(!is.na(thisDat$filter))>0){
      mod<-glm(n~offset(logDil),thisDat[!is.na(thisDat$filter),],family='poisson')
      fakeDat<-data.frame(logDil=-log(1:50000))
      preds<-predict(mod,fakeDat)
      lines(exp(-fakeDat$logDil),exp(preds)+1)
      iuPerUl<-exp(coef(mod)['(Intercept)'])
      mtext(sprintf('%0.1f IU/ul (%s)',iuPerUl,treat),3)
    }
  }
}
dev.off()



tit<-readCounts('ice/out/2018-11-18_infectionAids2/counts.csv')
tit$moat<-grepl('moat',tit$plate)
tit$red<-grepl('red',tit$plate)
tit$spin<-grepl('spin',tit$plate)
tit$dil<-as.numeric(sub('^([0-9]+)x.*','\\1',tit$plate))*rep(c(1,3),each=4)[sapply(tit$row,function(xx)which(xx==LETTERS))]
virus<-c('89.6','CH40TF','CH406mo','Media')
tit$virus<-structure(rep(virus,2),.Names=LETTERS[1:8])[tit$row]
if(any(!tit$virusPlate %in% names(virus)))stop('Unknown plate')
tit$treat<-rep(c('poly','media','dextran','media'),c(5,1,5,1))[tit$col]
tit$treatDil<-c(60/2^(4:0),0,160/2^(0:4),0)[tit$col]
tit$treatSpin<-paste(ifelse(tit$spin,'spin','noSpin'),tit$treat)
tit<-tit[order(tit$virus,tit$treatSpin,tit$treatDil,tit$moat,tit$red,tit$dil),]
tit$filter<-ave(tit$n,tit$virus,tit$treatSpin,tit$treatDil,tit$moat,tit$red,tit$col,FUN=function(xx)ifelse(1:length(xx)>=which.max(xx)&xx>10,xx,NA))

pdf('out/2018-11-18_infectionAid.pdf',width=10,height=8)
uniqTreats<-unique(tit$treatSpin)
uniqTreats<-uniqTreats[!uniqTreats %in% c('noSpin media','spin media')]
uniqVirus<-unique(tit$virus)
uniqVirus<-uniqVirus[order(uniqVirus=='Media')]
for(virus in uniqVirus){
  for(treat in uniqTreats){
    treatDils<-c(0,unique(tit$treatDil[tit$treatSpin==treat]))
    par(mfrow=c(ifelse(length(treatDils)>1,2,1),ceiling(length(treatDils)/2)))
    for(treatDil in treatDils){
      if(treatDil==0)thisDat<-tit[tit$virus==virus&tit$spin==!grepl('noSpin',treat)&tit$dil>0&!tit$red&tit$treatDil==treatDil,]
      else thisDat<-tit[tit$virus==virus&tit$treatSpin==treat&tit$dil>0&!tit$red&tit$treatDil==treatDil,]
      thisDat$logDil<--log(thisDat$dil)
      withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s %s',virus,treat,treatDil),ylim=range(tit$n+1),xlim=range(thisDat$dil),pch=21,bg=is.na(thisDat$filter)+1,cex=2))
      if(sum(!is.na(thisDat$filter))>0){
        mod<-glm(n~offset(logDil),thisDat[!is.na(thisDat$filter),],family='poisson')
        fakeDat<-data.frame(logDil=-log(1:50000))
        preds<-predict(mod,fakeDat)
        lines(exp(-fakeDat$logDil),exp(preds)+1)
        iuPerUl<-exp(coef(mod)['(Intercept)'])
        mtext(sprintf('%0.1f IU/ul (%s)',iuPerUl,treat),3)
      }
    }
  }
}
dev.off()
ius<-withAs(xx=tit[!is.na(tit$filter)&!tit$red,],by(xx[,c('n','dil')],list(paste(xx$treatSpin,sprintf('%03d',xx$treatDil)),xx$virus),function(thisDat){
    thisDat$logDil<--log(thisDat$dil)
    mod<-glm(n~offset(logDil),thisDat,family='poisson')
    exp(coef(mod)['(Intercept)'])
}))

tit<-readCounts('ice/out/2019-01-02_infectionAids3/counts.csv')
tit$red<-grepl('red',tit$plate)
tit$spin<-grepl('spin',tit$plate)
tit$dilNum<-as.numeric(sapply(strsplit(tit$plate,'_'),'[[',2))+floor((tit$col-1)/4)
#50ul in well + 50ul virus, initial 400ul neat + 500ul media, 200ul remainder + 300ul at each dilution
tit$dil<-1/(.5*(400/900)*2.5^-(tit$dilNum-1))
viruses<-read.csv('ice/2019-01-02_infectionAids3/virus.csv',stringsAsFactors=FALSE,header=FALSE)[,1]
tit$virus<-viruses[tit$rowNum+((tit$col-1)%%4)*8]
tit$treat<-ifelse(grepl('dex',tit$plate),'Dextran',ifelse(grepl('poly',tit$plate),'Polybrene','Media'))
tit$treatSpin<-paste(ifelse(tit$spin,'spin','noSpin'),tit$treat)
tit<-tit[order(tit$virus,tit$treatSpin,tit$red,tit$dil),]
tit$filter<-ave(tit$n,tit$virus,tit$treatSpin,tit$red,tit$row,FUN=function(xx)ifelse(1:length(xx)>=which.max(xx)&xx>10,xx,NA))

pdf('out/2019-01-02_infectionAid.pdf',width=10,height=8)
uniqTreats<-unique(tit$treatSpin)
uniqTreats<-uniqTreats[order(!grepl('noSpin',uniqTreats),!grepl('Media',uniqTreats),!grepl('Polybrene',uniqTreats))]
uniqVirus<-unique(tit$virus)
uniqVirus<-uniqVirus[order(uniqVirus=='Media')]
for(virus in viruses){
  par(mfrow=c(2,ceiling(length(uniqTreats)/2)))
  for(treat in uniqTreats){
    thisDat<-tit[tit$virus==virus&tit$treatSpin==treat&!tit$red,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s',virus,treat),ylim=range(tit$n+1),xlim=c(1,max(tit$dil)),pch=21,bg=is.na(thisDat$filter)+1,cex=2))
    if(sum(!is.na(thisDat$filter))>0){
      mod<-glm(n~offset(logDil),thisDat[!is.na(thisDat$filter),],family='poisson')
      fakeDat<-data.frame(logDil=-log(1:50000))
      preds<-predict(mod,fakeDat)
      lines(exp(-fakeDat$logDil),exp(preds)+1)
      iuPerUl<-exp(coef(mod)['(Intercept)'])/100
      mtext(sprintf('%0.1f IU/ul (%s)',iuPerUl,treat),3)
    }
  }
}
dev.off()

ius<-withAs(xx=tit[!is.na(tit$filter)&!tit$red,],by(xx[,c('n','dil')],list(xx$virus,xx$treatSpin),function(thisDat){
    thisDat$logDil<--log(thisDat$dil)
    mod<-glm(n~offset(logDil),thisDat,family='poisson')
    exp(coef(mod)['(Intercept)'])/100
}))
slopes<-withAs(xx=tit[!is.na(tit$filter)&!tit$red,],by(xx[,c('n','dil')],list(xx$virus,xx$treatSpin),function(thisDat){
    thisDat$logDil<--log(thisDat$dil)
    mod<-glm(n~logDil,thisDat,family='poisson')
    coef(mod)['logDil']
}))


tit<-readCounts('ice/out/2019-01-14_infectivity/counts.csv')
viruses<-read.csv('ice/20190109_virusSubset.csv',stringsAsFactors=FALSE)
tit<-tit[tit$plate!='media_4_850',]
tit$red<-grepl('red',tit$plate)
tit$dextran<-ifelse(grepl('dex',tit$plate),'Dextran','Media')
tit$dilPlate<-as.numeric(sapply(strsplit(tit$plate,'_'),'[[',2))
tit$dilNum<-c(1:9,NA,1:9,NA,1:9,NA,1:9,1:9)[(tit$dilPlate-1)*12+tit$col]
tit$plateNum<-rep(c(1,NA,2,NA,3,NA,4,5),c(9,1,9,1,9,1,9,9))[(tit$dilPlate-1)*12+tit$col]
#50ul in well + 50ul virus, initial 400ul neat + 500ul media, 200ul remainder + 300ul at each dilution
tit$dil<-1/(.5*(400/700)*2^-(tit$dilNum-1))
tit$virus<-viruses[(tit$plateNum-1)*8+tit$rowNum,'id']
tit[is.na(tit$dil),'virus']<-NA
tit<-tit[order(tit$virus,tit$dextran,tit$red,tit$dil),]
tit$filter<-ave(tit$n,tit$virus,tit$dextran,tit$red,tit$row,FUN=function(xx)ifelse(1:length(xx)>=which.max(xx)&xx>5,xx,NA))


ius2<-withAs(xx=tit[!is.na(tit$filter)&!tit$red&!is.na(tit$virus),],by(xx[,c('n','dil')],list(xx$virus,xx$dextran),function(thisDat){
    thisDat$logDil<--log(thisDat$dil)
    mod<-glm(n~offset(logDil),thisDat,family='poisson')
    exp(coef(mod)['(Intercept)'])/100
}))
round(ius2[viruses$id,c('Media','Dextran')],2)
ius2Stan<-withAs(xx=tit[!tit$red&!is.na(tit$virus),],by(xx[,c('n','dil')],list(xx$virus,xx$dextran),function(thisDat){
  fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,tit[!tit$red&is.na(tit$dil),'n'])
  mean(as.matrix(fit)[,'baseIU'])/100
}))
pdf('out/2019-01-14_infectivity.pdf',width=10,height=8)
for(virus in viruses$id){
  par(mfrow=c(1,2))
  for(treat in c('Media','Dextran')){
    thisDat<-tit[tit$virus==virus&tit$dextran==treat&!tit$red&!is.na(tit$virus),]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s',virus,treat),ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg=is.na(thisDat$filter)+1,cex=2))
    thisIu<-ius2[virus,treat]
    thisIuStan<-ius2Stan[virus,treat]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    predsStan<-thisIuStan/fakeDils*100
    lines(fakeDils,preds+1)
    lines(fakeDils,predsStan+1,col='red')
    mtext(sprintf('%0.1f IU/ul (%s)',thisIuStan,treat),3)
  }
}
dev.off()



readBigBatch<-function(countsFile,virusFile){
  viruses<-read.csv(virusFile,skip=1,stringsAsFactors=FALSE)[,-1]$id
  tit<-readCounts(countsFile)
  tit$red<-grepl('red',tit$plate)
  tit$dextran<-ifelse(grepl('dex',tit$plate),'Dextran','Media')
  tit$plateNum<-as.numeric(sapply(strsplit(tit$plate,'_'),'[[',2))
  #50ul in well + 50ul virus, initial 400ul neat + 500ul media, 200ul remainder + 300ul at each dilution
  tit$dil<-1/(.5*(400/750)*(280/160)^-(tit$col-1))
  tit$virus<-viruses[(tit$plateNum-1)*8+tit$rowNum]
  if(any(is.na(tit$virus))){
    warning('Virus not found')
    browser()
  }
  if(any(!viruses %in% tit$virus))stop('Missing virus')
  tit<-tit[order(tit$virus,tit$dextran,tit$red,tit$dil),]
  tit$filter<-ave(tit$n,tit$virus,tit$dextran,tit$red,tit$row,FUN=function(xx)ifelse(1:length(xx)>=which.max(xx)&xx>5,xx,NA))
  ius<-withAs(xx=tit[!tit$red&!is.na(tit$virus),],by(xx[,c('n','dil','filter')],list(xx$virus,xx$dextran),function(thisDat){
      thisDat$logDil<--log(thisDat$dil)
      if(all(is.na(thisDat$filter)))return(NA)
      mod<-glm(n~offset(logDil),thisDat[!is.na(thisDat$filter),],family='poisson')
      exp(coef(mod)['(Intercept)'])/100
  }))
  iusStan<-do.call(rbind,parallel::mclapply(structure(viruses,.Names=viruses),function(virus,tit){
    sapply(c('Media'='Media','Dextran'='Dextran'),function(treat){
      thisDat<-tit[!tit$red&!is.na(tit$virus)&tit$virus==virus&tit$dextran==treat,c('n','dil')]
      fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,tit[!tit$red&grepl('MEDIA|Media',tit$virus),'n'])
      return(mean(as.matrix(fit)[,'baseIU'])/100)
  })},tit,mc.cores=20))
  iusLessThanX<-withAs(xx=tit[!tit$red&!is.na(tit$virus),],by(xx[,c('n','dil','filter')],list(xx$virus,xx$dextran),function(thisDat){
      thisDat$logDil<--log(thisDat$dil)
      if(all(is.na(thisDat$filter)))return(NA)
      firstSmall<-max(c(0,which(thisDat$n>750)))+1
      mod<-glm(n~offset(logDil),thisDat[firstSmall:nrow(thisDat),],family='poisson')
      exp(coef(mod)['(Intercept)'])/100
  }))
  return(list('tit'=tit,'viruses'=viruses,'ius'=ius,'iusBayes'=iusStan,'iusX'=iusLessThanX))
}
plotBigBatch<-function(tit,viruses,ius=NULL,ius2=NULL,ius3=NULL){
  for(virus in unique(viruses)){
    par(mfrow=c(1,2))
    for(treat in c('Media','Dextran')){
      thisDat<-tit[tit$virus==virus&tit$dextran==treat&!is.na(tit$virus),]
      if(nrow(thisDat)==0){plot.new();next}
      if('red' %in% colnames(tit))tit<-tit[!tit$red,]
      thisDat$logDil<--log(thisDat$dil)
      withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s',virus,treat),ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg=if('filter' %in% colnames(thisDat$filter))is.na(thisDat$filter)+1 else 'grey',cex=2,yaxt='n'))
      logAxis(2,offset=1,axisMin=1,las=1)
      axis(2,1,0,las=1)
      fakeDils=2^seq(0,16,length.out=1000)
      if(!is.null(ius)){
        preds<-ius[virus,treat]/fakeDils*100
        lines(fakeDils,preds+1)
        mtext(sprintf('%.1f IU/ul',ius[virus,treat]),3,line=-1)
      }
      if(!is.null(ius2)){
        predsStan<-ius2[virus,treat]/fakeDils*100
        lines(fakeDils,predsStan+1,col='red')
        mtext(sprintf('%0.1f IU/ul (%s)',ius2[virus,treat],treat),3)
      }
      if(!is.null(ius3)){
        predsStan<-ius3[virus,treat]/fakeDils*100
        lines(fakeDils,predsStan+1,col='blue')
      }
    }
  }
}
makeDf<-function(tit){
  data.frame('id'=rownames(tit$iusBayes),'media'=tit$iusBayes[,'Media'],'dextran'=tit$iusBayes[,'Dextran'],'run'='2019-01-17',stringsAsFactors=FALSE)
}

tit<-dnar::cacheOperation('work/2019-01-17_titration.Rdat',readBigBatch,'ice/out/2019-01-17_titration2/counts.csv','ice/2019-01-14-titration.csv',OVERWRITE=TRUE)
out<-makeDf(tit)
pdf('out/2019-01-17_infectivity2.pdf',width=10,height=8)
  plotBigBatch(tit$tit,tit$viruses,tit$ius,tit$iusBayes,tit$iusX)
dev.off()
tit<-dnar::cacheOperation('work/2019-01-18_titration.Rdat',readBigBatch,'ice/out/2019-01-18_titration2/counts.csv','ice/2019-01-15-titration.csv',OVERWRITE=TRUE)
out<-rbind(out,makeDf(tit))
pdf('out/2019-01-18_infectivity2.pdf',width=10,height=8)
  plotBigBatch(tit$tit,tit$viruses,tit$ius,tit$iusBayes,tit$iusX)
dev.off()
tit<-dnar::cacheOperation('work/2019-01-25_titration.Rdat',readBigBatch,'ice/out/2019-01-25_titration2/counts.csv','ice/2019-01-22-titration.csv',OVERWRITE=TRUE)
out<-rbind(out,makeDf(tit))
pdf('out/2019-01-25_infectivity.pdf2',width=10,height=8)
  plotBigBatch(tit$tit,tit$viruses,tit$ius,tit$iusBayes,tit$iusX)
dev.off()
tit<-dnar::cacheOperation('work/2019-01-28_titration.Rdat',readBigBatch,'ice/out/2019-01-28_titration2/counts.csv','ice/2019-01-23-titration.csv',OVERWRITE=TRUE)
out<-rbind(out,makeDf(tit))
pdf('out/2019-01-28_infectivity2.pdf',width=10,height=8)
  plotBigBatch(tit$tit,tit$viruses,tit$ius,tit$iusBayes,tit$iusX)
dev.off()

ius<-out[grepl('^VOA|^MM|^WEAU',out$id),]
write.csv(ius,'out/ius.csv')


voa<-readCounts('ice/out/2019-06-20_voa/counts.csv')
voaVirus<-read.csv('ice/2019-06-20_voa.csv',stringsAsFactors=FALSE)$virus
voa$dilPlate<-sub('^.*_([0-9]+)$','\\1',voa$plate)
voa$baseDilution<-ifelse(voa$dilPlate==2,2,6)*2
voa$rowNum<-(8:1)[voa$rowNum]
voa$dil<-voa$baseDilution*ifelse(voa$dilPlate==2,2,2.5)^(voa$rowNum-1)
voa$virus<-paste(ifelse(grepl('new',voa$plate),'new',''),voaVirus[voa$col],sep='')
voa$virus[voa$rowNum<=5&voa$plate=='newDex_5'&voa$col==12]<-'NHG L mutant'
voa$dil[voa$rowNum<=5&voa$plate=='newDex_5'&voa$col==12]<-c(4*5^(0:3),4)[voa$rowNum[voa$rowNum<=5&voa$plate=='newDex_5'&voa$col==12]]
voa$dil[voa$rowNum>5&voa$plate=='newDex_5'&voa$col==12]<-6*2.5^(voa$rowNum[voa$rowNum>5&voa$plate=='newDex_5'&voa$col==12]-6)
voa$dextran<-ifelse(grepl('[Dd]ex',voa$plate),'Dextran','Media')

runStan<-function(viruses,tit){
  iusStan<-do.call(rbind,parallel::mclapply(structure(viruses,.Names=viruses),function(virus,tit,...){
    sapply(c('Media'='Media','Dextran'='Dextran'),function(treat){
      thisDat<-tit[!is.na(tit$virus)&tit$virus==virus&tit$dextran==treat,c('n','dil')]
      if(nrow(thisDat)==0)return(NA)
      fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,tit[grepl('MEDIA|Media',tit$virus),'n'])
      return(mean(as.matrix(fit)[,'baseIU'])/100)
  })},tit,mc.cores=20))
}
voaStan<-cacheOperation('work/voa_062019.Rdat',runStan,structure(unique(voa$virus),.Names=unique(voa$virus)),voa[,c('virus','n','dil','dextran')])

#uniqVirus<-unique(voa$virus)
uniqVirus<-c(voaVirus,'NHG L mutant')
uniqVirus<-uniqVirus[orderIn(sub('new','',uniqVirus),voaVirus)]
pdf('out/2019-06-20_voa.pdf',width=10,height=8)
par(mar=c(4,4,1,.5))
plotBigBatch(voa,uniqVirus,voaStan)
dev.off()



tit<-readCounts('ice/out/2019-07-19_reboundTropism/counts.csv')
vir<-read.csv('ice/2019-07-19_reboundTropism.csv',stringsAsFactors=FALSE)$Name
virus<-data.frame(
  'id'=c(vir[1:48],'YU2','SG3','89.6','SG3+WEAU','YU2','SG3','SG3+WEAU','Media',vir[1:40]),
  'well'=sprintf('%s%d',rep(LETTERS[1:8],12),rep(1:12,each=8)),
  'baseDil'=rep(c(1,10,4),c(52,4,40)),
  stringsAsFactors=FALSE
)
rownames(virus)<-virus$well
tit$drug<-c('a'='AMD','m'='Mar','am'='AMD+Mar','nd'='No Drug')[sub('[0-9]','',tit$plate)]
tit$dilPlate<-as.numeric(sub('[a-z]+','',tit$plate))
tit$virus<-virus[tit$well,'id']
tit$baseDil<-virus[tit$well,'baseDil']
dilutes<-cumprod(c(100/25,200/40,200/80,200/80,200/80))
tit$dil<-tit$baseDil*dilutes[tit$dilPlate]
iusStan<-do.call(rbind,cacheOperation('work/tropism_20190719.Rdat',parallel::mclapply,structure(unique(tit$virus),.Names=unique(tit$virus)),function(virus,tit,...){
  sapply(structure(unique(tit$drug),.Names=unique(tit$drug)),function(treat){
    thisDat<-tit[!is.na(tit$virus)&tit$virus==virus&tit$drug==treat,c('n','dil')]
    if(nrow(thisDat)==0)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,tit[grepl('MEDIA|Media',tit$virus),'n'])
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit,mc.cores=20))
write.csv(iusStan,'out/miniInputTitration.csv')


uniqDrugs<-unique(tit$drug)
drugColors<-structure(rainbow.lab(length(uniqDrugs)),.Names=uniqDrugs)
pdf('out/reboundTropism.pdf')
for(ii in unique(tit$virus)){
  thisDat<-tit[tit$virus==ii,]
  thisDat$logDil<--log(thisDat$dil)
  withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=ii,ylim=range(tit$n+1)*c(1,2),xlim=range(tit$dil),type='n'))
  for(dd in uniqDrugs){
    withAs(zz=thisDat[thisDat$drug==dd,],points(zz$dil,zz$n+1,pch=21,cex=2,bg=drugColors[dd],lwd=3))
  }
  legend('topright',names(drugColors),col=drugColors,lty=1,inset=.01)
  for(dd in uniqDrugs){
    fakeDil<-2^seq(0,16,.1)
    preds<-iusStan[ii,dd]*100/fakeDil
    lines(fakeDil,preds+1,col=drugColors[dd])
  }
  amd<-1-iusStan[ii,'AMD']/iusStan[ii,'No Drug']
  mar<-1-iusStan[ii,'Mar']/iusStan[ii,'No Drug']
  amdMar<-1-iusStan[ii,'AMD+Mar']/iusStan[ii,'No Drug']
  mtext(sprintf('Inhibition:\nAMD: %0.2f%%\nMar: %0.2f%%\nAMD+Mar: %0.2f%%',amd*100,mar*100,amdMar*100),3,at=exp(log(max(tit$dil))*.9))
}
dev.off()

tit<-readCounts('ice/out/2019-08-09_tropismMini/counts.csv')
vir<-read.csv('ice/2019-07-19_reboundTropism.csv',stringsAsFactors=FALSE)$Name
virus<-data.frame(
  'id'=c(vir[1:48],paste(vir[1:36],'retro'),'SG3','YU2','SG3+WEAU','89.6','SG3','YU2',paste(vir[43:48],'retro')),
  'well'=sprintf('%s%d',rep(LETTERS[1:8],each=12),rep(1:12,8)),
  'baseDil'=rep(c(1,150/40,1,5,150/40),c(48,36,4,2,6)),
  'control'=rep(c(FALSE,TRUE,FALSE),c(48+36,6,6)),
  stringsAsFactors=FALSE
)
rownames(virus)<-virus$well
tit$drug<-c('a'='AMD','m'='Mar','am'='AMD+Mar','nd'='No Drug')[sub('[0-9]','',tit$plate)]
tit$dilPlate<-as.numeric(sub('[a-z]+','',tit$plate))
tit$virus<-virus[tit$well,'id']
tit$baseDil<-virus[tit$well,'baseDil']
dilutes<-cumprod(c(100/25,rep(150/21.4,3)))
tit$dil<-tit$baseDil*dilutes[tit$dilPlate]
iusStan<-do.call(rbind,cacheOperation('work/tropism_20190813.Rdat',parallel::mclapply,structure(unique(tit$virus),.Names=unique(tit$virus)),function(virus,tit,...){
  sapply(structure(unique(tit$drug),.Names=unique(tit$drug)),function(treat){
    thisDat<-tit[!is.na(tit$virus)&tit$virus==virus&tit$drug==treat,c('n','dil')]
    if(nrow(thisDat)==0)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit,mc.cores=20))
write.csv(iusStan,'out/miniOutputTitration.csv')
uniqDrugs<-unique(tit$drug)
drugColors<-structure(rainbow.lab(length(uniqDrugs)),.Names=uniqDrugs)
pdf('out/reboundTropism2.pdf')
for(ii in unique(virus$id[order(virus$control,virus$id)])){
  thisDat<-tit[tit$virus==ii,]
  thisDat$logDil<--log(thisDat$dil)
  withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=ii,ylim=range(tit$n+1)*c(1,2),xlim=range(tit$dil),type='n'))
  for(dd in uniqDrugs){
    withAs(zz=thisDat[thisDat$drug==dd,],points(zz$dil,zz$n+1,pch=21,cex=2,bg=drugColors[dd],lwd=3))
  }
  legend('topright',names(drugColors),col=drugColors,lty=1,inset=.01)
  for(dd in uniqDrugs){
    fakeDil<-2^seq(0,16,.1)
    preds<-iusStan[ii,dd]*100/fakeDil
    lines(fakeDil,preds+1,col=drugColors[dd])
  }
  amd<-1-iusStan[ii,'AMD']/iusStan[ii,'No Drug']
  mar<-1-iusStan[ii,'Mar']/iusStan[ii,'No Drug']
  amdMar<-1-iusStan[ii,'AMD+Mar']/iusStan[ii,'No Drug']
  mtext(sprintf('Inhibition:\nAMD: %0.2f%%\nMar: %0.2f%%\nAMD+Mar: %0.2f%%',amd*100,mar*100,amdMar*100),3,at=exp(log(max(tit$dil))*.9))
}
dev.off()

tit<-readCounts('ice/out/2019-08-29_macaqueTests/counts.csv')
tit<-tit[!grepl('_red',tit$plate),]
tit<-tit[order(tit$plate,tit$row,tit$col),]

ic50<-tit[grepl('^[1-5]$',tit$plate),]
rowKey<-structure(c(LETTERS[1:4],LETTERS[1:4]),.Names=LETTERS[1:8])
ic50$dedupeRow<-rowKey[ic50$row]
ic50$dil<-4*4^as.numeric(ic50$plate)
doseKey<-c(1:6,1:6)
ic50$ifnDose<-doseKey[ic50$col]
iusStan<-do.call(rbind,cacheOperation('work/ic50_20190829.Rdat',parallel::mclapply,structure(unique(ic50$dedupeRow),.Names=unique(ic50$dedupeRow)),function(row,tit,...){
  sapply(structure(unique(ic50$col),.Names=unique(ic50$col)),function(col){
    thisDat<-tit[tit$dedupeRow==row&tit$col==col,c('n','dil')]
    if(nrow(thisDat)==0)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit=ic50,mc.cores=20))
pdf('out/ic50_20190830.pdf',width=16,height=8)
  cols<-structure(rainbow.lab(length(unique(ic50$ifnDose)),alpha=.5),.Names=unique(ic50$dil))
  for(row in unique(ic50$dedupeRow)){
    thisDat<-ic50[ic50$dedupeRow==row,]
    #plot(thisDat$n+1,thisDat$dil,pch=21,bg=cols[thisDat$ifnDose],log='xy',xaxt='n',yaxt='n',xlab='Dilution',ylab='TZMbl count',main=row)
    #logAxis(las=1)
    #logAxis(1)
    fakeDil<-1:10000
    #for(ii in unique(thisDat$col))lines(fakeDil,iusStan[row,ii]/fakeDil*100+1,col=cols[doseKey[ii]])
    par(mfrow=c(2,6))
    for(ii in unique(thisDat$col)){
      plot(thisDat[thisDat$col==ii,'dil'],thisDat[thisDat$col==ii,'n']+1,log='xy',xaxt='n',yaxt='n',xlab='Dilution',ylab='TZMbl count',main=sprintf('%s%d IU/ul=%0.1f',row,ii,iusStan[row,ii]),pch=21,bg=cols[doseKey[ii]],ylim=range(ic50$n+1))
      lines(fakeDil,iusStan[row,ii]/fakeDil*100+1,col=cols[doseKey[ii]])
      logAxis(las=1)
      logAxis(1)
    }
    par(mfrow=c(1,1))
    plot(doseKey,iusStan[row,],bg=cols[doseKey],pch=21,xlab='IFN dose',ylab='IUs/ul',cex=2)
  }
dev.off()

mac<-tit[grepl('^[a-e]$',tit$plate),]
mac[mac$n>10&grepl('[a-c]',mac$plate),]

tit<-readCounts('ice/out/2019-09-10_macaqueCheck/counts.csv')
tit<-tit[!grepl('_red',tit$plate),]
tit<-tit[order(tit$plate,tit$row,tit$col),]
tit[tit$plate==5,'n']



tit<-readCounts('ice/out/2019-09-16_macIsolation/counts.csv')
max(tit$n[grep('mac_',tit$plate)])
tit<-tit[!grepl('mac_',tit$plate),]
tit<-tit[order(tit$plate,tit$col,tit$row),]
tit[tit$n>5&tit$plate!='control',c('plate','well','n')]
tit$treat<-NA
tit$treat[tit$plate=='hum_untreated']<-ifelse(tit$row[tit$plate=='hum_untreated'] %in% LETTERS[1:4],'Untreated','EFC')
tit$treat[tit$plate=='hum_retro1']<-ifelse(tit$row[tit$plate=='hum_retro1'] %in% LETTERS[1:4],'Retro','Retro+EFC')
tit$treat[tit$plate=='hum_retro2']<-ifelse(tit$row[tit$plate=='hum_retro2'] %in% LETTERS[1:2],'Retro+Poly',ifelse(tit$row[tit$plate=='hum_retro2'] %in% LETTERS[3:5],'Retro+Dex40','Retro+Dex10'))
withAs(tit=tit[tit$plate!='control',],tapply(tit$n>5,list(tit$treat,tit$col),sum))



tit<-readCounts('ice/out/2019-09-16_infectivity/counts.csv')
tit$baseDil<-100/25
tit$dilPlate<-NA
tit$dilPlate[!grepl('control',tit$plate)]<-as.numeric(sapply(strsplit(tit$plate[!grepl('control',tit$plate)],'[_x]'),'[[',2))
tit$dilPlate[grepl('control',tit$plate)&tit$row %in% LETTERS[5:8]]<-2^tit$col[grepl('control',tit$plate)&tit$row %in% LETTERS[5:8]]
tit$dilPlate[grepl('control',tit$plate)&tit$row %in% LETTERS[1:4]]<-1
tit$dil<-tit$baseDil*tit$dilPlate
#
controlVir<-c('E'='89.6','F'='NL43','G'='YU2','H'='SG3')
vir<-read.csv('ice/miniExpansionsMAP.csv',stringsAsFactors=FALSE,header=FALSE)
vir<-vir[vir$V3!='',3:15]
if(any(vir$V3!=c(LETTERS[c(1:7,1:3,1:8)])))stop('Problem reading virus')
virStack<-data.frame('row'=rep(vir$V3,12),'col'=rep(1:12,each=nrow(vir)),'id'=unlist(vir[,2:13]),'plate'=rep(rep(c('1_reb','1_other','2_reb'),c(7,3,8)),12),'method'=rep(rep(c('bead','retro'),c(7+3+4,4)),12),stringsAsFactors=FALSE)
rownames(virStack)<-paste(virStack$plate,virStack$row,virStack$col,sep='__')
#
rowLookup<-data.frame('plate'=rep(c('rebound1','other','rebound2'),each=8),'row'=rep(LETTERS[1:8],3),'origPlate'=c(rep(c('1_reb',NA),each=4),rep(c('1_reb','1_other','2_reb'),c(3,3,2)),rep('2_reb',8)),'origRow'=c(LETTERS[1:4],rep(NA,4),LETTERS[c(5:7,1:3,1:2)],LETTERS[1:8]),stringsAsFactors=FALSE)
rownames(rowLookup)<-paste(rowLookup$plate,rowLookup$row,sep='__')
#
tit$rowLookup<-paste(sub('_.*','',tit$plate),tit$row,sep='__')
if(any(!tit$rowLookup %in% rownames(rowLookup)))stop('Unmatched row')
tit$origPlate<-rowLookup[tit$rowLookup,'origPlate']
tit$origRow<-rowLookup[tit$rowLookup,'origRow']
tit$virLookup<-paste(tit$origPlate,tit$origRow,tit$col,sep='__')
if(any(!tit$virLookup %in% rownames(virStack)&!is.na(tit$origPlate)))stop('Unmatched well')
tit$vir<-virStack[tit$virLookup,'id']
tit$method<-virStack[tit$virLookup,'method']
tit[tit$plate=='rebound1_control'&tit$row %in% LETTERS[5:8],'vir']<-controlVir[tit[tit$plate=='rebound1_control'&tit$row %in% LETTERS[5:8],'row']]
if(any(is.na(tit$vir)))stop('Missing virus')
tit[is.na(tit$method),'method']<-'control'
tit$plateMethod<-paste(tit$origPlate,tit$method)
tit$vir[tit$vir=='']<-'Media'
iusStan<-do.call(rbind,cacheOperation('work/infectivity_20190917.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
  sapply(structure(unique(tit$plateMethod),.Names=unique(tit$plateMethod)),function(method){
    thisDat<-tit[tit$vir==vir&tit$plateMethod==method,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    print(thisDat)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit=tit,mc.cores=20))
write.csv(iusStan,'out/2019-09-19_miniInfect.csv')

  
pdf('out/2019-09-18_infectivity.pdf',width=10,height=8)
for(virus in unique(tit$vir)){
  for(treat in unique(tit$plateMethod[tit$vir==virus])){
    thisDat<-tit[tit$vir==virus&tit$plateMethod==treat,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s',virus,treat),ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2))
    thisIu<-ifelse(virus %in% rownames(iusStan)&treat %in% colnames(iusStan),iusStan[rownames(iusStan)==virus,treat],NA)
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
  }
}
dev.off()

tit<-readCounts('ice/out/2019-09-19_various/counts.csv')
tit$dil<-4*2^tit$col
tit$basePlate<-sub('_.*','',tit$plate)
tit$vir<-paste(tit$basePlate,tit$row)
iusStan<-cacheOperation('work/infectivity_20190919.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
    thisDat<-tit[tit$vir==vir,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    print(thisDat)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
},tit=tit,mc.cores=20)

tit<-readCounts('ice/out/2019-09-25_ic50/counts.csv')
tit$treat[tit$plate=='macInHuman_untreated']<-ifelse(tit$row[tit$plate=='macInHuman_untreated'] %in% LETTERS[1:4],'Untreated','EFC')
tit$treat[tit$plate=='macInHuman_retro1']<-ifelse(tit$row[tit$plate=='macInHuman_retro1'] %in% LETTERS[1:4],'Retro','Retro+EFC')
tit$treat[tit$plate=='macInHuman_retro2']<-ifelse(tit$row[tit$plate=='macInHuman_retro2'] %in% LETTERS[1:2],'Retro+Poly',ifelse(tit$row[tit$plate=='macInHuman_retro2'] %in% LETTERS[3:5],'Retro+Dex40','Retro+Dex10'))
withAs(tit=tit[tit$plate!='control',],tapply(tit$n>5,list(tit$treat,tit$col),sum))
tapply(tit$n>5,list(tit$treat,tit$col),sum)



tit<-readCounts('ice/out/2019-09-26_infectivity/counts.csv')
flask<-tit[tit$plate=='flasks'&tit$row %in% LETTERS[2:7],]
flask<-flask[order(flask$row,flask$col),]
flask$reduceRow<-structure(rep(LETTERS[c(2,4,6)],each=2),names=LETTERS[2:7])[flask$row]
tapply(flask$n,list(flask$reduceRow,flask$col),max)
tit<-tit[tit$plate !='flasks',]
tit$plateDil<-as.numeric(sub('x','',sapply(strsplit(tit$plate,'_'),'[',2)))
tit$dil<-NA
tit$dil[grep('control|negative',tit$plate)]<-4*2^tit$col[grep('control|negative',tit$plate)]
tit$dil[grep('^p[12]_',tit$plate)]<-4*tit$plateDil[grep('^p[12]_',tit$plate)]
tit$basePlate<-sub('_.*','',tit$plate)
#
controlVir<-structure(c('Old 89.6','Old NL43','New YU2','New SG3','Old SG3+WEAU','New 89.6','New 89.6','Media'),.Names=LETTERS[1:8])
vir<-read.csv('ice/miniExpansionsMAP.csv',stringsAsFactors=FALSE,header=FALSE)
vir<-vir[vir$V3!='',3:15]
if(any(vir$V3!=c(LETTERS[c(1:7,1:3,1:8)])))stop('Problem reading virus')
virStack<-data.frame('row'=rep(vir$V3,12),'col'=rep(1:12,each=nrow(vir)),'id'=unlist(vir[,2:13]),'plate'=rep(rep(c('1_reb','1_other','2_reb'),c(7,3,8)),12),'method'=rep(rep(c('bead','retro'),c(7+3+4,4)),12),stringsAsFactors=FALSE)
rownames(virStack)<-paste(virStack$plate,virStack$row,virStack$col,sep='__')
#
rowLookup<-data.frame('plate'=rep(c('p1','p2'),each=8),'row'=rep(LETTERS[1:8],2),'origPlate'=c(rep(c('1_reb','1_other','2_reb'),c(3,3,2)),rep('2_reb',8)),'origRow'=c(LETTERS[c(5:7,1:3,3:4)],LETTERS[1:8]),stringsAsFactors=FALSE)
rownames(rowLookup)<-paste(rowLookup$plate,rowLookup$row,sep='__')
#
tit$rowLookup<-paste(sub('_.*','',tit$plate),tit$row,sep='__')
if(any(!tit$rowLookup[grepl('^p[12]_',tit$plate)] %in% rownames(rowLookup)))stop('Unmatched row')
tit$origPlate<-rowLookup[tit$rowLookup,'origPlate']
tit$origRow<-rowLookup[tit$rowLookup,'origRow']
tit$virLookup<-paste(tit$origPlate,tit$origRow,tit$col,sep='__')
if(any(!tit$virLookup %in% rownames(virStack)&!is.na(tit$origPlate)))stop('Unmatched well')
tit$vir<-virStack[tit$virLookup,'id']
tit$method<-virStack[tit$virLookup,'method']
tit[grepl('control',tit$plate),'vir']<-controlVir[tit[grepl('control',tit$plate),'row']]
tit[grepl('negative',tit$plate),'vir']<-'Media'
if(any(is.na(tit$vir)))stop('Missing virus')
tit[is.na(tit$method),'method']<-'control'
tit$plateMethod<-paste(tit$origPlate,tit$method)
tit$vir[tit$vir=='']<-'Media'
#
iusStan<-do.call(rbind,cacheOperation('work/infectivity_20190927.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
  sapply(structure(unique(tit$plateMethod),.Names=unique(tit$plateMethod)),function(method){
    thisDat<-tit[tit$vir==vir&tit$plateMethod==method,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    print(thisDat)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit=tit,mc.cores=20))
write.csv(iusStan,'out/2019-09-27_miniInfect.csv')
#
pdf('out/2019-09-27_infectivity.pdf',width=10,height=8)
for(virus in unique(tit$vir)){
  for(treat in unique(tit$plateMethod[tit$vir==virus])){
    thisDat<-tit[tit$vir==virus&tit$plateMethod==treat,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s',virus,treat),ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2))
    thisIu<-ifelse(virus %in% rownames(iusStan)&treat %in% colnames(iusStan),iusStan[rownames(iusStan)==virus,treat],NA)
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
  }
}
dev.off()



tit<-readCounts('ice/out/2019-10-04_macaqueInHuman/counts.csv')
tapply(tit$n,list(tit$row,tit$col,tit$plate),c)

tit<-readCounts('ice/out/2019-10-09_various/counts.csv')
tit<-tit[!grepl('red',tit$plate),]
tit[tit$plate=='isolation_flasks_2','plate']<-'flasks_2'
withAs(tit=tit[grepl('isolation',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
withAs(tit=tit[grepl('flask',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
tit<-tit[grep('^ic50',tit$plate),]
tit$plateDil<-as.numeric(sub('x','',sapply(strsplit(tit$plate,'_'),'[',2)))
tit$plateDil[tit$row=='H']<-2^(tit$col[tit$row=='H']-1)
tit$dil<-4*tit$plateDil
tit$vir<-tit$well
tit$vir[tit$row=='G']<-'Media'
tit$vir[tit$row=='H']<-'SG3'
iusStan<-unlist(cacheOperation('work/tzmbl_ic50_20191009.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
    thisDat<-tit[tit$vir==vir,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
},tit=tit,mc.cores=20))
write.csv(data.frame('infect'=iusStan),'out/2019-10-09_tzmblIc50.csv')
pdf('out/2019-10-09_tzmblIc50.pdf',width=10,height=8)
for(virus in unique(tit$vir[order(tit$row,tit$col)])){
    thisDat<-tit[tit$vir==virus,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=virus,ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2))
    thisIu<-iusStan[virus]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
}
dev.off()
round(tapply(iusStan,list(sub('^([A-Z])([0-9]+)$','\\1',names(iusStan)),as.numeric(ifelse(grepl('^([A-Z])([0-9]+)$',names(iusStan)),sub('^([A-Z])([0-9]+)$','\\2',names(iusStan)),1))),c))




tit<-readCounts('ice/out/2019-10-18_isolation/counts.csv')
withAs(tit=tit[!grepl('red',tit$plate)&grepl('isolation',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
tit<-tit[grep('^ic50',tit$plate),]
tit$plateDil<-as.numeric(sub('x','',sapply(strsplit(tit$plate,'_'),'[',2)))
tit$plateDil[tit$row %in% c('G','H')]<-2^(tit$col[tit$row %in% c('G','H')]-1)*tit$plateDil[tit$row %in% c('G','H')]/10
tit$dil<-2*tit$plateDil
tit$vir<-tit$well
tit$vir[tit$row=='H']<-'Media'
tit$vir[tit$row=='G']<-'YU2'
iusStan<-unlist(cacheOperation('work/tzmbl_ic50_20191018.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
    thisDat<-tit[tit$vir==vir,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
},tit=tit,mc.cores=20))
write.csv(data.frame('infect'=iusStan),'out/2019-10-18_tzmblIc50.csv')
pdf('out/2019-10-18_tzmblIc50.pdf',width=10,height=8)
for(virus in unique(tit$vir[order(tit$row,tit$col)])){
    thisDat<-tit[tit$vir==virus,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=virus,ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2,xaxt='n'))
    logAxis(1)
    thisIu<-iusStan[virus]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
}
dev.off()
icTable<-tapply(iusStan,list(sub('^([A-Z])([0-9]+)$','\\1',names(iusStan)),as.numeric(ifelse(grepl('^([A-Z])([0-9]+)$',names(iusStan)),sub('^([A-Z])([0-9]+)$','\\2',names(iusStan)),1))),c)
pdf('out/2019-10-19_ic50.pdf',height=10,width=4)
par(mfrow=c(6,2),mar=c(5,4,1,1))
for(ii in 1:6){
  ns<-icTable[ii,]
  doses<-c(concBeta6,concBeta6)
  zeroFill<-doses[2]/100
  plot(ifelse(doses==0,zeroFill,doses),ns,log='xy',ylab='IU/ul',xlab='IFNb dose',xaxt='n',yaxt='n',main=rownames(icTable)[ii],mgp=c(2.5,.7,0))
  logAxis(las=1)
  logAxis(1,axisMin=min(concBeta6[-1]))
  axis(1,zeroFill,0)
  plot(ifelse(doses==0,zeroFill,doses),ns/exp(mean(log(ns[c(1,7)]))),log='xy',ylab='Proportion max IU',xlab='IFNb dose',xaxt='n',yaxt='n',main=rownames(icTable)[ii],mgp=c(2.5,.7,0))
  logAxis(las=1,mgp=c(3,.8,0))
  logAxis(1,axisMin=min(concBeta6[-1]))
  abline(h=.5,lty=2)
  axis(1,zeroFill,0)
}
dev.off()


tit<-readCounts('ice/out/2019-10-25_isolation/counts.csv')
withAs(tit=tit[!grepl('red',tit$plate)&grepl('isolation',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))

tit<-readCounts('ice/out/2019-10-31_tropism/counts.csv')
withAs(tit=tit[!grepl('red',tit$plate)&grepl('macIso',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
tit<-tit[!grepl('macIso',tit$plate),]
tit$drug<-c('a'='AMD','m'='Mar','am'='AMD+Mar','nd'='No drug')[sapply(strsplit(tit$plate,'[-_]'),'[',1)]
tit$dilPlate<-sapply(strsplit(tit$plate,'[-_]'),'[',2)
vir<-read.csv('data/tropsimIds.csv',stringsAsFactors=FALSE) #from Tropisms and Infect. IDs email
rownames(vir)<-vir$Wells
tit$virus<-vir[tit$well,'IDs']
tit$virus[is.na(tit$virus)]<-'Media'
tit$virus[is.na(tit$virus)|grepl('media',tit$plate)]<-'Media'
tit$drug[grepl('media',tit$plate)]<-'No drug'
#tit$dil<-(275/50)^as.numeric(tit$dilPlate)
tit$dil<-4*(275/50)^as.numeric(tit$dilPlate) #added *4 12/31/2019

iusStan<-do.call(rbind,cacheOperation('work/tropism_20191031.Rdat',parallel::mclapply,structure(unique(tit$virus),.Names=unique(tit$virus)),function(virus,tit,...){
  sapply(structure(unique(tit$drug),.Names=unique(tit$drug)),function(treat){
    thisDat<-tit[!is.na(tit$virus)&tit$virus==virus&tit$drug==treat,c('n','dil')]
    if(nrow(thisDat)==0)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit,mc.cores=20))
write.csv(iusStan,'out/tropism_20191031.csv')
uniqDrugs<-unique(tit$drug)
drugColors<-structure(rainbow.lab(length(uniqDrugs)),.Names=uniqDrugs)
pdf('out/tropism_20191031.pdf')
for(ii in vir$IDs){
  thisDat<-tit[tit$virus==ii,]
  thisDat$logDil<--log(thisDat$dil)
  withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=ii,ylim=range(tit$n+1)*c(1,2),xlim=range(tit$dil),type='n'))
  for(dd in uniqDrugs){
    withAs(zz=thisDat[thisDat$drug==dd,],points(zz$dil,zz$n+1,pch=21,cex=2,bg=drugColors[dd],lwd=3))
  }
  legend('topright',names(drugColors),col=drugColors,lty=1,inset=.01)
  for(dd in uniqDrugs){
    fakeDil<-2^seq(0,16,.1)
    preds<-iusStan[ii,dd]*100/fakeDil
    lines(fakeDil,preds+1,col=drugColors[dd])
  }
  amd<-1-iusStan[ii,'AMD']/iusStan[ii,'No drug']
  mar<-1-iusStan[ii,'Mar']/iusStan[ii,'No drug']
  amdMar<-1-iusStan[ii,'AMD+Mar']/iusStan[ii,'No drug']
  mtext(sprintf('Inhibition:\nAMD: %0.2f%%\nMar: %0.2f%%\nAMD+Mar: %0.2f%%',amd*100,mar*100,amdMar*100),3,at=exp(log(max(tit$dil))*.9))
}
dev.off()
percInhib<-(1-iusStan[,colnames(iusStan)!='No drug']/iusStan[,'No drug'])*100
maxN<-tapply(tit$n,tit$virus,max)
write.csv(percInhib[names(maxN)[maxN>20],],'out/tropismPercentInhibition_20191108.csv')









tit<-readCounts('ice/out/2019-11-04_t12Infectivity/counts.csv')
withAs(tit=tit[grepl('control',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
tit$dilPlate<-sub('x','',sapply(strsplit(tit$plate,'_'),'[',2))
tit$isControl<-grepl('control',tit$plate)
tit$isT12<-grepl('t12',tit$plate)
tit[tit$isControl,'dil']<-4*2^tit$col[tit$isControl] #2 fold dilution series starting at 2x and 25ul into 75ul
tit[!tit$isControl,'dil']<-4*2*as.numeric(tit$dilPlate[!tit$isControl]) #50ul-50ul to start (low volume) and 25ul into 75ul
controlVirs<-c('89.6','SG3','YU2','A08.21-2F4','A09-IMC','601-r1','SL92b','SG3+WEAU')
tit$vir[tit$isControl]<-controlVirs[tit$rowNum[tit$isControl]]
vir<-read.csv('ice/miniExpansionsMAP.csv',stringsAsFactors=FALSE,header=FALSE)
vir<-vir[vir$V3!='',3:15]
if(any(vir$V3!=c(LETTERS[c(1:7,1:3,1:8)])))stop('Problem reading virus')
virStack<-data.frame('row'=rep(vir$V3,12),'col'=rep(1:12,each=nrow(vir)),'id'=unlist(vir[,2:13]),'plate'=rep(rep(c('1_reb','1_other','2_reb'),c(7,3,8)),12),'method'=rep(rep(c('bead','retro'),c(7+3+4,4)),12),stringsAsFactors=FALSE)
rownames(virStack)<-paste(virStack$plate,virStack$row,virStack$col,sep='__')
#
rowLookup<-data.frame('plate'=rep(c('p1','p2'),each=8),'row'=rep(LETTERS[1:8],2),'origPlate'=c(rep(c('1_reb','1_other','2_reb'),c(3,3,2)),rep('2_reb',8)),'origRow'=c(LETTERS[c(5:7,1:3,3:4)],LETTERS[1:8]),stringsAsFactors=FALSE)
rownames(rowLookup)<-paste(rowLookup$plate,rowLookup$row,sep='__')
#
tit$rowLookup<-paste(sub('_.*','',tit$plate),tit$row,sep='__')
if(any(!tit$rowLookup[grepl('^p[12]_',tit$plate)] %in% rownames(rowLookup)))stop('Unmatched row')
tit$origPlate<-rowLookup[tit$rowLookup,'origPlate']
tit$origRow<-rowLookup[tit$rowLookup,'origRow']
tit$virLookup<-paste(tit$origPlate,tit$origRow,tit$col,sep='__')
if(any(!tit$virLookup %in% rownames(virStack)&!is.na(tit$origPlate)))stop('Unmatched well')
tit$vir[!tit$isControl]<-virStack[tit$virLookup[!tit$isControl],'id']
tit$vir[tit$vir=='']<-'Media'
tit$method<-virStack[tit$virLookup,'method']
tit$method[tit$isControl]<-'control'
tit$plateMethod<-paste(tit$origPlate,tit$method,ifelse(tit$isT12,'+T12',''))
#
iusStan<-do.call(rbind,cacheOperation('work/infectivityT12_20191104.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
  sapply(structure(unique(tit$plateMethod),.Names=unique(tit$plateMethod)),function(method){
    thisDat<-tit[tit$vir==vir&tit$plateMethod==method,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit=tit,mc.cores=20))
write.csv(iusStan,'out/2019-11-04_miniInfectT12.csv')

pdf('out/2019-11-04_infectivityT12.pdf',width=10,height=8)
for(virus in unique(tit$vir)){
  for(treat in sort(unique(tit$plateMethod[tit$vir==virus]))){
    thisDat<-tit[tit$vir==virus&tit$plateMethod==treat,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s %s',virus,treat),ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2))
    thisIu<-ifelse(virus %in% rownames(iusStan)&treat %in% colnames(iusStan),iusStan[rownames(iusStan)==virus,treat],NA)
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
  }
}
dev.off()

tit<-tit[order(tit$plate,tit$row,tit$col),]
minusT12<-tit[!tit$isT12,]
plusT12<-tit[tit$isT12,]
selector<-minusT12$n>10&plusT12$n>10&minusT12$n<1000&plusT12$n<1000
minusT12<-minusT12[selector,]
plusT12<-plusT12[selector,]
if(!all(minusT12$vir==plusT12$vir)&all(minusT12$vir==plusT12$vir))stop('Mismatch in minus vs plus T12')
minusIU<-iusStan[,c('2_reb retro ','2_reb bead ','1_other bead ','NA control ','1_reb bead ')]
plusIU<-iusStan[,sprintf('%s+T12',colnames(minus))]
plusMinusIU<-data.frame('minus'=as.vector(minusIU),'plus'=as.vector(plusIU))
plusMinusIU<-plusMinusIU[!is.na(plusMinusIU$minus),]
if(any(is.na(plusMinusIU$plus)))stop('Mismatch in IU NAs')
pdf('out/t12Compare.pdf')
par(mar=c(4,4,1.1,.1),las=1)
plot(minusT12$n,plusT12$n,log='xy',xlab='IU with no T1249',ylab='IU with T1249')
abline(0,1,col='red')
thisFit<-lm(plus~offset(minus),data.frame(minus=log(minusT12$n),plus=log(plusT12$n)))
fakeDat<-data.frame('minus'=log(1:5000))
preds<-predict(thisFit,fakeDat,interval='conf')
lines(exp(fakeDat$minus),exp(preds[,'fit']),col='blue')
polygon(exp(c(fakeDat$minus,rev(fakeDat$minus))),exp(c(preds[,'lwr'],rev(preds[,'upr']))),col='#0000FF33',border=NA)
thisFit<-lm(plus~offset(minus),data.frame(minus=log(minusT12$area),plus=log(plusT12$area)))
fakeDat<-data.frame('minus'=seq(log(100),log(100000),length.out=1000))
preds<-predict(thisFit,fakeDat,interval='conf')
plot((minusT12$area),(plusT12$area),log='xy',xlab='Area with no T1249',ylab='Area with T1249',xaxt='n',yaxt='n')
lines(exp(fakeDat$minus),exp(preds[,'fit']),col='blue')
polygon(exp(c(fakeDat$minus,rev(fakeDat$minus))),exp(c(preds[,'lwr'],rev(preds[,'upr']))),col='#0000FF33',border=NA)
logAxis(las=1);logAxis(1)
abline(0,1,col='red')
#plot((minusT12$area)/(minusT12$n),(plusT12$area)/(plusT12$n),log='xy',xlab='Area/IU with no T1249',ylab='Area/IU with T1249')
#abline(0,1,col='red')
plot(plusMinusIU$minus,plusMinusIU$plus,log='xy',xlab='IU/ul with no T1249',ylab='IU/ul with T1249',xaxt='n',yaxt='n')
thisFit<-lm(plus~offset(minus),data.frame(minus=log(plusMinusIU$minus),plus=log(plusMinusIU$plus)))
fakeDat<-data.frame('minus'=seq(log(.01),log(1e4),length.out=1000))
preds<-predict(thisFit,fakeDat,interval='conf')
abline(0,1,col='red')
lines(exp(fakeDat$minus),exp(preds[,'fit']),col='blue')
polygon(exp(c(fakeDat$minus,rev(fakeDat$minus))),exp(c(preds[,'lwr'],rev(preds[,'upr']))),col='#0000FF33',border=NA)
logAxis(las=1);logAxis(1)
dev.off()


tit<-readCounts('ice/out/2019-11-08_ic50/counts.csv')
withAs(tit=tit[grepl('isolation',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
tit<-tit[grep('imc',tit$plate),]
tit$plateDil<-as.numeric(sub('x','',sapply(strsplit(tit$plate,'_'),'[',2)))
tit$dil<-4*tit$plateDil
#20 and 50 were stained at normal time, rest were stained a day early
tit$vir<-tit$well
tit$vir[tit$plateDil %in% c('20','50')]<-sprintf('%s full',tit$vir[tit$plateDil %in% c('20','50')])
iusStan<-unlist(cacheOperation('work/tzmbl_ic50_20191108.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
    thisDat<-tit[tit$vir==vir,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
},tit=tit,mc.cores=20))
#write.csv(data.frame('infect'=iusStan),'out/2019-10-18_tzmblIc50.csv')

pdf('out/2019-11-08_tzmblIc50.pdf',width=10,height=8)
for(virus in unique(tit$vir[order(tit$row,tit$col)])){
    thisDat<-tit[tit$vir==virus,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=virus,ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2,xaxt='n'))
    logAxis(1)
    thisIu<-iusStan[virus]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
}
dev.off()

vir<-c('A'='MM33 TF_IMC', 'B'='MM33.13_IMC', 'C'='MM33.17_IMC', 'D'='A08.21-2P4_IMC', 'E'='A09r-1A2_IMC', 'F'='601r1__IMC', 'G'='MM33.1.13C1', 'H'='MM33.13.2D6') 
icTable<-withAs(iusStan=iusStan[grepl('full',names(iusStan))],tapply(iusStan,list(sub('^([A-Z])([0-9]+) full$','\\1',names(iusStan)),as.numeric(sub('^([A-Z])([0-9]+) full$','\\2',names(iusStan)),1)),c))
icTable2<-withAs(iusStan=iusStan[!grepl('full',names(iusStan))],tapply(iusStan,list(sub('^([A-Z])([0-9]+)$','\\1',names(iusStan)),as.numeric(ifelse(grepl('^([A-Z])([0-9]+)$',names(iusStan)),sub('^([A-Z])([0-9]+)$','\\2',names(iusStan)),1))),c))
ic50<-t(apply(icTable[,unlist(lapply(1:6,function(xx)c(xx,xx+6)))],1,function(xx)calculateBasicIc50(concBeta6,matrix(xx,ncol=length(xx)),vocal=TRUE)))
rownames(ic50)<-vir[rownames(ic50)]
write.csv(ic50,'out/2019-10-18_tzmblIc50Calc.csv')
ic50_2<-t(apply(icTable2[,unlist(lapply(1:6,function(xx)c(xx,xx+6)))],1,function(xx)calculateBasicIc50(concBeta6,matrix(xx,ncol=length(xx)),vocal=TRUE)))
rownames(ic50_2)<-vir[rownames(ic50_2)]
ic50s<-cbind(ic50[,c('ic50','percVres')],ic50_2[,c('ic50','percVres')])
colnames(ic50s)<-paste(rep(c('TZMbl full','TZMbl short'),each=2),colnames(ic50s))
p24s<-read.csv('out/imc_20191105_ic50.csv',row.names=1)
p24s50<-p24s[grepl('50$',rownames(p24s)),]
p24s100<-p24s[grepl('100$',rownames(p24s)),]
rownames(p24s50)<-sub('_human.*','',rownames(p24s50))
rownames(p24s100)<-sub('_human.*','',rownames(p24s100))
tmp<-cbind(p24s50[rownames(ic50s),c('ic50','vres')],p24s100[rownames(ic50s),c('ic50','vres')])
colnames(tmp)<-paste(rep(c('p24 1:50','p24 1:100'),each=2),colnames(tmp))
out<-cbind(ic50s,tmp)
colnames(out)<-sub('percVres','vres',colnames(out))
round(out[,grep('ic50',colnames(out))],3)
round(out[,grep('vres',colnames(out))],1)



pdf('out/2019-11-08_ic50.pdf',height=16,width=4)
par(mfrow=c(8,2),mar=c(5,4,1,1))
for(jj in 1:2){
for(ii in 1:8){
  if(jj==1)ns<-icTable[ii,]
  else ns<-icTable2[ii,]
  doses<-c(concBeta6,concBeta6)
  zeroFill<-doses[2]/100
  plot(ifelse(doses==0,zeroFill,doses),ns,log='xy',ylab='IU/ul',xlab='IFNb dose',xaxt='n',yaxt='n',main=rownames(icTable)[ii],mgp=c(2.5,.7,0))
  logAxis(las=1)
  logAxis(1,axisMin=min(concBeta6[-1]))
  axis(1,zeroFill,0)
  plot(ifelse(doses==0,zeroFill,doses),ns/exp(mean(log(ns[c(1,7)]))),log='xy',ylab='Proportion max IU',xlab='IFNb dose',xaxt='n',yaxt='n',main=paste(vir[rownames(icTable)[ii]],ifelse(jj==1,'full','short')),mgp=c(2.5,.7,0))
  logAxis(las=1,mgp=c(3,.8,0))
  logAxis(1,axisMin=min(concBeta6[-1]))
  abline(h=.5,lty=2)
  thisIc<-if(jj==1)ic50 else ic50_2
  abline(v=thisIc[vir[rownames(icTable)[ii]],'ic50'],lty=2)
  axis(1,zeroFill,0)
}}
dev.off()

tit<-readCounts('ice/out/2019-11-08_ic50/counts.csv')
tit<-tit[grep('mac',tit$plate),]
tit$plateDil<-as.numeric(sub('x','',sapply(strsplit(tit$plate,'_'),'[',2)))
tit$dil<-4*tit$plateDil
tit$vir<-tit$well
iusStan<-unlist(cacheOperation('work/tzmbl_mac_20191108.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
    thisDat<-tit[tit$vir==vir,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
},tit=tit,mc.cores=20))
write.csv(data.frame('infect'=iusStan),'out/2019-10-18_tzmblMacIc50.csv')

pdf('out/2019-11-08_tzmblMacIc50.pdf',width=10,height=8)
for(virus in unique(tit$vir[order(tit$row,tit$col)])){
    thisDat<-tit[tit$vir==virus,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=virus,ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2,xaxt='n'))
    logAxis(1)
    thisIu<-iusStan[virus]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
}
dev.off()

vir<-c('A'='SHIV.D.375M', 'B'='SHIV.D.RRj17-H11', 'C'='SHIV.D.RIa17-G4', 'D'='SHIV.D.RRj17.15-R3', 'E'='SHIV.D.RRj17.15-R5', 'F'='SL92b', 'G'='SHIV.D.RRj17.15-R2', 'H'='SHIV.AD8EO')
icTable<-withAs(iusStan=iusStan,tapply(iusStan,list(sub('^([A-Z])([0-9]+)$','\\1',names(iusStan)),as.numeric(sub('^([A-Z])([0-9]+)$','\\2',names(iusStan)),1)),c))
ic50s<-t(apply(icTable[,unlist(lapply(1:6,function(xx)c(xx,xx+6)))],1,function(xx)calculateBasicIc50(concBeta6,matrix(xx,ncol=length(xx)),vocal=TRUE)))
pdf('out/2019-10-19_macIc50.pdf',height=10,width=4)
par(mfrow=c(8,2),mar=c(5,4,1,1))
for(ii in 1:8){
  ns<-icTable[ii,]
  doses<-c(concBeta6,concBeta6)
  zeroFill<-doses[2]/100
  plot(ifelse(doses==0,zeroFill,doses),ns,log='xy',ylab='IU/ul',xlab='IFNb dose',xaxt='n',yaxt='n',main=rownames(icTable)[ii],mgp=c(2.5,.7,0))
  logAxis(las=1)
  logAxis(1,axisMin=min(concBeta6[-1]))
  axis(1,zeroFill,0)
  plot(ifelse(doses==0,zeroFill,doses),ns/exp(mean(log(ns[c(1,7)]))),log='xy',ylab='Proportion max IU',xlab='IFNb dose',xaxt='n',yaxt='n',main=vir[rownames(icTable)[ii]],mgp=c(2.5,.7,0))
  logAxis(las=1,mgp=c(3,.8,0))
  logAxis(1,axisMin=min(concBeta6[-1]))
  abline(h=.5,lty=2)
  axis(1,zeroFill,0)
  abline(v=ic50s[ii,'ic50'],lty=2)
}
dev.off()
out<-ic50s[,c('ic50','percVres')]
colnames(out)[2]<-'vres'
rownames(out)<-vir[rownames(out)]






tit<-readCounts('ice/out/2019-11-11_isolation/counts.csv')
withAs(tit=tit[grepl('iso',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
tit<-tit[!grepl('iso',tit$plate),]
tit$vir<-tit$row
tit$method<-sub('control_','',sub('_50','_',tit$plate))
tit$dil<-2^tit$col
iusStan<-do.call(rbind,cacheOperation('work/controls_20191111.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
  sapply(structure(unique(tit$method),.Names=unique(tit$method)),function(method){
    thisDat<-tit[tit$vir==vir&tit$method==method,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit=tit,mc.cores=20))
round(iusStan[sort(rownames(iusStan)),])


tit<-readCounts('ice/out/2019-11-18_isolations/counts.csv'),'area'
withAs(tit=tit,tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
print(withAs(zz=tit[!grepl('fred|control',tit$plate)&tit$n>5&(tit$col!=1|grepl('pbmc',tit$plate)),],zz[order(zz$plate,zz$well),c('plate','well','n')]),row.names=FALSE)

tit<-readCounts('ice/out/2019-11-25_isolationsFred/counts.csv')
withAs(tit=tit,tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
print(withAs(zz=tit[!grepl('fred|control',tit$plate)&tit$n>5&(tit$col!=1|grepl('pbmc',tit$plate)),],zz[order(zz$plate,zz$well),c('plate','well','n')]),row.names=FALSE)
vir<-read.csv('ice/titration sample list_fred_20191122.csv',header=FALSE)[,-1]
vir<-rbind(vir,data.frame('V2'=c('YU2_1','SG3_1','896_1','Media','Media'),'V3'=''))
vir$id<-trimws(paste(vir$V2,vir$V3))
cols<-structure(c(1:5,1,1:5,1),.Names=12:1)
tit<-tit[grepl('fred|control',tit$plate),]
tit$plateDil<-suppressWarnings(as.numeric(sub('x','',sapply(strsplit(tit$plate,'_'),'[',2))))
tit$dil<-4*tit$plateDil*ifelse(tit$col<7,2,1)
tit$vir<-tit$well
tit$virId<-tit$rowNum+(cols[as.character(tit$col)]-1)*8
tit$vir<-vir[tit$virId,'id']
tit[grep('control',tit$plate),'dil']<-ifelse(tit[grep('control',tit$plate),'rowNum'] %in% 4:6,3,2)*ifelse(tit[grep('control',tit$plate),'rowNum'] %in% 4:6,240/100,2)^tit[grep('control',tit$plate),'col']
tit[grep('control',tit$plate),'vir']<-c(rep(c('YU2','SG3','896'),2),'Media','Media')[tit[grep('control',tit$plate),'rowNum']]
iusStan<-unlist(cacheOperation('work/fred_20191126.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
    thisDat<-tit[tit$vir==vir,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    #setting disperse to .1 to deal with low max counts
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1,disperse=.1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
},tit=tit,mc.cores=20))
pdf('out/2019-11-25_infectivity.pdf',width=10,height=8)
for(virus in unique(tit$vir[order(grepl('control',tit$plate),tit$virId)])){
    thisDat<-tit[tit$vir==virus,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=virus,ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2,xaxt='n'))
    logAxis(1)
    thisIu<-iusStan[virus]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
}
dev.off()
vir$iuPerUl<-iusStan[vir$id]
write.csv(vir,'out/20191126_infectivity.csv')


tit<-readCounts('ice/out/2019-12-02_isolationIc50Ifn/counts.csv')
withAs(tit=tit[!grepl('ifn',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
print(withAs(zz=tit[!grepl('ifn',tit$plate)&tit$n>5&(tit$col!=1|grepl('pbmc',tit$plate)),],zz[order(zz$plate,zz$well),c('plate','well','n')]),row.names=FALSE)
tit<-tit[grepl('^ifn_',tit$plate),]
tit$plate<-ifelse(tit$plate=='ifn_3x','ifn_15x',ifelse(tit$plate=='ifn_15x','ifn_3x',tit$plate))
tit$dil<-as.numeric(sub('x$','',sapply(strsplit(tit$plate,'_'),'[',2)))*4
vir<-paste(rep(c('Fresh','Incubated'),each=4),rep(c('YU2','SG3+WEAU','MM33.TF','MM33.13'),2))
ifnCols<-c(
  'Untreated 0',
  sprintf('IFNa2 %s',sapply(5.5/9^(0:4),format,scientific=FALSE,digits=2)),
  'Untreated 0',
  sprintf('IFNb %s',sapply(4400/100^(0:4),format,scientific=FALSE,digits=2))
)
tit$ifn<-ifnCols[tit$col]
tit$vir<-vir[tit$rowNum]
iusStan<-do.call(rbind,cacheOperation('work/ifn_20191202.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
  sapply(structure(unique(tit$ifn),.Names=unique(tit$ifn)),function(ifn){
    thisDat<-tit[tit$vir==vir&tit$ifn==ifn,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit=tit,mc.cores=20))
round(iusStan[sort(rownames(iusStan)),sort(colnames(iusStan))])
pdf('out/2019-12-02_ifn.pdf',width=10,height=8)
for(virus in unique(tit$vir)){
  for(ifn in unique(ifnCols)){
    thisDat<-tit[tit$vir==virus&tit$ifn==ifn,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=paste(virus,ifn),ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2,xaxt='n'))
    logAxis(1)
    thisIu<-iusStan[virus,ifn]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
  }
}
dev.off()
pdf('out/ifn_vs_tzmbl.pdf',width=8,height=12)
par(mfrow=c(4,2),mar=c(4.5,4,1,1))
for(ii in 1:2){
for(vir in sort(unique(sub('[^ ]+ ','',rownames(iusStan))))){
  for(ifn in c("IFNa2","IFNb")){
    thisAll<-iusStan[grepl(vir,rownames(iusStan),fixed=TRUE),]
    thisDat<-thisAll[,grepl(sprintf('%s|Untreated',ifn),colnames(iusStan))]
    doses<-as.numeric(sapply(strsplit(colnames(thisDat),' '),'[',2))
    plot(1,1,type='n',log=ifelse(ii==1,'xy','x'),xlim=range(doses[doses>0]),ylim=c(ifelse(ii==1,1,min(thisAll)),max(thisAll)),xlab=sprintf('%s dose',ifn),ylab='Infectivity (IU/ul)',xaxt='n',las=1,main=vir)
    dnar::logAxis(1)
    points(doses,thisDat[sprintf('Fresh %s',vir),],pch=21,bg='blue',cex=1.5)
    points(doses,thisDat[sprintf('Incubated %s',vir),],pch=21,bg='red',cex=1.5)
    abline(h=thisDat[paste(c('Fresh','Incubated'),vir),'Untreated 0'],lty=2,col=c('blue','red'))
    legend('bottomleft',c('Fresh IFN','Incubated IFN'),pch=21,pt.bg=c('blue','red'),inset=.01)
  }
}
}
dev.off()

tit<-readCounts('ice/out/2019-12-02_isolationIc50Ifn/counts.csv')
tit<-tit[grepl('^ifn[a2b]+_',tit$plate),]
tit$pNum<-sapply(strsplit(tit$plate,'_'),'[',2)
tit$dil<-as.numeric(sub('x$','',sapply(strsplit(tit$plate,'_'),'[',3)))*4
tit$ifn<-sapply(strsplit(tit$plate,'_'),'[',1)
virs<-read.csv('ice/2019-12-02_ic50.csv',stringsAsFactors=FALSE)
rownames(virs)<-sprintf('%s_p%d_%s',tolower(virs$ifnLabel),virs$Plate,LETTERS[virs$Row])
#remove row mixup
tit<-tit[!(grepl('ifnb_p3_512x',tit$plate)&tit$row=='D'),]
tit[grepl('ifnb_p3_512x',tit$plate)&tit$row=='F','row']<-'D'
tit$lookup<-sprintf('%s_%s_%s',tit$ifn,tit$pNum,tit$row)
if(any(!tit$lookup %in% rownames(virs)))stop('Missing virus info')
tit$vir<-paste(tit$pNum,virs[tit$lookup,'Samples'])
tit$ifn<-virs[tit$lookup,'IFN']
tit$ifnDose<-sprintf('%s__%s',tit$ifn,tit$col)

iusStan<-do.call(rbind,cacheOperation('work/tzmblIc50_20191202.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
  sapply(structure(unique(tit$ifnDose),.Names=unique(tit$ifnDose)),function(col){
    thisDat<-tit[tit$vir==vir&tit$ifnDose==col,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit=tit,mc.cores=20))
round(iusStan[sort(rownames(iusStan)),])


pdf('out/2019-12-02_tzmblIc50.pdf',width=10,height=8)
for(virus in unique(tit$vir[order(tit$plate,tit$row,tit$col)])){
  for(col in unique(tit$ifnDose[order(tit$plate,tit$row,tit$col)])){
    thisDat<-tit[tit$vir==virus&tit$ifnDose==col,]
    thisDat$logDil<--log(thisDat$dil)
    browser()
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=paste(virus,col),ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2,xaxt='n'))
    logAxis(1)
    thisIu<-iusStan[virus,col]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
  }
}
dev.off()



tit<-readCounts('ice/out/2019-12-09_isolation/counts.csv')
withAs(tit=tit[!grepl('control|infectivity',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
withAs(tit=tit[grepl('control',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
print(withAs(zz=tit[!grepl('control|infectivity',tit$plate)&tit$n>5&(tit$col!=1|grepl('pbmc',tit$plate)),],zz[order(zz$plate,zz$well),c('plate','well','n')]),row.names=FALSE)
withAs(tit=tit[grepl('control',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
zz<-withAs(tit=tit[grepl('control',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
pdf('out/rinse_vs_noRinse.pdf')
plot(data.frame('rinse'=as.vector(zz[,,'control']),'noRinse'=as.vector(zz[,,'control_noRinse'])),log='xy')
abline(0,1,lty=2)
dev.off()

tit<-readCounts('ice/out/2019-12-16_isolation/counts.csv')
withAs(tit=tit[!grepl('control|infectivity',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
print(withAs(zz=tit[!grepl('control|infectivity',tit$plate)&tit$n>5&(tit$col!=1|grepl('pbmc',tit$plate)),],zz[order(zz$plate,zz$well),c('plate','well','n')]),row.names=FALSE)


tit<-readCounts('ice/out/2019-12-21_tropism/counts.csv')
tit[tit$plate=='n_p2_2','plate']<-'p2_n_2'
tit[tit$plate=='p2_a2_2','plate']<-'p2_a_2'
vir<-read.csv('ice/2019-12-18_tropism.csv',stringsAsFactors=FALSE)
#vir[vir$Plate==1,3:10]<-rev(unlist(vir[vir$Plate==1,3:10])) #check for plate flip
virus<-data.frame('id'=unlist(vir[,3:10]),'plate'=rep(vir$Plate,8),'row'=rep(colnames(vir)[3:10],each=nrow(vir)),'col'=rep(vir$Row,8),stringsAsFactors=FALSE)
rownames(virus)<-paste(virus$plate,virus$row,virus$col,sep='')
tit$drug<-c('a'='AMD','m'='Mar','am'='AMD+Mar','n'='No Drug')[sapply(strsplit(tit$plate,'_'),'[[',2)]
tit$plateNum<-as.numeric(sub('p','',sapply(strsplit(tit$plate,'_'),'[[',1)))
tit$dilPlate<-as.numeric(sub('B','',sapply(strsplit(tit$plate,'_'),'[[',3)))
tit$dilPlate[tit$plate=='p2_n_3']<-2 #assuming plating mistake
dilutes<-cumprod(c(100/25,rep(275/50,4)))
tit$dil<-dilutes[tit$dilPlate]
tit$well<-paste(tit$plateNum,tit$row,tit$col,sep='')
tit$virus<-virus[tit$well,'id']
dilVir<-virus[virus$plate==2&!grepl('IMC|Media',virus$id)&virus$col %in% 9:10,'id']
tit$dil[tit$virus %in% dilVir]<-tit$dil[tit$virus %in% dilVir]*5
iusStan<-do.call(rbind,cacheOperation('work/tropism_20191223.Rdat',parallel::mclapply,structure(unique(tit$virus),.Names=unique(tit$virus)),function(virus,tit,...){
  sapply(structure(unique(tit$drug),.Names=unique(tit$drug)),function(treat){
    thisDat<-tit[!is.na(tit$virus)&tit$virus==virus&tit$drug==treat,c('n','dil')]
    if(nrow(thisDat)==0)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit,mc.cores=20))
write.csv(iusStan,'out/20191223_tropism.csv')
uniqDrugs<-unique(tit$drug)
drugColors<-structure(rainbow.lab(length(uniqDrugs)),.Names=uniqDrugs)
pdf('out/reboundTropism_20191223.pdf')
for(ii in unique(virus$id[order(grepl('IMC|Media',virus$id),virus$id)])){
  thisDat<-tit[tit$virus==ii,]
  thisDat$logDil<--log(thisDat$dil)
  withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=ii,ylim=range(tit$n+1)*c(1,2),xlim=range(tit$dil),type='n'))
  for(dd in uniqDrugs){
    withAs(zz=thisDat[thisDat$drug==dd,],points(zz$dil,zz$n+1,pch=21,cex=2,bg=drugColors[dd],lwd=3))
  }
  legend('topright',names(drugColors),col=drugColors,lty=1,inset=.01)
  for(dd in uniqDrugs){
    fakeDil<-2^seq(0,16,.1)
    preds<-iusStan[ii,dd]*100/fakeDil
    lines(fakeDil,preds+1,col=drugColors[dd])
  }
  amd<-1-iusStan[ii,'AMD']/iusStan[ii,'No Drug']
  mar<-1-iusStan[ii,'Mar']/iusStan[ii,'No Drug']
  amdMar<-1-iusStan[ii,'AMD+Mar']/iusStan[ii,'No Drug']
  mtext(sprintf('Inhibition:\nAMD: %0.2f%%\nMar: %0.2f%%\nAMD+Mar: %0.2f%%',amd*100,mar*100,amdMar*100),3,at=exp(log(max(tit$dil))*.9))
}
dev.off()

zz<-withAs(tit=tit[grepl('n_3',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))

tit<-readCounts('ice/out/2020-01-06_manolisIsolation/counts.csv')
withAs(tit=tit[!grepl('control|infectivity|red',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))

tit<-readCounts('ice/out/2020-01-13_miniExpansion_isolation/counts.csv')
withAs(tit=tit[!grepl('iso',tit$plate),],tapply(tit$n,list(tit$row,tit$col,tit$plate),c))
tit<-tit[!grepl('iso',tit$plate),]
tit$plateDil<-suppressWarnings(as.numeric(sub('x','',sapply(strsplit(tit$plate,'_'),'[',2))))
tit$dil<-2*tit$plateDil
tit$dil[grepl('control',tit$plate)]<-2*2^(tit$col[grepl('control',tit$plate)]-1)
tit$vir<-paste(sub('_.*','',tit$plate),tit$well)
tit$vir[grepl('control',tit$plate)]<-paste('control',tit$row[grepl('control',tit$plate)])
iusStan<-unlist(cacheOperation('work/mini_20200113.Rdat',parallel::mclapply,structure(unique(tit$vir),.Names=unique(tit$vir)),function(vir,tit,...){
    thisDat<-tit[tit$vir==vir,c('n','dil'),drop=FALSE]
    if(nrow(thisDat)<2)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
},tit=tit,mc.cores=20))
round(tapply(iusStan,list(sub('[0-9]*$','',names(iusStan)),as.numeric(sub('^$','1',sub('[^ ]+ [A-Z]([0-9]*)$','\\1',names(iusStan))))),c))

pdf('out/2020-01-13_miniExpansion.pdf',width=10,height=8)
for(virus in unique(tit$vir[order(grepl('control',tit$plate),tit$vir)])){
    thisDat<-tit[tit$vir==virus,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=virus,ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2,xaxt='n'))
    logAxis(1)
    thisIu<-iusStan[virus]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
}
dev.off()


round(tapply(iusStan,list(sub('[0-9]*$','',names(iusStan)),as.numeric(sub('^$','1',sub('[^ ]+ [A-Z]([0-9]*)$','\\1',names(iusStan))))),c))




yy<-read.csv('ice/2018-04-27_iceSpin/ice_spin/counted/ice_spin.csv',row.names=1)
xx<-readCounts('ice/out/2018-04-27_iceSpin/counts.csv')
xx<-xx[xx$plate=='ice_spin',]
xx<-xx[order(xx$col,xx$row),]
xx$count<-unlist(yy)
plot(xx$count,xx$n,log='xy')
abline(0,1)

yy<-read.csv('ice/2018-05-07_iceTropism/tropism/counted/tropism.csv',row.names=1)
xx<-readCounts('ice/out/2018-05-07_iceTropism/counts.csv')
xx<-xx[xx$plate=='tropism',]
xx<-xx[order(xx$col,xx$row),]
xx$count<-unlist(yy)
pdf('out/counter_vs_script.pdf',width=4,height=4)
par(mgp=c(2.5,.8,0),mar=c(3.4,3.4,.1,.1))
lims<-range(c(xx$count+1,xx$n+1))
plot(xx$count+1,xx$n+1,log='xy',xlab='Count on machine',ylab='Count by script',xaxt='n',yaxt='n',xlim=lims,ylim=lims)
logAxis(las=1,offset=1)
logAxis(1,offset=1)
axis(1,1,0)
axis(2,1,0,las=1)
abline(0,1)
dev.off()

tit<-readCounts('ice/out/2020-02-20_20200203_ifnDepth/counts.csv')
tit<-tit[!grepl('25k_',tit$plate),]
virus<-c('MM33.13 IMC','SG3+WEAU','SL92b','NL43-GFP','YU2','MM33.TF IMC','89.6','Media')
tit$virus<-virus[tit$rowNum]
tit$dil<-(2^(0:11))[tit$col]
iusStan<-do.call(rbind,cacheOperation('work/ifn_20200203.Rdat',parallel::mclapply,structure(unique(tit$virus),.Names=unique(tit$virus)),function(virus,tit,...){
  sapply(structure(unique(tit$plate),.Names=unique(tit$plate)),function(treat){
    thisDat<-tit[!is.na(tit$virus)&tit$virus==virus&tit$plate==treat,c('n','dil')]
    if(nrow(thisDat)==0)return(NA)
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dil,1)
    return(mean(as.matrix(fit)[,'baseIU'])/100)
})},tit,mc.cores=20))

pdf('out/ifnTzmbl.pdf')
par(mfrow=c(3,2))
for(ii in sort(unique(tit$virus))){
  for(jj in c(sort(unique(tit$plate)),'')){
    if(jj==''){plot(1,1,bty='n',xlab='',ylab='',type='n',xaxt='n',yaxt='n');next()}
    thisDat<-tit[tit$virus==ii&tit$plate==jj,]
    thisDat$logDil<--log(thisDat$dil)
    withAs(zz=thisDat,plot(zz$dil,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=paste(ii,jj),ylim=range(tit$n+1),xlim=c(1,max(tit$dil,na.rm=TRUE)),pch=21,bg='blue',cex=2,xaxt='n'))
    logAxis(1)
    thisIu<-iusStan[ii,jj]
    fakeDils=1:50000
    preds<-thisIu/fakeDils*100
    lines(fakeDils,preds+1)
    mtext(sprintf('%0.1f IU/ul',thisIu),3)
}
}
dev.off()


