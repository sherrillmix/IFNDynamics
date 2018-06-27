library(dnar)
if(!exists('dat'))source('readNewData.R')
readCounts<-function(countFile){
  counts<-read.csv(countFile,stringsAsFactors=FALSE)
  counts$plate<-basename(counts$dir)
  counts$well<-sub('\\.CTL$','',counts$file)
  counts$col<-as.numeric(sub('^[A-Z]','',counts$well))
  counts$row<-sub('[0-9]+$','',counts$well)
  return(counts)
}
readPlateViruses<-function(plateFile,virusFile){
  plate<-read.csv(plateFile,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
  viruses<-read.csv(virusFile,stringsAsFactors=FALSE)
  plate<-apply(plate,2,trimws)
  viruses<-rbind(viruses,c('Empty','Empty'))
  rownames(viruses)<-viruses$id
  if(any(!unlist(plate) %in% viruses$id))warning('Unknown virus ',unique(unlist(plate)[!unlist(plate) %in% viruses$id]))
  plateIds<-data.frame('row'=rep(rownames(plate),ncol(plate)),'col'=rep(colnames(plate),each=nrow(plate)),'vId'=as.vector(unlist(plate)),stringsAsFactors=FALSE)
  plateIds$well<-sprintf('%s%s',plateIds$row,plateIds$col)
  rownames(plateIds)<-plateIds$well
  plateIds$virus<-viruses[plateIds$vId,'sample']
  return(plateIds)
}

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



