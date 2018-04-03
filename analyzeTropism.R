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
iceProps<-t(apply(iceMeans,1,function(xx)xx/xx[1]))

icePats<-unique(na.omit(dat[rownames(iceMeans),'pat']))
times<-as.numeric(sub('h','',colnames(iceMeans)))
pdf('out/iceMeans.pdf')
  plot(1,1,type='n',xlab='Hours on ice',ylab='Proportion of 0h TZM-BL infectivity',xlim=range(times),ylim=range(iceProps))
  apply(iceProps,1,function(xx)lines(times,xx))
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
  plot(dat[rownames(iceMeans),'time'],iceProps[,'48h'],ylab='Proportion TZM-BL infectivity after 8h on ice',xlab='Days after onset of symptoms',pch=21,bg=patCols[dat[rownames(iceMeans),'pat']],las=1,cex=1.4,mgp=c(2.5,.7,0))
  axis(4,iceProps[!rownames(iceMeans) %in% rownames(dat),'48h'],sub('CoJRFL_','',rownames(iceProps)[!rownames(iceMeans) %in% rownames(dat)]),las=1)
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
tapply(iceTrop$n,list(iceTrop$row,iceTrop$col,iceTrop$plate),c)

iceLuc<-read.csv('ice/2018-04-02-iceTropism/2018-04-02_ice.csv',check.names=FALSE,row.names=1)[,-13]
tropLuc<-read.csv('ice/2018-04-02-iceTropism/2018-04-02_tropism.csv',check.names=FALSE,row.names=1)[,-13]
