source('readReboundData.R')

long<-read.csv('out/allLongitudinal.csv',row.names=1,stringsAsFactors=FALSE)
imcNames<-c('MM33.TF'='MM33.TF','MM33.13'='MM33.13','MM33.17'='MM33.17','601r1'='601.REB.r1','A09r-1A2'='A09.REB.1A2','A08 P1A5'='A08.REB.1A5','A08.21-2P4'='A08.REB.2F4','9201reb'='9242.REB.r1')

imcFiles<-c(
  'out/fred_20200224_ic50.csv',
  'out/imc_20191105_ic50.csv',
  'out/imc_20191115_ic50.csv',
  'out/zapRedo_20190712.csv',
  'out/a08_20200424_ic50.csv',
  'out/man_20200424_ic50.csv'
)

imcIc50s<-lapply(imcFiles,read.csv,stringsAsFactors=FALSE,row.names=1)
imcIc50s<-lapply(imcIc50s,function(xx)xx[!grepl('13C1|2A4|6B5|8B5|1E2',rownames(xx))&!grepl('^(Reb|QVOA)\\.',rownames(xx))&grepl('MM33|601|9201|A08|A09',rownames(xx)),])
desiredCols<-c('ic50.IFNa2','ic50.IFNb','repCap.IFNa2','repCap.IFNb','vres.IFNa2','vres.IFNb')
standardized<-do.call(rbind,mapply(function(xx,yy){
  if(yy=='out/imc_20191105_ic50.csv'){
    colnames(xx)<-sprintf('%s.IFNb',colnames(xx))
    xx<-xx[grepl('50$',rownames(xx)),]
  }
  if(yy=='out/zapRedo_20190712.csv')colnames(xx)[colnames(xx)=='IFNa2_repCap.1']<-'IFNb_repCap'
  colnames(xx)<-sub('Vres','vres',sub('IC50','ic50',sub('(IFN[^_.]+)[_.](.*)','\\2.\\1',colnames(xx))))
  xx<-xx[!grepl('MM33.17.2A4|MM33.1.13C1',rownames(xx)),]
  xx$file<-yy
  xx$virus<-imcNames[trimws(sub('[ _]+IMC.*','',sub('^.*MM33[ ._](TF|13|17).*$','MM33.\\1',rownames(xx))))]
  xx[,desiredCols[!desiredCols %in% colnames(xx)]]<-NA
  xx[,c(desiredCols,'virus','file')]
},imcIc50s,imcFiles,SIMPLIFY=FALSE))
standardized<-standardized[order(standardized$virus,standardized$file),]

isoNames<-c('MM33.TF'='MM33.01','MM33.13'='MM33.13.2D6','MM33.17'='MM33.17.2A4',"601.REB.r1"='601.*4(A7|A8|C1|B4)',"A09.REB.1A2"='A09.REB.1A2',"A08.REB.1A5"='A08.REB.1A5',"A08.REB.2F4"='A08.REB.2F4',"9242.REB.r1"='601.REB.(4B4|4A8|4A7|4C1)')[imcNames]
imcIsos<-lapply(isoNames,function(xx)combined[grep(xx,combined$virus),])
for(ii in names(isoNames)[grep('MM33',names(isoNames))]){
  thisDat<-long[grep(isoNames[ii],long$id),]
  ifna2<-na.omit(unlist(thisDat[,grep('IFNa2.*IC50',colnames(thisDat))]))
  ifnb<-na.omit(unlist(thisDat[,grep('IFNb.*IC50',colnames(thisDat))]))
  repCap<-na.omit(unlist(thisDat[,grep('replication',colnames(thisDat))]))
  n<-max(c(length(ifna2),length(ifnb),length(repCap)))
  imcIsos[[ii]]<-data.frame('ic50_IFNa2'=c(ifna2,rep(NA,n-length(ifna2))),'ic50_IFNb'=c(ifnb,rep(NA,n-length(ifnb))),'repCap'=c(repCap,rep(NA,n-length(repCap))))
}



pos<-structure(1:length(imcNames),.Names=imcNames)
pdf('out/imc_vs_isolates.pdf',width=8,height=4)
plotFunc<-function(imcVals,isoVals,ylab='IFNa2 IC50 (pg/ml)',pos,log='y'){
  plot(1,1,type='n',xlab='',ylab=ylab,xlim=range(pos)+c(-.5,.5),xaxs='i',ylim=range(c(unlist(imcVals),unlist(isoVals)),na.rm=TRUE),xaxt='n',log=log,yaxt=ifelse(log=='','s','n'),mgp=c(2.5,.5,0),las=1,tcl=-.3)
  if(log=='y')dnar::logAxis(las=1,mgp=c(3,.6,0))
  #dnar::slantAxis(1,pos,names(pos))
  for(ii in 1:0)axis(1,pos[1:length(pos)%%2==ii],names(pos)[1:length(pos)%%2==ii],mgp=c(3,1.9+ii*-1.5,0),tcl=-1.8+ii*1.5)
  imcX<-rep(pos,sapply(imcVals,function(xx)length(unlist(xx))))
  isoX<-rep(pos,sapply(isoVals,function(xx)length(unlist(xx))))
  #spreadImc<-ave(unlist(imcVals),imcX,FUN=function(xx)beeswarm::swarmx(rep(0,length(xx)),xx,cex=.45)$x)
  #spreadIso<-ave(unlist(isoVals),isoX,FUN=function(xx)beeswarm::swarmx(rep(0,length(xx)),xx,cex=.45)$x)
  spreadImc<-ave(unlist(imcVals),imcX,FUN=function(xx){zz<-rep(NA,length(xx));zz[!is.na(xx)]<-vipor::offsetX(xx[!is.na(xx)],width=.2);zz})
  spreadIso<-ave(unlist(isoVals),isoX,FUN=function(xx){
    if(sum(!is.na(xx))<10){
      beeswarm::swarmx(rep(0,length(xx)),xx)$x
    }else{
      zz<-rep(NA,length(xx));zz[!is.na(xx)]<-vipor::offsetX(xx[!is.na(xx)],width=.2);zz
  }})
  points(imcX+spreadImc+.23,unlist(imcVals),pch=21,bg='#FF000033')
  points(isoX+spreadIso-.23,unlist(isoVals),pch=21,bg='#0000FF33')
  abline(v=pos[-1]-.5,col='#00000033')
  legend('bottomleft',c('Isolate','IMC'),pch=21,pt.bg=c('#0000FF33','#FF000033'),inset=c(-.092,-.335),xpd=NA)
}
par(mar=c(5.1,3.5,.1,.1))
plotFunc(by(standardized[,'ic50.IFNa2'],standardized$virus,function(xx)xx)[imcNames],lapply(imcIsos,function(xx)xx[,'ic50_IFNa2'])[imcNames],pos=pos)
plotFunc(by(standardized[,'ic50.IFNb'],standardized$virus,function(xx)xx)[imcNames],lapply(imcIsos,function(xx)xx[,'ic50_IFNb'])[imcNames],pos=pos,ylab='IFNb IC50 (pg/ml)')
plotFunc(by(standardized[,c('repCap.IFNa2','repCap.IFNb')],standardized$virus,function(xx)xx)[imcNames],lapply(imcIsos,function(xx)xx[,'repCap'])[imcNames],pos=pos,log='',ylab='Replicative capacity (p24 ng/ml)')
dev.off()

#infectivity
#read.csv('out/20191223_tropism.csv')
#read.csv('out/macInf_20200226.csv')
#out/macIsoInf_20200210.csv:1
#out/macIsoInf_20200220.csv:1
#out/macInf_20200226.csv

