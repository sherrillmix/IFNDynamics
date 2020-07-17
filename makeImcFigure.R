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
pos<-structure(1:length(imcNames),.Names=imcNames)


plotFunc<-function(imcVals,isoVals,ylab='IFNa2 IC50 (pg/ml)',pos,log='y',showLegend=TRUE,showBar=FALSE){
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
  if(showBar){
    ranges<-do.call(rbind,tapply(unlist(imcVals),imcX,range,na.rm=TRUE))
    ranges<-ranges[ranges[,1]!=ranges[,2],]
    segments(as.numeric(rownames(ranges))+.23,ranges[,1],as.numeric(rownames(ranges))+.23,ranges[,2],col='#00000033',lwd=2)
    isSingle<-ave(!is.na(unlist(imcVals)),imcX,FUN=sum)==1
    points(imcX[isSingle]+spreadImc[isSingle]*.3+.23,unlist(imcVals)[isSingle],pch=21,bg='#FF000033')
    #segments(imcX[!isSingle]+.18,unlist(imcVals)[!isSingle],imcX[!isSingle]+.28,unlist(imcVals)[!isSingle],pch=21,col='#FF000033',lwd=2)
    means<-tapply(unlist(imcVals)[!isSingle],imcX[!isSingle],function(xx)exp(mean(log(xx),na.rm=TRUE)))
    print(means)
    points(as.numeric(names(means))+.23,means,pch=21,bg='#FF000033')
    points(imcX[!isSingle]+spreadImc[!isSingle]*.3+.23,unlist(imcVals)[!isSingle],pch=21,bg='#FF000033',col='#000000CC',cex=.4,lwd=.2)
    #segments(imcX[!isSingle]+.18,unlist(imcVals)[!isSingle],imcX[!isSingle]+.28,unlist(imcVals)[!isSingle],pch=21,col='#FF000033',lwd=2)
    #segments(imcX[!isSingle]+.18,unlist(imcVals)[!isSingle],imcX[!isSingle]+.28,unlist(imcVals)[!isSingle],pch=21,col='#FF000033',lwd=2)
  }else{
    points(imcX+spreadImc+.23,unlist(imcVals),pch=21,bg='#FF000033')
  }
  points(isoX+spreadIso-.23,unlist(isoVals),pch=21,bg='#0000FF33')
  abline(v=pos[-1]-.5,col='#00000033')
  if(showLegend)legend('bottomleft',c('Isolate','IMC'),pch=21,pt.bg=c('#0000FF33','#FF000033'),inset=c(-.092,-.39),xpd=NA)
}
imcIsos<-lapply(isoNames,function(xx)combined[grep(xx,combined$virus),])
pdf('out/Fig._S7.pdf',width=8,height=7)
par(mar=c(5.1,3.5,.1,.1),mfrow=c(2,1))
plotFunc(by(standardized[,'ic50.IFNa2'],standardized$virus,function(xx)xx)[imcNames],lapply(imcIsos,function(xx)xx[,'ic50_IFNa2'])[imcNames],pos=pos,showLegend=FALSE,showBar=TRUE)
text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.001,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'A',xpd=NA,adj=c(0,1),cex=2)
plotFunc(by(standardized[,'ic50.IFNb'],standardized$virus,function(xx)xx)[imcNames],lapply(imcIsos,function(xx)xx[,'ic50_IFNb'])[imcNames],pos=pos,ylab='IFNb IC50 (pg/ml)',showBar=TRUE)
text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.001,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'B',xpd=NA,adj=c(0,1),cex=2)
dev.off()


#infectivity
#read.csv('out/20191223_tropism.csv')
#read.csv('out/macInf_20200226.csv')
#out/macIsoInf_20200210.csv:1
#out/macIsoInf_20200220.csv:1
#out/macInf_20200226.csv


imcIsos<-lapply(isoNames,function(xx)combined[grep(xx,combined$virus),])
for(ii in names(isoNames)[grep('MM33',names(isoNames))]){
  thisDat<-long[grep(isoNames[ii],long$id),]
  ifna2<-na.omit(unlist(thisDat[,grep('IFNa2.*IC50',colnames(thisDat))]))
  ifnb<-na.omit(unlist(thisDat[,grep('IFNb.*IC50',colnames(thisDat))]))
  repCap<-na.omit(unlist(thisDat[,grep('replication',colnames(thisDat))]))
  n<-max(c(length(ifna2),length(ifnb),length(repCap)))
  imcIsos[[ii]]<-data.frame('ic50_IFNa2'=c(ifna2,rep(NA,n-length(ifna2))),'ic50_IFNb'=c(ifnb,rep(NA,n-length(ifnb))),'repCap'=c(repCap,rep(NA,n-length(repCap))))
}
pdf('out/imc_vs_isolates.pdf',width=8,height=4)
par(mar=c(5.1,3.5,.1,.1))
plotFunc(by(standardized[,'ic50.IFNa2'],standardized$virus,function(xx)xx)[imcNames],lapply(imcIsos,function(xx)xx[,'ic50_IFNa2'])[imcNames],pos=pos)
plotFunc(by(standardized[,'ic50.IFNb'],standardized$virus,function(xx)xx)[imcNames],lapply(imcIsos,function(xx)xx[,'ic50_IFNb'])[imcNames],pos=pos,ylab='IFNb IC50 (pg/ml)')
plotFunc(by(standardized[,c('repCap.IFNa2','repCap.IFNb')],standardized$virus,function(xx)xx)[imcNames],lapply(imcIsos,function(xx)xx[,'repCap'])[imcNames],pos=pos,log='',ylab='Replicative capacity (p24 ng/ml)')
dev.off()
#system('pdftk out/imc_vs_isolates.pdf cat 1 2 output tmp.pdf; pdfjam tmp.pdf --nup 1x2 --outfile out/Fig._S7.pdf')


tmp<-long
tmp$study<-'MM'
tmp$ic50_IFNa2<-tmp$ic50
tmp$ic50_IFNb<-tmp$beta
tmp$class<-'Long'
allIc50<-rbind(combined[!combined$study %in% c('MM','Transmission'),c('study','pat','ic50_IFNa2','ic50_IFNb','class')],tmp[,c('study','pat','ic50_IFNa2','ic50_IFNb','class')])
allIc50<-allIc50[!(is.na(combined$ic50_IFNb)&is.na(combined$ic50_IFNa2)),]
allIc50$display<-ifelse(allIc50$study=='MM',allIc50$pat,ifelse(allIc50$class=='Rebound','Rebound','Outgrowth'))
#cols<-structure(c('#8dd3c7','#ffffb3','#bebada','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#fb8072'),.Names=c(sort(unique(tmp$pat)),"Rebound","Outgrowth"))
cols<-structure(c(sprintf('%saa',c(
        "#af7890",
        "#7262cf",
        "#afb441",
        "#4eb491",
        "#d24f37",
        "#d59948",
        "#7181ca",
        "#62813c",
        "#3e9cb7",
        "#9c6631"
        )),'#9EC0E1ee','#E581A0ee'),.Names=c(sort(unique(tmp$pat)),"Outgrowth","Rebound"))
pdf('out/ifna2_vs_ifnb.pdf',width=5,height=6)
  par(mar=c(8,3.4,.1,.1))
  plot(1,1,type='n',xaxt='n',yaxt='n',xlim=range(allIc50$ic50_IFNa2,na.rm=TRUE),ylim=range(allIc50$ic50_IFNb,na.rm=TRUE),log='xy',ylab='IFNb IC50 (pg/ml)',mgp=c(2.5,1,0),xlab='',cex=1.3)
  title(xlab='IFNa2 IC50 (pg/ml)',mgp=c(2,1,0))
  dnar::logAxis(1,mgp=c(3,.7,0))
  dnar::logAxis(las=1,mgp=c(3,.7,0))
  dnar::withAs(allIc50=allIc50[allIc50$class=='Long',],points(allIc50$ic50_IFNa2,allIc50$ic50_IFNb,pch=21,bg=cols[allIc50$display]))
  dnar::withAs(allIc50=allIc50[allIc50$class!='Long',],points(allIc50$ic50_IFNa2,allIc50$ic50_IFNb,pch=21,bg=cols[allIc50$display]))
  legend('bottom',names(cols),pch=21,pt.bg=cols,pt.cex=1.3,ncol=4,inset=-.35,xpd=NA,x.intersp=.9,text.width=NULL)
dev.off()

