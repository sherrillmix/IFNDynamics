manOut<-read.csv('out/man_20200424_ic50.csv',stringsAsFactors=FALSE)
oldMan<-read.csv('out/reboundVoaData.csv',stringsAsFactors=FALSE)
oldMan<-oldMan[oldMan$pat=='S22'&oldMan$Type1=='Rebound',]
oldMan$day<-0
oldMan$ic50.IFNa2<-oldMan$ic50_IFNa2
oldMan$ic50.IFNb<-oldMan$ic50_IFNb
manOut<-manOut[!is.na(manOut$pat),c('pat','day','ic50.IFNa2','ic50.IFNb')]
manOut<-rbind(manOut,oldMan[,colnames(manOut)])

cols<-sprintf('%s99',c('1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')[-1])
patCols<-structure(cols[1:length(unique(manOut$pat))],.Names=unique(manOut$pat))
patZeroPos<-structure((1:length(unique(manOut$pat)))*4-8,.Names=unique(manOut$pat))
pdf('out/longAti_Ic50.pdf',width=8,height=3.5)
par(mfrow=c(1,2),mar=c(3.4,3.6,.1,.8))
for(ifn in c('ic50.IFNa2','ic50.IFNb')){
  plot(ifelse(manOut$day==0,patZeroPos[manOut$pat],manOut$day),manOut[,ifn],ylim=range(manOut[,ifn]),xlim=range(manOut$day[!is.na(manOut$day)]),log='y',yaxt='n',ylab=sub('ic50\\.','IC50 ',ifn),pch=21,xlab='Days after detectable rebound',bg=patCols[manOut$pat],cex=1.2,mgp=c(2.3,.8,0))
  dnar::logAxis(las=1)
  legend('bottomleft',names(patCols),inset=0,pch=21,pt.bg=patCols,bty='n')
}
dev.off()


