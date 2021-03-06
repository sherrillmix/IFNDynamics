
rt<-read.csv('data/RTs 01.05.2018.csv',skip=1,stringsAsFactors=FALSE)
rt$old<-suppressWarnings(as.numeric(rt$old.RT)
rt$new<-suppressWarnings(as.numeric(rt$New.RT))

rt2<-read.csv('data/20190109_virusSubset.csv',stringsAsFactors=FALSE)
rt2<-rt2[!is.na(rt2$number),]
rt2$old<-as.numeric(rt2$old.RT)
rt2$new<-as.numeric(rt2$New.RT)

lims<-range(c(rt$old,rt$new,rt2$old,rt2$new),na.rm=TRUE)
pdf('out/oldNewRT.pdf',width=4,height=4)
  par(mar=c(3.75,3.75,.2,.2))
  dnar::withAs(rt=rt[grep('MM',rt$X.2),],plot(rt$old,rt$new,xlab='Old RT',ylab='New RT',log='xy',xlim=lims,ylim=lims,yaxt='n',xaxt='n',mgp=c(2.75,1,0)))
  dnar::logAxis(las=1)
  dnar::logAxis(1)
  points(rt2$old,rt2$new,col='blue')
  abline(0,1,lty=2)
  fit<-lm(new~0+old,dat=log(rbind(rt[,c('old','new')],rt2[,c('old','new')])))
  fakeDat<-data.frame(old=seq(log(lims[1]),log(lims[2]),length.out=1000))
  pred<-predict(fit,fakeDat,interval='conf')
  pred2<-predict(fit,fakeDat,interval='pred')
  lines(exp(fakeDat[,'old']),exp(pred[,'fit']))
  polygon(exp(c(fakeDat$old,rev(fakeDat$old))),exp(c(pred[,'lwr'],rev(pred[,'upr']))),col='#00000022',border=NA)
  polygon(exp(c(fakeDat$old,rev(fakeDat$old))),exp(c(pred2[,'lwr'],rev(pred2[,'upr']))),col='#00000011',border=NA)
dev.off()
