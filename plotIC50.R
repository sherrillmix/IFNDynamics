library(dnar)
library(vipor)
ic50<-read.csv('IC50s for all subjects.csv',header=FALSE,stringsAsFactors=FALSE)[,-1:-2]
day<-as.numeric(as.vector(ic50[3,]))
patient<-unlist(ic50[1,])
ic50<-apply(ic50[-1:-3,],2,as.numeric)
ic50df<-data.frame('ic50'=as.vector(ic50),'patient'=rep(patient,each=nrow(ic50)),'day'=rep(day,each=nrow(ic50)),stringsAsFactors=FALSE)
ic50df<-ic50df[!is.na(ic50df$ic50),]
cols=rainbow.lab(length(unique(patient)),alpha=.6)
names(cols)<-unique(patient)
pdf('ic50_vs_time.pdf',height=4,width=6)
par(mar=c(3.5,3.5,.2,.2))
plot(ic50df$day,ic50df$ic50,log='xy',yaxt='n',pch=21,bg=cols[ic50df$patient],cex=1.4,col=cols[ic50df$patient],xaxt='n',xlab='Days after infection',ylab='IFNa2 IC50',mgp=c(2.5,.8,0))
logAxis(2,las=1)
axis(1,2^(1:10))
legend('bottomleft',names(cols),pch=21,pt.bg=cols,col=cols,pt.cex=1.4,inset=.02)
par(mar=c(3.5,3.5,1.2,.2))
for(ii in unique(patient)){
  plot(ic50df$day,ic50df$ic50,log='xy',yaxt='n',xaxt='n',xlab='Days after infection',ylab='IFNa2 IC50',mgp=c(2.5,.8,0),type='n',main=ii)
  withAs(xx=ic50df[ic50df$patient==ii,],points(xx$day,xx$ic50,pch=21,cex=1.4,col=cols[ii],bg=cols[ii]))
  logAxis(2,las=1)
  axis(1,2^(1:10))
}
dev.off()
