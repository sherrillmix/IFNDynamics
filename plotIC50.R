library(dnar)
library(vipor)
ic50<-read.csv('IC50s for all subjects.csv',header=FALSE,stringsAsFactors=FALSE)[,-1:-2]
day<-as.numeric(as.vector(ic50[3,]))
patient<-unlist(ic50[1,])
ic50<-apply(ic50[-1:-3,],2,as.numeric)
ic50df<-data.frame('ic50'=as.vector(ic50),'patient'=rep(patient,each=nrow(ic50)),'day'=rep(day,each=nrow(ic50)),stringsAsFactors=FALSE)
ic50df<-ic50df[!is.na(ic50df$ic50),]
ic50df$day2<-ic50df$day^2
cols<-rainbow.lab(length(unique(patient)),alpha=.6)
cols2<-rainbow.lab(length(unique(patient)),alpha=.2)
names(cols)<-names(cols2)<-unique(patient)
pdf('ic50_vs_time.pdf',height=4,width=6)
par(mar=c(3.3,3.3,1.2,.2))
plot(ic50df$day,ic50df$ic50,log='xy',yaxt='n',pch=21,bg=cols[ic50df$patient],cex=1.4,col=cols[ic50df$patient],xaxt='n',xlab='Days after infection',ylab='IFNa2 IC50',mgp=c(2.3,.8,0),xlim=range(c(ic50df$day),2048))
logAxis(2,las=1,mgp=c(3,.8,0))
axis(1,2^(1:12),mgp=c(3,.8,0))
legend('bottomleft',names(cols),pch=21,pt.bg=cols,col=cols,pt.cex=1.4,inset=.02)
for(ii in unique(patient)){
  fit<-lm(I(log(ic50))~day+day2,data=ic50df[ic50df$patient==ii,])
  plot(ic50df$day,ic50df$ic50,log='xy',yaxt='n',xaxt='n',xlab='Days after infection',ylab='IFNa2 IC50',mgp=c(2.3,.8,0),type='n',main=sprintf('%s r^2=%0.3f',ii,summary(fit)$adj.r.squared),xlim=range(c(ic50df$day),2048))
  logAxis(2,las=1,mgp=c(3,.8,0))
  axis(1,2^(1:12),mgp=c(3,.8,0))
  fakeDays<-2^2:2^12
  #pred<-fit$coef['(Intercept)']+fit$coef['day']*fakeDays+fit$coef['I(day^2)']*fakeDays^2
  pred<-exp(predict(fit,data.frame('day'=fakeDays,'day2'=fakeDays^2),interval='prediction'))
  lines(fakeDays,pred[,1],col=cols[ii])
  polygon(c(fakeDays,rev(fakeDays)),c(pred[,2],rev(pred[,3])),col=cols2[ii],border=NA)
  withAs(xx=ic50df[ic50df$patient==ii,],points(xx$day,xx$ic50,pch=21,cex=1.4,col=cols[ii],bg=cols[ii]))
}
dev.off()
