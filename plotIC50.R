library(dnar)
library(vipor)
readHeadedCsv<-function(fileName){
  ic50<-read.csv(fileName,header=FALSE,stringsAsFactors=FALSE)[,-1]
  labs<-ic50[,1]
  ic50<-ic50[,-1]
  day<-as.numeric(as.vector(ic50[3,]))
  patient<-unlist(ic50[1,])
  ic50<-apply(ic50[-1:-3,],2,as.numeric)
  ic50df<-data.frame('ic50'=as.vector(ic50),'patient'=rep(patient,each=nrow(ic50)),'day'=rep(day,each=nrow(ic50)),stringsAsFactors=FALSE)
  ic50df<-ic50df[!is.na(ic50df$ic50),]
  ic50df$day2<-ic50df$day^2
  return(ic50df) 
}

plotIc50<-function(ic50,ylab='IFNa2 IC50',outFile='out.pdf',fitAll=FALSE,legendPos='topright',isLog=TRUE){
  cols<-rainbow.lab(length(unique(ic50$patient)),alpha=.6)
  cols2<-rainbow.lab(length(unique(ic50$patient)),alpha=.2)
  names(cols)<-names(cols2)<-unique(ic50$patient)
  pdf(outFile,height=4,width=6)
    par(mar=c(3.3,3.3,1.2,.2))
    fit<-lm(I(log(ic50))~day+day2,data=ic50)
    fakeDays<-2^2:2^12
    pred<-exp(predict(fit,data.frame('day'=fakeDays,'day2'=fakeDays^2),interval='prediction'))
    plot(ic50$day,ic50$ic50,log=ifelse(isLog,'y',''),yaxt=ifelse(isLog,'n','s'),pch=21,bg=cols[ic50$patient],cex=1.4,col=cols[ic50$patient],xlab='Days after infection',ylab=ylab,mgp=c(2.3,.8,0))
    if(isLog)logAxis(2,las=1,mgp=c(3,.8,0))
    if(fitAll){
      lines(fakeDays,pred[,1],col='#00000022')
      polygon(c(fakeDays,rev(fakeDays)),c(pred[,2],rev(pred[,3])),col='#00000011',border=NA)
      title(main=sprintf('All r^2=%0.3f',summary(fit)$adj.r.squared))
    }
    legend(legendPos,names(cols),pch=21,pt.bg=cols,col=cols,pt.cex=1.4,inset=.02)
    for(ii in unique(ic50$patient)){
      fit<-lm(I(log(ic50))~day+day2,data=ic50[ic50$patient==ii,])
      plot(ic50$day,ic50$ic50,log=ifelse(isLog,'y',''),yaxt=ifelse(isLog,'n','s'),xlab='Days after infection',ylab=ylab,mgp=c(2.3,.8,0),type='n',main=sprintf('%s r^2=%0.3f',ii,summary(fit)$adj.r.squared))
      if(isLog)logAxis(2,las=1,mgp=c(3,.8,0))
      pred<-exp(predict(fit,data.frame('day'=fakeDays,'day2'=fakeDays^2),interval='prediction'))
      lines(fakeDays,pred[,1],col=cols[ii])
      polygon(c(fakeDays,rev(fakeDays)),c(pred[,2],rev(pred[,3])),col=cols2[ii],border=NA)
      withAs(xx=ic50[ic50$patient==ii,],points(xx$day,xx$ic50,pch=21,cex=1.4,col=cols[ii],bg=cols[ii]))
    }
  dev.off()
}

alphaIc50<-readHeadedCsv('data/IC50s for all subjects_2_alpha.csv')
alphaVres<-readHeadedCsv('data/IC50s for all subjects_2_alphaVres.csv')
betaIc50<-readHeadedCsv('data/IC50s for all subjects_2_beta.csv')
betaVres<-readHeadedCsv('data/IC50s for all subjects_2_betaVres.csv')

plotIc50(alphaIc50,'IFNa2 IC50','out/ifna2_ic50.pdf')
plotIc50(alphaVres,'IFNa2 Vres','out/ifna2_vres.pdf',isLog=FALSE)
plotIc50(betaIc50,'IFNb IC50','out/ifnb_ic50.pdf')
plotIc50(betaVres,'IFNb Vres','out/ifnb_vres.pdf',isLog=FALSE)

