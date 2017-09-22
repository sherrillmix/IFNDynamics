library(dnar)
library(vipor)
readHeadedCsv<-function(fileName,hasDays=TRUE){
  ic50<-read.csv(fileName,header=FALSE,stringsAsFactors=FALSE)[,-1]
  labs<-ic50[,1]
  ic50<-ic50[,-1]
  if(hasDays)day<-as.numeric(as.vector(ic50[3,]))
  else day<-as.numeric(as.vector(ic50[2,]))*7
  patient<-unlist(ic50[1,])
  ic50<-apply(ic50[-1:-ifelse(hasDays,3,2),],2,as.numeric)
  ic50df<-data.frame('ic50'=as.vector(ic50),'patient'=rep(patient,each=nrow(ic50)),'day'=rep(day,each=nrow(ic50)),stringsAsFactors=FALSE)
  ic50df<-ic50df[!is.na(ic50df$ic50),]
  ic50df$day2<-ic50df$day^2
  return(ic50df) 
}

plotIc50<-function(ic50,ylab='IFNa2 IC50',outFile='out.pdf',fitAll=TRUE,legendPos='topright',isLog=TRUE,excludeFromAll=NULL){
  cols<-rainbow.lab(length(unique(ic50$patient)),alpha=.6)
  cols2<-rainbow.lab(length(unique(ic50$patient)),alpha=.2)
  names(cols)<-names(cols2)<-sort(unique(ic50$patient))
  fits<-list()
  ic50All<-ic50[!ic50$patient %in% excludeFromAll,]
  pdf(outFile,height=4,width=6)
    par(mar=c(3.3,3.3,1.2,.2))
    fit<-lm(I(log(ic50))~day+day2,data=ic50All)
    fits[['all']]<-fit
    fakeDays<-2^2:2^12
    pred<-exp(predict(fit,data.frame('day'=fakeDays,'day2'=fakeDays^2),interval='prediction'))
    plot(ic50All$day,ic50All$ic50,log=ifelse(isLog,'y',''),yaxt=ifelse(isLog,'n','s'),pch=21,bg=cols[ic50All$patient],cex=1.4,col=cols[ic50All$patient],xlab='Days after infection',ylab=ylab,mgp=c(2.3,.8,0),las=1)
    if(isLog)logAxis(2,las=1,mgp=c(3,.8,0))
    if(fitAll){
      lines(fakeDays,pred[,1],col='#00000022')
      polygon(c(fakeDays,rev(fakeDays)),c(pred[,2],rev(pred[,3])),col='#00000011',border=NA)
      title(main=sprintf('All r^2=%0.3f',summary(fit)$adj.r.squared))
    }
    legend(legendPos,names(cols)[names(cols)%in%ic50All$patient],pch=21,pt.bg=cols,col=cols[names(cols)%in%ic50All$patient],pt.cex=1.4,inset=.02)
    for(ii in sort(unique(ic50$patient))){
      message(ii)
      fit<-lm(I(log(ic50))~day+day2,data=ic50[ic50$patient==ii,])
      fits[[ii]]<-fit
      plot(ic50$day,ic50$ic50,log=ifelse(isLog,'y',''),yaxt=ifelse(isLog,'n','s'),xlab='Days after infection',ylab=ylab,mgp=c(2.3,.8,0),type='n',main=sprintf('%s r^2=%0.3f',ii,summary(fit)$adj.r.squared),las=1)
      if(isLog)logAxis(2,las=1,mgp=c(3,.8,0))
      pred<-exp(predict(fit,data.frame('day'=fakeDays,'day2'=fakeDays^2),interval='prediction'))
      lines(fakeDays,pred[,1],col=cols[ii])
      polygon(c(fakeDays,rev(fakeDays)),c(pred[,2],rev(pred[,3])),col=cols2[ii],border=NA)
      withAs(xx=ic50[ic50$patient==ii,],points(xx$day,xx$ic50,pch=21,cex=1.4,col=cols[ii],bg=cols[ii]))
    }
  dev.off()
  return(invisible(fits))
}


alphaIc50<-readHeadedCsv('data/IC50s for all subjects_2_alpha.csv')
alphaIc50[alphaIc50$patient=='MM33','patient']<-'MM33.2'
alphaIc50<-rbind(
  alphaIc50,
  readHeadedCsv('data/IC50_patient5_alpha.csv',FALSE),
  readHeadedCsv('data/IC50_patient33_redo_alpha.csv',FALSE)
)

alphaVres<-readHeadedCsv('data/IC50s for all subjects_2_alphaVres.csv')
alphaVres[alphaVres$patient=='MM33','patient']<-'MM33.2'
alphaVres<-rbind(
  alphaVres,
  readHeadedCsv('data/IC50_patient5_alphaVres.csv',FALSE),
  readHeadedCsv('data/IC50_patient33_redo_alphaVres.csv',FALSE)
)
betaIc50<-readHeadedCsv('data/IC50s for all subjects_2_beta.csv')
betaVres<-readHeadedCsv('data/IC50s for all subjects_2_betaVres.csv')

plotIc50(alphaIc50,'IFNa2 IC50','out/ifna2_ic50.pdf',legendPos='top',excludeFromAll='MM33.2')
plotIc50(alphaVres,'IFNa2 Vres','out/ifna2_vres.pdf',isLog=FALSE,legendPos='top',excludeFromAll='MM33.2')
plotIc50(betaIc50,'IFNb IC50','out/ifnb_ic50.pdf')
plotIc50(betaVres,'IFNb Vres','out/ifnb_vres.pdf',isLog=FALSE)

