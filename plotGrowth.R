library(dnar)
library(vipor)
marvinQVOAs<-c('MM23\\.18|MM34\\.25')
grow<-read.csv('data/Plasma vs QVOA for Marvin.csv',stringsAsFactors=FALSE)
grow$study<-ifelse(grepl(marvinQVOAs,grow$SampleID1),'London Cohort after ART',
  ifelse(grepl('^MM',grow$SampleID1),'London Cohort',
  ifelse(grepl('^CH',grow$SampleID1),'Transmission Pairs',
  ifelse(grepl('^B[0-9]{3}_',grow$SampleID1),'QVOA Nussenzweig',
  ifelse(grepl('^[0-9]{3}[-_]',grow$SampleID1),'QVOA Seaman',
  ifelse(grepl('^[0-9]{6}',grow$SampleID1),'QVOA Bar/Siliciano',
      NA))))))
if(any(is.na(grow$study)))stop('Problem identifying study')
tmp<-colnames(grow)[grep('[dD]ay',colnames(grow))]
dayCols<-tapply(tmp,sub('\\.\\.ng.*','',sub('day','Day',sub('([0-9])\\.[0-9]$','\\1',tmp))),c)
growMean<-do.call(cbind,lapply(dayCols,function(xx)apply(grow[,xx],1,mean)))
growMin<-do.call(cbind,lapply(dayCols,function(xx)apply(grow[,xx],1,max)))
growMax<-do.call(cbind,lapply(dayCols,function(xx)apply(grow[,xx],1,min)))
days<-as.numeric(sub('Day\\.','',colnames(growMean)))
studyCol<-rainbow.lab(length(unique(grow$study)))
studyCol<-studyCol[order(vanDerCorput(length(studyCol)))]
names(studyCol)<-unique(grow$study)
pdf('out/isolates_vs_qvoa.pdf',width=9,height=6)
  par(mfrow=c(2,3))
  for(ii in unique(grow$study)){
    selector<-grow$study==ii
    plot(1,1,type='n',xlab='Day',xlim=c(1,7),ylim=range(growMean),ylab='p24 (ng/ul)',main=ii,las=1)
    apply(growMean[selector,],1,function(xx)lines(days,xx,col='#00000033',lwd=2))
  }
dev.off()
pdf('out/isolates_vs_qvoa_day7.pdf',width=10,height=6,useDingbats=FALSE)
  par(mar=c(7,3.5,.1,.1))
  plot(1,1,type='n',xlab='',xlim=c(1,nrow(grow))+c(-1.5,1.5),ylim=range(cbind(growMin,growMax)),ylab='p24 (ng/ul)',las=1,xaxt='n',mgp=c(2.5,.9,0),xaxs='i')
  axis(1,1:nrow(grow),grow$SampleID1,las=2,cex.axis=.6,mgp=c(3,.6,0))
  segments(1:nrow(grow),growMin[,'Day.7'],1:nrow(grow),growMax[,'Day.7'],col=studyCol[grow$study],lwd=2)
  points(rep(1:nrow(grow),length(dayCols[['Day.7']])),unlist(grow[,dayCols[['Day.7']]]),col=studyCol[grow$study],pch=c(1,1,2,2),cex=.7)
  points(1:nrow(grow),growMean[,'Day.7'],pch=21,bg=studyCol[grow$study],cex=1.5)
  legend('top',names(studyCol),pch=21,pt.bg=studyCol,ncol=3,pt.cex=1.7,inset=0,bty='n')
dev.off()
