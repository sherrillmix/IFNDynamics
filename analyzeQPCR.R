library(dnar)
multi<-read.csv('ice/2018-07-06_sybrRT/Mulicomponent.csv',skip=47,stringsAsFactors=FALSE)
multi$rn<-as.numeric(sub(',','',multi$FAM))/as.numeric(sub(',','',multi$ROX))
multi$well<-trimws(multi$Well.Position)
multi$standard<-ave(multi$rn,multi$well,FUN=function(xx)xx-xx[1])

#amp<-read.csv('ice/2018-07-06_sybrRT/Amplification.csv',stringsAsFactors=FALSE,skip=47)

results<-read.csv('ice/2018-07-06_sybrRT/Results.csv',stringsAsFactors=FALSE,skip=47)
results$well<-trimws(results$Well.Position)
results$col<-as.numeric(gsub('[^0-9]','',results$well))
results$ct<-as.numeric(ifelse(results$CT=='Undetermined',0,results$CT))

pdf('out/qpcr.pdf',width=10,height=8)
par(mfrow=c(2,4))
maxCycle<-50
melt<-multi[multi$Cycle>maxCycle,]
multi<-multi[multi$Cycle<=maxCycle,]
xlim<-c(1,maxCycle)
#ylim<-range(multi$rn)
ylim<-c(0,1.75)
rowDef<-c('A'='ProtoScript','B'='SuperScript','C'='89.6','D'='SG3','E'='H20','F'='GAPDH HeLa DNA','G'='GAPDH RNA no RT','H'='GAPDH H20')
for(ii in LETTERS[1:8]){
  thisDat<-multi[grepl(sprintf('^[ ]*%s',ii),multi$well),]
  plot(1,1,xlim=xlim,ylim=ylim,type='n',xlab='Cycle',ylab='Rn',las=1,main=rowDef[gsub('[^A-H]','',ii)])
  cols<-structure(rev(rainbow.lab(length(unique(thisDat$well)))),.Names=unique(thisDat$well))
  for(jj in unique(thisDat$well)){
    withAs(xx=thisDat[thisDat$well==jj,],lines(xx$Cycle,xx$standard,col=cols[jj]))
  }
}
for(ii in LETTERS[1:8]){
  thisDat<-results[grepl(sprintf('^[ ]*%s',ii),results$well),]
  divisor<-ifelse(ii %in% LETTERS[6:8],5,3)
  plot(divisor^(thisDat$col-1),thisDat$ct,xlab='Dilution',ylab='Ct',las=1,main=rowDef[gsub('[^A-H]','',ii)],log='x',ylim=c(0,50),xaxt='n')
  logAxis(1)
}
dev.off()

