wells<-colnames(trans)[grep('well',colnames(trans))]
measures<-sub('\\..*','',wells)
trans<-trans[order(trans$patient,trans$time),]
patCols<-rainbow.lab(length(unique(trans$patient)),alpha=.8)
names(patCols)<-unique(trans$patient)
measureCols<-rainbow.lab(length(unique(measures)),alpha=.8)
names(measureCols)<-unique(measures)
pdf('inflammation.pdf')
for(ii in unique(measures)){
  thisCols<-wells[measures==ii]
  ylims<-range(trans[,thisCols])*c(.8,1)
  plot(1,1,type='n',xlab='Days following onset of symptoms',ylab=ii,las=1,ylim=ylims,xlim=range(trans$time),log='y')
  meanMeasure<-apply(log10(trans[,thisCols]),1,mean)
  print(summary(lm(meanMeasure~time*patient,data=trans)))
  for(jj in thisCols){
    for(kk in unique(trans$patient)){
      lines(trans[trans$patient==kk,'time'],trans[trans$patient==kk,jj],col=patCols[kk],lwd=2)
    }
  }
  legend('bottomright',names(patCols),col=patCols,lwd=2,lty=1,bty='n')
}
for(kk in unique(trans$patient)){
  ylims<-range(apply(trans[,wells],2,function(xx)xx/max(xx)))
  plot(1,1,type='n',xlab='Days following onset of symptoms',ylab='Proportion of maximum measurement',las=1,ylim=ylims,xlim=range(trans$time),log='y',main=kk)
  for(ii in unique(measures)){
    thisCols<-wells[measures==ii]
    for(jj in thisCols){
        lines(trans[trans$patient==kk,'time'],trans[trans$patient==kk,jj]/max(trans[trans$patient==kk,jj]),col=measureCols[ii],lwd=2)
      }
    }
  legend('bottomright',names(measureCols),col=measureCols,lwd=2,lty=1,bty='n')
}
dev.off()

pdf('inflam_vs_cd4_and_vl.pdf')
  par(mar=c(4,4,4,4))
  for(ii in unique(measures)){
    thisCols<-wells[measures==ii]
    meanMeasure<-apply(log10(trans[,thisCols]),1,mean)
    lab<-sprintf('Regression of %s on VL p-value: %0.03f',ii,summary(lm(meanMeasure~trans$vl))$coef[2,'Pr(>|t|)'])
    plot(rep(trans$vl,2),unlist(trans[,thisCols]),bg=patCols[trans$patient],log='xy',ylab=ii,xlab='Viral load',las=1,pch=21,cex=1.2,xaxt='n',main=lab)
    logAxis(1,addExtra=TRUE)
    legend('bottomright',names(patCols),pt.bg=patCols,bty='n',pch=21,inset=c(-.15,0),xpd=NA,xjust=0)
  }
  for(ii in unique(measures)){
    thisCols<-wells[measures==ii]
    t.test(trans$cd4<500,meanMeasure)
    meanMeasure<-apply(log10(trans[,thisCols]),1,mean)
    lab<-sprintf('One sided t test p-value: %0.03f',t.test(meanMeasure~trans$cd4<500,alternative='less')$p.value)
    plot(rep(trans$cd4,2),unlist(trans[,thisCols]),bg=patCols[trans$patient],log='xy',ylab=ii,xlab='CD4 count',las=1,pch=21,cex=1.2,xaxt='n',main=lab)
    abline(v=500,lty=2)
    logAxis(1,addExtra=TRUE)
    legend('bottomright',names(patCols),pt.bg=patCols,bty='n',pch=21,inset=c(-.15,0),xpd=NA,xjust=0)
  }
dev.off()

