
env<-dnar::read.fa('work/alignHmm_2019-05-07/Env.nt.fasta.gz')
hxEnv<-env[1,]
env<-env[-1,]
env$mm<-sapply(strsplit(env$name,'\\.'),'[',1)
env$time<-as.numeric(sapply(strsplit(env$name,'\\.'),'[',4))
env$sgs<-(sapply(strsplit(env$name,'\\.'),'[',3))
env$source<-(sapply(strsplit(env$name,'\\.'),'[',2))
cols<-c(A = "green",T = "red", C = "blue", G = "yellow", `-` = "grey")
for(early in c(TRUE,FALSE)){
  pdf(sprintf('out/highlight_%s.pdf',ifelse(early,'first','all')),height=10,width=8)
    if(early)par(mfrow=c(3,3),mar=c(4,3.5,1,.6))
    else par(mar=c(3,3,1,.1))
    counter<-1
    for(ii in sort(unique(env$mm))){
      message(ii)
      thisDat<-env[env$mm==ii&nchar(dnar::degap(env$seq))>nchar(dnar::degap(hxEnv$seq))-500&env$sgs=='SGS'&env$source=='PLAS',]
      if(nrow(thisDat)==0)next()
      minTime<-min(thisDat$time)
      if(early)thisDat<-thisDat[thisDat$time<=minTime+30,]
      #determine ordering
      splits<-do.call(rbind,strsplit(thisDat[,'seq'],''))
      hxSplit<-strsplit(hxEnv$seq,'')[[1]]
      colSelect<-!(apply(splits=='-',2,mean)==1&hxSplit=='-')
      hxSeq<-paste(hxSplit[colSelect],collapse='')
      splits<-splits[,colSelect]
      consensus<-apply(splits[thisDat$time==minTime,],2,dnar::mostAbundant)
      dists<-outer(1:nrow(thisDat),1:nrow(thisDat),function(xx,yy)mapply(function(xxx,yyy)sum(splits[xxx,]!=splits[yyy,]),xx,yy))
      consDist<-apply(splits,1,function(xx)sum(xx!=consensus))
      clOrder<-sapply(1:nrow(thisDat),function(xx)which(hclust(as.dist(dists))$order==xx))
      newOrder<-order(thisDat$time,consDist!=0,-clOrder)
      thisDat<-thisDat[newOrder,]
      splits<-splits[newOrder,]
      plot(1,1,type='n',ylim=c(nrow(thisDat),1)+c(.6,-.6),xlim=c(1,ncol(splits)),bty='n',yaxt='n',xaxt='n',main=ii,ylab='',xlab='',yaxs='i',xaxs='i')
      abline(h=1:nrow(thisDat),col='#00000033')
      diffs<-t(apply(splits,1,function(xx)ifelse(xx!=consensus,xx,NA)))
      for(jj in 1:nrow(diffs)){
        if(any(selector<-!is.na(diffs[jj,]))) segments(which(selector),rep(jj,sum(selector))-.4,which(selector),rep(jj,sum(selector))+.4,col=cols[diffs[jj,selector]],pch='|')
      }
      timePos<-tapply(1:nrow(thisDat),thisDat$time,mean)
      axis(2,timePos,names(timePos),las=1,lwd=0,mgp=c(3,.2,0))
      xpos<-dnar::convertLineToUser(.1,2)
      segments(xpos,tapply(1:nrow(thisDat),thisDat$time,min),xpos,tapply(1:nrow(thisDat),thisDat$time,max),xpd=NA)
      xpos<-dnar::convertLineToUser(.1,4)
      segments(xpos,tapply(1:nrow(thisDat),thisDat$time,min),xpos,tapply(1:nrow(thisDat),thisDat$time,max),xpd=NA)
      xLabs<-c(1,seq(500,nchar(dnar::degap(hxSeq)),500))
      axis(1,dnar::noGap2Gap(hxSeq,xLabs),xLabs)
      if(!early|counter==8)mtext('HXB2 env position (nt)',1,line=2)
      if(!early|counter%%3==1)mtext('Time after onset of symptoms',2,line=2)
      counter<-counter+1
    }
  dev.off()
}
