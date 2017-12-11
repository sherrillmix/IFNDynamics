xx<-read.csv('data/EJs CD4 VL.csv',stringsAsFactors=FALSE)

ejs<-do.call(rbind,lapply(seq(1,ncol(xx),7),function(ii){
    out<-xx[,ii+0:5]
    colnames(out)<-c('id','date','dfosx','vl','cd4','notes')
    return(out)
}))

ejs<-ejs[!is.na(ejs$date)&ejs$date!=''&ejs$date!='Date'&ejs$date!='Sample date',]
ejs<-ejs[ejs$dfosx!='SEMEN '&!is.na(ejs$dfosx),]
#ejs$pat<-fillDown(ifelse(grepl('EJ',ejs$id),sub('\\..*$','',ejs$id),NA))
ejs$pat<-fillDown(gsub(' ','',sub('\\(.*','',sub('\\..*','',sub('EJ?','',ejs$id)))))
ejs$dfosx<-as.numeric(ejs$dfosx)
ejs$cd4<-as.numeric(sub('\\(.*','',sub(',','',ejs$cd4)))
ejs$vl<-as.numeric(sub('[,><] *','',ejs$vl))
#ejs$firstBelow500<-ave(ejs$cd4,ejs$pat,FUN=function(xx)(1:length(xx))==min(c(Inf,which(xx<500))))
perEj<-data.frame('maxCd4'=tapply(ejs$cd4,ejs$pat,function(xx)max(c(-Inf,xx),na.rm=TRUE)))
perEj$firstBelow400<-sapply(rownames(perEj),function(xx){id<-min(c(Inf,which(ejs[ejs$pat==xx,'cd4']<400)));ejs[ejs$pat==xx,][id,'dfosx']})
perEj$lastAbove400<-sapply(rownames(perEj),function(xx){id<-max(c(-Inf,which(ejs[ejs$pat==xx,'cd4']>400)));ejs[ejs$pat==xx,][id,'dfosx']})
perEj$lastAbove800<-sapply(rownames(perEj),function(xx){id<-max(c(-Inf,which(ejs[ejs$pat==xx,'cd4']>400)));ejs[ejs$pat==xx,][id,'dfosx']})


pdf('out/allPats.pdf',width=12,height=9)
  par(mfrow=c(5,5),mar=c(6,4,1,6))
  xlims<-range(ejs$dfosx,na.rm=TRUE)
  ylims<-range(ejs$cd4,na.rm=TRUE)
  ylims2<-range(ejs$vl,na.rm=TRUE)
  for(ii in unique(ejs$pat)){
    thisEj<-ejs[ejs$pat==ii,]
    plot(thisEj$dfosx,thisEj$cd4,xlim=xlims,ylim=ylims,main=sprintf('EJ%s',ii),type='l',xlab='DFoSx',ylab='CD4 count',las=1)
    par(new=TRUE)
    plot(thisEj$dfosx,thisEj$vl/1000,type='l',yaxt='n',xlab='',ylab='',xlim=xlims,ylim=ylims2/1000,xaxt='n',col='red')
    axis(4,pretty(ylims2/1000),col='red',las=1)
    mtext('Viral load (x1000)',4,line=3,cex=.7)
  }
dev.off()
