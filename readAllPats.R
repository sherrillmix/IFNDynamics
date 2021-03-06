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
ejs[ejs$pat=='6911','pat']<-'69'
ejs$dfosx<-as.numeric(ejs$dfosx)
ejs$cd4<-as.numeric(sub('[A-Za-z].*$','',sub('\\(.*','',sub('[?,]','',ejs$cd4))))
ejs$vl<-as.numeric(sub('[A-Za-z].*$','',sub('N/D|\\?','',gsub('[,><] *','',ejs$vl))))
ejs$ej<-sprintf('EJ%s',ejs$pat)
ejs$date[ejs$date==' 12/07/2001'&ejs$ej=='EJ50']<-'7/12/2001'
ejs$date[ejs$date=='6/9/2001'&ejs$ej=='EJ52']<-'9/6/2001'
#ejs$firstBelow500<-ave(ejs$cd4,ejs$pat,FUN=function(xx)(1:length(xx))==min(c(Inf,which(xx<500))))
perEj<-data.frame('maxCd4'=tapply(ejs$cd4,ejs$pat,function(xx)max(c(-Inf,xx),na.rm=TRUE)))
perEj$firstBelow400<-sapply(rownames(perEj),function(xx){id<-min(c(Inf,which(ejs[ejs$pat==xx,'cd4']<400)));ejs[ejs$pat==xx,][id,'dfosx']})
perEj$lastAbove400<-sapply(rownames(perEj),function(xx){id<-max(c(-Inf,which(ejs[ejs$pat==xx,'cd4']>400)));ejs[ejs$pat==xx,][id,'dfosx']})
perEj$lastAbove800<-sapply(rownames(perEj),function(xx){id<-max(c(-Inf,which(ejs[ejs$pat==xx,'cd4']>400)));ejs[ejs$pat==xx,][id,'dfosx']})


pats<-unique(ejs$pat)
desiredPats<-c(23,42,89,92,104,67)
#pats<-pats[order(!pats %in% desiredPats,as.numeric(pats))]
pats<-c(desiredPats,pats[!pats %in% desiredPats])
pdf('out/allPats.pdf',width=12,height=8)
  par(mfrow=c(5,5),mar=c(6,4,1,6))
  xlims<-range(ejs$dfosx,na.rm=TRUE)
  ylims<-range(ejs$cd4,na.rm=TRUE)
  ylims2<-range(ejs$vl,na.rm=TRUE)
  for(ii in pats){
    thisEj<-ejs[ejs$pat==ii&!is.na(ejs$dfosx),]
    withAs(thisEj=thisEj[!is.na(thisEj$cd4),],plot(thisEj$dfosx,thisEj$cd4,xlim=xlims,ylim=ylims,main=sprintf('EJ%s',ii),type='l',xlab='DFoSx',ylab='CD4 count',las=1))
    points(thisEj$dfosx,thisEj$cd4,cex=.7)
    par(new=TRUE)
    withAs(thisEj=thisEj[!is.na(thisEj$vl),],plot(thisEj$dfosx,thisEj$vl,type='l',yaxt='n',xlab='',ylab='',xlim=xlims,ylim=ylims2,xaxt='n',col='red',log='y'))
    logAxis(4,las=1,col.ticks='red',cex.axis=.9,col.axis='red')
    points(thisEj$dfosx,thisEj$vl,cex=.7,col='red')
    axis(1,thisEj$dfosx,rep('',nrow(thisEj)),tcl=-.25)
    #axis(4,pretty(ylims2/1000),col='red',las=1)
    mtext('Viral load',4,line=2.7,cex=.7,col='red')
  }
dev.off()
