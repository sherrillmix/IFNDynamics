pre<-read.csv('data/prelim data 032120.csv',stringsAsFactors=FALSE,skip=6)

mac<-list(
  zb722=read.csv('data/macrophage replication data 050820_1.csv',skip=5,nrow=7,stringsAsFactors=FALSE,check.names=FALSE),
  zb725=read.csv('data/macrophage replication data 050820_1.csv',skip=16,nrow=7,stringsAsFactors=FALSE,check.names=FALSE),
  zb668=read.csv('data/macrophage replication data 050820_1.csv',skip=27,nrow=7,stringsAsFactors=FALSE,check.names=FALSE),
  zb620=read.csv('data/macrophage replication data 050820_2.csv',skip=5,nrow=7,stringsAsFactors=FALSE,check.names=FALSE),
  zb624=read.csv('data/macrophage replication data 050820_2.csv',skip=16,nrow=7,stringsAsFactors=FALSE,check.names=FALSE),
  zb710=read.csv('data/macrophage replication data 050820_2.csv',skip=27,nrow=7,stringsAsFactors=FALSE,check.names=FALSE),
  zb28=read.csv('data/macrophage replication data 050820_3.csv',skip=7,nrow=7,stringsAsFactors=FALSE,check.names=FALSE),
  zb31=read.csv('data/macrophage replication data 050820_3.csv',skip=18,nrow=7,stringsAsFactors=FALSE,check.names=FALSE)
)

nameFixes<-c('QVOA.post.A08-P6D3_UT'='A08 UT P6D3','QVOA.post.A08-P6F6_UT'='A08 UT P6F6','QVOA.post.A08-P5E4_UT'='A08 UT P5E4','9201 293T'='Reb. 9201_293T stock','601 293T'='Reb. 601r1_293T stock','A08 2F5 293T'='Reb. A08-2F4_293T stock','A09 293T'='Reb. A09-1A2_293T stock','A08 1A5 293T'='Reb. A08-1A5_293T stock')
macro<-do.call(rbind,mapply(function(xx,yy){
  goodCol<-!(1:ncol(xx) %in% 1) & !is.na(xx[1,])
  #xx<-xx[!apply(is.na(xx),1,all),]
  print(yy)
  out<-data.frame('donor'=rep(yy,nrow(xx)*sum(goodCol)),'day'=rep(xx$`days post infection`,sum(goodCol)),'virus'=rep(colnames(xx)[goodCol],each=nrow(xx)),'p24'=unlist(xx[,goodCol]),stringsAsFactors=FALSE)
  out$virus[out$virus %in% names(nameFixes)]<-nameFixes[out$virus[out$virus %in% names(nameFixes)]]
  out
},mac,names(mac),SIMPLIFY=FALSE))

#dupes<-c('A08 UT P5E4','A08 UT P6D3','A08 UT P6F6','A08 BE P4C8','A08 BE P5D5','A08 BE P5E1')
#macro<-macro[!macro$virus %in% dupes,]

finalNames<-read.csv('data/Table S5_05.21.2020-1.csv',stringsAsFactors=FALSE,skip=4)
finalNames<-finalNames[finalNames$ID!=''&nchar(finalNames$ID)<80,]
macNames<-data.frame('orig'=unique(macro$virus),stringsAsFactors=FALSE)
select<-!grepl(' Q2 ',macNames$orig)&!grepl('-BE4C8| BE7E2',macNames$orig)&!grepl('Ao2 pre',macNames$orig)&!grepl('A14|A02|A13',macNames$orig)&!grepl('P51A3 L80',macNames$orig)&!grepl(' BE ',macNames$orig)&!grepl('_BE$',macNames$orig)
bak<-macNames[!select,]
macNames<-macNames[select,,drop=FALSE]
rownames(macNames)<-macNames$orig
macNames$final<-NA
macNames[macNames$orig%in% finalNames$ID,'final']<-macNames[macNames$orig%in% finalNames$ID,'orig']
select<-sub(' ?TF$','',macNames$orig)%in% finalNames$ID
macNames[select,'final']<-sub(' ?TF$','',macNames[select,'orig'])
select<-sub('CH58','CH058',sub('TF 293T$','',macNames$orig))%in% finalNames$ID
macNames[select,'final']<-sub('CH58','CH058',sub('TF 293T','',macNames[select,'orig']))
nameConversion<-c('MM33 TF'='MM33.TF','Reb. 601r1_293T stock'='601.REB.r1','Reb. 9201_293T stock'='9242.REB.r1','RHPA TF'='RHPA','UG21'='UG021','UG24'='UG024','CH42 TF'='CH042','CH77 TF'='CH077','CH58 TF'='CH058','MM33 17'='MM33.17','9203 preATI GI2'='9244.VOA.G12','CH492 P51A3 293T'='CH492','CH492 P51A3'='CH492')
macNames[macNames$orig %in% names(nameConversion),'final']<-nameConversion[macNames[macNames$orig %in% names(nameConversion),'orig']]
macNames[is.na(macNames$final),'orig']
macNames$fix<-sub('QVOA.post.A08-P([0-9][A-Z][0-9]).*','A08.VOA.\\1',sub('Reb. A0([98])-([0-9][A-Z][0-9])_.*','A0\\1.REB.\\2',sub('Reb. A08.21-','A08.REB.',sub(' UT P','.VOA.',sub('BEAT','',sub('9201','9242',sub('9202','9243',sub('9203','9244',sub('9207','9241',sub(' reb ','.REB.',sub(' (preATI|postATI|Q2) ','.VOA.',macNames$orig)))))))))))
macNames[is.na(macNames$final)&macNames$fix %in% finalNames$ID,'final']<-macNames[is.na(macNames$final)&macNames$fix %in% finalNames$ID,'fix']
if(nrow(macNames[is.na(macNames$final),])>0)stop('Missing assignment for mac data')
if(any(!finalNames$ID %in% macNames$final))stop('Missing assignment for final table')
macro$fix<-macNames[macro$virus,'final']

auc<-function(xx,yy,low=min(xx),high=max(xx)){
  integrate(approxfun(xx,yy),low,high)
}
area<-do.call(rbind,lapply(structure(unique(macro$virus),.Names=unique(macro$virus)),function(vir)sapply(structure(unique(macro$donor),.Names=unique(macro$donor)),function(don){thisDat<-macro[macro$virus==vir&macro$donor==don,];if(nrow(thisDat)==0)return(NA);auc(thisDat$day,thisDat$p24,low=2,high=20)$value})))

stackArea<-data.frame('virus'=rep(rownames(area),ncol(area)),'donor'=rep(colnames(area),each=nrow(area)),'auc'=as.vector(area),row.names=NULL)
summary(mod<-lm(I(log(auc))~virus+donor+0,data=stackArea[!is.na(stackArea$auc),]))
noNa<-stackArea[!is.na(stackArea$auc),]
plot(noNa[,'auc'],exp(predict(mod)),log='xy')
abline(0,1)
coef<-sapply(rownames(area),function(xx)mod$coefficients[sprintf('virus%s',xx)])
names(coef)<-rownames(area)
pred<-predict(mod)
predMat<-do.call(rbind,lapply(structure(rownames(area),.Names=rownames(area)),function(vir)sapply(structure(colnames(area),.Names=colnames(area)),function(don){if(any(tmp<-noNa$virus==vir&noNa$donor==don))pred[tmp]else NA})))


s3<-read.csv('out/S3Update_20200226.csv',stringsAsFactors=FALSE)
rownames(s3)<-s3$Isolate.ID

virDf<-data.frame('virus'=rownames(predMat),'rebound'=grepl('[rR]eb',rownames(predMat)),'beta'=grepl('[_. -]BE',rownames(predMat)),stringsAsFactors=FALSE)
virDf$tf<-grepl('TF',virDf$virus)
virDf$voa<-grepl('VOA|pre|post|A[0-9]+ UT',virDf$virus)&!virDf$rebound
virDf$chronic<-grepl('CH492|MM33[ .]1[37]',virDf$virus)
virDf$class<-ifelse(virDf$rebound,'Rebound',ifelse(virDf$beta,'IFNb-selected',ifelse(virDf$tf,'TF',ifelse(virDf$voa,'VOA',ifelse(virDf$chronic,'Chronic','')))))
write.csv(virDf[,c('virus','class')],'out/macroVirusDoublecheck.csv',row.names=FALSE)
classes<-read.csv('data/macroVirusDoublecheck-FBRupdate.csv',stringsAsFactors=FALSE,row.names=1)
if(any(!virDf$virus %in% rownames(classes)))stop('Missing virus class')
virDf$type2<-classes[virDf$virus,'class']
virDf$auc<-coef[virDf$virus]
virDf$tropism<-NA
#Reb..A08.2F5_293T.stock from email
virDf$tropism[virDf$virus %in% c('CH470TF 293T','CH58TF 293T','Reb. 9201_293T stock','Reb. 601r1_293T stock','Reb. A09-1A2_293T stock','Reb. A08-2F5_293T stock','Reb. A08-2F5_293T stock','CH492 P51A3 293T')]<-'R5'
virDf$tropism[virDf$virus %in% c('Reb. A08-1A5_293T stock','A08 UT P7C1','A08 BE P4E6','A08 UT P5E2')]<-'X4'
virDf$tropism[virDf$virus %in% c('A08 UT P7F8','A08 UT P8E8')]<-'Dual'
virDf$tropism[virDf$virus %in% s3$Isolate.ID]<-s3[virDf$virus[virDf$virus %in% s3$Isolate.ID],'tropism']


pdf('out/macroRep.pdf',width=20,height=3.5)
  par(mfrow=c(1,length(unique(macro$donor))),mar=c(4.4,4,1.1,6))
  ylim<-range(macro$p24)
  for(ii in unique(macro$donor)){
    plot(1,1,type='n',ylim=ylim,ylab='',yaxt='n',xlim=range(macro$day),xlab='Day',log='y',mgp=c(1.75,.5,0),tcl=-.3,main=ii)
    title(ylab='p24 concentration (ng/ml)',mgp=c(2.75,1,0))
    dnar::logAxis(las=1)
    thisDat<-macro[macro$donor==ii,]
    #thisDat<-thisDat[order(thisDat$virus,thisDat$day),]
    for(jj in unique(thisDat$virus)){
      if(any(thisDat$virus==jj))for(kk in seq(1,sum(thisDat$virus==jj),7))lines(thisDat[thisDat$virus==jj,'day'][kk+0:6],thisDat[thisDat$virus==jj,'p24'][kk+0:6])
    }
    isLast<-thisDat$day==max(thisDat$day)
    thisLast<-thisDat[isLast,]
    thisLast<-thisLast[order(thisLast$p24),]
    clusters<-cutree(hclust(dist(log(thisLast$p24))),h=.15)
    nClust<-ave(clusters,clusters,FUN=length)
    probs<-nClust>1
    probAdjust<-thisLast[probs,'p24']*1.15^ave(clusters[probs],clusters[probs],FUN=function(xx)1:length(xx)-length(xx)/2-.5)
    axis(4,thisLast$p24[!probs],thisLast$virus[!probs],las=1,cex.axis=.4,tcl=-.1,mgp=c(1,.2,0))
    axis(4,thisLast$p24[probs],rep('',sum(probs)),las=1,cex.axis=.4,tcl=-.1,mgp=c(1,.2,0))
    axis(4,probAdjust,thisDat[isLast,'virus'][probs],las=1,cex.axis=.4,tcl=-.1,mgp=c(1,.5,0),xpd=NA,tick=FALSE)
    segments(dnar::convertLineToUser(.1,4),thisLast[probs,'p24'],dnar::convertLineToUser(.5,4),probAdjust,xpd=NA,lwd=.5)
  }
  par(mfrow=c(4,ceiling(length(unique(macro$virus))/4)),mar=c(0,0,0,0))
  ylim<-range(macro$p24)
  donorCol<-structure(dnar::rainbow.lab(length(unique(macro$donor))),.Names=unique(macro$donor))
  donorOffset<-structure(2^seq(-.5,.5,length.out=length(unique(macro$donor))),.Names=unique(macro$donor))
  for(ii in unique(macro$virus)){
    plot(1,1,type='n',ylim=ylim,ylab='',yaxt='n',xlim=range(macro$day),xlab='Day',log='y',mgp=c(1.75,.5,0),tcl=-.3,xaxt='n')
    mtext(ii,3,line=-1,cex=.35)
    #title(ylab='p24 concentration (ng/ml)',mgp=c(2.75,1,0))
    #dnar::logAxis(las=1)
    thisDat<-macro[macro$virus==ii,]
    #thisDat<-thisDat[order(thisDat$donor,thisDat$day),]
    for(jj in unique(thisDat$donor)){
      #hard coding 7 days
      for(kk in seq(1,sum(thisDat$donor==jj),7))lines(thisDat[thisDat$donor==jj,'day'][kk+0:6],thisDat[thisDat$donor==jj,'p24'][kk+0:6]*ifelse(thisDat[thisDat$donor==jj,'p24'][kk+0:6]==.01,donorOffset[jj],1),col=donorCol[jj])
    }
  }
  plot(1,1,type='n',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  legend('right',inset=.07,names(donorCol),col=donorCol,lty=1,cex=.5)
dev.off()



pdf('out/macroAucHeat.pdf',height=12,width=7)
  par(mar=c(3.3,6,.1,3.1))
  breaks<-exp(seq(min(log(area),na.rm=TRUE)*.99,max(log(area),na.rm=TRUE)*1.01,length.out=101))
  cols<-rev(heat.colors(100))
  plotFunc<-function(area){
  image(1:ncol(area),1:nrow(area),t(area),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='')
  nas<-which(is.na(area),arr.ind=TRUE)
  rect(nas[,2]-.5,nas[,1]-.5,nas[,2]+.5,nas[,1]+.5,col='#CCCCCC',border=NA)
  abline(v=2:ncol(area)-.5,h=2:nrow(area)-.5,col='#777777')
  box()
  axis(1,1:ncol(area),colnames(area),tcl=-.3,mgp=c(1,.4,0))
  axis(2,1:nrow(area),rownames(area),las=1,cex.axis=.5)
  axis(4,1:nrow(area),sprintf('%.02f',exp(coef[rownames(area)])),las=1,cex.axis=.5)
  dnar::insetScale(log10(breaks),cols,insetPos=c(.0175,.01,.025,.24),at=0:3,labels=sapply(0:3,function(xx)as.expression(bquote(10^.(xx)))),main='p24 AUC')
  }
  #plotFunc(area)
  plotFunc(area[order(coef[rownames(area)]),])
dev.off()


pdf('out/macroCompare.pdf',width=10,height=5)
  xPos<-structure(1:length(unique(virDf$type2)),.Names=unique(virDf$type2))
  posOffset<-ave(virDf$auc,virDf$type2,FUN=function(xx)seq(-.05,.05,length.out=length(xx)))
  plot(xPos[as.character(virDf$type2)]+posOffset,exp(virDf$auc),yaxt='n',xaxt='n',log='y',ylab='Inferred p24 AUC',xlab='',pch=21,bg=ifelse(is.na(virDf$tropism),NA,ifelse(virDf$tropism=='R5','#FF000033','#0000FF33')))
  abline(h=.2,lty=2)
  dnar::logAxis(las=1)
  axis(1,xPos,names(xPos))
  legend('topright',inset=.01,c('R5','X4/Dual'),pch=21,pt.bg=c('#FF000033','#0000FF33'))
  tmp<-virDf[virDf$tropism=='R5'&!is.na(virDf$tropism),]
  xPos<-structure(1:length(levels(tmp$type2)),.Names=levels(tmp$type2))
  posOffset<-ave(tmp$auc,tmp$type2,FUN=function(xx)seq(-.05,.05,length.out=length(xx)))
  #plot(xPos[as.character(tmp$type2)]+posOffset,exp(tmp$auc),yaxt='n',xaxt='n',log='y',ylab='Base p24 AUC',xlab='',main='R5 viruses')
  #abline(h=.2,lty=2)
  #dnar::logAxis(las=1)
  #axis(1,xPos,names(xPos))
dev.off()


pdf('out/macPredHeat.pdf')
  image(1:ncol(area),1:nrow(area),exp(t(predMat)),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='')
  axis(1,1:ncol(area),colnames(area),tcl=-.3,mgp=c(1,.4,0))
  axis(2,1:nrow(area),rownames(area),las=1,cex.axis=.5)
  nas<-which(is.na(area),arr.ind=TRUE)
  rect(nas[,2]-.5,nas[,1]-.5,nas[,2]+.5,nas[,1]+.5,col='#CCCCCC',border=NA)
  abline(v=2:ncol(area)-.5,h=2:nrow(area)-.5,col='#777777')
  box()
  dnar::insetScale(log10(breaks),cols,insetPos=c(.0175,.01,.025,.3),at=0:3,labels=sapply(0:3,function(xx)as.expression(bquote(10^.(xx)))),main='p24 AUC')
dev.off()

#virDf$virus %in% seqs$name
seqs<-dnar::read.fa('a08Seqs/AlignSeq.nt.fasta')
seqs$trim<-substring(seqs$seq,6500,8500)
pdf('test.pdf')
dnaplotr::plotDNA(seqs$seq)
dnaplotr::plotDNA(seqs$trim)
dev.off()
dists<-levenR::leven(seqs$trim[-1],nThreads=3)
seqSplit<-do.call(rbind,strsplit(seqs$trim,''))
dists2<-do.call(rbind,lapply(2:nrow(seqSplit),function(xx)sapply(2:nrow(seqSplit),function(yy)sum(seqSplit[xx,]!=seqSplit[yy,]))))
rownames(dists2)<-seqs$name[-1]
tree<-phangorn::NJ(dists2)

a08<-virDf[grep('A08',virDf$virus),]
a08$id<-sapply(strsplit(sub('_(BE|UT)$','',sub('[ _]293T.*','',a08$virus)),'[ _.-]'),tail,1)
hits<-sapply(a08$id,grep,seqs$name)
names(hits)<-a08$virus
hits[['Reb. A08.21-7C1']]<-grep('21[._-]7C1',seqs$name)
hits[['Reb. A08.21-7D3']]<-grep('21[._-]7D3',seqs$name)
hits[sapply(hits,length)!=1]
if(any(sapply(hits,length)>1))stop('Multi hit')
a08$seqId<-seqs$name[sapply(hits,function(xx)if(length(xx)==0) NA else xx)]
probs<-is.na(a08$ic50)
hits2<-sapply(a08[probs,'id'],function(xx)which(grepl('A08',s3$Isolate.ID)&grepl(xx,s3$Isolate.ID)))
if(any(sapply(hits2,length)>1))stop('Multi hit')
a08[probs,'ic50']<-sapply(hits2,function(xx)if(length(xx)==1)s3[xx,'IFN.b.IC50..pg.ml.'] else NA)
a08$isRebound<-grepl('Reb',a08$virus)
aucs<-sapply(tree$tip.label,function(xx)if(any(tmp<-a08$seqId==xx&!is.na(a08$seqId)))a08[tmp,'auc'] else NA)
trops<-sapply(tree$tip.label,function(xx)if(any(tmp<-a08$seqId==xx&!is.na(a08$seqId)))a08[tmp,'tropism'] else NA)
ics<-sapply(tree$tip.label,function(xx)if(any(tmp<-a08$seqId==xx&!is.na(a08$seqId)))a08[tmp,'ic50'] else NA)
rebound<-sapply(tree$tip.label,function(xx)if(any(tmp<-a08$seqId==xx&!is.na(a08$seqId)))a08[tmp,'isRebound'] else NA)



breaks<-exp(seq(min(aucs,na.rm=TRUE)-.01,max(aucs,na.rm=TRUE)+.01,length.out=101))
aucCut<-cut(exp(aucs),breaks)
aucCols<-dnar::rainbow.lab(100)
breaks2<-exp(seq(min(log(s3$IFN.b.IC50..pg.ml.),na.rm=TRUE)-.01,max(log(s3$IFN.b.IC50..pg.ml.),na.rm=TRUE)+.01,length.out=101))
icCut<-cut(ics,breaks2)
icCols<-rev(heat.colors(120))[-1:-20]

tropCols<-c('X4'='#0000FF','R5'='#FF0000','Dual'='#880088')
library(ggtree)
pdf('out/A08MacroTree.pdf')
  out<-ggtree(tree)+#,size=.2,ladderize=FALSE)
    geom_tiplab(color=tropCols[trops],size=min(7.5,1500/length(tree$tip.label)),label='-',vjust=.35)+
    geom_tiplab(color=aucCols[as.numeric(aucCut)],size=min(7.5,1500/length(tree$tip.label)),label='  -',vjust=.35)+
    geom_tiplab(color=icCols[as.numeric(icCut)],size=min(7.5,1500/length(tree$tip.label)),label='    -',vjust=.35)+
    geom_tiplab(size=5,label=ifelse(!is.na(rebound)&rebound,'         *',''),vjust=.8)+
    geom_tiplab(size=2,label=ifelse(grepl('M[0-9]+$',tree$tip.label),sprintf('%s',sub('^.*(M[0-9]+)$','\\1',tree$tip.label)),''),vjust=.4)+
    geom_treescale(offset=-5,fontsize=2.5,x=.01,y=-.05*length(tree$tip.label))
  #print(out)
  vp<-grid::viewport()
  par(mar=c(0,0,0,0))
  plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  print(out,vp=vp)
  ticks<-log10(c(.2,.5,1,2,4))
  ticks<-ticks[10^ticks>min(breaks)&10^ticks<max(breaks)]
  legend(grconvertX(.175,'npc','user'),grconvertY(.05,'npc','user'),names(tropCols),fill=tropCols,cex=.65,yjust=.4,bty='n')
  #dnar::insetScale(log(breaks),aucCols,main='Macrophage growth (AUC)',insetPos =c(0.05, 0.3, 0.07, 0.5),at=log(10^ticks),lab=sapply(ticks,function(xx)as.expression(bquote(10^.(xx)))),cex=.7)
  dnar::insetScale(log(breaks),aucCols,main='Macrophage growth (AUC)',insetPos =c(0.05, 0.3, 0.07, 0.5),at=log(10^ticks),lab=10^ticks,cex=.7)
  ticks2<-pretty(log10(breaks2))
  ticks2<-ticks2[10^ticks2>min(breaks2)&10^ticks2<max(breaks2)]
  dnar::insetScale(log(breaks2),icCols,main='IFNb IC50 (pg/ml)',insetPos =c(0.05, 0.55, 0.07, 0.75),at=log(10^ticks2),lab=sapply(ticks2,function(xx)as.expression(bquote(10^.(xx)))),cex=.7)
  legend(grconvertX(.8,'npc','user'),grconvertY(.05,'npc','user'),'Rebound',pch='*',pt.cex=1,cex=.7,yjust=.4,bty='n')
dev.off()

pdf('out/A08_ic50_vs_macro.pdf',height=5,width=5)
par(mar=c(4,4,.3,.5))
plot(a08$ic50,exp(ifelse(a08$auc<log(.15),log(.15),a08$auc)),bg=sprintf('%s66',tropCols[a08$tropism]),pch=21,ylab='Macrophage replication (AUC)',log='xy',las=1,xaxt='n',xlab='IFNb IC50 (pg/ml)')
dnar::logAxis(1)
dev.off()
