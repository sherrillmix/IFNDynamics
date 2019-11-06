
miniIn<-read.csv('out/miniInputTitration.csv',row.names=1)
miniIn$id<-rownames(miniIn)
miniOut<-read.csv('out/miniOutputTitration.csv',row.names=1)
miniOut$id<-sub(' retro','',rownames(miniOut))
miniOut$retro<-grepl('retro',rownames(miniOut))

combo<-merge(miniOut[miniOut$retro,c('id','No.Drug')],miniOut[!miniOut$retro,c('id','No.Drug')],by.y='id',by.x='id',all.x=TRUE,all.y=TRUE,suffixes=c('.retro',''))
plot(combo$No.Drug,combo$No.Drug.retro);abline(0,1)
combo<-merge(combo,miniIn[,c('No.Drug','id')],by.y='id',by.x='id',suffixes=c('','.in'),all.x=TRUE,all.y=TRUE)
rt<-read.csv('data/Mini_expansion_RT.csv',row.names=1)
rt<-rt[,!grepl('^Input',colnames(rt))]
colnames(rt)<-sub('\\.+$','',sub('RT\\.ng\\.ml\\.\\.','',colnames(rt)))
#beads and retronectin were flipped on plate
colnames(rt)[c(which(colnames(rt)=='Beads'),which(colnames(rt)=='Retron'))]<-c('Retron','Beads')
rownames(rt)<-sub('^Rb','Rebound',rownames(rt))
colnames(rt)<-sprintf('RT_%s',colnames(rt))
#the two batches were switched in plate (photo of plate lid cleared up)
if(any(!combo$id %in% c(rownames(rt),'SG3','YU2','SG3+WEAU','89.6','Media')))stop('Missing RT')
combo<-cbind(combo,rt[combo$id,])
combo$infectivity_Retro<-combo$No.Drug.retro/combo$RT_Retron/1000
combo$infectivity_Bead<-combo$No.Drug/combo$RT_Beads/1000
combo$infectivity_Old<-combo$No.Drug.in/combo$RT_Old/1000
pdf('out/infectivityCompare.pdf')
  par(mar=c(3,3.2,.1,.1))
  lims<-dnar::withAs(xx=c(combo$infectivity_Old,combo$infectivity_Bead,combo$infectivity_Retro),range(xx[xx<Inf],na.rm=TRUE))
  lims<-NULL
  dnar::withAs(combo=combo[!is.na(combo$infectivity_Retro)&!is.na(combo$infectivity_Bead)&combo$infectivity_Retro+combo$infectivity_Bead<Inf,],plot(combo$infectivity_Retro,combo$infectivity_Bead,ylab='Bead isolate infectivity (IU/pg RT)',pch=21,bg='#00000033',cex=1.3,mgp=c(2.3,.8,0),las=1,xlab='',xlim=lims,ylim=lims,log='xy'))
  title(xlab='Retronectin isolate infectivity (IU/pg RT)',mgp=c(1.9,1,0))
  abline(0,1,lty=1)
  dnar::withAs(combo=combo[!is.na(combo$infectivity_Old)&!is.na(combo$infectivity_Bead)&combo$infectivity_Old+combo$infectivity_Bead<Inf,],plot(combo$infectivity_Old,combo$infectivity_Bead,ylab='Bead isolate infectivity (IU/pg RT)',pch=21,bg='#00000033',cex=1.3,mgp=c(2.3,.8,0),las=1,xlab='',xlim=lims,ylim=lims))
  title(xlab='Input isolate infectivity (IU/pg RT)',mgp=c(1.9,1,0))
  abline(0,1,lty=1)
  dnar::withAs(combo=combo[!is.na(combo$infectivity_Old)&!is.na(combo$infectivity_Retro)&combo$infectivity_Old+combo$infectivity_Retro<Inf,],plot(combo$infectivity_Old,combo$infectivity_Retro,ylab='Retro isolate infectivity (IU/pg RT)',pch=21,bg='#00000033',cex=1.3,mgp=c(2.3,.8,0),las=1,xlab='',xlim=lims,ylim=lims,log='xy'))
  title(xlab='Input isolate infectivity (IU/pg RT)',mgp=c(1.9,1,0))
  abline(0,1,lty=1)
dev.off()
newNames<-c('No.Drug.retro'='IU_Retro','No.Drug'='IU_Bead','No.Drug.in'='IU_Old','RT_Retron'='RT_Retro','RT_Beads'='RT_Bead')
for(ii in names(newNames))colnames(combo)[colnames(combo)==ii]<-newNames[ii]
write.csv(combo[!is.na(combo$RT_Old),],'out/miniInfectivity_20190820.csv',row.names=FALSE)

rep<-read.csv("data/miniRepCap.csv",stringsAsFactors=FALSE,header=FALSE)
firstExp<-read.csv('data/Marvin infectivity RT titers 07042019-1.csv',stringsAsFactors=FALSE)
firstExp<-firstExp[firstExp$X!='CAM13 control',]
firstExp$inf<-as.numeric(ifelse(grepl('DIV',firstExp$IU.ng.RT),NA,firstExp$IU.ng.RT))/1000
rownames(firstExp)<-firstExp$X
inFirst<-combo$id %in% firstExp$X
firstExp[combo[inFirst,'id'],'inf']<-combo[inFirst,'infectivity_Bead']
firstExp$class<-rep(c('Rebound','VOA','MM acute','TP chronic','MM chronic','TP acute'),c(48,33,8,6,8,8))
pdf('out/infectivityByClass.pdf',width=8,height=4)
  for(repCap in c(FALSE,TRUE)){
    tmp<-firstExp[!is.na(firstExp$inf)&firstExp$inf<Inf&firstExp$inf>0,]
    if(repCap)tmp<-tmp[tmp$X %in% rep$V2,]
    uniqC<-unique(tmp$class)
    pos<-structure(1:length(unique(tmp$class)),.Names=uniqC[order(!grepl('MM',uniqC),!grepl('TP',uniqC),grepl('Rebound|acute',uniqC))])
    par(mar=c(3,3.5,.2,.2),lheight=.8)
    plot(1,1,type='n',las=1,ylab='Infectivity (IU/pg RT)',xlim=c(1,length(unique(firstExp$class))+.2),ylim=range(tmp$inf),xaxt='n',log='y',yaxt='n',mgp=c(2.3,1,.3),xlab='')
    dnar::logAxis(las=1)
    axis(1,pos,sub(' ','\n',names(pos)),padj=1,mgp=c(1,.1,0))
    offset<-ave(tmp$inf,tmp$class,FUN=function(xx)beeswarm::swarmx(0,xx,cex=1.2)[[1]])
    points(pos[(tmp$class)]+offset,tmp$inf,cex=1.2,bg=ifelse(tmp$X %in% combo$id,'#0000FF','#FF0000'),pch=21)
  }
dev.off()

newInfect<-read.csv('out/2019-09-19_miniInfect.csv',row.names=1)
newInfectCont<-newInfect[!is.na(newInfect$NA.control),c('NA.control'),drop=FALSE]
newInfect<-newInfect[rownames(newInfect)!='Media',c('X2_reb.bead','X1_reb.bead','X1_other.bead','X2_reb.retro')]
newInfect<-newInfect[apply(!is.na(newInfect),1,sum)>0,]
if(any(!is.na(newInfect$X1_other.bead)&!is.na(newInfect$X2_reb.retro)))stop('Overlapping values')
if(any(!is.na(newInfect$X1_reb.bead)&!is.na(newInfect$X2_reb.retro)))stop('Overlapping values')
if(any(!is.na(newInfect$X1_reb.bead)&!is.na(newInfect$X1_other.bead)))stop('Overlapping values')
newInfect$redoBead<-newInfect$X2_reb.bead
newInfect$redoRetro<-newInfect$X2_reb.retro
newInfect$redoFirst<-ifelse(is.na(newInfect$X1_reb.bead),newInfect$X1_other.bead,newInfect$X1_reb.bead)

combo$id %in% rownames(newInfect)
rownames(newInfect)[!rownames(newInfect) %in% combo$id]
newInfect$id<-rownames(newInfect)
rt$id<-rownames(rt)
firstExp$id<-rownames(firstExp)
all<-merge(combo[,c('id','IU_Retro','IU_Bead','IU_Old','infectivity_Retro','infectivity_Bead','infectivity_Old','RT_Bead','RT_Retro','RT_Old')],newInfect[,c('id','redoBead','redoRetro','redoFirst')],all.x=TRUE,all.y=TRUE)
all<-merge(firstExp[,c('id','ng.RT.uL','class')],all,all.y=TRUE)
colnames(all)[colnames(all)=='ng.RT.uL']<-'RT_First'
all$redoBeadInf<-all$redoBead/all$RT_Bead/1000
all$redoFirstInf<-all$redoFirst/all$RT_First/1000
all$inf<-ifelse(is.na(all$redoBeadInf),all$redoFirstInf,all$redoBeadInf)
all$inf[all$inf==Inf]<-NA

tmp<-merge(firstExp[,c('id','avg.IU.uL')],newInfect[,c('id','redoBead','redoRetro')])

pdf('out/infectivityByClass_20190919.pdf',width=6,height=4)
  for(repCap in c(FALSE,TRUE)){
    tmp<-all[!is.na(all$inf)&!is.na(all$class),]
    if(repCap)tmp<-tmp[tmp$id %in% rep$V2,]
    uniqC<-unique(tmp$class)
    pos<-structure(1:length(unique(tmp$class)),.Names=uniqC[order(!grepl('MM',uniqC),!grepl('TP',uniqC),grepl('Rebound|acute',uniqC))])
    par(mar=c(3,3.5,.2,.2),lheight=.8)
    plot(1,1,type='n',las=1,ylab='Infectivity (IU/pg RT)',xlim=c(1,length(unique(tmp$class))+.2)+c(-.3,.3),ylim=range(tmp$inf),xaxt='n',log='y',yaxt='n',mgp=c(2.3,1,.3),xlab='')
    dnar::logAxis(las=1)
    axis(1,pos,sub(' ','\n',names(pos)),padj=1,mgp=c(1,.1,0))
    offset<-ave(tmp$inf,tmp$class,FUN=function(xx)beeswarm::swarmx(0,xx,cex=1.2)[[1]])
    points(pos[(tmp$class)]+offset,tmp$inf,cex=1.2,pch=21,bg='#00000055')
  }
dev.off()


redoInfect<-read.csv('out/2019-09-27_miniInfect.csv',row.names=1)
redoInfectCont<-redoInfect[!is.na(redoInfect$NA.control),c('NA.control'),drop=FALSE]
redoInfect<-redoInfect[rownames(redoInfect)!='Media',c('X2_reb.bead','X1_reb.bead','X1_other.bead','X2_reb.retro')]
redoInfect<-redoInfect[apply(!is.na(redoInfect),1,sum)>0,]
if(any(!is.na(redoInfect$X1_other.bead)&!is.na(redoInfect$X2_reb.retro)))stop('Overlapping values')
if(any(!is.na(redoInfect$X1_reb.bead)&!is.na(redoInfect$X2_reb.retro)))stop('Overlapping values')
if(any(!is.na(redoInfect$X1_reb.bead)&!is.na(redoInfect$X1_other.bead)))stop('Overlapping values')
redoInfect$redoBead<-redoInfect$X2_reb.bead
redoInfect$redoRetro<-redoInfect$X2_reb.retro
redoInfect$redoFirst<-ifelse(is.na(redoInfect$X1_reb.bead),redoInfect$X1_other.bead,redoInfect$X1_reb.bead)
rownames(redoInfect)[!rownames(redoInfect) %in% combo$id]
redoInfect$id<-rownames(redoInfect)
#
rt$id<-rownames(rt)
firstExp$id<-rownames(firstExp)
all<-merge(combo[,c('id','IU_Retro','IU_Bead','IU_Old','infectivity_Retro','infectivity_Bead','infectivity_Old','RT_Bead','RT_Retro','RT_Old')],redoInfect[,c('id','redoBead','redoRetro','redoFirst')],all.x=TRUE,all.y=TRUE)
all<-merge(firstExp[,c('id','ng.RT.uL','class')],all,all.y=TRUE)
colnames(all)[colnames(all)=='ng.RT.uL']<-'RT_First'
all$redoBeadInf<-all$redoBead/all$RT_Bead/1000
all$redoFirstInf<-all$redoFirst/all$RT_First/1000
all$inf<-ifelse(is.na(all$redoBeadInf),all$redoFirstInf,all$redoBeadInf)
all$inf[all$inf==Inf]<-NA
#
pdf('out/infectivityByClass_20190927.pdf',width=6,height=4)
  for(repCap in c(FALSE,TRUE)){
    tmp<-all[!is.na(all$inf)&!is.na(all$class),]
    if(repCap)tmp<-tmp[tmp$id %in% rep$V2,]
    uniqC<-unique(tmp$class)
    pos<-structure(1:length(unique(tmp$class)),.Names=uniqC[order(!grepl('MM',uniqC),!grepl('TP',uniqC),grepl('Rebound|acute',uniqC))])
    par(mar=c(3,3.5,.2,.2),lheight=.8)
    plot(1,1,type='n',las=1,ylab='Infectivity (IU/pg RT)',xlim=c(1,length(unique(tmp$class))+.2)+c(-.3,.3),ylim=range(tmp$inf),xaxt='n',log='y',yaxt='n',mgp=c(2.3,1,.3),xlab='')
    dnar::logAxis(las=1)
    axis(1,pos,sub(' ','\n',names(pos)),padj=1,mgp=c(1,.1,0))
    offset<-ave(tmp$inf,tmp$class,FUN=function(xx)beeswarm::swarmx(0,xx,cex=1.2)[[1]])
    points(pos[(tmp$class)]+offset,tmp$inf,cex=1.2,pch=21,bg='#00000055')
  }
dev.off()


repCap<-read.csv('data/Data Master 2019 _repCap.csv',stringsAsFactors=FALSE)
sapply(sapply(repCap$ID,grep,all$id),length)
infRow<-sapply(sub('  .*','',sub('Rebound |QVOA ','',trimws(repCap$ID))),function(xx){out<-grep(xx,all$id,fixed=TRUE);if(length(out)!=1)NA else out})
repCap$infectivity<-all[infRow,'inf']
tmp<-all[! 1:nrow(all) %in% infRow&!all$id %in% c('SG3','SG3+WEAU','YU2')&!is.na(all$inf),c('id','inf')]
colnames(tmp)<-c('ID','infectivity')
tmp[,colnames(repCap)[!colnames(repCap) %in% colnames(tmp)]]<-NA
out<-rbind(repCap,tmp)
write.csv(out,'out/mini_repInf.csv',row.names=FALSE)
inTrop<-read.csv('out/miniInputTitration.csv',row.names=1)
inTropProp<-inTrop[,colnames(inTrop)!='No.Drug']/inTrop$No.Drug*100
#outTrop<-read.csv('out/miniOutputTitration.csv',row.names=1)
#outTropProp<-outTrop[,colnames(outTrop)!='No.Drug']/outTrop$No.Drug
colnames(inTropProp)[colnames(inTropProp)=='AMD.Mar']<-'AMD-Mar'
inTropProp<-inTropProp[!rownames(inTropProp) %in% c('YU2','Media','SG3','89.6','SG3+WEAU'),]
out[,colnames(inTropProp)]<-inTropProp[out$ID,]
leftovers<-inTropProp[!rownames(inTropProp) %in% out$ID,]
leftovers$ID<-rownames(leftovers)
leftovers[,colnames(out)[!colnames(out) %in% colnames(leftovers)]]<-NA
out<-rbind(out,leftovers[,colnames(out)])
xx<-read.csv('rebound/out/rebQvoaMaster.csv',stringsAsFactors=FALSE)
rownames(xx)<-xx$ID
out[is.na(out[,1]),1]<-xx[out[is.na(out[,1]),"ID"],'Type']
write.csv(out,'out/mini_repInfTrop.csv',row.names=FALSE)



check<-all[,c('id','class','redoBead','RT_Bead','redoFirst','RT_First','inf')]
check$iu<-ifelse(is.na(check$redoFirst),check$redoBead,check$redoFirst)
check$rt<-ifelse(is.na(check$redoFirst),check$RT_Bead,check$RT_First)
check[!is.na(check$inf),c('id','rt','iu','inf')]
write.csv(check[!is.na(check$inf),c('id','rt','iu','inf')],'out/infectivity_check.csv',row.names=FALSE)
pdf('out/infectivity_reboundVoa.pdf',width=4,height=4)
par(mar=c(3.1,3.1,.1,.1))
withAs(xx=check[!is.na(check$inf)&check$class %in% c("VOA","Rebound"),],plot(xx$rt,xx$iu,log='xy',xlab='RT (ng/ul)',ylab='IU/ul',pch=21,bg=ifelse(xx$class=='VOA','blue','red'),xaxt='n',yaxt='n',mgp=c(2.1,1,0)))
dnar::logAxis(las=1)
dnar::logAxis(1)
dev.off()

redoInfect2<-read.csv('out/2019-11-04_miniInfectT12.csv',row.names=1)
colnames(redoInfect2)<-sub('\\.$','',colnames(redoInfect2))
#redoInfect2<-redoInfect2[grepl('T12',colnames(redoInfect2)),]
#colnames(redoInfect2)<-sub('\\..T12$','',colnames(redoInfect2))
redoInfect2<-redoInfect2[rownames(redoInfect2)!='Media',c('X2_reb.bead','X1_reb.bead','X1_other.bead','X2_reb.retro','X2_reb.bead..T12','X1_reb.bead..T12','X1_other.bead..T12','X2_reb.retro..T12')]
redoInfect2<-redoInfect2[apply(!is.na(redoInfect2),1,sum)>0,]
if(any(!is.na(redoInfect2$X1_other.bead)&!is.na(redoInfect2$X2_reb.retro)))stop('Overlapping values')
if(any(!is.na(redoInfect2$X1_reb.bead)&!is.na(redoInfect2$X2_reb.retro)))stop('Overlapping values')
if(any(!is.na(redoInfect2$X1_reb.bead)&!is.na(redoInfect2$X1_other.bead)))stop('Overlapping values')
redoInfect2$redoBead<-redoInfect2$X2_reb.bead
redoInfect2$redoRetro<-redoInfect2$X2_reb.retro
redoInfect2$redoFirst<-ifelse(is.na(redoInfect2$X1_reb.bead),redoInfect2$X1_other.bead,redoInfect2$X1_reb.bead)
redoInfect2$redoBeadT12<-redoInfect2$X2_reb.bead..T12
redoInfect2$redoRetroT12<-redoInfect2$X2_reb.retro..T12
redoInfect2$redoFirstT12<-ifelse(is.na(redoInfect2$X1_reb.bead..T12),redoInfect2$X1_other.bead..T12,redoInfect2$X1_reb.bead..T12)
rownames(redoInfect2)[!rownames(redoInfect2) %in% combo$id]
redoInfect2$id<-rownames(redoInfect2)
#
rt$id<-rownames(rt)
firstExp$id<-rownames(firstExp)
all<-merge(combo[,c('id','IU_Retro','IU_Bead','IU_Old','infectivity_Retro','infectivity_Bead','infectivity_Old','RT_Bead','RT_Retro','RT_Old')],redoInfect2[,c('id','redoBead','redoRetro','redoFirst','redoBeadT12','redoRetroT12','redoFirstT12')],all.x=TRUE,all.y=TRUE)
all<-merge(firstExp[,c('id','ng.RT.uL','class')],all,all.y=TRUE)
colnames(all)[colnames(all)=='ng.RT.uL']<-'RT_First'
all$redoBeadInf<-all$redoBead/all$RT_Bead/1000
all$redoFirstInf<-all$redoFirst/all$RT_First/1000
all$inf<-ifelse(is.na(all$redoBeadInf),all$redoFirstInf,all$redoBeadInf)
all$inf[all$inf==Inf]<-NA
all$redoBeadInfT12<-all$redoBeadT12/all$RT_Bead/1000
all$redoFirstInfT12<-all$redoFirstT12/all$RT_First/1000
all$infT12<-ifelse(is.na(all$redoBeadInfT12),all$redoFirstInfT12,all$redoBeadInfT12)
all$infT12[all$infT12==Inf]<-NA
#
pdf('out/infectivityByClass_20191104.pdf',width=6,height=4)
  for(t12 in c(FALSE,TRUE)){
    col<-ifelse(t12,'infT12','inf')
    tmp<-all[!is.na(all[,col])&!is.na(all$class),]
    if(repCap)tmp<-tmp[tmp$id %in% rep$V2,]
    uniqC<-unique(tmp$class)
    pos<-structure(1:length(unique(tmp$class)),.Names=uniqC[order(!grepl('MM',uniqC),!grepl('TP',uniqC),grepl('Rebound|acute',uniqC))])
    par(mar=c(3,3.5,1.2,.2),lheight=.8)
    plot(1,1,type='n',las=1,ylab='Infectivity (IU/pg RT)',xlim=c(1,length(unique(tmp$class))+.2)+c(-.3,.3),ylim=range(tmp[,col]),xaxt='n',log='y',yaxt='n',mgp=c(2.3,1,.3),xlab='',main=sprintf('%s T12',ifelse(t12,'With','Without')))
    dnar::logAxis(las=1)
    axis(1,pos,sub(' ','\n',names(pos)),padj=1,mgp=c(1,.1,0))
    offset<-ave(tmp[,col],tmp$class,FUN=function(xx)beeswarm::swarmx(0,xx,cex=1.2)[[1]])
    points(pos[(tmp$class)]+offset,tmp[,col],cex=1.2,pch=21,bg='#00000055')
  }
dev.off()

rawRT<-read.csv('data/Mini expansion RT math-19.31.10.csv',stringsAsFactors=FALSE)
rawIds<-read.csv('data/Mini expansion IDs-19.31.10.csv',stringsAsFactors=FALSE)
rawIds<-rawIds[rawIds[,2]!='H',]
rawIds[c(7,12:14),3:14]<-paste(unlist(rawIds[c(7,12:14),3:14]),'retro')
newRt<-data.frame('id'=unlist(rawIds[,3:14]),'rawRT'=unlist(rawRT[,4:15]),stringsAsFactors=FALSE)
#adjust for 40x standard dil and 25x sample dil
newRt$rt<-as.numeric(sub('[<>]','',newRt$rawRT))/40*25
newRt$isRetro<-grepl(' retro$',newRt$id)
newRt$id<-sub(' retro$','',newRt$id)
newRt<-newRt[!grepl('SG3 \\+ WEAU',newRt$id)&newRt$id!='',]
all$redoRTBead<-sapply(all$id,function(xx)if(any(selector<-newRt$id==xx&!newRt$isRetro))newRt[selector,'rt'] else NA)
all$redoRTRetro<-sapply(all$id,function(xx)if(any(selector<-newRt$id==xx&newRt$isRetro))newRt[selector,'rt'] else NA)

all$inf2<-ifelse(is.na(all$redoBead),all$redoFirst/all$redoRTBead,all$redoBead/all$redoRTBead)/1000
all$inf2[all$inf2==Inf]<-NA
all$infT122<-ifelse(is.na(all$redoBead),all$redoFirstT12/all$redoRTBead,all$redoBeadT12/all$redoRTBead)/1000
all$infT122[all$infT122==Inf]<-NA

pdf('out/infectivityByClass_20191104_newRT.pdf',width=6,height=4)
  for(t12 in c(FALSE,TRUE)){
    col<-ifelse(t12,'infT122','inf2')
    tmp<-all[!is.na(all[,col])&!is.na(all$class),]
    if(repCap)tmp<-tmp[tmp$id %in% rep$V2,]
    uniqC<-unique(tmp$class)
    pos<-structure(1:length(unique(tmp$class)),.Names=uniqC[order(!grepl('MM',uniqC),!grepl('TP',uniqC),grepl('Rebound|acute',uniqC))])
    par(mar=c(3,3.5,1.2,.2),lheight=.8)
    plot(1,1,type='n',las=1,ylab='Infectivity (IU/pg RT)',xlim=c(1,length(unique(tmp$class))+.2)+c(-.3,.3),ylim=range(tmp[,col]),xaxt='n',log='y',yaxt='n',mgp=c(2.3,1,.3),xlab='',main=sprintf('%s T12',ifelse(t12,'With','Without')))
    dnar::logAxis(las=1)
    axis(1,pos,sub(' ','\n',names(pos)),padj=1,mgp=c(1,.1,0))
    offset<-ave(tmp[,col],tmp$class,FUN=function(xx)beeswarm::swarmx(0,xx,cex=1.2)[[1]])
    points(pos[(tmp$class)]+offset,tmp[,col],cex=1.2,pch=21,bg='#00000055')
  }
dev.off()

tmp<-all[!is.na(all$redoRTBead),]
tmp$oldRt<-ifelse(is.na(tmp$RT_Bead),tmp$RT_Old,tmp$RT_Bead)
tmp$newRt<-tmp$redoRTBead
tmp[,c('oldRt','newRt')]
pdf('out/rtCompare.pdf')
plot(tmp$oldRt+.0001,tmp$newRt+.0001,xlab='Old RT measurement',ylab='New RT measurement',log='xy',xaxt='n',yaxt='n')
logAxis(las=1)
logAxis(1)
abline(0,1)
dev.off()





