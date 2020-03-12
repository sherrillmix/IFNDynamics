library(dnar)
miniIn<-read.csv('out/miniInputTitration.csv',row.names=1)
miniIn$id<-rownames(miniIn)
miniOut<-read.csv('out/miniOutputTitration.csv',row.names=1)
miniOut$id<-sub(' retro','',rownames(miniOut))
miniOut$retro<-grepl('retro',rownames(miniOut))

combo<-merge(miniOut[miniOut$retro,c('id','No.Drug')],miniOut[!miniOut$retro,c('id','No.Drug')],by.y='id',by.x='id',all.x=TRUE,all.y=TRUE,suffixes=c('.retro',''))
plot(combo$No.Drug,combo$No.Drug.retro);abline(0,1)
combo<-merge(combo,miniIn[,c('No.Drug','id')],by.y='id',by.x='id',suffixes=c('','.in'),all.x=TRUE,all.y=TRUE)
rt<-read.csv('data/Mini_expansion_RT_20191106.csv',row.names=1)
rt<-rt[,!grepl('^Input',colnames(rt))]
colnames(rt)<-sub('\\.+$','',sub('RT\\.ng\\.ml\\.\\.','',colnames(rt)))
rownames(rt)<-sub('^Rb','Rebound',rownames(rt))
colnames(rt)<-sprintf('RT_%s',colnames(rt))
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
otherTrop<-read.csv('out/tropism_20191031.csv',row.names=1)
colnames(otherTrop)[colnames(otherTrop)=='No.drug']<-'No.Drug'
vir<-read.csv('data/tropsimIds.csv',stringsAsFactors=FALSE) #from Tropisms and Infect. IDs email
otherTrop$isIsolate<-!rownames(otherTrop) %in% vir[grepl('[FGH]',vir$Wells),'IDs']
inTrop<-read.csv('out/miniInputTitration.csv',row.names=1)
inTrop$isIsolate<-TRUE
inTrop<-rbind(inTrop,otherTrop[,colnames(inTrop)])
inTrop<-inTrop[!rownames(inTrop) %in% c('YU2','Media','Media1','SG3','89.6','SG3+WEAU'),]
inTropAll<-inTrop
inTrop<-inTrop[inTrop$No.Drug>2,]
inTropProp<-inTrop[,colnames(inTrop)!='No.Drug']/inTrop$No.Drug*100
#outTrop<-read.csv('out/miniOutputTitration.csv',row.names=1)
#outTropProp<-outTrop[,colnames(outTrop)!='No.Drug']/outTrop$No.Drug
colnames(inTropProp)[colnames(inTropProp)=='AMD.Mar']<-'AMD-Mar'
#inTropProp<-inTropProp[!rownames(inTropProp) %in% c('YU2','Media','Media1','SG3','89.6','SG3+WEAU'),]
out[,colnames(inTropProp)]<-inTropProp[out$ID,]
leftovers<-inTropProp[!rownames(inTropProp) %in% out$ID,]
leftovers$ID<-rownames(leftovers)
leftovers[,colnames(out)[!colnames(out) %in% colnames(leftovers)]]<-NA
out<-rbind(out,leftovers[,colnames(out)])
xx<-read.csv('rebound/out/rebQvoaMaster.csv',stringsAsFactors=FALSE)
rownames(xx)<-xx$ID
out[is.na(out[,1]),1]<-xx[out[is.na(out[,1]),"ID"],'Type']
tmp<-out[out$ID %in% xx$ID,]
tmp<-tmp[order(tmp$ID),]
write.csv(out,'out/mini_repInfTrop.csv',row.names=FALSE)
pdf('test.pdf',height=20)
par(mar=c(2,14,1,1))
image(1:3,1:nrow(tmp),t(as.matrix(tmp[,c('AMD','Mar','AMD-Mar')])),col=rev(heat.colors(100)),xaxt='n',xlab='',yaxt='n',ylab='')
axis(1,1:3,c('AMD','Mar','AMD-Mar'))
axis(2,1:nrow(tmp),tmp$ID,las=1,cex=.1)
box()
dev.off()


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





xx<-read.csv('data/Data Master 2019 _ReboundAndQVOA_20191219.csv',stringsAsFactors=FALSE)
rownames(inTropAll)[!rownames(inTropAll) %in% xx$ID]
xx[!xx$ID %in% rownames(inTropAll),'ID']
out2<-cbind(xx,inTropAll[xx$ID,])
out2$id<-sub('^Rebound P','Rebound A09 P',out2$ID)
rownames(out2)<-NULL
write.csv(out2,'out/ReboundAndQVOA_20191219.csv')
bigTrop<-read.csv('out/20191223_tropism.csv',stringsAsFactors=FALSE,row.names=1)
bigTrop<-bigTrop[!grepl('Media|IMC',rownames(bigTrop))&!grepl('^RM',rownames(bigTrop)),]
rownames(bigTrop)<-sub('Reb. ','',sub('Rebound 50','A09',rownames(bigTrop)))
rownames(bigTrop)<-sub(' D-14 | W12 ',' ',rownames(bigTrop))
matches<-lapply(structure(rownames(bigTrop),.Names=rownames(bigTrop)),function(xx){
  regex<-sprintf('%s[ _-]',gsub('[ -]+','[ ._-].*',trimws(sub('\\([^)]+\\)','',xx))))
  hits<-out2$ID[grep(regex,sprintf('%s ',out2$id))]
  hits
})
if(any(sapply(matches,length)>1))stop('Multiple ID matches')
rownames(bigTrop)[(sapply(matches,length)==0)&!grepl('- (UT|BE) -',rownames(bigTrop))]
message('Missing IDs:')
dummy<-sapply(out2[!out2$ID %in% unlist(matches),'ID'],message)
rownames(bigTrop)[sapply(matches,length)==1]<-unlist(matches)
colnames(bigTrop)<-sprintf('%s_1218',colnames(bigTrop))
out2<-cbind(out2,bigTrop[out2$ID,])
virus<-read.csv('ice/2019-12-18_tropism.csv',stringsAsFactors=FALSE)
virus[,3:10]<-sub('Reb. ','',sub('Rebound 50','A09',unlist(virus[,3:10])))
virus[,3:10]<-sub(' D-14 | W12 ',' ',unlist(virus[,3:10]))
p1<-unlist(virus[virus$Plate==1&virus$Row %in% 1:12,3:10])
p2_10<-unlist(virus[virus$Plate==2&virus$Row %in% 9:9,3:10])
cols<-ifelse(out2$ID %in% unlist(matches)[names(unlist(matches)) %in% p1],'red',
  ifelse(out2$ID %in% unlist(matches)[names(unlist(matches)) %in% p2_10],'blue','black'))
cols<-ifelse(out2$ID %in% rownames(otherTrop),'red','black')
pdf('test.pdf');plot(out2$No.Drug+1,out2$No.Drug_1218+1,log='xy',col=out2$isIsolate+1);abline(0,1);dev.off()
non31<-out2[out2$Type!='Rebound','ID'][!out2[out2$Type!='Rebound','ID'] %in% rownames(otherTrop)]
plot(out2$AMD+1,out2$AMD_1218+1,log='xy',col=cols,cex=1+(out2$ID %in% non31));abline(0,1)

odd<-out2[out2$No.Drug_1218/out2$No.Drug>10&!is.na(out2$No.Drug+out2$No.Drug_1218),'ID']
tmpName<-names(unlist(matches))[unlist(matches) %in% odd]
write.csv(tmpName,'tmp.csv')

out2$inf_1218<-out2$No.Drug_1218/out2$RT/1000




# collect all tzmbl into one file. Need: master file, RT isolate, ND, A, M, AM isolate, RT_MiniBead, ND... MiniBead, RT_MiniRetro, ND... MiniRetro
xx<-read.csv('data/Data Master 2019 _ReboundAndQVOA_20191219.csv',stringsAsFactors=FALSE)
yy<-read.csv("data/RT for Fred's Isolates_20200109.csv",stringsAsFactors=FALSE)
colnames(yy)<-c('ID','RT')
yy[,colnames(xx)[!colnames(xx) %in% colnames(yy)]]<-NA
collect<-rbind(xx,yy[,colnames(xx)])
colnames(collect)[colnames(collect)=='RT']<-'RT_Isolate'
#include large tropism run of original isolates
bigTrop2<-bigTrop
#fix order on "- BE -" virus and remove QVOA_
rownames(bigTrop2)<-sub('^QVOA_','',sub('(.*) - ([A-Z]+) - (.*)','\\1 \\3 \\2',rownames(bigTrop2)))
matches<-lapply(structure(rownames(bigTrop2),.Names=rownames(bigTrop2)),function(xx){
  regex<-sprintf('%s[ _-]',gsub('[ -]+','[ ._-].*',trimws(sub('\\([^)]+\\)','',xx))))
  collect$ID[grep(regex,sprintf('%s ',collect$ID))]
})
if(any(sapply(matches,length)>1))stop('Duplicate ID found')
tmp<-collect[1,]
tmp[,]<-NA
tmp<-tmp[rep(1,sum(sapply(matches,length)==0)),]
tmp$ID<-names(matches)[sapply(matches,length)==0]

collect<-rbind(collect,tmp)
matches<-lapply(structure(rownames(bigTrop2),.Names=rownames(bigTrop2)),function(xx){
  regex<-sprintf('%s[ _-]',gsub('[ -]+','[ ._-].*',trimws(sub('\\([^)]+\\)','',xx))))
  collect$ID[grep(regex,sprintf('%s ',collect$ID))]
})
if(any(sapply(matches,length)!=1))stop('Problem matching big tropism')
if(any(table(unlist(matches))!=1))stop('Collision in big tropism')
rownames(bigTrop2)<-unlist(matches)
colnames(bigTrop2)<-sub('_.*','_Isolate',colnames(bigTrop2))
collect<-cbind(collect,bigTrop2[collect$ID,])
#include additional tropism from "all" collected
tmp<-all[all$id %in% collect$ID,]
tmp$RT_Bead_Combine<-apply(tmp[,c(ifelse(is.na(tmp$RT_Bead),'RT_First','RT_Bead'),'redoRTBead')],1,function(xx)mean(xx,na.rm=TRUE))
tmp$RT_Retro_Combine<-apply(tmp[,c('RT_Retro','redoRTRetro')],1,function(xx)mean(xx,na.rm=TRUE))
tmp$IU_Bead_Combine<-apply(tmp[,c('IU_Bead',ifelse(is.na(tmp$RT_Bead),'redoFirst','redoBead'))],1,function(xx)mean(xx,na.rm=TRUE))
tmp$isFirstBead<-is.na(tmp$RT_Bead)
tmp$IU_Retro_Combine<-apply(tmp[,c('IU_Retro','redoRetro')],1,function(xx)mean(xx,na.rm=TRUE))
tmp$IU_Retro_T122<-tmp$redoRetroT12
tmp$IU_Bead_T122<-ifelse(is.na(tmp$RT_Bead),tmp$redoFirstT12,tmp$redoBeadT12)
tmp$IU_Isolate_old<-tmp[,'IU_Old']
rownames(tmp)<-tmp$id
tmp2<-tmp[,c('RT_Bead_Combine','RT_Retro_Combine','IU_Bead_Combine','IU_Retro_Combine','IU_Retro_T122','IU_Bead_T122','IU_Isolate_old','isFirstBead')]
colnames(tmp2)<-sub('_Combine','',colnames(tmp2))
collect<-cbind(collect,tmp2[collect$ID,])
collect$IU_Isolate<-collect$No.Drug_Isolate
#note these are not all isolates (!isIsolate = first miniExpansion)
tmpIn<-inTropAll[,1:5]
colnames(tmpIn)<-sprintf('%s_old',colnames(tmpIn))
collect<-cbind(collect,tmpIn[collect$ID,])
collect$Type[is.na(collect$Type)&grepl('Reb\\.',collect$ID)]<-'Rebound'
collect$Type[is.na(collect$Type)&grepl('QVOA\\.post',collect$ID)]<-'VOA - Post ATI'
collect$Type[is.na(collect$Type)&grepl('voa.*\\(pre\\)',collect$ID)]<-'VOA - Pre ATI'
collect$Type[is.na(collect$Type)&grepl('A08.*\\(pre\\)',collect$ID)]<-'VOA - Pre ATI'
collect$Type[is.na(collect$Type)&grepl('A08.*\\(post\\)',collect$ID)]<-'VOA - Post ATI'
collect$Type[is.na(collect$Type)&grepl('818B.*\\(pre\\)',collect$ID)]<-'VOA - Pre ATI'
reboundIsolates<-c("PIP-017 (IFN-a2b) - 8F", "S-22 7D4", "BEAT-030 (IFN-a2b) - 5E", "PIP-008 (IFN-a2b) - 2A", "PIP-017 (IFN-a2b) - 8B", "BEAT-030 (IFN-a2b) - 7E", "S-07 (R-3) - 2A")
collect$Type[is.na(collect$Type)&collect$ID %in% reboundIsolates]<-'Rebound'
collect$Type[collect$ID=='9201 12F3']<-'VOA - Pre ATI'
collect$ID[is.na(collect$Type)]
collect$Patient[is.na(collect$Patient)&grepl('A08[. -]',collect$ID)]<-'A08'
collect$Patient[is.na(collect$Patient)&grepl('A01[. -]',collect$ID)]<-'A01'
collect$Patient[is.na(collect$Patient)&grepl('A09[. -]',collect$ID)]<-'A09'
collect$Patient[is.na(collect$Patient)&grepl('601[r. -]',collect$ID)]<-'601'
collect$Patient[is.na(collect$Patient)&grepl('9201[. _-]',collect$ID)]<-'9201'
collect$Patient[is.na(collect$Patient)&grepl('818B[. _-]',collect$ID)]<-'A09'
collect$Patient[is.na(collect$Patient)]<-sub('[_ ].*','',collect$ID[is.na(collect$Patient)])
collect$Study<-trimws(collect$Study)
collect$Study[is.na(collect$Study)&grepl('^A0[189]$',collect$Patient)]<-'VRC01'
collect$Study[is.na(collect$Study)&collect$Patient=='9201']<-'3BNC117 / 10-1074'
collect$Study[is.na(collect$Study)&collect$Patient=='601']<-'3BNC117'
collect$Study[is.na(collect$Study)&grepl('^S-[0-9]+$',collect$Patient)]<-'INTERRUPT'
collect$Study[is.na(collect$Study)&grepl('^BEAT-[0-9]+$',collect$Patient)]<-'BEAT'
collect$Study[is.na(collect$Study)&grepl('^PIP-[0-9]+$',collect$Patient)]<-'PIP'
collect$ID[is.na(collect$Study)]
collect<-collect[!grepl('293T',collect$ID),]
collect<-collect[collect$ID!='Reb. 9201-r1_UT_macro',]
collect$infectivity<-collect$IU_Isolate/collect$RT_Isolate/1000
write.csv(collect[,!colnames(collect) %in% c('simple','id')],'out/qvoaReboundCollected_2020110.csv',row.names=FALSE)
write.csv(collect[,c('Study','Type','Patient','ID','RT_Isolate','IU_Isolate')],'out/qvoaReboundIU_2020110.csv',row.names=FALSE)
picks<-collect[collect$RT_Isolate>.01 & collect$IU_Isolate>200&!is.na(collect$RT_Isolate)&!is.na(collect$IU_Isolate)&!is.na(collect$Type),c('ID','Patient','Type','IU_Isolate','RT_Isolate')]
#print(picks,row.names=FALSE)
#print(nrow(picks))
write.csv(picks,'highIuPicks.csv',row.names=FALSE)

newIc50<-read.csv('out/fred_20200224_ic50.csv',stringsAsFactors=FALSE)
newIc50<-newIc50[!grepl('MM33|PIP',newIc50$sample),]
newIc50$repCap<-(newIc50$repCap.IFNa2+newIc50$repCap.IFNb)/2
newIc50[newIc50$max.IFNa2<300,c('ic50.IFNa2','ic50.IFNb','vres.IFNa2','vres.IFNb')]<-NA
newIc50[newIc50$max.IFNa2<500,c('vres.IFNa2','vres.IFNb')]<-NA
rownames(newIc50)<-newIc50$sample
if(any(!rownames(newIc50) %in% collect$ID))stop('Mismatched ID')
#newIc50[!newIc50$sample %in% collect$ID,'sample']
collect[collect$ID %in% rownames(newIc50),c('ic50_IFNa2','ic50_IFNb','vres_IFNa2','vres_IFNb','repCap')]<-newIc50[collect[collect$ID %in% rownames(newIc50),'ID'],c('ic50.IFNa2','ic50.IFNb','vres.IFNa2','vres.IFNb','repCap')] 
neut<-read.csv('data/neuts infect 01142020 sCD4 rebound results table.csv',stringsAsFactors=FALSE)
neut<-neut[neut$ID!='',!grepl('^X\\.[0-9]+$',colnames(neut))]
if(any(!neut$ID %in% collect$ID & neut$Patient !='293T stock'))stop('Missing ID')
collect<-merge(collect,neut[,c('ID','sCD4.IC50','CD4.Ig.IC50')],by.x='ID',by.y='ID',all.x=TRUE)
collect$tropism<-apply(collect[,c('No.Drug_Isolate','AMD_Isolate','Mar_Isolate','AMD.Mar_Isolate')],1,function(xx){
  minNoDrug<-5
  threshold<-20
  if(any(is.na(xx)))return(NA)
  if(xx[1]<minNoDrug)return(NA)
  perc<-xx[2:4]/xx[1]*100
  if(perc[3]>threshold)warning('Problem AMD+Mar virus')
  if(all(perc[1:2]<threshold))stop('Problem virus')
  if(all(perc[1:2]>threshold))return('Dual')
  else return(c('R5','X4')[perc[1:2]>threshold])
})

write.csv(collect[,!colnames(collect) %in% c('simple','id')],'out/qvoaReboundCollected_20200224.csv',row.names=FALSE)


s3<-read.csv('data/Table S3.csv',stringsAsFactors=FALSE)
s3<-s3[!s3$Virus.ID %in% c('','Patient ID'),]
if(any(!s3$Virus.ID %in% collect$ID))stop('Missing ID')
s3out<-merge(s3,collect[,c('ID','RT_Isolate','IU_Isolate','AMD_Isolate','Mar_Isolate','AMD.Mar_Isolate','ic50_IFNa2','vres_IFNa2','ic50_IFNb','vres_IFNb','repCap')],by.x='Virus.ID',by.y='ID',all.x=TRUE,all.y=TRUE)
neut<-read.csv('data/neuts infect 01142020 sCD4 rebound results table.csv',stringsAsFactors=FALSE)
neut<-neut[neut$ID!='',!grepl('^X\\.[0-9]+$',colnames(neut))]
if(any(!neut$ID %in% s3out$Virus.ID & neut$Patient !='293T stock'))stop('Missing ID')
s3out<-merge(s3out,neut[,c('ID','sCD4.IC50','CD4.Ig.IC50')],by.x='Virus.ID',by.y='ID',all.x=TRUE)
s3out<-s3out[orderIn(s3out$Virus.ID,s3$Virus.ID),]
s3out$infectivity<-s3out$IU_Isolate/s3out$RT_Isolate/1000
s3out$AMD_percent<-s3out$AMD_Isolate/s3out$IU_Isolate*100
s3out$Mar_percent<-s3out$Mar_Isolate/s3out$IU_Isolate*100
s3out$AMDMar_percent<-s3out$AMD.Mar_Isolate/s3out$IU_Isolate*100
s3out[s3out$IU_Isolate<1&!is.na(s3out$IU_Isolate),c('AMD_percent','Mar_percent','AMDMar_percent')]<-NA
bak<-s3out
s3out<-s3out[,!colnames(s3out) %in% c('IU_Isolate','RT_Isolate','AMD_Isolate','Mar_Isolate','AMD.Mar_Isolate')]
if(any(!s3$Virus.ID %in% s3out$Virus.ID))stop('Lost virus')
if(any(s3out$Virus.ID[1:nrow(s3)]!=s3$Virus.ID))stop('Messed up order')
s3out$tropism<-apply(s3out[,c('AMD_percent','Mar_percent','AMDMar_percent')],1,function(xx){
  threshold<-20
  if(any(is.na(xx)))return(NA)
  if(xx[3]>threshold)warning('Problem AMD+Mar virus')
  if(all(xx[1:2]>threshold))return('Dual')
  else return(c('R5','X4')[xx[1:2]>threshold])
})
write.csv(s3out,'out/qvoaReboundIUNeut_2020124.csv',row.names=FALSE)

pdf('out/qvoaReboundTropism.pdf',height=12,width=3.5)
par(mar=c(2.5,6,.1,1.5))
cols<-rev(heat.colors(100))
breaks<-seq(0,100,1)
thisMat<-t(as.matrix(s3out[,c('AMD_percent','Mar_percent','AMDMar_percent')]))
thisMat[thisMat>100]<-100
image(1:3,1:nrow(s3out),thisMat,xaxt='n',yaxt='n',xlab='',ylab='',col=cols,breaks=)
axis(2,1:nrow(s3out),s3out$Virus.ID,las=1,cex.axis=.5,mgp=c(3,.4,0),tcl=-.3)
axis(1,1:3,c('AMD','Mar','AMD+Mar'),mgp=c(3,.4,0),tcl=-.3)
nas<-which(is.na(t(as.matrix(s3out[,c('AMD_percent','Mar_percent','AMDMar_percent')]))),arr.ind=TRUE)
rect(nas[,1]-.5,nas[,2]-.5,nas[,1]+.5,nas[,2]+.5,col='grey',border=NA)
axis(4,1:nrow(s3out),s3out$tropism,las=1,cex.axis=.5,mgp=c(3,.4,0),tcl=-.3)
box()
par(lheight=.7)
dnar::insetScale(breaks=breaks,col=cols,at=c(0,50,100),labels=c(0,50,'>100'),insetPos=c(.008,.035,.015,.25),main='Percent of untreated\ninfectivity',cex=.7)
dev.off()



#inTropAll and miniOut have additional tropism if needed

pdf('out/mini_vs_isolate_infectivity.pdf',width=4,height=4)
par(mar=c(2,4,1,0.1))
thisDat<-collect[!is.na(collect$Type)&!is.na(collect$RT_Bead)&!is.na(collect$IU_Bead),]
thisDat$type<-sub(' .*','',thisDat$Type)
pos<-structure(1:length(unique(thisDat$type)),.Names=unique(thisDat$type))
xSpread<-vipor::offsetX(log10(thisDat$IU_Bead/thisDat$RT_Bead/1000),thisDat$type)
tmp<-c(collect$IU_Bead/collect$RT_Bead/1000,collect$IU_Isolate/collect$RT_Isolate/1000)
ylim<-range(tmp[tmp<Inf],na.rm=TRUE)
plot(pos[thisDat$type]+xSpread,thisDat$IU_Bead/thisDat$RT_Bead/1000,log='y',main='Miniexpansion',ylab='Infectivity',yaxt='n',xaxt='n',xlab='',ylim=ylim)
axis(1,pos,names(pos))
logAxis(las=1)
thisDat<-collect[!is.na(collect$Type)&!is.na(collect$RT_Isolate)&!is.na(collect$IU_Isolate),]
thisDat$type<-sub(' .*','',thisDat$Type)
thisDat$inf<-thisDat$IU_Isolate/thisDat$RT_Isolate/1000
pos<-structure(1:length(unique(thisDat$type)),.Names=unique(thisDat$type))
xSpread<-vipor::offsetX(log10(thisDat$inf),thisDat$type)
plot(pos[thisDat$type]+xSpread,thisDat$inf,log='y',main='Isolates',ylab='Infectivity',yaxt='n',xaxt='n',xlab='',ylim=ylim)
logAxis(las=1)
axis(1,pos,names(pos))
dev.off()
t.test(log(thisDat$inf[thisDat$type=='Rebound'&!is.na(thisDat$type)&thisDat$inf<Inf]),log(thisDat$inf[thisDat$type=='VOA'&!is.na(thisDat$type)&thisDat$inf<Inf]))

pdf('out/IU_problem.pdf')
plot(collect$No.Drug_old,collect$IU_Isolate,log='xy',col=(!collect$isIsolate_old)+1,ylab='12-24 Isolate IU/ul',xlab='Previous "isolate" IU/ul');abline(0,1)
legend('bottomright',inset=.01,c('Rows FGH','Other'),col=c('red','black'),pch=1)
plot(collect$No.Drug_old,ifelse(collect$isIsolate_old,collect$IU_Isolate,collect$IU_Bead),log='xy',ylab='12-24 Isolate (or previously measured miniexpansion measurement) IU/ul',xlab='Previous "isolate" (or miniexpansion) IU/ul',col=(!collect$isIsolate_old)+1);abline(0,1)
legend('bottomright',inset=.01,c('Compared to bead mini','Compared to isolate'),col=c('red','black'),pch=1)
#calculate infectivity for the miniexpansion and isolates here
plot(collect$IU_Isolate,collect$IU_Bead,log='xy',ylab='Miniexpansion IU/ul',xlab='Original isolate IU/ul',col=ifelse(collect$isFirstBead,'blue','black'));abline(0,1)
points(collect$IU_Isolate,collect$IU_Retro,col='red')
legend('topleft',inset=.01,c('Bead first','Bead second','Retro'),col=c('blue','black','red'),pch=1)
plot(collect$IU_Isolate/collect$RT_Isolate/1000,collect$IU_Bead/collect$RT_Bead/1000,log='xy',ylab='Miniexpansion infectivity',xlab='Original isolate infectivity',col=ifelse(collect$isFirstBead,'blue','black'));abline(0,1)
points(collect$IU_Isolate/collect$RT_Isolate/1000,collect$IU_Retro/collect$RT_Retro/1000,col='red')
legend('topleft',inset=.01,c('Bead first','Bead second','Retro'),col=c('blue','black','red'),pch=1)
dev.off()

pdf('out/Inf_problem.pdf')
#calculate infectivity for the miniexpansion and isolates here
plot(collect$IU_Isolate,collect$IU_Bead,log='xy',ylab='Miniexpansion IU/ul',xlab='Original isolate IU/ul',bg=ifelse(collect$isFirstBead,'blue','red'),pch=21);abline(0,1)
legend('topleft',inset=.01,c('Bead first','Bead second'),pt.bg=c('blue','red'),pch=21)
plot(collect$IU_Isolate/collect$RT_Isolate/1000,collect$IU_Bead/collect$RT_Bead/1000,log='xy',ylab='Miniexpansion infectivity',xlab='Original isolate infectivity',bg=ifelse(collect$isFirstBead,'blue','red'),pch=21);abline(0,1)
legend('topleft',inset=.01,c('Bead first','Bead second'),pt.bg=c('blue','red'),pch=21)
dev.off()


pdf('test.pdf');plot(collect$IU_Isolate,collect$IU_Bead,log='xy');abline(0,1);dev.off()
pdf('test.pdf');plot(collect$No.Drug_Isolate,collect$IU_Isolate_old);abline(0,1);dev.off()
#pdf('test.pdf');plot(tmp$redoRetro,tmp$IU_Retro);abline(0,1);dev.off()
#pdf('test.pdf');plot(tmp$redoBead,tmp$IU_Bead);abline(0,1);dev.off()

s3<-read.csv('data/Table S3.nearfinal_20200226.csv',skip=1,stringsAsFactors=FALSE)
s3<-s3[s3$Isolate.ID!='',colnames(s3)!='X']
s3out<-merge(s3,collect[!is.na(collect$Type)&!is.na(collect$Patient),c('ID','Study','Patient','Type','repCap','ic50_IFNa2','ic50_IFNb','infectivity','sCD4.IC50','CD4.Ig.IC50','tropism')],by.x='Isolate.ID',by.y='ID',all.x=TRUE,all.y=TRUE,suffixes=c('','.TMP'))
s3out[is.na(s3out$Type),'Type']<-s3out[is.na(s3out$Type),'Type.TMP']
rownames(s3out)<-s3out$Isolate.ID
s3out<-s3out[c(s3$Isolate.ID,rownames(s3out)[!rownames(s3out) %in% s3$Isolate.ID]),]
s3out$Patient<-sub(' \\(.*','',s3out$Patient)
s3out$simpleType<-c('Rebound'='REBOUND','VOA'='VOA','VOA - Pre ATI'='VOA_PRE','VOA - Post ATI'='VOA_POST')[s3out$Type]
s3out$identifier<-sub('^[PM]([0-9])','\\1',sub('.*[ _-]','',trimws(sub('_BE$|_UT$','',sub(' \\(?pre.*| post.*','',sub(' bulk','',sub(' *pre-ATI\\(wk.*','',s3out$Isolate.ID)))))))
s3out$pubID<-sprintf('%s.%s.%s',s3out$Patient,s3out$simpleType,s3out$identifier)
if(any(table(s3out$pubID)>1))stop('Duplicate name created')
fillFunc<-function(df,col1,col2){selector<-is.na(df[,col1])&!is.na(df[,col2]);df[selector,col1]<-df[selector,col2];df}
s3out<-fillFunc(s3out,'Replicative.Capacity..ng.p24.ml.','repCap')
s3out<-fillFunc(s3out,'IFN.a2.IC50..pg.ml.','ic50_IFNa2')
s3out<-fillFunc(s3out,'IFN.b.IC50..pg.ml.','ic50_IFNb')
s3out<-fillFunc(s3out,'Infectivity..IU.pg.RT.','infectivity')
s3out<-fillFunc(s3out,'sCD4.IC50..ug.ml.','sCD4.IC50')
s3out<-fillFunc(s3out,'CD4.Ig.IC50..ug.ml.','CD4.Ig.IC50')
s3out<-fillFunc(s3out,'Isolate.Tropism','tropism')
s3out$tropism==s3out$Isolate.Tropism
write.csv(s3out[,!colnames(s3out) %in% c('simpleType','identifier')],'out/S3Update_20200226.csv',row.names=FALSE)


alias92<-c('9241'='9207', '9244'='9203', '9243'='9202', '9242'='9201', '601'='601')
alias92Rev<-structure(names(alias92),.Names=alias92)
collect$alias<-ifelse(collect$Patient %in% names(alias92Rev),alias92Rev[collect$Patient],collect$Patient)
reg<-'^.*wk([0-9-]+).*$'
collect$type2<-ifelse(collect$Patient %in% alias92&grepl(reg,collect$ID),sprintf('%s (week %s)',collect$Type,sub(reg,'\\1',collect$ID)),collect$Type)
lapply(unique(collect$Study),function(xx)withAs(collect=collect[collect$Study==xx&!is.na(collect$ic50_IFNb),],sum(table(collect$alias,collect$type2,collect$Study))))
lapply(unique(collect$Study),function(xx)withAs(this=collect[collect$Study==xx&!is.na(collect$ic50_IFNb),],(table(this$alias,factor(this$type2,unique(collect$type2)),this$Study)))) # 
levs<-sort(unique(collect$type2))
out<-do.call(rbind,lapply(unique(collect$alias[!is.na(collect$ic50_IFNb)]),function(xx)withAs(this=collect[collect$alias==xx&!is.na(collect$ic50_IFNb),],do.call(data.frame,c('study'=this$Study[1],'pat'=this$alias[1],as.list(table(factor(this$type2,levs))),stringsAsFactors=FALSE)))))
colnames(out)[-1:-2]<-levs
out<-out[order(out$study,out$pat),]
write.csv(out,'out/isolatesWithIFNb.csv',row.names=FALSE)


p24<-read.csv('data/p24_01.17.2020.csv',row.names=1,stringsAsFactors=FALSE)
p24[,1:24]<-as.numeric(sub('[<>]','',unlist(p24[,1:24])))
p24[9:16,c(22,24)]<-NA
means<-t(apply(p24[,1:24],1,tapply,rep(1:12,each=2),mean,na.rm=TRUE))
p24_2<-read.csv('data/p24_01.22.2020.csv',row.names=1,stringsAsFactors=FALSE)
p24_2[,1:24]<-as.numeric(sub('[<>]','',unlist(p24_2[,1:24])))
p24_2<-p24_2[1:16,1:24]
means2<-t(apply(p24_2[,1:24],1,tapply,rep(1:12,each=2),mean,na.rm=TRUE))
virus<-read.csv('ice/2019-12-18_tropism.csv',stringsAsFactors=FALSE)
vir<-rbind(t(virus[1:12,10:3]),t(virus[13:24,10:3])) 
#miniP24<-data.frame('virus'=as.vector(vir),'p24'=as.vector(means),row.names=1:192,stringsAsFactors=FALSE)
miniP24<-data.frame('virus'=as.vector(vir),'p24'=as.vector(means),'p24_2'=as.vector(means2),row.names=1:192,stringsAsFactors=FALSE)
rownames(collect)<-collect$ID
miniP24$virus<-sub('Reb. ','',sub('Rebound 50','A09',miniP24$virus))
miniP24$virus<-sub(' D-14 | W12 ',' ',miniP24$virus)
matches<-lapply(structure(miniP24$virus,.Names=miniP24$virus),function(xx){
  regex<-sprintf('%s[ _-]',gsub('[ _-]+','[ ._-].*',sub('(UT|BE) -(.*)','\\2 \\1',trimws(sub('\\([^)]+\\)','',xx)))))
  hits<-collect$ID[grep(regex,sprintf('%s ',collect$ID))]
  hits<-hits[!grepl('macro',hits)]
  hits
})
head(matches[sapply(matches,length)!=1&!grepl('^RM',names(matches))&!grepl('IMC',names(matches))])
miniP24$convert<-unlist(matches[sapply(matches,length)==1&!grepl('^RM',names(matches))&!grepl('IMC',names(matches))])[miniP24$virus]
miniP24[,'rt']<-collect[miniP24$convert,'RT_Isolate']
miniP24[,'iu']<-collect[miniP24$convert,'IU_Isolate']

pdf('out/miniP24.pdf',width=5,height=5)
  par(mar=c(3.3,3.3,1,.1),mgp=c(2.3,.8,0))
  plot(miniP24$rt,miniP24$p24,ylab='Output miniexpansion 1/17 p24 (pg/ml)',xlab='Input isolate RT',log='y',yaxt='n',pch=21,bg='#00000022',cex=1.2)
  #logAxis(1)
  logAxis(las=1)
  plot(miniP24$iu,miniP24$p24,ylab='Output miniexpansion 1/17 p24 (pg/ml)',xlab='Input isolate IU/ul',log='xy',xaxt='n',yaxt='n',pch=21,bg='#00000022',cex=1.2)
  logAxis(1)
  logAxis(las=1)
  plot(miniP24$p24,miniP24$p24_2,log='xy',xlab='Mini expansion 1/17 p24 (pg/ml)',ylab='Mini expansion 1/22 p24 (pg/ml)',xaxt='n',yaxt='n',pch=21,bg='#00000022',cex=1.2)
  logAxis(1)
  logAxis(las=1)
  abline(0,1,lty=2)
  plot(miniP24$iu,miniP24$p24_2,ylab='Output miniexpansion 1/22 p24 (pg/ml)',xlab='Input isolate IU/ul',log='xy',xaxt='n',yaxt='n',pch=21,bg='#00000022',cex=1.2,ylim=range(miniP24$p24))
  logAxis(1)
  logAxis(las=1)
dev.off()
