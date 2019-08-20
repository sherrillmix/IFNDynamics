
if(!exists('dat'))source('readNewData.R')
new<-read.csv('out/20190725_outliers.csv',row.names=1)
new<-new[,grep('old',colnames(new))]
new$beta<-grepl('Beta$',rownames(new))
new$id<-sub(' (Alpha|Beta)$','',rownames(new))
if(any(new$id[1:(nrow(new)/2)]!=new$id[(nrow(new)/2+1):nrow(new)]))stop('Id mismatch')
combo<-data.frame('id'=new$id[1:(nrow(new)/2)],
  'ic50'=new$old_ic50[1:(nrow(new)/2)],'beta'=new$old_ic50[(nrow(new)/2+1):nrow(new)],
  'repCapAlpha'=new$old_repCap[1:(nrow(new)/2)],'repCapBeta'=new$old_repCap[(nrow(new)/2+1):nrow(new)],
  'vres'=new$old_vres[1:(nrow(new)/2)],'vresBeta'=new$old_vres[(nrow(new)/2+1):nrow(new)],
  stringsAsFactors=FALSE)
combo$id[combo$id=='MM55.03.2A1']<-'MM55.03.2A1.bulk'
combo$id[combo$id=='WEAU.14.1C1']<-'WEAU.14.1C1.bulk'
if(any(!combo$id %in% rownames(dat)))stop('Unrecognized ID')
combo$origIc50<-dat[combo$id,'ic50']
combo$origBeta<-dat[combo$id,'beta']
combo$origRepCap<-dat[combo$id,'replication']
combo$time<-dat[combo$id,'time']
combo$rawOrigIc50<-apply(dat[combo$id,ifna2_ic50],1,function(xx)paste(round(xx[!is.na(xx)],4),collapse=','))
combo$rawOrigBeta<-apply(dat[combo$id,ifnb_ic50],1,function(xx)paste(round(xx[!is.na(xx)],4),collapse=','))
pdf('out/compareRedo.pdf')
plot(combo$ic50,combo$origIc50,log='xy',xaxt='n',yaxt='n',xlab='IFNa2 IC50 (previous)',ylab='IFNa2 IC50 (redo)')
logAxis(1)
logAxis(las=1)
abline(0,1)
plot(combo$beta,combo$origBeta,log='xy',xaxt='n',yaxt='n',xlab='IFNb IC50 (previous)',ylab='IFNb IC50 (redo)')
logAxis(1)
logAxis(las=1)
abline(0,1)
dev.off()

pdf('out/compareRedo_longitudinal.pdf',height=8,width=6)
#par(mfrow=c(3,4))
par(mar=c(5.1,3.5,1,.3),mfrow=c(2,1))
for(ii in unique(dat$pat)){
  thisDat<-dat[dat$pat==ii,]
  plot(thisDat$time/7,thisDat$ic50,main=sprintf('%s alpha',ii),xlab='Weeks after symptoms',ylab='IFNa2 IC50',log='y',yaxt='n',ylim=range(dat$ic50,na.rm=TRUE),xlim=range(dat$time/7),pch=21,bg='grey',mgp=c(2.5,.8,0))
  logAxis(las=1)
  thisCombo<-combo[grep(ii,combo$id),]
  points(thisCombo$time/7+10,thisCombo$ic50,bg='red',pch=21)
  segments(thisCombo$time/7+10,thisCombo$ic50,thisCombo$time/7,thisCombo$origIc50,col='red')
  plot(thisDat$time/7,thisDat$beta,main=sprintf('%s beta',ii),xlab='Weeks after symptoms',ylab='IFNb IC50',log='y',yaxt='n',ylim=range(dat$beta,na.rm=TRUE),xlim=range(dat$time/7),pch=21,bg='grey',mgp=c(2.5,.8,0))
  logAxis(las=1)
  thisCombo<-combo[grep(ii,combo$id),]
  points(thisCombo$time/7+10,thisCombo$beta,bg='red',pch=21)
  segments(thisCombo$time/7+10,thisCombo$beta,thisCombo$time/7,thisCombo$origBeta,col='red')
}
dev.off()

pdf('out/compareRedo_repCap.pdf',width=4,height=4)
par(mar=c(3.1,4,.2,.2))
lims<-range(c(combo$repCapAlpha,combo$repCapBeta,combo$origRepCap))
plot(combo$repCapAlpha,combo$repCapBeta,ylab='New replicative capacity (IFNb)',mgp=c(2.8,.6,0),las=1,xlab='')
title(xlab='New replicative capacity (IFNa2)',mgp=c(2,0,0))
abline(0,1)
plot(combo$origRepCap,combo$repCapAlpha,ylab='New replicative capacity (IFNa2)',xlab='',xlim=lims,ylim=lims,mgp=c(2.8,.6,0),las=1)
title(xlab='Previous replicative capacity',mgp=c(2,0,0))
segments(combo$origRepCap,combo$origRepCap,combo$origRepCap,combo$repCapAlpha,col='#00000011')
abline(0,1)
plot(combo$origRepCap,combo$repCapBeta,ylab='New replicative capacity (IFNb)',xlab='',xlim=lims,ylim=lims,mgp=c(2.8,.6,0),las=1)
segments(combo$origRepCap,combo$origRepCap,combo$origRepCap,combo$repCapBeta,col='#00000011')
title(xlab='Previous replicative capacity',mgp=c(2,0,0))
abline(0,1)
dev.off()

rownames(combo)<-combo$id
out<-combo[rownames(dat),]
newNames<-c('ic50'='IFNa2 August2019 Pooled Donor IC50 (pg/ml)','vres'='IFNa2 August2019 Pooled Donor Vres at 5.5 pg/ml','beta'='IFNb August2019 Pooled Donor IC50 (pg/ml)','vresBeta'='IFNb August2019 Pooled Donor Vres at 4,400 pg/ml')
for(ii in names(newNames))colnames(out)[colnames(out)==ii]<-newNames[ii]
rownames(out)<-rownames(dat)
write.csv(out[,newNames],'out/redoIfn.csv')

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
dnar::withAs(combo=combo[!is.na(combo$infectivity_Retro)&!is.na(combo$infectivity_Bead)&combo$infectivity_Retro+combo$infectivity_Bead<Inf,],plot(combo$infectivity_Retro,combo$infectivity_Bead,ylab='Bead isolate infectivity (IU/pg RT)',pch=21,bg='#00000033',cex=1.3,mgp=c(2.3,.8,0),las=1,xlab=''))
title(xlab='Retronectin isolate infectivity (IU/pg RT)',mgp=c(1.9,1,0))
abline(0,1,lty=1)
dnar::withAs(combo=combo[!is.na(combo$infectivity_Old)&!is.na(combo$infectivity_Bead)&combo$infectivity_Old+combo$infectivity_Bead<Inf,],plot(combo$infectivity_Old,combo$infectivity_Bead,ylab='Bead isolate infectivity (IU/pg RT)',pch=21,bg='#00000033',cex=1.3,mgp=c(2.3,.8,0),las=1,xlab=''))
title(xlab='Input isolate infectivity (IU/pg RT)',mgp=c(1.9,1,0))
abline(0,1,lty=1)
dnar::withAs(combo=combo[!is.na(combo$infectivity_Old)&!is.na(combo$infectivity_Retro)&combo$infectivity_Old+combo$infectivity_Retro<Inf,],plot(combo$infectivity_Old,combo$infectivity_Retro,ylab='Retro isolate infectivity (IU/pg RT)',pch=21,bg='#00000033',cex=1.3,mgp=c(2.3,.8,0),las=1,xlab=''))
title(xlab='Input isolate infectivity (IU/pg RT)',mgp=c(1.9,1,0))
abline(0,1,lty=1)
dev.off()
newNames<-c('No.Drug.retro'='IU_Retro','No.Drug'='IU_Bead','No.Drug.in'='IU_Old','RT_Retron'='RT_Retro','RT_Beads'='RT_Bead')
for(ii in names(newNames))colnames(combo)[colnames(combo)==ii]<-newNames[ii]
write.csv(combo[!is.na(combo$RT_Old),],'out/miniInfectivity_20190820.csv',row.names=FALSE)

head(combo)

