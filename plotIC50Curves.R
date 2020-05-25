source('functions.R')

source('readMeta.R')
source('readNewData.R')


if(!exists('dilutes')){
  dilutions<-read.csv('VOA MM-cohort/dilutions.csv',stringsAsFactors=FALSE)
  dilutes<-readIfns("VOA MM-cohort/IC50 Alpha for New patients BULK isol. .xls")
  newDilutes<-readIfns('data/For Scott 03.06.2018.xlsx',c(1:2,5:20))
  newDilutes$sample<-sprintf('%s new',newDilutes$sample)
  dilutes<-rbind(dilutes,newDilutes)
  dilutes[,1:20][dilutes[,1:20]<lowerLimit]<-lowerLimit
  dilutes$dilution<-sapply(dilutes$sample,function(xx)if(xx %in% dilutions$sample[dilutions$ifn=='alpha']) dilutions[dilutions$ifn=='alpha'&dilutions$sample==xx,'dilution'] else 50)
  dilutes[,1:20]<-dilutes[,1:20]*dilutes$dilution/1000
  dilutesBeta<-readIfns("VOA MM-cohort/IC50 Beta for all BULK isolates .xlsx")
  dilutesBeta[,1:20][dilutesBeta[,1:20]<lowerLimit]<-lowerLimit
  dilutesBeta$dilution<-sapply(dilutesBeta$sample,function(xx)if(xx %in% dilutions$sample[dilutions$ifn=='beta']) dilutions[dilutions$ifn=='beta'&dilutions$sample==xx,'dilution'] else 50)
  dilutesBeta[,1:20]<-dilutesBeta[,1:20]*dilutesBeta$dilution/1000
} 


pdf('out/allCurves.pdf',width=6,height=4)
  plotIfns(dilutes,concAlpha,'IFNa2 concentration (pg/ml)')
  plotIfns(dilutes,concAlpha,'IFNa2 concentration (pg/ml)',log='x')
dev.off()
pdf('out/allCurvesBeta.pdf',width=6,height=4)
  plotIfns(dilutesBeta,concBeta,'IFNb concentration (pg/ml)')
  plotIfns(dilutesBeta,concBeta,'IFNb concentration (pg/ml)',log='x')
dev.off()
pdf('out/newCurves.pdf',width=6,height=4)
  plotIfns(newDilutes,concAlpha,'IFNa2 concentration (pg/ml)')
  plotIfns(newDilutes,concAlpha,'IFNa2 concentration (pg/ml)',log='x')
dev.off()
#system('pdftk A=out/allCurves.pdf B=out/allCurvesBeta.pdf cat A164 A180 B208 B218 output out/example_ic50.pdf')

pdf('out/allCurvesSimple.pdf',width=5,height=4)
  plotIfns(dilutesBeta,concBeta,'IFNb concentration (pg/ml)',log='x',showVres=FALSE,showMax=TRUE,showMain=FALSE,showLegend=FALSE,showPercent=TRUE)
dev.off()


pdf('out/allCurvesCombo.pdf',width=10,height=5)
  plotDualIfns(dilutes,dilutesBeta,concAlpha,concBeta)
dev.off()

ic50Fits<-calcIc50s(dilutes,concAlpha)
ic50FitsBeta<-calcIc50s(dilutesBeta,concBeta)



oldFiles<-list.files('ic50','xlsx?$',full.names=TRUE)
oldDilutes<-lapply(oldFiles,function(xx){message(xx);readIfns(xx,dilCol=1)})
names(oldDilutes)<-basename(oldFiles)
oldDilutes<-mapply(function(xx,yy){xx$sample<-sprintf('%s (%s)',xx$sample,yy);xx},oldDilutes,names(oldDilutes),SIMPLIFY=FALSE)
oldAll<-do.call(rbind,oldDilutes)
oldAll$dil<-sub('.*\\(dil = ([0-9:]+)\\)','\\1',oldAll$dil)
oldAll$file<-rep(names(oldDilutes),sapply(oldDilutes,nrow))
oldAll$id<-sub('[ _][Bb]ulk','.bulk',sub('19P4B5','19.P4B5',sub(' .*$','',sub('Bulk- ([^ ]+) ','\\1.bulk ',sub('(UK|MM) ','\\1',sub(' - ','.',oldAll$sample))))))
#oldAll$id[grepl('^[0-9]',oldAll$id)]<-sprintf('UK%s',oldAll$id[grepl('^[0-9]',oldAll$id)])
oldAll$id<-sub('[.-]P([0-9])','.\\1',convertUKMM(oldAll$id,mmLookup))
oldAll$id[grepl('^MM',oldAll$id)]<-sapply(oldAll$id[grepl('^MM',oldAll$id)],fixDecimals)
oldAll<-oldAll[nchar(oldAll$id)>2,]
rowIds<-lapply(oldAll$id,grep,dat$id)
probs<-which(sapply(rowIds,length)>1)
#assuming "bulk" was not left off
rowIds[probs]<-lapply(oldAll$id[probs],function(xx)grep(sprintf('%s$',xx),dat$id))
#if(any(probs))rowIds[probs]<-mapply(function(bulk,rows)rows[grepl('bulk',dat[rows,'id'],ignore.case=TRUE)],grepl('bulk',oldAll$sample[probs],ignore.case=TRUE),rowIds[probs],SIMPLIFY=FALSE)
if(any(sapply(rowIds,length)>1))stop('Found ambiguous ID')
oldAll$row<-NA
oldAll$row[sapply(rowIds,length)==1]<-unlist(rowIds[sapply(rowIds,length)==1])
#oldAll$correctedId<-NA
#oldAll$correctedId[grep('UK',oldAll$id)]<-sapply(strsplit(oldAll$id[grep('UK',oldAll$id)],'[-._]'),function(xx){sprintf('%s.%s',fixDecimals(sprintf('%s.%s',names(ejLookup)[ejLookup==sub('UK','EJ',xx[1])],xx[2])),sub('^P','',xx[3]))})
#rowIds<-lapply(oldAll$correctedId[grep('UK',oldAll$id)],grep,dat$id)
#oldAll[grep('UK',oldAll$id),'row'][sapply(rowIds,length)==1]<-unlist(rowIds[sapply(rowIds,length)==1])
rowIds<-lapply(gsub(' ','',sub('[. _-][Bb][Uu][Ll][Kk]','_bulk',oldAll$id[is.na(oldAll$row)])),grep,sub('[.-]P([0-9])','.\\1',convertUKMM(sub('[ ._-][Bb][Uu][Ll][Kk]','.bulk',dat$Patient.Original.ID),mmLookup)))
oldAll[is.na(oldAll$row),'row'][sapply(rowIds,length)==1]<-unlist(rowIds[sapply(rowIds,length)==1])
oldAll$isAlpha<-grepl('[Aa]lpha|IC50-a2',oldAll$sample)
oldAll$isBeta<-grepl('[Bb]eta',oldAll$sample)
probs<-oldAll$isAlpha+oldAll$isBeta!=1
if(any(probs)){
  oldAll$isAlpha[probs]<-grepl('\\(alpha2\\)',rownames(oldAll)[probs])
  oldAll$isBeta[probs]<-grepl('\\(beta\\)',rownames(oldAll)[probs])
}
if(any(oldAll$isAlpha+oldAll$isBeta!=1))stop('Unclear if sample is alpha or beta')
oldAll$marvinIc50<-ifelse(oldAll$isAlpha,dat[oldAll$row,'ic50'],dat[oldAll$row,'beta'])
oldAll$marvinVres<-ifelse(oldAll$isAlpha,dat[oldAll$row,'vres'],dat[oldAll$row,'betaVres'])
oldAll$newId<-dat[oldAll$row,'id']

write.csv(oldAll[grepl('and|amd|   |x',oldAll$dil)&!is.na(oldAll$dil),c('V1','V2','file','sample','dil')],'out/dilProblems.csv',row.names=FALSE)
fixes<-read.csv('data/dilProblems.csv',stringsAsFactors=FALSE)
fixes<-unique(fixes[,-1:-2])
rownames(fixes)<-fixes$sample
oldAll[grepl('and|amd|   |x',oldAll$dil)&!is.na(oldAll$dil),'dil']<-fixes[oldAll[grepl('and|amd|   |x',oldAll$dil)&!is.na(oldAll$dil),'sample'],'dil']
tmp<-unique(oldAll[grepl('and|amd|   |x',oldAll$dil)|is.na(oldAll$dil),c('file','sample','dil')])
tmp$dil<-''
write.csv(tmp,'out/dilProblems2.csv',row.names=FALSE)
oldAll$finalDil<-as.numeric(trimws(sub('1 to ?|1:','',oldAll$dil)))

oldAlpha<-oldAll[oldAll$isAlpha,]
oldBeta<-oldAll[oldAll$isBeta,]
pdf('out/oldCurves.pdf',width=6,height=4)
  oldAlphaFits<-plotIfns(oldAlpha,concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE)
  oldBetaFits<-plotIfns(oldBeta,concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE)
  plotIfns(oldAlpha,concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=TRUE,log='x')
  plotIfns(oldBeta,concBeta,'IFNb concentration (pg/ml)',condenseTechs=TRUE,log='x')
dev.off()
tmp<-tempfile()
pdf(tmp,width=6,height=4)
  selectIds<-c('MM14.02.2A1.bulk','MM14.10.2C1.bulk')
  #02.2A1 alpha 1/50
  #02.2A1 beta 1/50
  #10.2C1 alpha 1/50
  #10.2C1 beta 1/50
  for(ii in selectIds){
    thisDat<-oldAlpha[oldAlpha$newId==ii&!is.na(oldAlpha$newId),,drop=FALSE]
    thisDat$sample<-ii
    plotIfns(thisDat,concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=TRUE,log='x',multiple=50/1000)
    thisDat<-oldBeta[oldBeta$newId==ii&!is.na(oldBeta$newId),,drop=FALSE]
    thisDat$sample<-ii
    plotIfns(thisDat,concBeta,'IFNb concentration (pg/ml)',condenseTechs=TRUE,log='x',multiple=50/1000)
  }
dev.off()
system(sprintf('pdftk %s cat 2 4 6 8 output out/example_ic50.pdf',tmp))

#find all mm33 curves
oldMm33<-oldAll[grepl('MM33',oldAll$newId),]
oldMm33<-oldMm33[apply(is.na(oldMm33[,1:20]),1,sum)<5,]
#mislabeled
oldMm33[oldMm33$newId=='MM33.19.2D2'&round(oldMm33$V1) %in% c(21618,22813),'newId']<-'MM33.19.3D2'
mm33<-readIfns('data/Missing IC50 for curves.xlsx',exclude=0,minRows=2,dilCol=1,maxRows=80)
mm33$isBeta<-grepl('[Bb]eta',mm33$sample)
mm33$isAlpha<-grepl('[Aa]lpha',mm33$sample)
if(any(mm33$isAlpha+mm33$isBeta!=1))stop('Ambiguous alpha/beta')
mm33$newId<-sub(' .*','',mm33$sample)
mm33IdsAlpha<-dat[dat$pat=='MM33'&!is.na(dat$ic50),'id'] 
mm33IdsBeta<-dat[dat$pat=='MM33'&!is.na(dat$beta),'id'] 
#dat[dat$pat=='MM33'&!is.na(dat$beta),'id'] %in% c(oldMm33$newId[oldMm33$isBeta],mm33$newId[mm33$isBeta])

targetCols<-c(sprintf('V%d',1:20),'isBeta','newId','sample')
allMm33<-rbind(oldMm33[,targetCols],mm33[,targetCols])
allMm33$sum<-apply(round(allMm33[,1:20]),1,sum,na.rm=TRUE)
allMm33<-allMm33[!duplicated(paste(allMm33$newId,allMm33$sum)),]
allMm33<-allMm33[apply(allMm33[,1:2],1,mean,na.rm=TRUE)>200&!is.na(allMm33[,'V1']),]
dummy<-sapply(mm33IdsAlpha[!mm33IdsAlpha %in% allMm33$newId[!allMm33$isBeta]],message,' missing alpha')
dummy<-sapply(mm33IdsBeta[!mm33IdsBeta %in% allMm33$newId[allMm33$isBeta]],message,' missing beta')
if(any(duplicated(allMm33$sum)))stop('Duplicate IC50 found')
allMm33[,c('ic50','vres','max')]<-t(apply(allMm33[,1:20],1,function(xx){calculateBasicIc50(concBeta,matrix(xx,nrow=1))})[c(1,3,4),])
propMm33<-t(apply(t(do.call(rbind,lapply(1:10,function(ii)apply(allMm33[,(ii-1)*2+1:2],1,mean,na.rm=TRUE)))),1,function(xx)xx/xx[1]))
propById<-by(propMm33,list(allMm33$isBeta,allMm33$newId),apply,2,median,na.rm=TRUE)
propAlpha<-do.call(rbind,propById['FALSE',])
propBeta<-do.call(rbind,propById['TRUE',])
alphaDays<-sapply(rownames(propAlpha),function(xx)dat[dat$id==xx,'time'])
betaDays<-sapply(rownames(propBeta),function(xx)dat[dat$id==xx,'time'])
if(any(is.na(alphaDays))|any(is.na(betaDays)))stop('Missing day')
fakeDays<-0:max(dat[dat$pat=='MM33'&!dat$qvoa,'time'])
cols<-rainbow.lab(length(fakeDays),alpha=.4)
names(cols)<-fakeDays
pdf('out/mm33_curves.pdf',width=5,height=5)
  plotFunc<-function(concBeta,propBeta,xlab='IFNb concentration (pg/ml)',col='#00000033',...){
    if(length(col)==1)col<-rep(col,nrow(propBeta))
    else if(length(col)!=nrow(propBeta))stop('Color length mismatch with proportions')
    par(mar=c(3.5,3,1.5,.5))
    origConcs<-concs<-concBeta
    zeroOffset<-.01
    fakeZero<-min(concs[concs>0])*zeroOffset
    concs[concs==0]<-fakeZero
    plot(1,1,type='n',xlab=xlab,ylab='',log='x',las=1,xaxt='n',mgp=c(1.8,1,0),ylim=c(0,1),yaxt='n',xlim=range(concs),bty='l',...)
    title(ylab='Proportion of untreated p24',mgp=c(2.1,1,0))
    axis(2,pretty(c(0,1.2)),las=1,tcl=-.2,mgp=c(3,.5,0))
    abline(h=.5,lty=2)
    sapply(1:nrow(propBeta),function(xx)lines(concs,propBeta[xx,],col=col[xx]))
    dnar::logAxis(1,axisMin=min(origConcs[origConcs>0]),mgp=c(2.5,.7,0))
    axis(1,min(origConcs[origConcs>0])*zeroOffset,0,mgp=c(2.5,.7,0))
  }
  plotFunc(concAlpha,propAlpha,xlab='IFNa2 concentration (pg/ml)')
  plotFunc(concBeta,propBeta)
dev.off()
pdf('out/mm33_curves_by_time.pdf',width=10,height=5)
  par(mfrow=c(2,4))
  for(ii in unique(alphaDays)){
    theseIds<-dat[dat$time==ii&dat$pat=='MM33','id']
    plotFunc(concAlpha,propAlpha[rownames(propAlpha) %in% theseIds,],xlab='IFNa2 concentration (pg/ml)',main=sprintf('%d DFOSx',ii),col=cols[as.character(ii)])
  }
  plotFunc(concAlpha,propAlpha,xlab='IFNa2 concentration (pg/ml)',col=cols[as.character(alphaDays)],main='Combined')
  par(mfrow=c(2,4))
  for(ii in unique(betaDays)){
    theseIds<-dat[dat$time==ii&dat$pat=='MM33','id']
    plotFunc(concBeta,propBeta[rownames(propBeta) %in% theseIds,],main=sprintf('%d DFOSx',ii),col=cols[as.character(ii)])
  }
  plotFunc(concBeta,propBeta,col=cols[as.character(betaDays)],main='Combined')
dev.off()



if(FALSE){
  oldAlphaVresIc50<-as.data.frame(do.call(rbind,by(oldAlpha,oldAlpha$sample,function(xx)findVresIc50(concAlpha,xx[,1:20]))))
  #oldAlphaVresIc50$id<-oldAlpha$id
  notFoundAlpha<-oldAlphaVresIc50[is.na(oldAlphaVresIc50$row),]
  oldAlphaVresIc50<-oldAlphaVresIc50[!is.na(oldAlphaVresIc50$row),]
  oldAlphaVresIc50$newId<-dat[oldAlphaVresIc50$row,'id']
  oldAlphaVresIc50$marvinIc50<-dat[oldAlphaVresIc50$row,'ic50']
  oldAlphaVresIc50$marvinVres<-dat[oldAlphaVresIc50$row,'vres']
  oldAlphaVresIc50$ic50Diff<-log2(oldAlphaVresIc50$ic50/oldAlphaVresIc50$marvinIc50)

  oldBetaVresIc50<-as.data.frame(do.call(rbind,by(oldBeta,oldBeta$sample,function(xx)findVresIc50(concBeta,xx[,1:20]))))
  #oldBetaVresIc50$id<-oldBeta$id
  notFoundBeta<-oldBetaVresIc50[is.na(oldBetaVresIc50$row),]
  oldBetaVresIc50<-oldBetaVresIc50[!is.na(oldBetaVresIc50$row),]
  oldBetaVresIc50$newId<-dat[oldBetaVresIc50$row,'id']
  oldBetaVresIc50$marvinIc50<-dat[oldBetaVresIc50$row,'beta']
  oldBetaVresIc50$marvinVres<-dat[oldBetaVresIc50$row,'betaVres']
  oldBetaVresIc50<-as.data.frame(do.call(rbind,by(oldBeta,oldBeta$sample,function(xx)findVresIc50(concBeta,xx[,1:20]))))
  oldAlphaFits<-do.call(rbind,oldAlphaFits)
  oldBetaFits<-do.call(rbind,oldBetaFits)
}


oldBasics<-rbind(
  as.data.frame(do.call(rbind, by(oldAll[oldAll$isAlpha,],oldAll$sample[oldAll$isAlpha],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE],dil=if(any(is.na(xx$finalDil)))-1 else xx$finalDil)}))), as.data.frame(do.call(rbind, by(oldAll[oldAll$isBeta,],oldAll$sample[oldAll$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE],dil=if(any(is.na(xx$finalDil)))-1 else xx$finalDil)})))
)
oldBasics[oldBasics$repCap<0,'repCap']<-NA
oldBasics$marvinIc50<-sapply(rownames(oldBasics),function(xx)oldAll[oldAll$sample==xx,'marvinIc50'][1])
oldBasics$marvinVres<-sapply(rownames(oldBasics),function(xx)oldAll[oldAll$sample==xx,'marvinVres'][1])
oldBasics$newId<-sapply(rownames(oldBasics),function(xx)oldAll[oldAll$sample==xx,'newId'][1])
oldBasics$isAlpha<-sapply(rownames(oldBasics),function(xx)oldAll[oldAll$sample==xx,'isAlpha'][1])
oldBasics$row<-sapply(rownames(oldBasics),function(xx)oldAll[oldAll$sample==xx,'row'][1])
newIc50<-mapply(function(isAlpha,ic50,id){
  cols<-if(isAlpha)ifna2_ic50 else ifnb_ic50
  out<-unlist(dat[dat$id==id,cols])
  diffs<-abs(out-ic50)
  if(all(is.na(diffs)))return(NA)
  out[diffs==min(diffs,na.rm=TRUE)&!is.na(diffs)]
},oldBasics$isAlpha,oldBasics$ic50,oldBasics$newId)
newVres<-mapply(function(isAlpha,ic50,id){
  out<-dat[dat$id==id,if(isAlpha)ifna2_vres else ifnb_vres]
  diffs<-abs(out-ic50)
  if(all(is.na(diffs)))return(NA)
  out[diffs==min(diffs,na.rm=TRUE)&!is.na(diffs)]
},oldBasics$isAlpha,oldBasics$vres,oldBasics$newId)
oldBasics[!is.na(newIc50),'marvinIc50']<-newIc50[!is.na(newIc50)]
oldBasics[!is.na(newVres),'marvinVres']<-newVres[!is.na(newVres)]
oldBasics$ic50Disagree<-(abs(log2(oldBasics$marvinIc50/oldBasics$ic50))>1&apply(!is.na(oldBasics[,c('ic50','marvinIc50')]),1,all))|(is.na(oldBasics$ic50)!=is.na(oldBasics$marvinIc50))
oldBasics$vresDisagree<-(abs(oldBasics$marvinVres-oldBasics$vres)>1&apply(!is.na(oldBasics[,c('vres','marvinVres')]),1,all))|(is.na(oldBasics$marvinVres)!=is.na(oldBasics$vres))
oldBasics$notFound<-is.na(oldBasics$newId)
probs<-oldBasics$ic50Disagree|oldBasics$vresDisagree|oldBasics$notFound
problems<-oldBasics[probs,]
withAs(zz=problems,write.csv(zz[order(is.na(zz$newId),is.na(zz$marvinVres),is.na(zz$ic50)),],'out/oldIc50Check.csv'))
tmp<-tapply(oldBasics[!probs,'repCap'],oldBasics[!probs,'newId'],c)
tmp2<-tapply(rownames(oldBasics[!probs,]),oldBasics[!probs,'newId'],c)
out<-data.frame('id'=rep(names(tmp),sapply(tmp,'length')),'repCap'=unlist(tmp),'sample'=unlist(tmp2),stringsAsFactors=FALSE)
write.csv(out,'out/collectedRepCaps.csv',row.names=FALSE)


newP24<-readIfns('data/p24 +vertical std. 5.31.2018 controls  test (plate 1 and plate 2_ alpha2).xls',exclude=0,firstCol=4,nameCol=2,ifnCol=1)
newP24$isAlpha<-grepl('alpha',newP24$ifn,ignore.case=TRUE)
newP24$isBeta<-grepl('beta',newP24$ifn,ignore.case=TRUE)
newP24$sample<-sprintf('%s(%s)',sub('\\(.*','',newP24$sample),newP24$ifn)
out<-rbind(as.data.frame(do.call(rbind, by(newP24[newP24$isAlpha,],newP24$sample[newP24$isAlpha],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))), as.data.frame(do.call(rbind, by(newP24[newP24$isBeta,],newP24$sample[newP24$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE])}))))
write.csv(out,'out/newP24_ic50.csv')

pdf('out/newIc50.pdf',width=6,height=4)
  plotIfns(newP24[newP24$isAlpha,],concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(newP24[newP24$isBeta,],concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()



newP24Rerun<-readIfns('data/p24 +vertical std. 6.12.2018 controls  test .xls',exclude=0,firstCol=4,nameCol=2,ifnCol=1)
newP24Rerun$isAlpha<-grepl('alpha',newP24Rerun$ifn,ignore.case=TRUE)
newP24Rerun$isBeta<-grepl('beta',newP24Rerun$ifn,ignore.case=TRUE)
newP24Rerun$sample<-sprintf('%s(%s)',sub('\\(.*','',newP24Rerun$sample),newP24Rerun$ifn)
out2<-rbind(as.data.frame(do.call(rbind, by(newP24Rerun[newP24Rerun$isAlpha,],newP24Rerun$sample[newP24Rerun$isAlpha],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))), as.data.frame(do.call(rbind, by(newP24Rerun[newP24Rerun$isBeta,],newP24Rerun$sample[newP24Rerun$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE])}))))
write.csv(out2,'out/newP24Rerun_ic50.csv')

pdf('out/newIc50Redo.pdf',width=6,height=4)
  plotIfns(newP24Rerun[newP24Rerun$isAlpha,],concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(newP24Rerun[newP24Rerun$isBeta,],concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()


weauRaw<-readIfns('data/IC50 alpha and beta for WEAU and few EJs BULK isolates 3.xlsx',ifnCol=1)
weau<-do.call(rbind, by(weauRaw,paste(weauRaw$sample,weauRaw$ifn),function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))
tmp<-weauRaw
tmp$sample<-paste(weauRaw$sample,weauRaw$ifn)
pdf('out/weauCheck.pdf',width=6,height=4)
  plotIfns(tmp,concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

filterWeau<-weauRaw
filterWeau[,1:20]<- t(apply(filterWeau[,1:20],1,function(xx){
      xx<-ifelse(xx>30000,NA,xx)
      first<-max(which.max(c(xx,-Inf))-1,which(is.na(xx)))
      if(first>1)xx[1:(first)]<-NA
      xx
}))
filterWeau[,1:20]<-filterWeau[,1:20]*as.numeric(filterWeau$ifn)
#weauFilt<-do.call(rbind, by(filterWeau,paste(filterWeau$sample),function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))
dilCol<-c('20'='red','2'='blue')
pdf('out/weauFilter.pdf',width=6,height=4)
  for(ii in seq(1,8,2)){
    thisDat<-filterWeau[ii+c(0,1,8,9),]
    plot(1,1,type='n',ylim=range(filterWeau[,1:20],na.rm=TRUE),xlim=c(min(concAlpha[-1])/100,max(concAlpha)),log='xy',xlab='IFNa2',ylab='p24',xaxt='n',yaxt='n')
    logAxis(las=1)
    logAxis(1)
    tmp<-concAlpha
    concAlpha[1]<-concAlpha[2]/100
    for(ii in 1:nrow(thisDat)){
      points(tmp,thisDat[ii,seq(1,20,2)],bg=dilCol[thisDat$ifn[ii]],pch=20+ii)
      points(tmp,thisDat[ii,seq(2,20,2)],bg=dilCol[thisDat$ifn[ii]],pch=20+ii)
    }
  }
dev.off()

weauRaw2<-readIfns('data/IC50 alpha and beta for WEAU and few EJs BULK isolates 2.xlsx',ifnCol=1)
weau2<-do.call(rbind, by(weauRaw2,paste(weauRaw2$sample,weauRaw2$ifn),function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))
write.csv(weau2[unique(paste(weauRaw2$sample,weauRaw2$ifn)),],'out/weauAlphaIc50.csv')
tmp<-weauRaw2
tmp$sample<-paste(weauRaw2$sample,weauRaw2$ifn)
pdf('out/weau2Check.pdf',width=6,height=4)
  plotIfns(tmp,concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()
justWeau<-as.data.frame(weau2[grep('WEAU',rownames(weau2)),])
justWeau$visit<-as.numeric(sapply(strsplit(rownames(justWeau),'\\.'),'[[',2))
justWeau<-justWeau[grepl('50$',rownames(justWeau)),]
conts<-as.data.frame(weau2[grep('Controls.*20$',rownames(weau2)),])
#assuming all visits are in table and id starts at 1
weauMeta<-read.csv('meta/weau.csv',stringsAsFactors=FALSE)
justWeau$time<-weauMeta[justWeau$visit,'Time']
if(!exists('dat'))source('readNewData.R')
pdf('out/prelimWeau.pdf',width=11,height=6)
par(mar=c(3.7,4,1.2,5.5),mfrow=c(1,2))
ylim<-range(c(conts$ic50,justWeau$ic50,conts2$ic50),na.rm=TRUE)
plot(justWeau$time/7,justWeau$ic50,xlab='Time (weeks) ',ylab='IFNa2 IC50',las=1,log='y',yaxt='n',ylim=ylim,main="WEAU",xlim=c(0,150),cex=1.4,pch=21,bg='#00000033',mgp=c(2.5,.7,0),tcl=-.4)
#text(justWeau$time/7,justWeau$ic50,rownames(justWeau),cex=.4)
logAxis(las=1)
slantAxis(4,conts$ic50,sub(' .*','',rownames(conts)),las=1,location=.8,axisArgs=list(tcl=-.3))
withAs(just15=dat[dat$pat=='MM15',],plot(just15$time/7,just15$ic50,xlab='Time (weeks) ',ylab='IFNa2 IC50',las=1,log='y',yaxt='n',ylim=ylim,main="MM15",xlim=c(0,150),cex=1.4,pch=21,bg='#00000033',tcl=-.4,mgp=c(2.5,.7,0)))
logAxis(las=1)
withAs(conts2=dat[dat$Patient.Original.ID %in% sub(' .*','',rownames(conts)),],slantAxis(4,conts2$ic50,conts2$id,las=1,location=.8,axisArgs=list(tcl=-.3)))
dev.off()

weauRaw3<-readIfns('data/IC50 alpha and beta for WEAU and few EJs BULK isolates .xlsx',ifnCol=1)
weauRaw3$isBeta<-grepl('beta|Beta',weauRaw3$sheet)
weauRaw3$sample<-paste(weauRaw3$sample,weauRaw3$ifn)
weau3<-rbind(as.data.frame(do.call(rbind, by(weauRaw3[!weauRaw3$isBeta,],weauRaw3$sample[!weauRaw3$isBeta],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))), as.data.frame(do.call(rbind, by(weauRaw3[weauRaw3$isBeta,],weauRaw3$sample[weauRaw3$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE])}))))
weau3$replication<-weau3$max*as.numeric(sapply(rownames(weau3),function(xx)weauRaw3[weauRaw3$sample==xx,'ifn'][1]))/1000
write.csv(weau3,'out/weau3_ic50.csv')
pdf('out/weau3Check.pdf',width=6,height=4)
  plotIfns(weauRaw3[!weauRaw3$isBeta,],concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(weauRaw3[weauRaw3$isBeta,],concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()


montaner<-readIfns('data/IC50 - Rebounds from Montaner_Dilutions fixed.xls',ifnCol=1)
montaner$dilution<-montaner$ifn
montaner$isBeta<-grepl('beta|Beta',montaner$sheet)
montaner$sample<-paste(montaner$sample,montaner$ifn)
montanerIc<-rbind(as.data.frame(do.call(rbind, by(montaner[!montaner$isBeta,],montaner$sample[!montaner$isBeta],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dilution']))}))), as.data.frame(do.call(rbind, by(montaner[montaner$isBeta,],montaner$sample[montaner$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dilution']))}))))
write.csv(montanerIc,'out/montaner_ic50.csv')
pdf('out/montanerCheck.pdf',width=6,height=4)
  plotIfns(montaner[!montaner$isBeta,],concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(montaner[montaner$isBeta,],concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

subsetConc<-c(1,3,6,10)
minusCols<--(1:20)[-sapply(subsetConc,function(xx)(xx-1)*2+1:2)]
pdf('out/montanerCheckSubset.pdf',width=6,height=4)
  plotIfns(montaner[!montaner$isBeta,minusCols],concAlpha[subsetConc],'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(montaner[montaner$isBeta,minusCols],concBeta[subsetConc],'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()
system(sprintf('pdftk A=out/montanerCheck.pdf B=out/montanerCheckSubset.pdf shuffle A B output out/montanerCompare.pdf'))



weauRaw3<-readIfns('data/IC50 alpha and beta for WEAU and few EJs BULK isolates .xlsx',ifnCol=1)
weauRaw3$isBeta<-grepl('beta|Beta',weauRaw3$sheet)
weauRaw3$sample<-paste(weauRaw3$sample,weauRaw3$ifn)
weau3<-rbind(as.data.frame(do.call(rbind, by(weauRaw3[!weauRaw3$isBeta,],weauRaw3$sample[!weauRaw3$isBeta],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))), as.data.frame(do.call(rbind, by(weauRaw3[weauRaw3$isBeta,],weauRaw3$sample[weauRaw3$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE])}))))
select<-c(1,seq(2,10,2))
selectCol<-unlist(lapply(select,function(xx)(xx-1)*2+1:2))
weau3Less<-rbind(as.data.frame(do.call(rbind, by(weauRaw3[!weauRaw3$isBeta,],weauRaw3$sample[!weauRaw3$isBeta],function(xx){calcBasicIc50(concAlpha[select],xx[,selectCol,drop=FALSE])}))), as.data.frame(do.call(rbind, by(weauRaw3[weauRaw3$isBeta,],weauRaw3$sample[weauRaw3$isBeta],function(xx){calcBasicIc50(concBeta[select],xx[,selectCol,drop=FALSE])}))))
colnames(weau3Less)<-sprintf('6Conc %s',colnames(weau3Less))
round(cbind(weau3,weau3Less),3)


reboundRaw<-read6ConcIfns('data/p24 repeat IC50 rebounds 12.14.2018.xls',dilCol=1)
reboundRaw$isBeta<-TRUE
reboundRedo<-as.data.frame(do.call(rbind, by(reboundRaw[reboundRaw$isBeta,],reboundRaw$sample[reboundRaw$isBeta],function(xx){
    calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))
})))
#write.csv(reboundRedo,'out/reboundRedoIc50.csv')
newNames<-sub(' \\(.+$','',rownames(reboundRedo))
oldNames<-sub(' \\(.+$','',rownames(rebound))
for(ii in  1:nrow(reboundRedo)){
  oldRow<-which(oldNames==newNames[ii]&grepl('Beta|beta',rownames(rebound)))
  if(length(oldRow)!=1)stop('Ambiguous redo')
  rebound[oldRow,]<-reboundRedo[ii,]
}
write.csv(rebound,'out/reboundIc50.csv')

reboundRedoConvert<-convertIfn6(reboundRaw[,1:24])
reboundRedoConvert$sample<-rep(reboundRaw$sample,each=2)
pdf('out/reboundRedoIc50.pdf',width=6,height=4)
  plotIfns(reboundRedoConvert,concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

voaRaw<-readIfns('data/IC50s - new VOAs 01.10.2019.xlsx',minRows=8,ifnCol=1)
voaRaw$dilution<-voaRaw$ifn
voaRaw$isBeta<-grepl('beta|Beta',voaRaw$sheet)
voa<-rbind(as.data.frame(do.call(rbind, by(voaRaw[!voaRaw$isBeta,],voaRaw$sample[!voaRaw$isBeta],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dilution']))}))), as.data.frame(do.call(rbind, by(voaRaw[voaRaw$isBeta,],voaRaw$sample[voaRaw$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dilution']))}))))
write.csv(voa,'out/voa_2019-01-10.csv')

pdf('out/voaCheck.pdf',width=6,height=4)
  plotIfns(voaRaw[!voaRaw$isBeta,],concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(voaRaw[voaRaw$isBeta,],concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

rebound130Raw<-read6ConcIfns('data/Rebound IC50s a2-b 6 doses.xlsx',dilCol=1)
rebound130Raw$isBeta<-grepl('beta|Beta',rebound130Raw$sheet)
rebound130<-rbind(
  as.data.frame(do.call(rbind, by(rebound130Raw[!rebound130Raw$isBeta,],rebound130Raw$sample[!rebound130Raw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(rebound130Raw[rebound130Raw$isBeta,],rebound130Raw$sample[rebound130Raw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
write.csv(rebound130,'out/rebound_2019-01-30_Ic50.csv')

rebound130Convert<-convertIfn6(rebound130Raw[,1:24])
rebound130Convert$sample<-rep(rebound130Raw$sample,each=2)
rebound130Convert$isBeta<-grepl('Beta|beta',rebound130Convert$sample)
pdf('out/rebound_2019-01-30_Ic50.pdf',width=6,height=4)
  plotIfns(rebound130Convert[!rebound130Convert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(rebound130Convert[rebound130Convert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

tmp1<-readIfns('tmpMM55/IC50 Alpha for New patients BULK isol. .xls',ifnCol=1)
tmp2<-readIfns('tmpMM55/IC50 alpha and beta for BULK isolates .xlsx',ifnCol=1)
tmp2$sheet<-paste(tmp2$sheet,ifelse(grepl('[4-8]',tmp2$sheet),'beta','alpha'))
tmp3<-readIfns('tmpMM55/IC50 alpha and beta for WEAU and few EJs BULK isolates .xlsx',ifnCol=1)
tmp4<-readIfns('tmpMM55/IC50 Beta for all BULK isolates .xlsx',ifnCol=1)
tmp4$sheet<-paste(tmp4$sheet,'beta')
mm55<-rbind(tmp1,tmp2,tmp3,tmp4)
mm55<-mm55[grep('EJ101|MM55|UK101',mm55$sample),]
mm55$isBeta<-grepl('beta|Beta',mm55$sheet)
#mm55$sample<-paste(mm55$sample,mm55$ifn)
mm55<-mm55[order(sub('EJ101|MM55','',mm55$sample)),]
mm55Ic<-rbind(as.data.frame(do.call(rbind, by(mm55[!mm55$isBeta,],mm55$sample[!mm55$isBeta],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))), as.data.frame(do.call(rbind, by(mm55[mm55$isBeta,],mm55$sample[mm55$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE])}))))
pdf('out/mm55Check.pdf',width=6,height=4)
  plotIfns(mm55[!mm55$isBeta,],concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(mm55[mm55$isBeta,],concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

imcRaw<-read6ConcIfns('ice/p24_imc_ic50_20190705.xlsx',dilCol=1,exclude=1:2)
#tech rep 2,4,6,8... have 2x dilution
imcRaw[,seq(2,24,2)]<-imcRaw[,seq(2,24,2)]*2
imcRaw[2,2]<-NA
imcRaw$dil2<-imcRaw$dilution
imcRaw$isBeta<-grepl('beta|Beta',imcRaw$sheet)
imc<-rbind(
  as.data.frame(do.call(rbind, by(imcRaw[!imcRaw$isBeta,],imcRaw$sample[!imcRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(imcRaw[imcRaw$isBeta,],imcRaw$sample[imcRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
imcConvert<-convertIfn6(imcRaw[,1:24])
imcConvert$sample<-rep(sprintf('%02d %s',1:nrow(imcRaw),imcRaw$sample),each=2)
imcConvert$isBeta<-grepl('Beta|beta',imcConvert$sample)
pdf('out/imc_ic50_20190705_ic50.pdf',width=12,height=8)
  plotIfns(imcConvert[!imcConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(imcConvert[imcConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()


oldRaw<-readIfns('data/Outliers p24 raw values from Old Machine.xlsx',exclude=0,dilCol=1,ifnCol=27)
oldRaw<-oldRaw[!grepl("^ID ",oldRaw$sample),]
oldRaw<-do.call(rbind,by(oldRaw,oldRaw$sample,function(xx)xx[as.numeric(xx$dil)==min(as.numeric(xx$dil)),]))
oldRaw$isBeta<-grepl('beta|Beta',oldRaw$ifn)
oldRaw$sample<-paste(sub(' .*','',oldRaw$sample),oldRaw$ifn)
old<-rbind(as.data.frame(do.call(rbind, by(oldRaw[!oldRaw$isBeta,],oldRaw$sample[!oldRaw$isBeta],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dil']))}))), as.data.frame(do.call(rbind, by(oldRaw[oldRaw$isBeta,],oldRaw$sample[oldRaw$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dil']))}))))
newRaw<-readIfns('data/Outliers p24 raw values from New Machine.xlsx',exclude=0,dilCol=1,ifnCol=27)
newRaw<-newRaw[!grepl("^ID ",newRaw$sample),]
newRaw<-do.call(rbind,by(newRaw,newRaw$sample,function(xx)xx[as.numeric(xx$dil)==min(as.numeric(xx$dil)),]))
newRaw$isBeta<-grepl('beta|Beta',newRaw$ifn)
newRaw$sample<-paste(sub(' .*','',newRaw$sample),newRaw$ifn)
new<-rbind(as.data.frame(do.call(rbind, by(newRaw[!newRaw$isBeta,],newRaw$sample[!newRaw$isBeta],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dil']))}))), as.data.frame(do.call(rbind, by(newRaw[newRaw$isBeta,],newRaw$sample[newRaw$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dil']))}))))

pdf('out/old_201907.pdf',width=6,height=4)
  plotIfns(oldRaw[!oldRaw$isBeta,],concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(oldRaw[oldRaw$isBeta,],concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()
pdf('out/new_201907.pdf',width=6,height=4)
  plotIfns(newRaw[!newRaw$isBeta,],concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(newRaw[newRaw$isBeta,],concBeta,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

if(any(rownames(new)!=rownames(old)))stop('Misordered rownames')
tmp<-cbind(new[,1,drop=FALSE],old[,1,drop=FALSE])
colnames(tmp)<-c('new_ic50','old_ic50')
write.csv(tmp,'out/old_vs_new_ic50.csv')

colnames(new)<-sprintf('new_%s',colnames(new))
colnames(old)<-sprintf('old_%s',colnames(old))
tmp<-cbind(old,new)
tmp$foldDiff<-2^abs(log2(tmp$new_ic50/tmp$old_ic50))
write.csv(tmp,'out/20190725_outliers.csv')

newPendRaw<-read6ConcIfns('data/Pending p24 raw values from New Machine.xlsx',dilCol=1,exclude=0,ifnCol=27)
newPendRaw$isBeta<-grepl('beta|Beta',newPendRaw$ifn)
newPendRaw<-newPendRaw[!grepl('^[OP] ',newPendRaw$sample),]
newPendRaw$sample<-paste(newPendRaw$sample,newPendRaw$ifn)#sub('\\(.*','',newPendRaw$sample)
newPendRaw$sample<-sub('alpha and beta','',newPendRaw$sample)
newPend<-rbind(
  as.data.frame(do.call(rbind, by(newPendRaw[!newPendRaw$isBeta,],newPendRaw$sample[!newPendRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(newPendRaw[newPendRaw$isBeta,],newPendRaw$sample[newPendRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
newPendConvert<-convertIfn6(newPendRaw[,1:24])
newPendConvert$sample<-rep(sprintf('%02d %s',1:nrow(newPendRaw),newPendRaw$sample),each=2)
newPendConvert$isBeta<-grepl('Beta|beta',newPendConvert$sample)
pdf('out/newPending_20190725_ic50.pdf',width=12,height=8)
  plotIfns(newPendConvert[!newPendConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(newPendConvert[newPendConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

oldPendRaw<-read6ConcIfns('data/Pending p24 raw values from Old Machine.xlsx',dilCol=1,exclude=0,ifnCol=27)
oldPendRaw$sample<-sub('alpha and beta','',oldPendRaw$sample)
oldPendRaw$isBeta<-grepl('beta|Beta',oldPendRaw$ifn)
oldPendRaw<-oldPendRaw[!grepl('^[OP] ',oldPendRaw$sample),]
oldPendRaw$sample<-paste(oldPendRaw$sample,oldPendRaw$ifn)#paste(sub('\\(.*','',oldPendRaw$sample),oldPendRaw$ifn)
oldPend<-rbind(
  as.data.frame(do.call(rbind, by(oldPendRaw[!oldPendRaw$isBeta,],oldPendRaw$sample[!oldPendRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(oldPendRaw[oldPendRaw$isBeta,],oldPendRaw$sample[oldPendRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
prob<-'QVOA A08 Pre ATI - M5 (6.3.19_plate 1-  ) Alpha'
oldPend[rownames(oldPend)==prob,]<-withAs(xx=oldPendRaw[oldPendRaw$sample==prob,],calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE])[1,],dil=as.numeric(xx[,'dilution'])[1]))
oldPendConvert<-convertIfn6(oldPendRaw[,1:24])
oldPendConvert$sample<-rep(sprintf('%02d %s',1:nrow(oldPendRaw),oldPendRaw$sample),each=2)
oldPendConvert$isBeta<-grepl('Beta|beta',oldPendConvert$sample)
pdf('out/oldPending_20190725_ic50.pdf',width=12,height=8)
  plotIfns(oldPendConvert[!oldPendConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(oldPendConvert[oldPendConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()


if(any(rownames(newPend)!=rownames(oldPend)))stop('Misordered rownames')
tmpNew<-newPend
tmpOld<-oldPend
colnames(tmpNew)<-sprintf('new_%s',colnames(tmpNew))
colnames(tmpOld)<-sprintf('old_%s',colnames(tmpOld))
tmp<-cbind(tmpOld,tmpNew)
tmp$foldDiff<-2^abs(log2(tmp$new_ic50/tmp$old_ic50))
write.csv(tmp,'out/20190725_pending.csv')

oldPendRaw2<-read6ConcIfns('data/Reb and QVOA p24  re-run 07.29.2019 - Old Machine.xlsx',dilCol=1,exclude=0,ifnCol=27)
oldPendRaw2<-oldPendRaw2[!grepl('^[OP] ',oldPendRaw2$sample),]
oldPendRaw2$isBeta<-grepl('beta|Beta',oldPendRaw2$ifn)
oldPendRaw2$sample<-paste(oldPendRaw2$sample,oldPendRaw2$ifn)
oldPendRaw2$sample<-sub('\\([^(]+\\) (Alpha|Beta)','\\1',oldPendRaw2$sample)
oldPendRaw2$dilution<-as.numeric(sub('1 to ','',oldPendRaw2$dil))
oldPend2<-rbind(
  as.data.frame(do.call(rbind, by(oldPendRaw2[!oldPendRaw2$isBeta,],oldPendRaw2$sample[!oldPendRaw2$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(oldPendRaw2[oldPendRaw2$isBeta,],oldPendRaw2$sample[oldPendRaw2$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
oldPend2$ic50[is.na(oldPend2$ic50)&grep('Alpha',rownames(oldPend2))]<-max(concAlpha6)
write.csv(oldPend2,'out/oldPending2_20190801.csv')
oldPendConvert2<-convertIfn6(oldPendRaw2[,1:24])
oldPendConvert2$sample<-rep(sprintf('%02d %s',1:nrow(oldPendRaw2),oldPendRaw2$sample),each=2)
oldPendConvert2$isBeta<-grepl('Beta|beta',oldPendConvert2$sample)
pdf('out/oldPending2_20190801_ic50.pdf',width=12,height=8)
  plotIfns(oldPendConvert2[!oldPendConvert2$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(oldPendConvert2[oldPendConvert2$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()


#repCapRaw<-readIfns('data/Resending raw data for Scot (plate7 and 14).xlsx',dilCol=1,exclude=0,ifnCol=27)
#repCapRaw$isBeta<-grepl('beta|Beta',repCapRaw$ifn)
#repCapRaw$sample<-paste(sub(' .*','',repCapRaw$sample),repCapRaw$ifn)
#repCap<-rbind(as.data.frame(do.call(rbind, by(repCapRaw[!repCapRaw$isBeta,],repCapRaw$sample[!repCapRaw$isBeta],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dil']))}))), as.data.frame(do.call(rbind, by(repCapRaw[repCapRaw$isBeta,],repCapRaw$sample[repCapRaw$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE],dil=as.numeric(xx[,'dil']))}))))



tmp<-data.frame('id'=oldPend2$id,'high'=oldPend[sapply(oldPend2$id,function(xx)min(which(oldPend$id==xx))),'ic50'],'low'=oldPend2[,'ic50'])
tmp$high[is.na(tmp$high)]<-max(concAlpha)
pdf('out/high_vs_low.pdf')
  withAs(tmp=tmp[grep('Alpha',tmp$id),],plot(tmp$high,tmp$low,log='xy',xaxt='n',yaxt='n',xlab='High base measurement (first)',ylab='Low base measurement (second)',main='IFNa2',xlim=range(c(tmp$low,tmp$high)),ylim=range(c(tmp$low,tmp$high))))
  withAs(tmp=tmp[grep('Alpha',tmp$id),],text(tmp$high,tmp$low*1.05,tmp$id,cex=.2))
  logAxis(las=1)
  logAxis(1)
  abline(0,1)
  withAs(tmp=tmp[grep('Beta',tmp$id),],plot(tmp$high,tmp$low,log='xy',xaxt='n',yaxt='n',xlab='High base measurement (first)',ylab='Low base measurement (second)',main='IFNb',xlim=range(c(tmp$low,tmp$high)),ylim=range(c(tmp$low,tmp$high))))
  withAs(tmp=tmp[grep('Beta',tmp$id),],text(tmp$high,tmp$low*1.1,tmp$id,cex=.2))
  logAxis(las=1)
  logAxis(1)
  abline(0,1)
dev.off()

rebRedoRaw<-read6ConcIfns('data/Rebound Redos.xlsx',dilCol=1,exclude=0,ifnCol=27)
rebRedoRaw<-rebRedoRaw[!grepl('Standard',rebRedoRaw$sample),]
rebRedoRaw$dilution<-as.numeric(sub('1:','',rebRedoRaw$dilution))
rebRedoRaw$isBeta<-grepl('beta|Beta',rebRedoRaw$ifn)
rebRedoRaw$sample<-paste(rebRedoRaw$sample,rebRedoRaw$ifn,rebRedoRaw$dil)#paste(sub('\\(.*','',rebRedoRaw$sample),rebRedoRaw$ifn)
rebRedo<-rbind(
  as.data.frame(do.call(rbind, by(rebRedoRaw[!rebRedoRaw$isBeta,],rebRedoRaw$sample[!rebRedoRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(rebRedoRaw[rebRedoRaw$isBeta,],rebRedoRaw$sample[rebRedoRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
redoOut<-rebRedo[grepl('500$',rownames(rebRedo)),]
rownames(redoOut)<-sub('  +',' ',sub(' \\(P.*\\)','',sub('Rebound ','',rownames(redoOut))))
rownames(redoOut)<-sub(' 500$','',sub('alpha','Alpha',sub('beta','Beta',rownames(redoOut))))
write.csv(redoOut,'out/reboundRedo_20190909_ic50.csv')

rebConvert<-convertIfn6(rebRedoRaw[,1:24])
#rebConvert$sample<-rep(sprintf('%02d %s',1:nrow(rebRedoRaw),rebRedoRaw$sample),each=2)
rebConvert$sample<-rep(sprintf('%s',rebRedoRaw$sample),each=2)
rebConvert$isBeta<-grepl('Beta|beta',rebConvert$sample)
pdf('out/reboundRedo_20190909_ic50.pdf',width=12,height=12)
  par(mfrow=c(2,1))
  plotIfns(rebConvert[!rebConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(rebConvert[rebConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

oldPend$id<-sub('\\([^)]+\\) (Alpha|Beta)$','\\1',rownames(oldPend))
oldPend$source<-'pend1'
oldPend2$id<-rownames(oldPend2)
oldPend2$source<-'pend2'
redoOut$source<-'redo'
rownames(redoOut)<-sub('^A08','Rb A08',rownames(redoOut))
redoOut$id<-rownames(redoOut)
if(any(!oldPend2$id %in% oldPend$id))stop('Mismatch between pending')
out<-rbind(
  oldPend[!oldPend$id %in% oldPend2$id[!grepl('UK61|A08\\.21_P2F4',oldPend2$id)]&!grepl('UK61|A08\\.21_P2F4',oldPend$id),],
  oldPend2[!grepl('UK61|A08\\.21_P2F4',oldPend2$id),]
)
rownames(out)<-sub(' \\(6\\..*\\) ',' ',rownames(out))
if(!any(rownames(redoOut) %in% rownames(out)))stop('Redos not in out')
out<-rbind(out[!rownames(out) %in% rownames(redoOut),],redoOut)
if(nrow(oldPend[!grepl('UK61|A08\\.21_P2F4',oldPend$id),])!=nrow(out))stop('Combined pending rows does not match')
if(any(!oldPend$id[!grepl('UK61|A08\\.21_P2F4',oldPend$id)] %in% out$id))stop('Missing combined ID')
rownames(out)<-sub('A09 Rebound','Rebound A09',sub('^S-','Rebound S-',sub('^BEAT','Rebound BEAT',sub('^Rb','Rebound',rownames(out)))))
rownames(out)<-sub('^92','QVOA 92',rownames(out))
print(table(out$source))
write.csv(out,'out/pendingCombined_20190815_ic50.csv')




concBetaMac<-c(0,4400*10^(-4:0))
macRaw<-read6ConcIfns('ice/p27 ResMac beta IC50 08.16.2019.xlsx',dilCol=2,exclude=0,minRows=3,nameCol=1)
macRaw$sample<-sub(' \\(Sheet2\\)','',macRaw$sample)
#p27 protocol includes a 1:5 dilution
macRaw$dilution<-as.numeric(macRaw$dilution)*5
macRaw<-macRaw[,c(colnames(macRaw)[1:12],'sample','dilution','sheet')]
macRaw$isBeta<-TRUE
macConvert<-macRaw
macConvert[,1:12]<-macConvert[,rep(1:6,each=2)+rep(c(0,6),6)]
mac<-as.data.frame(do.call(rbind, by(macConvert,macRaw$sample,function(xx)calcBasicIc50(concBetaMac,xx[,1:12,drop=FALSE],dil=as.numeric(xx[,'dilution'])))))
pdf('out/macaqueIsolates_20190816.pdf',width=6,height=8)
  plotIfns(macConvert,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50,ylab='p27 concentration (pg/ml)')
dev.off()

tmp<-oldPendRaw2[grepl('(BEAT-004_P4E3|BEAT-004_P4D1|9207_C19|9207_A10).*Beta',oldPendRaw2$sample),-(rep((which(concBeta %in% concBetaMac)-1)*2,each=2)+1:2)]
colnames(tmp)[1:12]<-3:14
tmp$sample<-sprintf('Human %s %d',ifelse(grepl('BEAT',tmp$sample),'Rebound','Outgrowth'),ave(1:nrow(tmp),grepl('BEAT',tmp$sample),FUN=function(xx)1:length(xx)))
macConvert2<-rbind(macConvert[,1:13],tmp[,1:13])
macConvert2[macConvert2$sample=='SHIV.D.RRj17-H11 (plasma isolate)','sample']<-'SHIV.D #1 (plasma isolate)'
macConvert2[macConvert2$sample=='SHIV.D.RIa17-G4 (plasma isolate)','sample']<-'SHIV.D #2 (plasma isolate)'
macConvert2[macConvert2$sample=='SHIV.D.191859.375M.dCT (293T viral stock)','sample']<-'SHIV.D (293T stock)'
macConvert2$sample<-sub(' viral ',' ',macConvert2$sample)
nonControl<-macConvert2[!grepl('Rebound|Outgrowth',macConvert2$sample),'sample']
col<-c('#92B2D1FF','#e5c494','#66c2a5','#fc8d62','#a6d854','#D47794FF')
col<-rep(col,c(2,rep(1,nrow(macConvert2)-4),2))
names(col)<-macConvert2[order(!grepl('Rebound',macConvert2$sample),grepl('Outgrowth',macConvert2$sample),!grepl('AD8',macConvert2$sample),grepl('1',macConvert2$sample)),'sample']
pch<-rep(c(22,21,22),c(2,nrow(macConvert2)-4,2))
names(pch)<-names(col)
pdf('out/macaqueStack.pdf',height=4.5,,width=6.5)
plotStackedIfns(macConvert2,concBetaMac,'IFNb concentration (pg/ml)',col=rev(col),bty='n',pch=pch)
dev.off()


imcRaw<-read6ConcIfns('data/IC50 Beta 10.21.2019.xlsx',dilCol=1,exclude=1,ifnCol=27)
imcRaw<-imcRaw[!grepl('Standard',imcRaw$sample),]
imcRaw$dilution<-as.numeric(sub('1 to |1:','',imcRaw$dilution))
imcRaw$isBeta<-TRUE #grepl('beta|Beta',imcRaw$ifn)
imcRaw$sample<-paste(imcRaw$sample,imcRaw$dil)#paste(imcRaw$sample,imcRaw$ifn,imcRaw$dil)#paste(sub('\\(.*','',imcRaw$sample),imcRaw$ifn)
imc<-rbind(
  #as.data.frame(do.call(rbind, by(imcRaw[!imcRaw$isBeta,],imcRaw$sample[!imcRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(imcRaw[imcRaw$isBeta,],imcRaw$sample[imcRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
imcOut<-imc[,]
rownames(imcOut)<-sub('  +',' ',sub(' \\(p.*\\)','',rownames(imcOut)))
write.csv(imcOut,'out/imc_20191105_ic50.csv')

imcConvert<-convertIfn6(imcRaw[,1:24])
#imcConvert$sample<-rep(sprintf('%02d %s',1:nrow(imcRaw),imcRaw$sample),each=2)
imcConvert$sample<-rep(sprintf('%s',imcRaw$sample),each=2)
imcConvert$isBeta<-TRUE#grepl('Beta|beta',imcConvert$sample)
pdf('out/imc_20191105_ic50.pdf',width=6,height=6)
  plotIfns(imcConvert[imcConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()


imcRaw<-read6ConcIfns('data/IC50 11.14.2019.xlsx',dilCol=1,exclude=1,ifnCol=27)
imcRaw<-imcRaw[!grepl('Standard|EMPTY|FRED',imcRaw$sample),]
imcRaw$dilution<-as.numeric(sub('1 to |1:','',imcRaw$dilution))
imcRaw$isBeta<-grepl('beta|Beta',imcRaw$sheet)
imcRaw$orig<-imcRaw$sample
imcRaw$sample<-paste(imcRaw$sample,ifelse(imcRaw$isBeta,'IFNb','IFNa2'))#paste(imcRaw$sample,imcRaw$ifn,imcRaw$dil)#paste(sub('\\(.*','',imcRaw$sample),imcRaw$ifn)
imcRaw$plate<-sub('(Plate [0-9]).*','\\1',imcRaw$sheet)
imcRaw$sample<-sprintf('%s%s',imcRaw$sample,ifelse(grepl('MM33.*IMC',imcRaw$sample),sprintf(' %s',imcRaw$plate),''))
imcRaw$sample<-sub('  +',' ',sub(' \\([Pp].*\\)','',imcRaw$sample))
imcRaw$sampleIfn<-imcRaw$sample
imcRaw$sample<-paste(imcRaw$sample,imcRaw$dil)
selector<-apply(imcRaw[,c(1:2,13:14)],1,mean)>10000
if(any(imcRaw[selector,'dilution']>50))stop('High dilution sample blow out')
origImc<-imcRaw
imcRaw<-imcRaw[!selector,]
imc<-rbind(
  as.data.frame(do.call(rbind, by(imcRaw[!imcRaw$isBeta,],imcRaw$sample[!imcRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(imcRaw[imcRaw$isBeta,],imcRaw$sample[imcRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
imcOut<-imc[,]
#rownames(imcOut)<-sub('  +',' ',sub(' \\([Pp].*\\)','',rownames(imcOut)))
write.csv(imcOut,'out/imc_20191114_ic50.csv')

imcRaw$minDil<-ave(imcRaw$dilution,imcRaw$sampleIfn,FUN=min)
selectRaw<-imcRaw[imcRaw$minDil==imcRaw$dilution,]
selectRaw$sample<-selectRaw$sampleIfn
if(any(apply(selectRaw[,c(11:12,23:24)],1,mean)<30&selectRaw$dilution==105))stop('Low vres')
imcSelect<-withAs(imcRaw=selectRaw,rbind(
  as.data.frame(do.call(rbind, by(imcRaw[!imcRaw$isBeta,],imcRaw$sample[!imcRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(imcRaw[imcRaw$isBeta,],imcRaw$sample[imcRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
))
selectRaw$id<-sub('IFN[ab]2? ?','',selectRaw$sample)
out<-data.frame('IFNa2 IC50'=rep(NA,length(unique(selectRaw$id))),'IFNa2 Vres'=NA,'IFNb IC50'=NA,'IFNb Vres'=NA,'repCap'=NA,'IFNa2 repCap'=NA,'IFNb repCap'=NA,stringsAsFactors=FALSE,row.names=unique(selectRaw$id),check.names=FALSE)
for(ii in rownames(out)){
  if(any(selector<-grepl(gsub(' Plate','.*Plate',ii),rownames(imcSelect))&grepl('IFNa2',rownames(imcSelect))))out[ii,c('IFNa2 IC50','IFNa2 Vres','IFNa2 repCap')]<-imcSelect[selector,c(1:2,4)]
  if(sum(selector)>1)stop('What')
  if(any(selector<-grepl(gsub(' Plate','.*Plate',ii),rownames(imcSelect))&grepl('IFNb',rownames(imcSelect))))out[ii,c('IFNb IC50','IFNb Vres','IFNb repCap')]<-imcSelect[selector,c(1:2,4)]
  if(any(selector<-grepl(gsub(' Plate','.*Plate',ii),rownames(imcSelect))))out[ii,c('repCap')]<-mean(imcSelect[selector,'repCap'])
}

write.csv(out,'out/imc_20191115_ic50.csv')


imcConvert<-convertIfn6(imcRaw[,1:24])
#imcConvert$sample<-rep(sprintf('%02d %s',1:nrow(imcRaw),imcRaw$sample),each=2)
imcConvert$sample<-rep(sprintf('%s',imcRaw$sample),each=2)
imcConvert$isBeta<-TRUE#grepl('Beta|beta',imcConvert$sample)
pdf('out/imc_20191114_ic50.pdf',width=6,height=6)
  plotIfns(imcConvert[!imcConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(imcConvert[imcConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

macRaw<-read6ConcIfns('data/10.21.2019 - IC50 Beta Mac.Isol. - Human CD4 cells.xlsx',dilCol=1,nameCol=2,nCol=12)
macRaw$sample<-sprintf('%s (human cells)',sub(' \\(human.*','',sub(' \\(p27.*\\)','',macRaw$sample)))
macRaw2<-read6ConcIfns('data/12.24.2019 - IC50 Beta Mac.Isol. - Macaque CD4 cells.xlsx',dilCol=1,nameCol=2,nCol=12)
macRaw2$sample<-sprintf('%s (macaque cells)',sub(' \\(human.*','',sub(' \\(p27.*\\)','',macRaw2$sample)))
macRaw<-rbind(macRaw,macRaw2)
macRaw<-macRaw[!grepl('Standard',macRaw$sample),]
#p27 protocol includes a 1:5 dilution
macRaw$dilution<-as.numeric(sub('1 to ','',macRaw$dilution))*5
macRaw$isBeta<-TRUE
macConvert<-macRaw
macConvert[,1:12]<-macConvert[,rep(1:6,each=2)+rep(c(0,6),6)]
mac<-as.data.frame(do.call(rbind, by(macConvert,macRaw$sample,function(xx)calcBasicIc50(concBeta6,xx[,1:12,drop=FALSE],dil=as.numeric(xx[,'dilution'])))))
macConvert[,1:12]<-macConvert[,1:12]+100
pdf('out/macaqueIc50_20191021.pdf',width=5,height=4)
  plotIfns(macConvert,concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50,ylab='p27 concentration (pg/ml)',showFit=FALSE)
dev.off()
pdf('out/macaqueIc50_20191021_pretty.pdf',width=10,height=8)
  par(mfrow=c(2,2))
  plotIfns(macConvert[apply(macConvert[,1:2],1,mean)>600&grepl('maca',macConvert$sample),],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50,ylab='p27 concentration (pg/ml)',showFit=FALSE)
dev.off()


repP24<-readIfns('data/For Replicative  Capacity in supplemental table.xlsx',exclude=0,minRows=2,dilCol=1)
repP24$isAlpha<-TRUE
out<-as.data.frame(do.call(rbind, by(repP24,repP24$sample,function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE],dil=as.numeric(xx$dil))})))
write.csv(out,'out/repCapFind_20200213_ic50.csv')
pdf('out/repCapFind_20200213.pdf',width=6,height=4)
  plotIfns(repP24,concAlpha,'IFNa2 concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

fredRaw<-read6ConcIfns('data/IC50 12.24.2019.xlsx',dilCol=1,ifnCol=27)
fredRaw$sample<-trimws(sub('\\(.*$','',sub('\\(IMC\\).*\\(PLATE ([0-9])','IMC Plate \\1(',fredRaw$sample)))
fredRaw<-fredRaw[!grepl('Standard',fredRaw$sample),]
fredRaw$dilution<-as.numeric(sub('1 to ','',fredRaw$dilution))
fredRaw$isBeta<-fredRaw$ifn=='Beta'
fredRaw$sampleIfn<-sprintf('%s %s',fredRaw$sample,ifelse(fredRaw$isBeta,'IFNa2','IFNb'))
switch<-c('MM33 TF IMC Plate 7'='MM33.13  IMC Plate 7','MM33.13  IMC Plate 7'='MM33 TF IMC Plate 7')
fredRaw$sample[fredRaw$sample %in% names(switch)]<-switch[fredRaw$sample[fredRaw$sample %in% names(switch)]]
fred<-list(
  as.data.frame(do.call(rbind, by(fredRaw[!fredRaw$isBeta,],fredRaw$sample[!fredRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(fredRaw[fredRaw$isBeta,],fredRaw$sample[fredRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
fred<-lapply(fred,function(xx){xx$sample<-rownames(xx);xx})
fredOut<-merge(fred[[1]],fred[[2]],by.x='sample',by.y='sample',all.x=TRUE,all.y=TRUE,suffixes=c('.IFNa2','.IFNb'))
write.csv(fredOut,'out/fred_20200224_ic50.csv',row.names=FALSE)


fredConvert<-convertIfn6(fredRaw[,1:24])
fredConvert$sample<-rep(sprintf('%s',fredRaw$sample),each=2)
fredConvert$isBeta<-rep(fredRaw$isBeta,each=2)
pdf('out/fred_20200224_ic50.pdf',width=6,height=6)
  plotIfns(fredConvert[!fredConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(fredConvert[fredConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

manRaw<-read6ConcIfns('data/IC50 - A08 and Manolis - 03.23.2020.xlsx',dilCol=1)
manRaw$sample<-trimws(sub('\\(Alpha\\)|\\(Beta\\)','',manRaw$sample))
manRaw<-manRaw[!grepl('Standard',manRaw$sample),]
manRaw$dilution<-as.numeric(sub('1:','',manRaw$dilution))
manRaw$isBeta<-grepl('Beta',manRaw$sheet)
manRaw$sampleIfn<-sprintf('%s %s',manRaw$sample,ifelse(manRaw$isBeta,'IFNa2','IFNb'))
man<-list(
  as.data.frame(do.call(rbind, by(manRaw[!manRaw$isBeta,],manRaw$sample[!manRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(manRaw[manRaw$isBeta,],manRaw$sample[manRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
man<-lapply(man,function(xx){xx$sample<-rownames(xx);xx})
manOut<-merge(man[[1]],man[[2]],by.x='sample',by.y='sample',all.x=TRUE,all.y=TRUE,suffixes=c('.IFNa2','.IFNb'))
manInfo<-read.csv('out/manolis_allIds_20200330.csv')
rownames(manInfo)<-manInfo$flaskId
manOut$flaskId<-sub('.*(M[0-9][0-9]).*','\\1',manOut$sample)
manOut$flaskId[grepl('MM[0-9][0-9]',manOut$sample)]<-NA
manOut$pat<-manInfo[manOut$flaskId,'sPat']
manOut$day<-manInfo[manOut$flaskId,'time']
write.csv(manOut,'out/man_20200424_ic50.csv',row.names=FALSE)
pdf('out/manIc50.pdf')
par(mfrow=c(2,2))
for(ifn in c('ic50.IFNa2','ic50.IFNb')){
for(ii in unique(manOut$pat[!is.na(manOut$pat)])){
  withAs(xx=manOut[manOut$pat==ii&!is.na(manOut$pat),],plot(xx$day,xx[,ifn],ylim=range(manOut[,ifn]),xlim=range(manOut$day[!is.na(manOut$day)]),log='y',yaxt='n',main=ii,ylab=sub('ic50\\.','IC50 ',ifn),pch=21,bg='red',xlab='Days after detectable rebound'))
  dnar::logAxis(las=1)
}
}
dev.off()

manConvert<-convertIfn6(manRaw[,1:24])
manConvert$sample<-rep(sprintf('%s',manRaw$sample),each=2)
manConvert$isBeta<-rep(manRaw$isBeta,each=2)
pdf('out/man_20200424_ic50.pdf',width=6,height=6)
  plotIfns(manConvert[!manConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(manConvert[manConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

a08Raw<-read6ConcIfns('data/p24 04.01.2020 - A08 only.xlsx',dilCol=1,exclude=0,ifnCol=27)
a08Raw$sample<-trimws(sub('\\(Sheet1\\)','',a08Raw$sample))
a08Raw<-a08Raw[!grepl('^[OP]$',a08Raw$sample),]
a08Raw$dilution<-as.numeric(sub('1 to ','',a08Raw$dilution))
a08Raw$isBeta<-grepl('[Bb]eta',a08Raw$ifn)
a08Raw$select<-(a08Raw$dilution==20&grepl('1E2|MM33.13',a08Raw$sample))|a08Raw$dilution==100&!grepl('1E2',a08Raw$sample)
bak<-a08Raw
a08Raw<-a08Raw[a08Raw$select,]

a08<-list(
  as.data.frame(do.call(rbind, by(a08Raw[!a08Raw$isBeta,],a08Raw$sample[!a08Raw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(a08Raw[a08Raw$isBeta,],a08Raw$sample[a08Raw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
a08<-lapply(a08,function(xx){xx$sample<-rownames(xx);xx})
a08Out<-merge(a08[[1]],a08[[2]],by.x='sample',by.y='sample',all.x=TRUE,all.y=TRUE,suffixes=c('.IFNa2','.IFNb'))
write.csv(a08Out,'out/a08_20200424_ic50.csv',row.names=FALSE)

a08Convert<-convertIfn6(a08Raw[,1:24])
a08Convert$sample<-rep(sprintf('%s',a08Raw$sample),each=2)
a08Convert$isBeta<-rep(a08Raw$isBeta,each=2)
pdf('out/a08_20200424_ic50.pdf',width=6,height=6)
  plotIfns(a08Convert[!a08Convert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(a08Convert[a08Convert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

#ic50s<-lapply(c('out/imc_20191105_ic50.csv','out/imc_20191115_ic50.csv'),read.csv,stringsAsFactors=FALSE,row.names=1)
#for(ii in c('2[PFpf]4','A09r','601','9201','1A5')){
  #message('#################',ii)
  #print(lapply(ic50s,function(xx)xx[grep(ii,rownames(xx)),]))
#}
