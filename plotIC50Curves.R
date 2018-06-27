library(xlsx)
library(dnar)
source('functions.R')

#concAlpha<-c(0, 0.001, 0.003, 0.008, 0.023, 0.068, 0.204, 0.611, 1.833, 5.5)
concAlpha<-c(0,5.55*3^(-8:0))
#concBeta<-c(0, 4.4E-05, 0.00044, 0.0044, 0.044, 0.44, 4.4, 44, 440, 4400)
concBeta<-c(0,4400*10^(-8:0))
lowerLimit<-50



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

pdf('out/allCurvesCombo.pdf',width=10,height=5)
  plotDualIfns(dilutes,dilutesBeta,concAlpha,concBeta)
dev.off()

ic50Fits<-calcIc50s(dilutes,concAlpha)
ic50FitsBeta<-calcIc50s(dilutesBeta,concBeta)



oldFiles<-list.files('ic50','xlsx?$',full.names=TRUE)
oldDilutes<-lapply(oldFiles,function(xx){message(xx);readIfns(xx)})
names(oldDilutes)<-basename(oldFiles)
oldDilutes<-mapply(function(xx,yy){xx$sample<-sprintf('%s (%s)',xx$sample,yy);xx},oldDilutes,names(oldDilutes),SIMPLIFY=FALSE)
oldAll<-do.call(rbind,oldDilutes)
oldAll$id<-sub('[ _][Bb]ulk','.bulk',sub('19P4B5','19.P4B5',sub(' .*$','',sub('Bulk- ([^ ]+) ','\\1.bulk ',sub('(UK|MM) ','\\1',sub(' - ','.',oldAll$sample))))))
#oldAll$id[grepl('^[0-9]',oldAll$id)]<-sprintf('UK%s',oldAll$id[grepl('^[0-9]',oldAll$id)])
oldAll$id<-sub('[.-]P([0-9])','.\\1',convertUKMM(oldAll$id,mmLookup))
oldAll$id[grepl('^MM',oldAll$id)]<-sapply(oldAll$id[grepl('^MM',oldAll$id)],fixDecimals)
rowIds<-lapply(oldAll$id,grep,dat$id)
probs<-which(sapply(rowIds,length)>1)
if(any(probs))rowIds[probs]<-mapply(function(bulk,rows)rows[grepl('bulk',dat[rows,'id'],ignore.case=TRUE)],grepl('bulk',oldAll$sample[probs],ignore.case=TRUE),rowIds[probs],SIMPLIFY=FALSE)
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
  as.data.frame(do.call(rbind, by(oldAll[oldAll$isAlpha,],oldAll$sample[oldAll$isAlpha],function(xx){calcBasicIc50(concAlpha,xx[,1:20,drop=FALSE])}))), as.data.frame(do.call(rbind, by(oldAll[oldAll$isBeta,],oldAll$sample[oldAll$isBeta],function(xx){calcBasicIc50(concBeta,xx[,1:20,drop=FALSE])})))
)
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
problems<-oldBasics[oldBasics$ic50Disagree|oldBasics$vresDisagree|oldBasics$notFound,]
withAs(zz=problems,write.csv(zz[order(is.na(zz$newId),is.na(zz$marvinVres),is.na(zz$ic50)),],'out/oldIc50Check.csv'))


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



