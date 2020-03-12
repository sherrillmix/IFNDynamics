library(xlsx)
source('functions.R')

zapRaw<-read6ConcIfns('data/IC50 Alpha and Beta ZAP - For Scott.xlsx',dilCol=1)
zapRaw2<-read6ConcIfns('data/IC50 Alpha for ZAP..1 to 100 dil.xlsx',dilCol=1,exclude=0,minRows=1)
#50 fold dilution seems to have gone a bit crazy
zapRaw2$sample<-sub('\\(Alpha2\\)','(Alpha)',zapRaw2$sample)
zapRaw[zapRaw$sample=='pHIV-1 M subtype B RHGA (chronic virus) (Alpha)',]<-zapRaw2[zapRaw2$sample=='pHIV-1 M subtype B RHGA (chronic virus) (Alpha)',]
zapRaw3<-read6ConcIfns('data/p24 ZAP 06.21.2019.xlsx',dilCol=1,minRows=1)
zapRaw[zapRaw$sample %in% zapRaw3$sample,]<-zapRaw3
zapRaw$isBeta<-grepl('beta|Beta',zapRaw$sheet)
zapRaw4<-read6ConcIfns('data/p24 ZAP 06.20.2019.xlsx',dilCol=1,minRows=1,exclude=1:2)
#zapRaw[zapRaw$sample=='pCR-XL-TOPO_HIV-1 M subtype B CH106 (TF) (Alpha)',]<-zapRaw4[zapRaw4$sample=='pCR-XL-TOPO_HIV-1 M subtype B CH106 (TF) (Alpha)'&zapRaw4$dilution==100,]
#zapRaw[zapRaw$sample=='Rebound A08.21_P2F4 (Alpha)',]<-zapRaw4[zapRaw4$sample=='Rebound A08.21_P2F4 (Alpha)'&zapRaw4$dilution==50,]
zap<-rbind(
  as.data.frame(do.call(rbind, by(zapRaw[!zapRaw$isBeta,],zapRaw$sample[!zapRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(zapRaw[zapRaw$isBeta,],zapRaw$sample[zapRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
write.csv(zap[!grepl('NHG L mutant',rownames(zap)),],'out/zap_2019-06-21.csv')
zapConvert<-convertIfn6(zapRaw[,1:24])
zapConvert$sample<-rep(zapRaw$sample,each=2)
zapConvert$isBeta<-grepl('Beta|beta',zapConvert$sample)
zapConvert<-zapConvert[!grepl('NHG L mutant',zapConvert$sample),] 
pdf('out/zap_20190619_ic50.pdf',width=12,height=8)
  par(mfrow=c(3,2))
  plotIfns(zapConvert[!zapConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  #plotIfns(zapConvert[!zapConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50,log='x',scaleMax=TRUE,ylims=c(0,1.4))
  plotIfns(zapConvert[zapConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  #plotIfns(zapConvert[zapConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50,log='x',scaleMax=TRUE,ylims=c(0,1.1))
dev.off()
zapConvert$sample<-sub('pCR-XL-TOPO_HIV-1 M subtype B ','',sub(' *\\((Alpha|Beta)\\)$','',sub('pHIV-1 M subtype B ','',zapConvert$sample)))
cols<-structure(rainbow.lab(length(unique(zapConvert$sample))+2)[-c(2,6)],.Names=unique(zapConvert$sample))
pdf('out/stackedZapIC50.pdf')
plotStackedIfns(zapConvert[!zapConvert$isBeta,],concAlpha6,'IFNa2 concentration (pg/ml)',cols)
plotStackedIfns(zapConvert[zapConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',cols)
dev.off()


tmp<-zap[!grepl('NHG L mutant',rownames(zap)),1:2]
tmpA<-tmp[grep('Alpha',rownames(tmp)),]
tmpB<-tmp[grep('Beta',rownames(tmp)),]
rownames(tmpB)<-sub(' (Beta)','',rownames(tmpB),fixed=TRUE)
tmp<-cbind(tmpA,tmpB[sub(' (Alpha)','',rownames(tmpA),fixed=TRUE),])
colnames(tmp)<-c('IFNa2_ic50','IFNa2_vres','IFNb_ic50','IFNb_vres')
print(tmp,digits=2)

tmp<-zap[!grepl('TOPO',rownames(zap)),4,drop=FALSE]
tmpA<-tmp[grep('Alpha',rownames(tmp)),,drop=FALSE]
tmpB<-tmp[grep('Beta',rownames(tmp)),,drop=FALSE]
tmp<-cbind(tmpA,tmpB[sub(' (Alpha)','',rownames(tmpA),fixed=TRUE),])
colnames(tmp)<-c('IFNa2_repCap','IFNb_repCap')


tmp<-zap[,4,drop=FALSE]
tmpA<-tmp[grep('Alpha',rownames(tmp)),,drop=FALSE]
tmpB<-tmp[grep('Beta',rownames(tmp)),,drop=FALSE]
rownames(tmpB)<-sub(' (Beta)','',rownames(tmpB),fixed=TRUE)
tmp<-cbind(tmpA,tmpB[sub(' (Alpha)','',rownames(tmpA),fixed=TRUE),])
colnames(tmp)<-c('repCap_alpha','repCap_beta')

source('rebound/readData.R',chdir=TRUE)
zapCombo<-zap[!grepl('NHG L mutant',rownames(zap))&grepl('Alpha|alpha',rownames(zap)),]
zapCombo$virus<-sub(' in LM region| \\(chronic virus\\)|_P2F4','',sub('pCR-XL-TOPO_HIV-1 M subtype B ','',sub(' *\\((Alpha|Beta)\\)$','',sub('pHIV-1 M subtype B ','',rownames(zapCombo)))))

#pos<-c(structure(1:length(zapCombo),.Names=zapCombo$virus),pos+nrow(zapCombo)+1.5)
pdf('out/zap_alpha.pdf',width=12,height=5.5)
  layout(t(1:2),width=c(8,4))
  pars<-plotQvoa2(combo$ic50,combo$label,pos,combo$class,combo$study,combo$speed,ylab='IFNa2 IC50 (pg/ml)',mar=c(6.9,4,.1,0))
  par(mar=c(6.9,0,.1,4))
  plot(1,1,type='n',ylim=pars$ylim,log='y',yaxt='n',xaxt='n',xlim=c(1,nrow(zapCombo)),xlab='',ylab='')
  slantAxis(1,1:nrow(zapCombo),zapCombo$virus)
  points(1:nrow(zapCombo),zapCombo$ic50,pch=21,bg='grey',cex=1.5)
  #abline(h=10^(-2:0),xpd=NA)
dev.off()

source('rebound/readBeta.R',chdir=TRUE)
zapCombo<-zap[!grepl('NHG L mutant',rownames(zap))&grepl('Beta|beta',rownames(zap)),]
zapCombo$virus<-sub(' in LM region| \\(chronic virus\\)|_P2F4','',sub('pCR-XL-TOPO_HIV-1 M subtype B ','',sub(' *\\((Alpha|Beta)\\)$','',sub('pHIV-1 M subtype B ','',rownames(zapCombo)))))
pdf('out/zap_beta.pdf',width=12,height=5.5)
  layout(t(1:2),width=c(8,4))
  pars<-plotQvoa2(combo$beta,combo$label,pos,combo$class,combo$study,combo$speed,mar=c(6.9,4,.1,0),ylab='IFNb IC50 (pg/ml)')
  par(mar=c(6.9,0,.1,4))
  plot(1,1,type='n',ylim=pars$ylim,log='y',yaxt='n',xaxt='n',xlim=c(1,nrow(zapCombo)),xlab='',ylab='')
  slantAxis(1,1:nrow(zapCombo),zapCombo$virus)
  points(1:nrow(zapCombo),zapCombo$ic50,pch=21,bg='grey',cex=1.5)
  #abline(h=10^(-2:0),xpd=NA)
dev.off()


zapRedoRaw<-read6ConcIfns('data/IC50 Alpha and Beta ZAP 06.29.2019 (for Scott).xlsx',dilCol=1)
zapRedoRaw$isBeta<-grepl('beta|Beta',zapRedoRaw$sheet)
zapRedoRaw$dilution<-as.numeric(sub('1 to ','',zapRedoRaw$dilution))
zapRedo<-rbind(
  as.data.frame(do.call(rbind, by(zapRedoRaw[!zapRedoRaw$isBeta,],zapRedoRaw$sample[!zapRedoRaw$isBeta],function(xx) calcBasicIc50(concAlpha6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution']))))),
  as.data.frame(do.call(rbind, by(zapRedoRaw[zapRedoRaw$isBeta,],zapRedoRaw$sample[zapRedoRaw$isBeta],function(xx) calcBasicIc50(concBeta6,convertIfn6(xx[,1:24,drop=FALSE]),dil=as.numeric(xx[,'dilution'])))))
)
zapRedoConvert<-convertIfn6(zapRedoRaw[,1:24])
zapRedoConvert$sample<-rep(zapRedoRaw$sample,each=2)
zapRedoConvert$isBeta<-grepl('Beta|beta',zapRedoConvert$sample)
#write.csv(zap[!grepl('NHG L mutant',rownames(zap)),],'out/zap_2019-06-21.csv')
pdf('out/zap_20190712_ic50.pdf',width=6,height=4)
  #par(mfrow=c(1,2))
  plotIfns(zapRedoConvert[!zapRedoConvert$isBeta,],concAlpha6,'IFNa concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
  plotIfns(zapRedoConvert[zapRedoConvert$isBeta,],concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50)
dev.off()

tmp<-zapRedo[!grepl('NHG L mutant',rownames(zapRedo)),c(1:2,4)]
tmpA<-tmp[grep('Alpha',rownames(tmp)),]
tmpB<-tmp[grep('Beta',rownames(tmp)),]
rownames(tmpB)<-sub(' (Beta)','',rownames(tmpB),fixed=TRUE)
tmp<-cbind(tmpA,tmpB[sub(' (Alpha)','',rownames(tmpA),fixed=TRUE),])
colnames(tmp)<-c('IFNa2_ic50','IFNa2_vres','IFNa2_repCap','IFNb_ic50','IFNb_vres','IFNa2_repCap')
write.csv(tmp,'out/zapRedo_20190712.csv')
rownames(tmp)<-sub(' in LM region| \\(chronic virus\\)|_P2F4','',sub('pCR-XL-TOPO_HIV-1 M subtype B ','',sub(' *\\((Alpha|Beta)\\)$','',sub('pHIV-1 M subtype B ','',rownames(tmp)))))
print(tmp[,],digits=2)

tmp<-zapRedoConvert[!grepl('MM33|UK61|CH106 CG low',zapRedoConvert$sample),]
tmp$sample<-trimws(gsub('[^A-Za-z0-9._ ()]','',sub('\\(IMC\\)|\\(chronic virus\\)|in LM region','',sub('pCR-XL-TOPO_HIV-1 M subtype B ','',sub(' *\\((Alpha|Beta)\\)$','',sub('pHIV-1 M subtype B ','',tmp$sample))))))
uniq<-unique(tmp$sample)
base<-sub('[_. ].*','',uniq)
#cols<-structure(rainbow.lab(length(uniq)),.Names=uniq)
colBrew<-c('#8dd3c7','#eeee55','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69')
cols<-structure(structure(colBrew[1:length(uniq)],.Names=unique(base))[base],.Names=uniq)
pch<-structure(ave(1:length(base),base,FUN=function(xx)20+1:length(xx)),.Names=uniq)
pdf('out/stackedZapIC50Redo.pdf')
plotStackedIfns(tmp[!tmp$isBeta,],concAlpha6,'IFNa2 concentration (pg/ml)',cols,pch=pch)
plotStackedIfns(tmp[tmp$isBeta,],concBeta6,'IFNb concentration (pg/ml)',cols,pch=pch)
dev.off()

source('rebound/readBeta.R',chdir=TRUE)
long<-read.csv('out/allLongitudinal.csv',row.names=1,stringsAsFactors=FALSE)
zapCombo<-zapRedo[grep('MM33.*Beta',rownames(zapRedo)),]
zapCombo$virus<-sub(' \\(.*$','',rownames(zapCombo))
zapCombo<-zapCombo[order(grepl('TF',zapCombo$virus),grepl('13',zapCombo$virus),decreasing=TRUE),]
#MM33.17.2A4 MM33.13.2D6
additional<-long[grepl('MM33.01|MM33.13|MM33.17',long$id)&!is.na(long$beta),]
additional$ic50<-additional$beta
additional$pos<-ifelse(grepl('MM33.01',additional$id),1,ifelse(additional$id=='MM33.13',2,3))-.15
combo$ic50<-combo$beta
compareIC50<-function(zapCombo,combo,ylab='IFNb IC50 (pg/ml)',additional=NULL){
  layout(t(1:2),width=c(8,4))
  pars<-plotQvoa2(combo$ic50,combo$label,pos,combo$class,combo$study,combo$speed,mar=c(6.9,4,.1,0),ylab=ylab)
  par(mar=c(6.9,0,.1,4))
  plot(1,1,type='n',ylim=pars$ylim,log='y',yaxt='n',xaxt='n',xlim=c(0.5,nrow(zapCombo)+.5),xlab='',ylab='')
  slantAxis(1,1:nrow(zapCombo),zapCombo$virus)
  points(1:nrow(zapCombo),zapCombo$ic50,pch=21,bg='grey',cex=1.5)
  if(!is.null(additional)){
    points(additional$pos,additional$ic50,bg='#FF000033',col=NA,pch=21)
  }
}
pdf('out/imc_beta.pdf',width=12,height=5.5)
  compareIC50(zapCombo,combo,additional=additional)
dev.off()
source('rebound/readData.R',chdir=TRUE)
zapCombo<-zapRedo[grep('MM33.*Alpha',rownames(zapRedo)),]
zapCombo$virus<-sub(' \\(.*$','',rownames(zapCombo))
zapCombo<-zapCombo[order(grepl('TF',zapCombo$virus),grepl('13',zapCombo$virus),decreasing=TRUE),]
additional<-long[grepl('MM33.01|MM33.13.2D6|MM33.17.2A4',long$id),]
additional$pos<-ifelse(grepl('MM33.01',additional$id),1,ifelse(grepl('MM33.13',additional$id),2,3))-.15
pdf('out/imc_alpha.pdf',width=12,height=5.5)
  compareIC50(zapCombo,combo,'IFNa2 IC50 (pg/ml)',additional)
dev.off()


source('rebound/readBeta.R',chdir=TRUE)
long<-read.csv('out/allLongitudinal.csv',row.names=1,stringsAsFactors=FALSE)
newImc<-read.csv('out/imc_20191105_ic50.csv',row.names=1)
newImc<-newImc[grepl('50$',rownames(newImc))&!grepl('MM33.1.13C1',rownames(newImc)),]
rownames(newImc)<-sub('_+IMC','',sub('_?(human CD4 expansion)? 50$','',rownames(newImc)))
newImc$virus<-rownames(newImc)
newImc<-newImc[order(grepl('MM33',newImc$virus),grepl('TF',newImc$virus),grepl('13',newImc$virus),decreasing=TRUE),]
zapCombo<-zapRedo[grep('MM33.*Beta',rownames(zapRedo)),]
zapCombo$virus<-sub(' +\\(.*$','',rownames(zapCombo))
zapCombo<-zapCombo[order(grepl('TF',zapCombo$virus),grepl('13',zapCombo$virus),decreasing=TRUE),]
tzmbl<-read.csv('out/2019-10-18_tzmblIc50Calc.csv',row.names=1)
tzmbl<-tzmbl[!grepl('MM33.13.2D6|MM33.1.13C1',rownames(tzmbl)),]
tzmbl$virus<-paste('TZMBL',rownames(tzmbl))
tzmbl$repCap<-NA
#zapCombo<-rbind(zapCombo,newImc,tzmbl[,colnames(newImc)])
zapCombo<-rbind(zapCombo,newImc)



#601 is a rake
imcAdditional<-combo[grepl('601_P4(A7|A8|C1|B4)|A09-1A2|A08.21_P2F4',combo$virus),]
imcAdditional$ic50<-imcAdditional$beta
imcAdditional$pos<-ifelse(grepl('601',imcAdditional$virus),4,ifelse(grepl('A08',imcAdditional$virus),5,6))-.15
#MM33.17.2A4 MM33.13.2D6
additional<-long[grepl('MM33.01|MM33.13.2D6|MM33.17',long$id)&!is.na(long$beta),]
additional$ic50<-additional$beta
additional$pos<-ifelse(grepl('MM33.01',additional$id),1,ifelse(grepl('MM33.13',additional$id),2,3))-.15
additional<-rbind(additional[,c('pos','ic50')],imcAdditional[,c('pos','ic50')])
combo$ic50<-combo$beta
compareIC50<-function(zapCombo,combo,ylab='IFNb IC50 (pg/ml)',additional=NULL){
  layout(t(1:2),width=c(8,4))
  pars<-plotQvoa2(combo$ic50,combo$label,pos,combo$class,combo$study,combo$speed,mar=c(6.9,4,.1,0),ylab=ylab)
  par(mar=c(6.9,0,.1,4))
  pos<-structure(1:length(unique(zapCombo$virus)),.Names=unique(zapCombo$virus))
  plot(1,1,type='n',ylim=pars$ylim,log='y',yaxt='n',xaxt='n',xlim=c(0.5,max(pos)+.5),xlab='',ylab='')
  slantAxis(1,pos,names(pos))
  points(pos[zapCombo$virus],zapCombo$ic50,pch=21,bg='grey',cex=1.5)
  if(!is.null(additional)){
    points(additional$pos,additional$ic50,bg='#FF000033',col=NA,pch=21)
  }
}
pdf('out/imc_beta2.pdf',width=12,height=5.5)
  compareIC50(zapCombo,combo,additional=additional)
dev.off()
