source('functions.R')
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


finalNames<-read.csv('data/Table S5_05.22.2020.csv',stringsAsFactors=FALSE,skip=4)
finalNames<-finalNames[finalNames$ID!=''&nchar(finalNames$ID)<80,]
finalNames$isPos<-grepl('\\*\\*\\*',finalNames$ID)
finalNames$ID<-sub('\\*\\*\\*','',finalNames$ID)
macNames<-data.frame('orig'=unique(macro$virus),stringsAsFactors=FALSE)
select<-!grepl('UG24',macNames$orig)&!grepl(' Q2 ',macNames$orig)&!grepl('-BE4C8| BE7E2',macNames$orig)&!grepl('Ao2 pre',macNames$orig)&!grepl('A14|A02|A13',macNames$orig)&!grepl('P51A3 L80',macNames$orig)&!grepl(' BE ',macNames$orig)&!grepl('_BE$',macNames$orig)
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
macro<-macro[!is.na(macro$fix),]
macro$bak<-macro$virus
macro$virus<-macro$fix

auc<-function(xx,yy,low=min(xx),high=max(xx)){
  condense<-tapply(yy,xx,function(xx)exp(mean(log(xx))))
  integrate(approxfun(as.numeric(names(condense)),condense),low,high)
}
area<-do.call(rbind,lapply(structure(unique(macro$virus),.Names=unique(macro$virus)),function(vir)sapply(structure(unique(macro$donor),.Names=unique(macro$donor)),function(don){thisDat<-macro[macro$virus==vir&macro$donor==don,];if(nrow(thisDat)==0)return(NA);auc(thisDat$day,thisDat$p24,low=2,high=20)$value})))

stackArea<-data.frame('virus'=rep(rownames(area),ncol(area)),'donor'=rep(colnames(area),each=nrow(area)),'auc'=as.vector(area),row.names=NULL,stringsAsFactors=FALSE)
baseDonor<-'zb725'
stackArea$donor<- factor(stackArea$donor,levels=c(baseDonor,unique(stackArea$donor[stackArea$donor!=baseDonor])))
mod<-lm(I(log(auc))~virus+donor+0,data=stackArea[!is.na(stackArea$auc),])
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
virDf$bak<-sapply(virDf$virus,function(xx)macro[macro$virus==xx,'bak'][1])
virDf$tf<-grepl('TF',virDf$virus)
virDf$voa<-grepl('VOA|pre|post|A[0-9]+ UT',virDf$virus)&!virDf$rebound
virDf$chronic<-grepl('CH492|MM33[ .]1[37]',virDf$virus)
virDf$class<-ifelse(virDf$rebound,'Rebound',ifelse(virDf$beta,'IFNb-selected',ifelse(virDf$tf,'TF',ifelse(virDf$voa,'VOA',ifelse(virDf$chronic,'Chronic','')))))
write.csv(virDf[,c('virus','class')],'out/macroVirusDoublecheck.csv',row.names=FALSE)
classes<-read.csv('data/macroVirusDoublecheck-FBRupdate.csv',stringsAsFactors=FALSE,row.names=1)
if(any(!virDf$bak %in% rownames(classes)))stop('Missing virus class')
virDf$type2<-classes[virDf$bak,'class']
virDf$auc<-coef[virDf$virus]
virDf$tropism<-sapply(virDf$virus,function(xx)finalNames[finalNames$ID==xx,'Virus.Tropism'])
virDf$isPos<-sapply(virDf$virus,function(xx)finalNames[finalNames$ID==xx,'isPos'])
rownames(virDf)<-virDf$virus
#Reb..A08.2F5_293T.stock from email
#virDf$tropism[virDf$virus %in% c('CH470TF 293T','CH58TF 293T','Reb. 9201_293T stock','Reb. 601r1_293T stock','Reb. A09-1A2_293T stock','Reb. A08-2F5_293T stock','Reb. A08-2F5_293T stock','CH492 P51A3 293T')]<-'R5'
#virDf$tropism[virDf$virus %in% c('Reb. A08-1A5_293T stock','A08 UT P7C1','A08 BE P4E6','A08 UT P5E2')]<-'X4'
#virDf$tropism[virDf$virus %in% c('A08 UT P7F8','A08 UT P8E8')]<-'Dual'
#virDf$tropism[virDf$virus %in% s3$Isolate.ID]<-s3[virDf$virus[virDf$virus %in% s3$Isolate.ID],'tropism']
source('readReboundData.R')
rownames(combined)<-combined$virus
nameConvert<-c("9242r1.REB13"='9242.REB.r1', "A08.REB.2F4"='A08.REB.2F4', "A08.REB.1A5"='A08.REB.1A5',  "A09.REB.1A2"='A09.REB.1A2',  "601r1.REB13"='601.REB.r1')
rownames(imc)<- nameConvert[imc$Isolate.ID2]
imc$ic50_IFNa2<-imc$IFNa2.IC50..pg.ml.5
imc$ic50_IFNb<-imc$IFNb.IC50..pg.ml.5
virDf$ic50_IFNa2<-rbind(combined[,'ic50_IFNa2',drop=FALSE],imc[,'ic50_IFNa2',drop=FALSE])[virDf$virus,'ic50_IFNa2']
virDf$ic50_IFNb<-rbind(combined[,'ic50_IFNb',drop=FALSE],imc[,'ic50_IFNb',drop=FALSE])[virDf$virus,'ic50_IFNb']
pdf('out/ic50_vs_macro.pdf')
plot(exp(virDf$auc),virDf$ic50_IFNa2,log='xy',xlab='Macrophage replication (p24 AUC)',ylab='IFNa2 IC50',yaxt='n',pch=21,bg=c('VOA'='#0000FF55','Rebound'='#FF000055')[virDf$type2])
dnar::logAxis(las=1)
plot(exp(virDf$auc),virDf$ic50_IFNb,log='xy',xlab='Macrophage replication (p24 AUC)',ylab='IFNb IC50',yaxt='n',pch=21,bg=c('VOA'='#0000FF55','Rebound'='#FF000055')[virDf$type2])
dnar::logAxis(las=1)
dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],plot(exp(virDf$auc),virDf$ic50_IFNa2,log='yx',xlab='Macrophage replication (p24 AUC)',ylab='IFNa2 IC50',pch=21,bg=c('VOA'='#0000FF55','Rebound'='#FF000055')[virDf$type2],las=1))
dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],plot(exp(virDf$auc),virDf$ic50_IFNb,log='yx',xlab='Macrophage replication (p24 AUC)',ylab='IFNb IC50',pch=21,bg=c('VOA'='#0000FF55','Rebound'='#FF000055')[virDf$type2],las=1))
dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],plot(rank(virDf$auc),rank(virDf$ic50_IFNb,na.last='keep'),log='y',xlab='Macrophage replication (p24 AUC)',ylab='IFNb IC50',pch=21,bg=c('VOA'='#0000FF55','Rebound'='#FF000055')[virDf$type2],las=1))
dev.off()

#Bayesian fit
library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
mac<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nDonor;
    int<lower=0> nObs;
    real auc[nObs];
    int<lower=0,upper=nDonor> donor[nObs];
    int<lower=0,upper=nVirus> virus[nObs];
    real censorBelow;
  }
  parameters {
    vector[nVirus] virusBeta;
    vector[nDonor] donorBeta;
    real<lower=0> obsSd;
    real<lower=0> donorSd;
    //matrix[nVirus,nDonor] virusDonorBeta;
  }
  transformed parameters{
    vector[nObs] predAuc;
    for(ii in 1:nObs)predAuc[ii]=virusBeta[virus[ii]]+donorBeta[donor[ii]];//+virusDonorBeta[virus[ii],donor[ii]];
  }
  model {
    for(ii in 1:nObs){
      if(auc[ii]>censorBelow)auc[ii]~normal(predAuc[ii],obsSd);
      else target += normal_lcdf(censorBelow|predAuc[ii],obsSd);
    }
    donorBeta~normal(0,donorSd);
    virusBeta~normal(0,10);
    obsSd~gamma(1,.1);
    donorSd~gamma(1,.1);
    //for(ii in 1:nVirus)virusDonorBeta[ii,]~double_exponential(0,1);
  }\n
'
macMod <- stan_model(model_code = mac)
areaDf<-do.call(rbind,lapply(structure(unique(macro$virus),.Names=unique(macro$virus)),function(vir){
  do.call(rbind,lapply(structure(unique(macro$donor),.Names=unique(macro$donor)),function(don){
      thisDat<-macro[macro$virus==vir&macro$donor==don,]
      if(nrow(thisDat)==0)return(NULL)
      auc<-sapply(1:(nrow(thisDat)/7),function(xx)auc(thisDat$day[1:7+7*(xx-1)],thisDat$p24[1:7+7*(xx-1)],low=2,high=20)$value)
      out<-data.frame('virus'=vir,'donor'=don,'auc'=auc,stringsAsFactors=FALSE)
  }))
})) 
donorId<-structure(1:length(unique(areaDf$donor)),.Names=sort(unique(areaDf$donor)))
virusId<-structure(1:length(unique(areaDf$virus)),.Names=unique(areaDf$virus))
dat=list(
  nVirus=length(unique(areaDf$virus)),
  nDonor=length(unique(areaDf$donor)),
  nObs=nrow(areaDf),
  auc=log(areaDf$auc),
  donor=donorId[areaDf$donor],
  virus=virusId[areaDf$virus],
  censorBelow=log(9)
)
fit <- sampling(macMod, data = dat, iter=8000, chains=40)
print(fit,c('predAuc','virusDonorBeta'),include=FALSE)

mat<-as.matrix(fit)
zz<-apply(mat[,grep('virusBeta',colnames(mat))],2,mean)
names(zz)<-names(virusId)
virDf$bayes<-zz[virDf$virus]
preds<-apply(mat[,grep('predAuc',colnames(mat))],2,mean)
pdf('test.pdf');plot(virDf$auc,ifelse(virDf$bayes<2.197,2.197,virDf$bayes));abline(0,1);plot(dat$auc,preds);abline(0,1);hist(dat$auc[dat$auc>log(9)]-preds[dat$auc>log(9)],breaks=20);plot(log(stackArea$auc[!is.na(stackArea$auc)]),pred);abline(0,1);dev.off()




check<-read.csv('data/Table S5..csv',stringsAsFactors=FALSE,skip=4)
bak<-check
check$virus<-sub('QVOA','VOA',sub('[*]+$','',check$ID))
check<-check[check$Type!='',]
if(any(!check$virus %in% macro$virus))stop('Missing virus')
if(any(!check$virus %in% rownames(virDf)))stop('Missing virus')
endPoint<-dnar::withAs(macro=macro[macro$day==20,],tapply(macro$p24,list(macro$virus,toupper(macro$donor)),mean))
zb<-colnames(check)[grep('ZB',colnames(check))]
endPoint[,]<-ifelse(is.na(endPoint),'',ifelse(endPoint>.5,'yes','no'))
which(!check[,zb]==endPoint[check$virus,zb],arr.ind=TRUE)
aucCheck<-area[check$virus,tolower(zb)]
min(aucCheck[check[,zb]=='yes'],na.rm=TRUE)
aucCheck[,]<-ifelse(is.na(aucCheck),'',ifelse(aucCheck>9,'yes','no'))
which(!check[,zb]==aucCheck,arr.ind=TRUE)
tmp<-area[,tolower(zb)]
tmp<-cbind(tmp,'average'=exp(virDf[rownames(tmp),'bayes']))
fixLabs<-sub('QVOA','VOA',sub('[*]+$','',bak$ID))
miss<-c(unique(fixLabs[!fixLabs %in% rownames(tmp)]),'FILL__')
fill<-cbind(area[1:length(miss),,drop=FALSE],'average'=NA);fill[,]<-NA
rownames(fill)<-miss
tmp<-rbind(tmp,fill)
tmp<-round(tmp,2)
tmp[tmp<9&!is.na(tmp)]<-'<9'
tmp[is.na(tmp)]<-''
colnames(tmp)<-sprintf('AUC_%s',toupper(colnames(tmp)))
expandTmp<-tmp[sub('^$','FILL__',sub('QVOA','VOA',sub('[*]+$','',bak$ID))),]
rownames(expandTmp)<-NULL
bak<-cbind(bak,expandTmp)
if(any((bak[,!grepl("AVERAGE",colnames(bak))&grepl("AUC",colnames(bak))]=='')!=(bak[,zb]=='')))stop('Missing AUC')
if(any(bak$ID %in% rownames(virDf)))
write.csv(bak,'out/TableS5_update.csv',row.names=FALSE)

pdf('out/ic50_vs_macro_logLog.pdf')
#dnar::withAs(virDf=virDf[virDf$type2 %in% c('Rebound','VOA'),],plot(exp(virDf$auc),virDf$ic50_IFNb,log='xy',xlab='Macrophage replication (p24 AUC)',ylab='IFNb IC50',pch=21,bg=c('VOA'='#0000FF55','Rebound'='#FF000055')[virDf$type2],las=1))
dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],plot(exp(virDf$auc),virDf$ic50_IFNb,log='xy',xlab='Macrophage replication (p24 AUC)',ylab='IFNb IC50',pch=21,bg=c('VOA'='#0000FF55','Rebound'='#FF000055')[virDf$type2],las=1))
dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],text(exp(virDf$auc),virDf$ic50_IFNb,virDf$virus,cex=.5,col='#00000099'))
dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],plot(exp(virDf$auc),virDf$ic50_IFNa2,log='xy',xlab='Macrophage replication (p24 AUC)',ylab='IFNa2 IC50',pch=21,bg=c('VOA'='#0000FF55','Rebound'='#FF000055')[virDf$type2],las=1))
dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],text(exp(virDf$auc),virDf$ic50_IFNa2,virDf$virus,cex=.5,col='#00000099'))
dev.off()

if(FALSE){
summary(dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],lm(log2(exp(virDf$auc))~log2(virDf$ic50_IFNb))))
(dnar::withAs(virDf=virDf[,],cor.test(exp(virDf$auc),log(virDf$ic50_IFNb),method='spearman')))
(dnar::withAs(virDf=virDf[,],cor.test(exp(virDf$auc),log(virDf$ic50_IFNa2),method='spearman')))
(dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],cor.test(exp(virDf$auc),log(virDf$ic50_IFNa2),method='spearman')))
(dnar::withAs(virDf=virDf[virDf$type2=='Rebound',],cor.test(exp(virDf$auc),log(virDf$ic50_IFNb),method='spearman')))
dnar::withAs(virDf=virDf[!virDf$isPos&virDf$type2 %in% c('TF','Chronic'),],fisher.test(table(virDf$type2,exp(virDf$auc)>500)))
dnar::withAs(virDf=virDf[!virDf$isPos&virDf$type2 %in% c('TF','Chronic'),],t.test(exp(virDf$auc)~virDf$type2))
dnar::withAs(virDf=virDf[!virDf$isPos&virDf$type2 %in% c('VOA','Rebound'),],fisher.test(table(virDf$type2,exp(virDf$auc)>500)))
dnar::withAs(virDf=virDf[!virDf$isPos&virDf$type2 %in% c('VOA','Rebound'),],t.test(exp(virDf$auc)~virDf$type2))
(dnar::withAs(virDf=virDf[virDf$type2=='VOA',],cor.test(exp(virDf$auc),log(virDf$ic50_IFNb),method='spearman')))
dnar::withAs(virDf=virDf[!virDf$isPos&virDf$type2 %in% c('TF','Chronic'),],wilcox.test(exp(virDf$auc)~virDf$type2))
dnar::withAs(virDf=virDf[!virDf$isPos&virDf$type2 %in% c('Rebound','VOA'),],wilcox.test(exp(virDf$auc)~virDf$type2))
}

pdf('out/macroRep.pdf',width=20,height=3.5)
  par(mfrow=c(1,length(unique(macro$donor))),mar=c(4.4,4,9.3,4))
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
    #clusters<-cutree(hclust(dist(log(thisLast$p24))),h=.1)
    #nClust<-ave(clusters,clusters,FUN=length)
    #probs<-nClust>1
    thisLast$newPos<-thisLast$p24
    for(kk in 2:nrow(thisLast))if(thisLast$newPos[kk]/thisLast$newPos[kk-1]<1.2)thisLast$newPos[kk]<-thisLast$newPos[kk-1]*1.3
    if(max(thisLast$newPos)>10^par('usr')[4])thisLast$newPos<-thisLast$newPos
    #probAdjust<-thisLast[probs,'p24']*1.15^ave(clusters[probs],clusters[probs],FUN=function(xx)1:length(xx)-length(xx)/2-.5)
    #axis(4,thisLast$p24[!probs],thisLast$virus[!probs],las=1,cex.axis=.4,tcl=-.1,mgp=c(1,.2,0))
    #axis(4,thisLast$p24[probs],rep('',sum(probs)),las=1,cex.axis=.4,tcl=-.1,mgp=c(1,.2,0))
    axis(4,thisLast$p24,rep('',nrow(thisLast)),las=1,cex.axis=.4,tcl=-.1,mgp=c(1,.2,0))
    for(kk in 1:nrow(thisLast))axis(4,thisLast$newPos[kk],thisLast[,'virus'][kk],las=1,cex.axis=.4,tcl=-.1,mgp=c(1,.5,0),xpd=NA,tick=FALSE)
    segments(dnar::convertLineToUser(.1,4),thisLast[,'p24'],dnar::convertLineToUser(.5,4),thisLast$newPos,xpd=NA,lwd=.5)
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
  #plot(1,1,type='n',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  legend('right',inset=.07,names(donorCol),col=donorCol,lty=1,cex=.5)
dev.off()

pdf('out/macroRep2.pdf',width=7,height=40)
  nVir<-length(unique(macro$virus))
  nDonor<-length(unique(macro$donor))
  layout(rbind(rep(0,nDonor+2),cbind(rep(0,nVir),matrix(1:c(nVir*nDonor),nrow=nVir,byrow=TRUE),rep(0,nVir)),rep(0,nDonor+2)),height=c(.3,rep(1,nVir),.75),width=c(.82,rep(1,nDonor),.95))
  par(mar=c(0,0,0,0))
  ylim<-range(macro$p24)*c(.5,2)
  virOrder<-sort(unique(macro$virus))
  donorOrder<-sort(unique(macro$donor))
  virOrder<-virOrder[order(coef[virOrder],decreasing=TRUE)]
  for(ii in virOrder){
    for(jj in donorOrder){
    #title(ylab='p24 concentration (ng/ml)',mgp=c(2.75,1,0))
    #dnar::logAxis(las=1)
    #thisDat<-thisDat[order(thisDat$donor,thisDat$day),]
      plot(1,1,type='n',ylim=ylim,ylab='',yaxt='n',xlim=c(0,22),xlab='Day',log='y',mgp=c(1.75,.5,0),tcl=-.3,xaxt='n')
      if(jj==donorOrder[1])dnar::logAxis(las=1,axisMax=max(macro$p24),axisMin=min(macro$p24),axisVals=c(-2,0,2))
      if(ii==virOrder[nVir])axis(1,c(2,10,20))
      if(ii==virOrder[ceiling(nVir/2)]&jj==donorOrder[1])mtext('p24 concentration (ng/ml)',2,3)
      if(jj==donorOrder[floor(nDonor/2)]&ii==virOrder[nVir])mtext('Days after infection',1,2.6)
      if(ii==virOrder[1])mtext(toupper(jj),3,cex=.7)
      if(jj==donorOrder[nDonor])text(par('usr')[2]+1,10^mean(par('usr')[3:4]),sub('VOA','QVOA',ii),cex=.7,xpd=NA,adj=0)
      mtext(paste(sub('VOA','QVOA',ii),toupper(jj)),3,line=-1,cex=.35)
      abline(h=.500,lty=2)
      thisDat<-macro[macro$virus==ii&macro$donor==jj,]
      if(nrow(thisDat)==0)next()
      #hard coding 7 days
      #for(kk in seq(1,sum(thisDat$donor==jj),7))lines(thisDat[thisDat$donor==jj,'day'][kk+0:6],thisDat[thisDat$donor==jj,'p24'][kk+0:6]*ifelse(thisDat[thisDat$donor==jj,'p24'][kk+0:6]==.01,donorOffset[jj],1),col=ifelse(thisDat[thisDat$donor==jj,'p24'][kk+6]>.5,'red','blue'))
      means<-tapply(thisDat[thisDat$donor==jj,'p24'],thisDat[thisDat$donor==jj,'day'],function(xx)exp(mean(log(xx))))
      #ifelse(tail(means,1)>.5,'red','blue')
      lines(as.numeric(names(means)),means,col=ifelse(area[ii,jj]>9,'red','blue'))
    }
  }
dev.off()
file.copy('out/macroRep2.pdf','out/Fig._S8.pdf',overwrite=TRUE)


pdf('out/macroAucHeat.pdf',height=12,width=7)
  par(mar=c(3.8,4,.1,.1))
  plotFunc<-function(area,log=TRUE,cutoff=9){
    #if(log)breaks<-exp(seq(min(log(area),na.rm=TRUE)*.99,max(log(area),na.rm=TRUE)*1.01,length.out=101))
    if(log)breaks<-c(min(area,na.rm=TRUE)*.99,exp(seq(log(cutoff),max(log(area),na.rm=TRUE)*1.01,length.out=100)))
    else breaks<-seq(min(area,na.rm=TRUE)*.99,max(area,na.rm=TRUE)*1.01,length.out=101)
    #if(log)cols<-c('white',rev(heat.colors(109))[-1:-10])
    #else cols<-rev(heat.colors(100))
    cols<-rev(heat.colors(100))
    image(1:ncol(area),1:nrow(area),t(area),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='')
    nas<-which(is.na(area),arr.ind=TRUE)
    rect(nas[,2]-.5,nas[,1]-.5,nas[,2]+.5,nas[,1]+.5,col='#CCCCCC',border=NA)
    abline(v=2:ncol(area)-.5,h=2:nrow(area)-.5,col='#777777')
    box()
    axis(1,1:ncol(area),toupper(colnames(area)),tcl=-.25,mgp=c(1,.3,0))
    axis(2,1:nrow(area),sub('VOA','QVOA',rownames(area)),las=1,cex.axis=.5,mgp=c(1,.5,0),tcl=-.3)
    #axis(4,1:nrow(area),sprintf('%.02f',exp(coef[rownames(area)])),las=1,cex.axis=.5)
    if(log){
      #addCol<-c(cols[1],'black',cols[-1])
      #addBreak<-c(breaks[1],8.5,breaks[-1])
      #dnar::insetScale(ifelse(1:length(breaks)==1,log(cutoff*.9),log10(breaks)),cols,insetPos=c(.015,.06,.0225,.35),at=c(log10(cutoff),1:3),labels=c(sprintf('<%d',cutoff),sapply(1:3,function(xx)as.expression(bquote(10^.(xx))))),main='Macrophage replication (p24 AUC)',labYOffset=c(-.3,0,0,0))
      dnar::insetScale(ifelse(1:length(breaks)==1,log(cutoff*.9),log10(breaks)),cols,insetPos=c(.015,.06,.0225,.35),at=c(1:3),labels=c(sapply(1:3,function(xx)as.expression(bquote(10^.(xx))))),main='Macrophage replication (p24 AUC)')
    }else{
      dnar::insetScale(breaks,cols,insetPos=c(.015,.06,.0225,.35),main='Macrophage replication (p24 AUC)')
    }
  }
  #plotFunc(area)
  plotFunc(area[order(coef[rownames(area)]),sort(colnames(area))])
  plotFunc(area[order(coef[rownames(area)]),sort(colnames(area))],log=FALSE)
dev.off()
system('pdftk out/macroAucHeat.pdf cat 2 output out/macroAucHeat_linear.pdf')
system('pdftk out/macroAucHeat.pdf cat 1 output out/Fig._S8.pdf')


pdf('out/macroCompare.pdf',width=5,height=3)
  par(mar=c(2.1,4,.1,.1))
  #xPos<-structure(1:length(unique(virDf$type2)),.Names=unique(virDf$type2))
  xPos<-c('TF'=3,'Chronic'=4,'Rebound'=1,'VOA'=2)
  if(any(!virDf$type2 %in% names(xPos)))stop('Extra type')
  #posOffset<-ave(virDf$auc,virDf$type2,FUN=function(xx)seq(-.1,.1,length.out=length(xx)))
  #plot(xPos[as.character(virDf$type2)]+posOffset,exp(virDf$auc),yaxt='n',xaxt='n',log='y',ylab='Inferred p24 AUC',xlab='',pch=21,bg=ifelse(is.na(virDf$tropism),NA,ifelse(virDf$tropism=='R5','#FF000033','#0000FF33')),type='n')
  plot(xPos[as.character(virDf$type2)],exp(virDf$auc),yaxt='n',xaxt='n',log='y',ylab='Macrophage replication (p24 AUC)',xlab='',type='n',xlim=range(xPos)+c(-.5,.5))
  posOffset<-ave(exp(virDf$auc),virDf$type2,FUN=function(xx)beeswarm::swarmx(rep(0,length(xx)),xx)$x)
  #points(xPos[as.character(virDf$type2)]+posOffset,exp(virDf$auc),pch=ifelse(virDf$isPos,8,21),col=ifelse(virDf$isPos,ifelse(virDf$tropism=='R5','#FF000077','#0000FF77'),'black'),bg=ifelse(is.na(virDf$tropism),NA,ifelse(virDf$tropism=='R5','#FF000033','#0000FF33')))
  points(xPos[as.character(virDf$type2[!virDf$isPos])]+posOffset[!virDf$isPos],exp(virDf$auc[!virDf$isPos]),pch=21,bg=ifelse(is.na(virDf$tropism[!virDf$isPos]),NA,ifelse(virDf$tropism[!virDf$isPos]=='R5','#FF000033','#0000FF33')))
  points(xPos[as.character(virDf$type2[virDf$isPos])]+posOffset[virDf$isPos],exp(virDf$auc[virDf$isPos]),pch='*',col=ifelse(virDf$tropism[virDf$isPos]=='R5','#FF000077','#0000FF77'),cex=1.5)
  abline(h=.2,lty=2)
  dnar::logAxis(las=1)
  axis(1,xPos,names(xPos),mgp=c(2,.7,0))
  legend('bottomleft',inset=c(-.16,-.203),c('R5','X4/Dual'),pch=21,pt.bg=c('#FF000033','#0000FF33'),xpd=NA,bty='n',y.intersp=.8)
  #abline(h=9,lty=2)
  #
  #
  cols<-c('Rebound'=unname(classCol['rebound']),'VOA'=unname(classCol['qvoa']),'TF'='#CCCCCC','Chronic'='#CCCCCC','Pos'='black')
  plot(xPos[as.character(virDf$type2)],exp(virDf$auc),xaxt='n',ylab='Macrophage replication (p24 AUC)',xlab='',type='n',xlim=range(xPos)+c(-.5,.5),las=1,mgp=c(3,.4,0),tcl=-.2,bty='l')
  posOffset<-ave(exp(virDf$auc),virDf$type2,FUN=function(xx)beeswarm::swarmx(rep(0,length(xx)),xx,cex=1.1)$x)
  points(xPos[as.character(virDf$type2)]+posOffset,exp(virDf$auc),pch=ifelse(virDf$tropism=='R5',21,22),bg=cols[ifelse(virDf$isPos,'Pos',virDf$type2)],cex=1.2)
  #points(xPos[as.character(virDf$type2[!virDf$isPos])]+posOffset[!virDf$isPos],exp(virDf$auc[!virDf$isPos]),pch=21,bg=ifelse(is.na(virDf$tropism[!virDf$isPos]),NA,ifelse(virDf$tropism[!virDf$isPos]=='R5','#FF000033','#0000FF33')))
  #points(xPos[as.character(virDf$type2[virDf$isPos])]+posOffset[virDf$isPos],exp(virDf$auc[virDf$isPos]),pch='*',col=ifelse(virDf$tropism[virDf$isPos]=='R5','#FF000077','#0000FF77'),cex=1.5)
  axis(1,xPos,sub('VOA','QVOA',names(xPos)),mgp=c(2,.7,0))
  legend('bottomleft',inset=c(-.18,-.203),c('R5','X4/Dual'),pch=c(21,22),xpd=NA,bty='n',y.intersp=.8,pt.cex=1.2)
  #abline(h=9,lty=2)
  #tmp<-virDf[virDf$tropism=='R5'&!is.na(virDf$tropism),]
  #xPos<-structure(1:length(levels(tmp$type2)),.Names=levels(tmp$type2))
  #posOffset<-ave(tmp$auc,tmp$type2,FUN=function(xx)seq(-.05,.05,length.out=length(xx)))
  #plot(xPos[as.character(tmp$type2)]+posOffset,exp(tmp$auc),yaxt='n',xaxt='n',log='y',ylab='Base p24 AUC',xlab='',main='R5 viruses')
  #abline(h=.2,lty=2)
  #dnar::logAxis(las=1)
  #axis(1,xPos,names(xPos))
dev.off()
system('pdftk out/macroCompare.pdf cat 2 output out/macroCompare_linear.pdf')

virDf$type3<-ifelse(virDf$isPos,'Pos',virDf$type2)
virDf$censor<-ifelse(exp(virDf$bayes)<9,9,exp(virDf$bayes))
virDf$uncensor<-exp(virDf$bayes)
pdf('out/macroComparePlusExample.pdf',width=3.5,height=5)
  for(logY in c('','y')){
  layout(matrix(1:2,nrow=2),height=c(4,3))
  cols<-c('Rebound'=unname(classCol['rebound']),'VOA'=unname(classCol['qvoa']),'TF'='#CCCCCC','Chronic'='#CCCCCC','Pos'='#333333')
  #
  selectExamples<-c('YU2','UG021','TYBE','601.REB.r1','A09.REB.1A2','9244.REB.9E6','9244.VOA.K2','9244.VOA.P11','9244.VOA.12J17','A08.REB.6D6','A08.REB.7C1','A08.VOA.1B5')
  thisDat<-macro[macro$virus %in% selectExamples&macro$donor=='zb31',]
  thisAvg<-tapply(thisDat$p24,list(thisDat$virus,thisDat$day),function(xx)exp(mean(log(xx))))
  par(mar=c(3.5,3.,.1,3.4))
  plot(thisDat$day,thisDat$p24,type='n',log='y',yaxt='n',xlab='',ylab='p24 (ng/ml)',bty='l',mgp=c(2.15,.3,0),tcl=-.3,xlim=c(0,max(thisDat$day)))
  title(xlab='Days after infection',mgp=c(1.3,1,0))
  dnar::logAxis(las=1,mgp=c(1,.6,0))
  thisTypes<-sapply(structure(rownames(thisAvg),.Names=rownames(thisAvg)),function(xx)virDf[virDf$virus==xx,'type3'])
  thisTrop<-sapply(structure(rownames(thisAvg),.Names=rownames(thisAvg)),function(xx)virDf[virDf$virus==xx,'tropism'])
  for(ii in selectExamples){
    lines(as.numeric(colnames(thisAvg)),thisAvg[ii,],col=cols[thisTypes[ii]],lwd=1.3)
    points(as.numeric(colnames(thisAvg)),thisAvg[ii,],pch=ifelse(thisTrop[ii]=='R5',21,22),bg=cols[thisTypes[ii]])
  }
  text(grconvertX(0.001,from='nfc'),grconvertY(.99,from='nfc'),'A',xpd=NA,adj=c(0,1),cex=2)
  abline(h=.5,lty=2)
  thisLast<-sort(thisAvg[,ncol(thisAvg)])
  lastPos<-10^seq(.9,2.75,length.out=length(thisLast))
  lastPos[thisLast<30]<-10^seq(-1,log10(.25),length.out=sum(thisLast<30))
  for(kk in 1:length(thisLast))axis(4,lastPos[kk],names(thisLast)[kk],las=1,cex.axis=.4,mgp=c(1,.2,0),xpd=NA,tick=FALSE,col.axis=cols[thisTypes[names(thisLast)[kk]]],cex.axis=.5)
  #,col=cols[thisTypes[names(thisLast)]])
  segments(20.54,thisLast,dnar::convertLineToUser(.15,4),lastPos,xpd=NA,lwd=1,col='#999999')
  par(mar=c(1.05,3.6,.8,0.6),lheight=.8)
  ###
  xPos<-c('TF'=3.1,'Chronic'=3.9,'Rebound'=1,'VOA'=2.2)
  if(any(!virDf$type2 %in% names(xPos)))stop('Extra type')
  plot(xPos[as.character(virDf$type2)],if(logY=='y')virDf$censor else virDf$uncensor,xaxt='n',ylab='Macrophage replication\n(average p24 AUC)',xlab='',type='n',xlim=range(xPos)+c(-.4,.25),las=1,mgp=c(1.95,.4,0),tcl=-.2,bty='l',log=logY,yaxt='n')
  if(logY=='y'){
    logAxis(las=1,mgp=c(1,.8,0))
  }else{
    axis(2,las=1,mgp=c(1,.25,0),tcl=-.15)
  }
  #posOffset<-ave(virDf$censor,virDf$type2,FUN=function(xx)beeswarm::swarmx(rep(0,length(xx)),xx,cex=ifelse(logY=='y',.85,.78))$x)
  posOffset<-ave(if(logY=='y')virDf$censor else virDf$uncensor,virDf$type2,FUN=function(xx)beeswarm::swarmx(rep(0,length(xx)),xx,cex=ifelse(logY=='y',.7,.8))$x)
  tmp<-cbind(virDf,posOffset)
  dnar::withAs(xx=tmp[order(posOffset),],points(xPos[as.character(xx$type2)]+xx$posOffset,if(logY=='y')xx$censor else xx$uncensor,pch=ifelse(xx$tropism=='R5',21,22),bg=cols[xx$type3],cex=1,lwd=ifelse(xx$virus %in% selectExamples,2.5,1)))
  for(ii in 1:length(xPos))axis(1,xPos[ii],sub('VOA','QVOA',names(xPos[ii])),mgp=c(2,.18,0),tcl=-.2)
  #legend('bottomleft',inset=c(-.25,-.215),c('R5','X4/Dual'),pch=c(21,22),xpd=NA,bty='n',y.intersp=.8,pt.cex=1.2)
  legend('topright',inset=c(-0.03,-.27),c('R5','X4/Dual'),pch=c(21,22),xpd=NA,bty='',y.intersp=.8,pt.cex=1.05,cex=.9)
  abline(h=9,lty=2)
  text(grconvertX(0.001,from='nfc'),grconvertY(1.09,from='nfc'),'B',xpd=NA,adj=c(0,1),cex=2)
  }
dev.off()
system('pdftk out/macroComparePlusExample.pdf cat 1 output out/Fig._5.pdf')




pdf('out/macPredHeat.pdf')
  breaks<-seq(min(area,na.rm=TRUE)*.99,max(area,na.rm=TRUE)*1.01,length.out=101)
  cols<-rev(heat.colors(100))
  image(1:ncol(area),1:nrow(area),exp(t(predMat)),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='')
  axis(1,1:ncol(area),colnames(area),tcl=-.3,mgp=c(1,.4,0))
  axis(2,1:nrow(area),rownames(area),las=1,cex.axis=.5)
  nas<-which(is.na(area),arr.ind=TRUE)
  rect(nas[,2]-.5,nas[,1]-.5,nas[,2]+.5,nas[,1]+.5,col='#CCCCCC',border=NA)
  abline(v=2:ncol(area)-.5,h=2:nrow(area)-.5,col='#777777')
  box()
  dnar::insetScale(log10(breaks),cols,insetPos=c(.0175,.01,.025,.3),at=0:3,labels=sapply(0:3,function(xx)as.expression(bquote(10^.(xx)))),main='p24 AUC')
dev.off()

if(FALSE){
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
}



