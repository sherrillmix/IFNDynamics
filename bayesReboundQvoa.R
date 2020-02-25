
library('rstan')
library(dnar)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('functions.R')

dat<-read.csv('out/allLongitudinal.csv',stringsAsFactors=FALSE)
acuteChronicMM<-dat[!dat$qvoa&(dat$time<60|dat$time>365),]
acuteChronicMM$class<-ifelse(acuteChronicMM$time<60,'Acute','Chronic')
acuteChronicMM$study<-'MM'
acuteChronicMM$virus<-acuteChronicMM$id
#distinct "patient" for each time point (doesn't work because acute)
#acuteChronicMM$pat<-acuteChronicMM$sample

stanCode<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nState;
    real ic50[nVirus];
    int<lower=1,upper=nPatient> patients[nVirus];
    int<lower=1,upper=nState> state[nVirus];
  }
  parameters {
    matrix[nPatient,nState] baseIc50Raw;
    vector[nState] stateMeans;
    vector[nState] stateSds;
    vector[nState] stateIsoSds;
  }
  transformed parameters{
    real expectedIC50[nVirus];
    for(ii in 1:nVirus){
      expectedIC50[ii]=stateMeans[state[ii]]+baseIc50Raw[patients[ii],state[ii]]*stateSds[state[ii]];
    }
  }
  model {
    ic50~normal(expectedIC50,stateIsoSds[patients]);
    stateIsoSds~gamma(1,.1);
    stateSds~gamma(1,.1);
    for(ii in 1:nState){
      baseIc50Raw[,ii]~normal(0,1);
    }
  }
'
mod <- stan_model(model_code = stanCode)


stanCode2<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nState;
    real ic50[nVirus];
    int<lower=1,upper=nPatient> patients[nVirus];
    int<lower=1,upper=nState> states[nVirus];
  }
  parameters {
    matrix[nPatient,nState] baseIc50Raw;
    vector[nState] stateMeans;
    vector[nState] stateSds;
    vector[nState] stateIsoSds;
  }
  transformed parameters{
    real expectedIC50[nVirus];
    for(ii in 1:nVirus){
      expectedIC50[ii]=stateMeans[1]+baseIc50Raw[patients[ii],1]*stateSds[1];
      if(states[ii]>1)expectedIC50[ii]=expectedIC50[ii]+stateMeans[states[ii]]+baseIc50Raw[patients[ii],states[ii]]*stateSds[states[ii]];
    }
  }
  model {
    ic50~normal(expectedIC50,stateIsoSds[patients]);
    stateIsoSds~gamma(1,.1);
    stateSds~gamma(1,.1);
    for(ii in 1:nState){
      baseIc50Raw[,ii]~normal(0,1);
    }
  }
'
mod2 <- stan_model(model_code = stanCode2)

stanCode3<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nState;
    int<lower=0> nStudy;
    real ic50[nVirus];
    int<lower=1,upper=nPatient> patients[nVirus];
    int<lower=1,upper=nState> states[nVirus];
    int<lower=1,upper=nStudy> studies[nVirus];
  }
  parameters {
    matrix[nPatient,nState] baseIc50Raw;
    vector[nState] stateMeans;
    vector[nStudy-1] studyMeansRaw;
    vector<lower=0>[nState] stateSds;
    vector<lower=0>[nState] stateIsoSds;
  }
  transformed parameters{
    real expectedIC50[nVirus];
    vector[nStudy] studyMeans;
    studyMeans[1]=0;
    studyMeans[2:nStudy]=studyMeansRaw;
    for(ii in 1:nVirus){
      expectedIC50[ii]=studyMeans[studies[ii]]+stateMeans[1]+baseIc50Raw[patients[ii],1]*stateSds[1];
      if(states[ii]>1)expectedIC50[ii]=expectedIC50[ii]+stateMeans[states[ii]]+baseIc50Raw[patients[ii],states[ii]]*stateSds[states[ii]];
    }
  }
  model {
    ic50~normal(expectedIC50,stateIsoSds[states]);
    stateIsoSds~gamma(1,.1);
    stateSds~gamma(1,.1);
    for(ii in 1:nState){
      baseIc50Raw[,ii]~normal(0,1);
    }
  }
'
mod3 <- stan_model(model_code = stanCode3)


source('rebound/readData.R',chdir=TRUE)
comboA<-combo[combo$study!='MM'|combo$class=='QVOA',c('virus','pat','class','ic50','study')]
comboA<-rbind(comboA,acuteChronicMM[!is.na(acuteChronicMM$ic50),colnames(comboA)])
fitBayes<-function(model,patient,states,studies,ic50,baseState='Acute',chains=50,baseStudy='Other',logFunc=log,...){
  patientId<-structure(1:length(unique(patient)),.Names=sort(unique(patient)))
  stateId<-structure(1:length(unique(states)),.Names=unique(states[order(states!=baseState)]))
  studyId<-structure(1:length(unique(studies)),.Names=sort(unique(studies[order(studies!=baseStudy)])))
  dat=list(
    nVirus=length(ic50),
    nPatient=max(patientId),
    nState=max(stateId),
    nStudy=max(studyId),
    ic50=logFunc(ic50),
    patients=patientId[patient],
    states=stateId[states],
    studies=studyId[as.character(studies)]
  )
  fit <- sampling(model, data = dat, iter=6000, chains=chains,thin=2,control=list(adapt_delta=.99,max_treedepth=15),...)
  return(list('fit'=fit,pats=patientId,states=stateId,studies=studyId,dat=dat))
}
#comboA$simpleClass<-ifelse(comboA$class %in% c('Donor','6 Month','Nadir','Last'),'Chronic',ifelse(grepl('Post-',comboA$label),'Post QVOA',comboA$class))
comboA$simpleClass<-ifelse(comboA$class %in% c('Chronic','Donor','6 Month','Nadir','Last'),'Chronic',comboA$class)
#class,pat,study,ic50,beta

#fit<-fitBayes(mod3,comboA$pat,comboA$simpleClass,comboA$study=='Transmission',comboA$ic50,chains=50)
fitA<-fitBayes(mod3,comboA$pat,comboA$simpleClass,ifelse(comboA$study %in% c('Transmission','BEAT'),comboA$study,'Other'),comboA$ic50,chains=50)

comboA$simpleClass<-ifelse(comboA$class %in% c('Chronic','Donor','6 Month','Nadir','Last'),'Chronic',ifelse(grepl('Post-',comboA$label),'Post QVOA',comboA$class))
fit2<-fitBayes(mod3,comboA$pat,comboA$simpleClass,comboA$study=='Transmission',comboA$ic50,chains=50)
#intentionally poor fit
comboA$simpleClass<-ifelse(comboA$class %in% c('Donor','6 Month','Nadir','Last','Rebound'),'Chronic',comboA$class)
fit3<-fitBayes(mod3,comboA$pat,comboA$simpleClass,comboA$study=='Transmission',comboA$ic50,chains=50)
bf1<-bridgesampling::bridge_sampler(fitA$fit)
bf2<-bridgesampling::bridge_sampler(fit2$fit)
bf3<-bridgesampling::bridge_sampler(fit3$fit)

bridgesampling::bayes_factor(bf2,bf1)

source('rebound/readBeta.R',chdir=TRUE)
comboB<-combo[combo$study!='MM'|combo$class=='QVOA',c('virus','pat','class','beta','study')]
comboB<-rbind(comboB,acuteChronicMM[!is.na(acuteChronicMM$beta),colnames(comboB)])
#remove adjustment
comboB$beta[comboB$study=='Transmission']<-comboB$beta[comboB$study=='Transmission']/6386
comboB$simpleClass<-ifelse(comboB$class %in% c('Donor','6 Month','Nadir','Last'),'Chronic',comboB$class)
#fitB<-fitBayes(mod3,comboB$pat,comboB$simpleClass,comboB$study=='Transmission',comboB$beta,chains=50)

fitB<-fitBayes(mod3,comboB$pat,comboB$simpleClass,ifelse(comboB$study %in% c('Transmission','BEAT'),comboB$study,'Other'),comboB$beta,chains=50)
save(fitA,fitB,file='work/qvoaBayes.Rdat')


comparePreds<-function(fit,ic50,ids,xWidth=.3,cols='red',log='y',expFunc=ifelse(log=='y',exp,c),...){
  mat<-as.matrix(fit$fit)
  n<-length(ic50)
  expects<-mat[,sprintf('expectedIC50[%d]',1:n)]
  sds<-mat[,sprintf('stateIsoSds[%d]',fit$states)]
  makeDense<-function(xx)density(xx,from=quantile(xx,.025),to=quantile(xx,.975))
  denses<-apply(expects,2,makeDense)
  denses2<-lapply(1:ncol(expects),function(ii)makeDense(rnorm(n*10,expects[,ii],sds[,fit$dat$states[ii]])))
  plot(1,1,type='n',xlim=c(0,n+1),ylim=range(ic50),log=log,xaxs='i',yaxt='n',xaxt='n',xlab='',...)
  if(log=='y'){
    logAxis(las=1)
  } else{
    axis(2,las=1)
  }
  for(ii in 1:n){#length(denses)){
    polygon(ii+c(denses[[ii]]$y,-rev(denses[[ii]]$y))/max(denses[[ii]]$y)*xWidth,expFunc(c(denses[[ii]]$x,rev(denses[[ii]]$x))),col='#00000033')
    polygon(ii+c(denses2[[ii]]$y,-rev(denses2[[ii]]$y))/max(denses2[[ii]]$y)*xWidth,expFunc(c(denses2[[ii]]$x,rev(denses2[[ii]]$x))),col='#00000033')
  }
  axis(1,1:n,ids,las=2,cex.axis=.5)
  points(1:n,ic50,pch='-',lwd=4,col=cols)
}

pdf('out/voaRebound_bayesCheck.pdf',width=40)
  comparePreds(fitA,comboA$ic50,comboA$virus)
  comparePreds(fitB,comboB$beta,comboB$virus)
dev.off()



plotSummary<-function(fit,ylab='IFNa2 IC50 (pg/ml)',mar=c(6.9,4,.1,.1),xWidth=.4,addAcute=TRUE,xaxis=TRUE,logYAxis=TRUE,reps=2,cols=structure(rep('#00000033',nStates),.Names=stateNames),cols2=cols){
  mat<-as.matrix(fit$fit)
  states<-mat[,sprintf('stateMeans[%d]',fit$states)]
  stateSds<-mat[,sprintf('stateSds[%d]',fit$states)]
  stateIsoSds<-mat[,sprintf('stateIsoSds[%d]',fit$states)]
  nStates<-ncol(states)
  if(addAcute){
    states[,2:nStates]<-states[,2:nStates]+states[,1]
    stateSds[,2:nStates]<-sqrt(stateSds[,2:nStates]^2+stateSds[,1]^2)
    stateIsoSds[,2:nStates]<-sqrt(stateSds[,2:nStates]^2+stateSds[,1]^2)
    ylim<-range(exp(fit$dat$ic50))
    stateNames<-names(fit$states)
  }else{
    states<-states[,-1,drop=FALSE]
    stateSds<-stateSds[,-1,drop=FALSE]
    stateIsoSds<-stateIsoSds[,-1,drop=FALSE]
    nStates<-nStates-1
    stateNames<-names(fit$states[-1,drop=FALSE])
  }
  stateNames[stateNames=='QVOA']<-'Outgrowth'
  statesPat<-do.call(cbind,lapply(1:nStates,function(ii)rnorm(nrow(states)*reps,states[,ii],stateSds[,ii])))
  statesIso<-do.call(cbind,lapply(1:nStates,function(ii)rnorm(nrow(states)*reps,states[,ii],sqrt(stateSds[,ii]^2+stateIsoSds[,ii]^2))))
  if(!addAcute)ylim<-exp(range(c(-1,1,apply(statesIso,2,quantile,c(.025,.975)))))
  makeDense<-function(xx)density(xx,from=quantile(xx,.025),to=quantile(xx,.975))
  denses<-apply(states,2,makeDense)
  denses2<-apply(statesPat,2,makeDense)
  denses3<-apply(statesIso,2,makeDense)
  #denses2<-lapply(1:ncol(expects),function(ii)makeDense(rnorm(n*10,expects[,ii],sds[,fit$dat$states[ii]])))
  par(mar=mar)
  plot(1,1,log='y',yaxt='n',ylab=ylab,xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim,xlim=c(.5,nStates+.5),las=1,mgp=c(2.45,1,0))
  if(logYAxis)logAxis(las=1,mgp=c(3,.7,0))
  else axis(2,c(.2,.5,1,2,5),c('0.2','0.5','1','2','5'),las=1)
  for(ii in 1:nStates){
    polygon(ii+c(denses[[ii]]$y,-rev(denses[[ii]]$y))/max(denses[[ii]]$y)*xWidth,exp(c(denses[[ii]]$x,rev(denses[[ii]]$x))),col=cols[stateNames[ii]])
    #polygon(ii+c(denses2[[ii]]$y,-rev(denses2[[ii]]$y))/max(denses2[[ii]]$y)*xWidth,exp(c(denses2[[ii]]$x,rev(denses2[[ii]]$x))),col=cols[stateNames[ii]])
    polygon(ii+c(denses3[[ii]]$y,-rev(denses3[[ii]]$y))/max(denses3[[ii]]$y)*xWidth,exp(c(denses3[[ii]]$x,rev(denses3[[ii]]$x))),col=cols2[stateNames[ii]])
  }
  if(xaxis)for(ii in 1:nStates)axis(1,ii,stateNames[ii])
  if(!addAcute)abline(h=1,lty=2)
  return(stateNames)
}
pdf('out/voaRebound_bayesSummary.pdf',height=3.5,width=8)
  layout(matrix(1:2,ncol=2),width=c(8,2.1))
  classCols<-structure(c("#9EC0E1E6", "#77BCA9B3", "#84C47DB3", "#9FB755B3", "#B99A4BB3", "#C77C62B3", "#E581A0E6"), .Names = c("Outgrowth", "Acute", "6 Month", "Donor", "Nadir", "Last", "Rebound")) 
  stCols<-sprintf('%s99',substring(classCols,1,7))
  stCols2<-sprintf('%s33',substring(classCols,1,7))
  names(stCols)<-names(stCols2)<-names(classCols)
  names(stCols2)[names(stCols2)=='Donor']<-names(stCols)[names(stCols)=='Donor']<-'Chronic'
  source('rebound/readData.R',chdir=TRUE)
  plotQvoa2(combo$ic50,combo$label,pos,combo$class,combo$study,combo$speed,ylab='IFNa2 IC50 (pg/ml)',mar=c(5.8,3.5,.1,.1),cex.axis=1.1,startDown=TRUE,pats=ifelse(combo$study %in% c('Transmission','MM'),NA,combo$pat))
  text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.0025,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'A',xpd=NA,adj=c(0,1),cex=2)
  #plotSummary(fitA)
  states<-plotSummary(fitA,addAcute=FALSE,ylab='Fold change from acute',mar=c(4.5,4,.1,1.6),cols=stCols,cols2=stCols2,xaxis=FALSE)
  slantAxis(1,1:length(states),states,textOffsets=c(-.6,-.4,-.2),location=.7,axisArgs=list(tcl=-.4))
  #Add B
  text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.03,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'B',xpd=NA,adj=c(0,1),cex=2)
  source('rebound/readBeta.R',chdir=TRUE)
  #remove adjustment based on single IC50 and use bayesian estimated
  combo$beta[combo$study=='Transmission']<-combo$beta[combo$study=='Transmission']/6386*2230
  #plotSummary(fitB,ylab='IFNb IC50 (pg/ml)')
  plotQvoa2(combo$beta,combo$label,pos,combo$class,combo$study,combo$speed,ylab='IFNb IC50 (pg/ml)',mar=c(5.8,3.5,.1,.1),cex.axis=1.1,startDown=TRUE,pats=ifelse(combo$study %in% c('Transmission','MM'),NA,combo$pat))
  text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.0025,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'C',xpd=NA,adj=c(0,1),cex=2)
  states<-plotSummary(fitB,addAcute=FALSE,ylab='Fold change from acute',mar=c(4.5,4,.1,1.6),cols=stCols,cols2=stCols2,xaxis=FALSE)
  slantAxis(1,1:length(states),states,textOffsets=c(-.6,-.4,-.2),location=.7,axisArgs=list(tcl=-.4))
  text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.03,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'D',xpd=NA,adj=c(0,1),cex=2)
dev.off()


#note using log10 instead of log to make plotting easier
logit<-function(p)log10(p)-log10(1-p)
invLogit<-function(x)10^(x)/(10^(x)+1)
mini<-read.csv('out/mini_repInfTrop.csv',stringsAsFactors=FALSE)
master<-read.csv('data/Rebounds and QVOAs2.csv',stringsAsFactors=FALSE)
rownames(master)<-master$ID
mini$repCap<-mini$ReCap..ng.of.p24.per.ml.
mini$p24Release<-as.numeric(sub('%$','',mini$p24.release....))/100
mini$pat<-master[mini$ID,'Patient']
mini$study<-master[mini$ID,'Study']
mini$class<-mini$X
mini$class[mini$class=="VOA"]<-"QVOA"
minorFits<-lapply(structure(c('repCap','p24Release','infectivity'),.Names=c('repCap','p24Release','infectivity')),function(ii){
  thisDat<-mini[!is.na(mini[,ii])&mini$class %in% c('Rebound','QVOA'),]
  fit<-fitBayes(mod3,thisDat$pat,thisDat$class,rep('A',nrow(thisDat)),thisDat[,ii],chains=50,baseState='QVOA',baseStudy='A',logFunc=if(ii %in% c('infectivity','repCap'))log else logit)
})

pdf('out/voaRebound_minor_bayesCheck.pdf',width=20)
  for(ii in 1:length(minorFits)){
    thisDat<-mini[!is.na(mini[,names(minorFits)[ii]])&mini$class %in% c('Rebound','QVOA'),]
    comparePreds(minorFits[[ii]],thisDat[,names(minorFits)[ii]],thisDat$ID,cols=ifelse(thisDat$class=='Rebound','red','blue'),log=ifelse(names(minorFits)[ii]=='p24Release','','y'),expFunc=ifelse(names(minorFits)[ii]=='p24Release',invLogit,exp),main=names(minorFits)[ii],ylab=names(minorFits)[ii])
  }
dev.off()

pdf('out/voaRebound_minor_bayesSummary.pdf',height=3.5)
  par(lheight=.73)
  for(ii in 1:length(minorFits)){
    par(mgp=c(2,1,0))
    plotSummary(minorFits[[ii]],addAcute=FALSE,ylab=ifelse(names(minorFits)[ii]=='p24Release','Fold change from VOA in odds ratio\nof release',sprintf('Fold change of %s\nfrom VOA',names(minorFits)[ii])),mar=c(2,4,.1,.1))
  }
dev.off()



ylabs<-c('infectivity'='Infectivity (IU/pg RT)','repCap'='Replicative capacity (ng p24/ml)','p24Release'='Proportion of p24 released')
plotMini<-function(ic50,label,pat,pos,class,study,ylab='IFNa2 IC50 (pg/ml)',mar=c(3,4,.1,.1),log='y',ylim=range(ic50)){
  spread<-offsetX(log10(ic50),label,width=.25)
  expFunc<-ifelse(log=='y',function(xx)10^xx,c)
  marSpace<-0
  par(mar=mar,lheight=.8)
  plot(pos[label]+spread,ic50,log=log,yaxt='n',ylab='',xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim)
  studyPos<-tapply(pos[label],study,mean)
  mtext(ylab,2,2.8,adj=.5,at=expFunc(par('usr')[3]+.4*diff(par('usr')[3:4])))
  #dnar::slantAxis(1,studyPos,names(studyPos))
  for(counter in 1:length(studyPos))axis(1,studyPos[counter],names(studyPos)[counter],padj=1,mgp=c(3,0+1.2*(1+counter)%%2,0),tcl=-.5+-1.2*(1+counter)%%2,cex.axis=1.1)
  if(log=='y'){
    if(diff(par('usr')[3:4])>2.5)logAxis(las=1)
    else axis(2,las=1,mgp=c(3,.8,0))
  }else{
    axis(2,las=1,mgp=c(3,.8,0))
  }
  ranges<-tapply(pos[label],pat,range)
  lapply(ranges[sapply(ranges,diff)>0],function(xx)rect(min(xx)-.4,expFunc(par('usr')[3]),max(xx)+.4,expFunc(par('usr')[4]),col='#00000033',border=NA))
  abline(v=pos,col='#00000055',lty=3)
  points(pos[label]+spread,ic50,pch=21,bg=classCols[class],col='#000000DD',lwd=1.5,cex=1.5)
  studyPos<-tapply(pos[label],study,mean)
  data.frame(patMin=tapply(pos[label],pat,min),patMax=tapply(pos[label],pat,max) )
  abline(v=tapply(pos[label],study,max)+.5+studySpace/2,lty=1)
}
mini<-mini[order(mini$study,mini$pat,mini$class),]
mini$label<-paste(mini$class,mini$pat)
mini$study<-sub('INTERRUPT','Interrupt',sub('OUTGROWTH','Longitudinal\noutgrowth',sub(' / ','/',mini$study)))
mini<-mini[mini$study!='Longitudinal\noutgrowth',] #check this isolate
pdf('out/repInfRel_qvoa_compare2.pdf',width=9,height=3.5)
  for(ii in 1:length(minorFits)){
    layout(matrix(1:2,ncol=2),width=c(8,2.2))
    thisDat<-mini[!is.na(mini[,names(minorFits)[ii]])&mini$class %in% c('Rebound','QVOA'),]
    thisDat<-thisDat[order(thisDat$study=='Longitudinal\noutgrowth'),]
    pos<-structure(1:length(unique(thisDat$label)),.Names=unique(thisDat$label))
    posStudy<-sapply(names(pos),function(xx)thisDat[thisDat$label==xx,'study'][1])
    studySpace<-1.5
    pos<-pos+cumsum(c(0,posStudy[-length(posStudy)]!=posStudy[-1]))*studySpace
    plotMini(thisDat[,names(minorFits)[ii]],thisDat$label,thisDat$pat,pos,thisDat$class,thisDat$study,ylab=ylabs[names(minorFits)[ii]],log=ifelse(names(minorFits)[ii]=='p24Release','','y'),ylim=if(names(minorFits)[ii]=='p24Release')c(0,1) else NULL)
    text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.001,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),LETTERS[(ii-1)*2+1],xpd=NA,adj=c(0,1),cex=2)
    par(mgp=c(2.2,1,0))
    plotSummary(minorFits[[ii]],addAcute=FALSE,ylab='Rebound fold change from VOA',mar=c(2.1,5,.1,.1),xaxis=FALSE,logYAxis=names(minorFits)[ii]=='infectivity')
    text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.02,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),LETTERS[(ii)*2],xpd=NA,adj=c(0,1),cex=2)
  }
dev.off()


meanCI<-function(xx,ci=c(.025,.975))c(mean(xx),quantile(xx,ci))
mat<-as.matrix(fitA$fit)
matB<-as.matrix(fitB$fit)
voaVsChron<-apply(mat[,c('stateMeans[2]','stateMeans[4]')],1,diff)
voaVsChronB<-apply(matB[,c('stateMeans[2]','stateMeans[4]')],1,diff)

voaVsReb<-apply(mat[,c('stateMeans[2]','stateMeans[3]')],1,diff)
voaVsRebB<-apply(matB[,c('stateMeans[2]','stateMeans[3]')],1,diff)
randomRebound<-rnorm(nrow(mat)*3,mat[,'stateMeans[1]']+mat[,'stateMeans[3]'],sqrt(mat[,'stateSds[1]']^2+mat[,'stateSds[3]']^2+mat[,'stateIsoSds[3]']^2))
randomVoaB<-rnorm(nrow(matB)*3,matB[,'stateMeans[1]']+matB[,'stateMeans[2]'],sqrt(matB[,'stateSds[1]']^2+matB[,'stateSds[2]']^2+matB[,'stateIsoSds[2]']^2))
randomVoa<-rnorm(nrow(mat)*3,mat[,'stateMeans[1]']+mat[,'stateMeans[2]'],sqrt(mat[,'stateSds[1]']^2+mat[,'stateSds[2]']^2+mat[,'stateIsoSds[2]']^2))
randomReboundB<-rnorm(nrow(matB)*3,matB[,'stateMeans[1]']+matB[,'stateMeans[3]'],sqrt(matB[,'stateSds[1]']^2+matB[,'stateSds[3]']^2+matB[,'stateIsoSds[3]']^2))
message('Voa vs chronic')
exp(meanCI(voaVsChron))
exp(meanCI(voaVsChronB))
message('Voa vs reb')
exp(meanCI(voaVsReb))
exp(meanCI(voaVsRebB))
stateCrI<-apply(mat[,grep('stateMeans',colnames(mat))],2,meanCI)
colnames(stateCrI)<-names(fitA$states)
message('alpha')
print(exp(stateCrI))
print(exp(-stateCrI))
stateCrIB<-apply(matB[,grep('stateMeans',colnames(mat))],2,function(xx)c(mean(xx),quantile(xx,c(.025,.975))))
colnames(stateCrIB)<-names(fitB$states)
message('beta')
print(exp(stateCrIB))
print(exp(-stateCrIB))

message('Mean IFNa2 acute:',format(exp(mean(mat[,'acuteMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(mat[,'acuteMean'],c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNb acute:',format(exp(mean(matB[,'acuteMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(matB[,'acuteMean'],c(.025,.975))),digits=4),collapse=' - '))
randomAcute<-rnorm(nrow(mat)*3,mat[,'acuteMean'],sqrt(mat[,'acuteSD']^2+mat[,'sigma']^2))


lapply(minorFits,function(xx)exp(meanCI(as.matrix(xx$fit)[,'stateMeans[2]'])))
lapply(minorFits,function(xx)exp(meanCI(as.matrix(xx$fit)[,'stateMeans[2]'])))
