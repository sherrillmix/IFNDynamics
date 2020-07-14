
library('rstan')
library(dnar)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('functions.R')
source('rebound/functions.R')

stanCode4_withMixture<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nState;
    int<lower=0> nStudy;
    real ic50[nVirus];
    int<lower=1,upper=nPatient> patients[nVirus];
    int<lower=1,upper=nState+1> states[nVirus];
    int<lower=1,upper=nStudy> studies[nPatient];
  }
  parameters {
    matrix[nPatient,nState] baseIc50Raw;
    vector[nState] stateMeans;
    vector[nStudy-1] studyMeansRaw;
    vector<lower=0>[nState] stateSds;
    vector<lower=0>[nState] stateIsoSds;
    real<lower=0,upper=1> postProp;
    real<lower=0,upper=1> preProp;
  }
  transformed parameters{
    matrix[nPatient,nState] expectedIC50;
    vector[nStudy] studyMeans;
    studyMeans[1]=0;
    studyMeans[2:nStudy]=studyMeansRaw;
    for(ii in 1:nPatient){
      for(jj in 1:(nState)){
        expectedIC50[ii,jj]=studyMeans[studies[ii]]+stateMeans[1]+baseIc50Raw[ii,1]*stateSds[1];
        if(jj>1)expectedIC50[ii,jj]=expectedIC50[ii,jj]+stateMeans[jj]+baseIc50Raw[ii,jj]*stateSds[jj];
      }
    }
  }
  model {
    for(ii in 1:nVirus){
      if(states[ii]<nState)ic50[ii]~normal(expectedIC50[patients[ii],states[ii]],stateIsoSds[states[ii]]);
      if(states[ii]==nState+1){
        target += log_sum_exp(
          log(postProp)+normal_lpdf(ic50[ii]|expectedIC50[patients[ii],2],stateIsoSds[2]),
          log(1-postProp)+normal_lpdf(ic50[ii]|expectedIC50[patients[ii],nState],stateIsoSds[nState])
        );
      }
      if(states[ii]==nState){
        target += log_sum_exp(
          log(preProp)+normal_lpdf(ic50[ii]|expectedIC50[patients[ii],2],stateIsoSds[2]),
          log(1-preProp)+normal_lpdf(ic50[ii]|expectedIC50[patients[ii],nState],stateIsoSds[nState])
        );
      }
    }
    postProp~beta(1,1);
    preProp~beta(1,1);
    stateIsoSds~gamma(1,1);
    stateSds~gamma(1,1);
    for(ii in 1:nState){
      baseIc50Raw[,ii]~normal(0,1);
    }
  }
'
mod4 <- stan_model(model_code = stanCode4_withMixture)

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

stanCodeNoDiff<-'
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
      expectedIC50[ii]=studyMeans[studies[ii]]+stateMeans[states[ii]]+baseIc50Raw[patients[ii],states[ii]]*stateSds[states[ii]];
      //if(states[ii]>1)expectedIC50[ii]=expectedIC50[ii]+stateMeans[states[ii]]+baseIc50Raw[patients[ii],states[ii]]*stateSds[states[ii]];
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
mod5 <- stan_model(model_code = stanCodeNoDiff)



source('readReboundData.R')

#note using log10 instead of log to make plotting easier
logit<-function(p)log10(p)-log10(1-p)
invLogit<-function(x)10^(x)/(10^(x)+1)

fitBayes<-function(model,patient,states,studies,ic50,baseState='Acute',chains=50,baseStudy='Other',logFunc=log,iter=6000,...){
  patientId<-structure(1:length(unique(patient)),.Names=sort(unique(patient)))
  stateId<-structure(1:length(unique(states)),.Names=unique(states[order(states!=baseState)]))
  studyId<-structure(1:length(unique(studies)),.Names=unique(studies[order(studies!=baseStudy)]))
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
  fit <- sampling(model, data = dat, iter=iter, chains=chains,thin=2,control=list(adapt_delta=.99,max_treedepth=15),...)
  return(list('fit'=fit,pats=patientId,states=stateId,studies=studyId,dat=dat))
}

minorFits<-lapply(structure(c('infectivity','repCap','p24Release'),.Names=c('infectivity','repCap','p24Release')),function(ii){
  thisDat<-combined[!is.na(combined[,ii])&!is.na(combined$voaVsRebound),]
  thisDat$preVsPost<-ifelse(thisDat$class=='Outgrowth','Pre-ATI',thisDat$class)
  fit<-fitBayes(mod5,thisDat$pat,thisDat$voaVsRebound,rep('A',nrow(thisDat)),thisDat[,ii],chains=50,baseState='Outgrowth',baseStudy='A',logFunc=if(ii %in% c('infectivity','repCap'))log else logit,iter=10000)
})

#pdf('out/voaRebound_minor_bayesCheck.pdf',width=20)
  #for(ii in 1:length(minorFits)){
    #thisDat<-combined[!is.na(combined[,names(combined)[ii]])&!is.na(combined$voaVsRebound),]
    #comparePreds(minorFits[[ii]],thisDat[,names(minorFits)[ii]],thisDat$ID,cols=ifelse(thisDat$class=='Rebound','red','blue'),log=ifelse(names(minorFits)[ii]=='p24Release','','y'),expFunc=ifelse(names(minorFits)[ii]=='p24Release',invLogit,exp),main=names(minorFits)[ii],ylab=names(minorFits)[ii])
  #}
#dev.off()

pdf('out/voaRebound_minor_bayesSummary.pdf',height=3.5)
  par(lheight=.73)
  for(ii in 1:length(minorFits)){
    par(mgp=c(2,1,0))
    plotSummary(minorFits[[ii]],addAcute=FALSE,ylab=ifelse(names(minorFits)[ii]=='p24Release','Fold change from VOA in odds ratio\nof release',sprintf('Fold change of %s\nfrom VOA',names(minorFits)[ii])),mar=c(2,4,.1,.1))
  }
dev.off()



ylabs<-c('infectivity'='Infectivity (IU/pg RT)','repCap'='Replicative capacity (ng p24/ml)','p24Release'='Proportion of p24 released')
plotMini<-function(ic50,label,pat,pos,class,study,ylab='IFNa2 IC50 (pg/ml)',mar=c(3,4,.1,.1),log='y',ylim=range(ic50)){
  spread<-vipor::offsetX(log10(ic50),label,width=.25)
  expFunc<-ifelse(log=='y',function(xx)10^xx,c)
  marSpace<-0
  par(mar=mar,lheight=.8)
  plot(pos[label]+spread,ic50,log=log,yaxt='n',ylab='',xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim)
  studyPos<-tapply(pos[label],study,mean)
  mtext(ylab,2,2.8,adj=.5,at=expFunc(par('usr')[3]+.4*diff(par('usr')[3:4])))
  #dnar::slantAxis(1,studyPos,names(studyPos))
  for(counter in 1:length(studyPos))axis(1,studyPos[counter],names(studyPos)[counter],padj=1,mgp=c(3,0+1.2*(1+counter)%%2,0),tcl=-.5+-1.2*(1+counter)%%2,cex.axis=1.1)
  if(log=='y'){
    if(diff(par('usr')[3:4])>1)logAxis(las=1)
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
pdf('out/repInfRel_qvoa_compare2.pdf',width=9,height=3.5)
  for(ii in 1:length(minorFits)){
    layout(matrix(1:2,ncol=2),width=c(8,2.2))
    thisDat<-combined[!is.na(combined[,names(minorFits)[ii]])&!is.na(combined$voaVsRebound),]
    thisDat<-thisDat[order(thisDat$study,thisDat$pat,thisDat$class=='Post-ATI',thisDat$class=='Rebound',thisDat$class=='Pre-ATI'),]
    pos<-structure(1:length(unique(thisDat$label)),.Names=unique(thisDat$label))
    posStudy<-sapply(names(pos),function(xx)thisDat[thisDat$label==xx,'study'][1])
    studySpace<-1.5
    pos<-pos+cumsum(c(0,posStudy[-length(posStudy)]!=posStudy[-1]))*studySpace
    #plotMini(thisDat[,names(minorFits)[ii]],thisDat$label,thisDat$pat,pos,thisDat$class,thisDat$study,ylab=ylabs[names(minorFits)[ii]],log=ifelse(names(minorFits)[ii]=='p24Release','','y'),ylim=if(names(minorFits)[ii]=='p24Release')c(0,1) else NULL)
    out<-plotQvoa2(thisDat[,names(minorFits)[ii]],thisDat$label,pos,thisDat$displayClass,thisDat$study,thisDat$speed,ylab=ylabs[names(minorFits)[ii]],mar=c(2.5,4.9,.1,.1),cex.axis=.9,startDown=TRUE,pats=thisDat$pat,classCols=classCols,labelXAxis=FALSE,log=ifelse(names(minorFits)[ii]=='p24Release','','y'),if(names(minorFits)[ii]=='p24Release')c(0,1) else range(thisDat[,names(minorFits)[ii]]))
    patPos<-tapply(out$pos,sub('BEAT-','',sub('.* ','',names(out$pos))),mean)
    patPos<-patPos[!grepl('^Acute|^Month|^Recipient|^Donor|^Nadir|^Year|^Last',names(patPos))]
    slantAxis(1,patPos,names(patPos),cex=.8,location=.8)
    text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.001,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),LETTERS[(ii-1)*2+1],xpd=NA,adj=c(0,1),cex=2)
    par(mgp=c(2.2,1,0))
    plotSummary(minorFits[[ii]],addAcute=FALSE,ylab='Fold change from outgrowth',mar=c(.5,5,.1,.1),xaxis=FALSE,logYAxis=names(minorFits)[ii]=='infectivity',subtract12=TRUE)
    text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.02,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),LETTERS[(ii)*2],xpd=NA,adj=c(0,1),cex=2)
  }
dev.off()
system('pdfjam out/repInfRel_qvoa_compare2.pdf --nup 1x3 --outfile tmp.pdf;pdfcrop tmp.pdf out/Fig._5.pdf')


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
pdf('out/voaRebound_minor_bayesCheck.pdf',width=20)
  for(ii in 1:length(minorFits)){
    thisDat<-combined[!is.na(combined[,names(minorFits)[ii]])&!is.na(combined$voaVsRebound),]
    comparePreds(minorFits[[ii]],thisDat[,names(minorFits)[ii]],thisDat$ID,cols=ifelse(thisDat$class=='Rebound','red','blue'),log=ifelse(names(minorFits)[ii]=='p24Release','','y'),expFunc=ifelse(names(minorFits)[ii]=='p24Release',invLogit,exp),main=names(minorFits)[ii],ylab=names(minorFits)[ii])
  }
dev.off()


zz<-as.matrix(minorFits$p24Release$fit)
exp(-meanCrI(zz[,'stateMeans[1]']-zz[,'stateMeans[2]']))
