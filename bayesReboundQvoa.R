#change to   
library('rstan')
library(dnar)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('functions.R')
source('rebound/functions.R')

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

source('readReboundData.R')
fitBayes2<-function(model,patient,states,studies,ic50,baseState='Acute',mixState='Post-ATI',chains=50,baseStudy='Other',logFunc=log,stateId=NULL,iter=6000,...){
  patientId<-structure(1:length(unique(patient)),.Names=sort(unique(patient)))
  studyId<-structure(1:length(unique(studies)),.Names=unique(studies[order(studies!=baseStudy)]))
  if(is.null(stateId))stateId=structure(1:length(unique(states)),.Names=unique(states[order(states!=baseState,states==mixState)]))
  if(any(!states %in% names(stateId)))stop('Missing state IDs')
  patientStudy<-tapply(studies,patient,unique)
  if(any(sapply(patientStudy,length)!=1))stop('Multiple studies for a single patient')
  dat=list(
    nVirus=length(ic50),
    nPatient=max(patientId),
    nState=max(stateId)-1,
    nStudy=max(studyId),
    ic50=logFunc(ic50),
    patients=patientId[patient],
    states=stateId[states],
    studies=studyId[patientStudy[names(patientId)]]
  )
  fit <- sampling(model, data = dat, iter=iter, chains=chains,thin=2,control=list(adapt_delta=.99,max_treedepth=15),...)
  return(list('fit'=fit,pats=patientId,states=stateId,studies=studyId,patientStudy=patientStudy,dat=dat))
}
fitA_withMix<-withAs(combined=combined[!is.na(combined$ic50_IFNa2),],fitBayes2(mod4,combined$pat,combined$simpleClass,ifelse(combined$study %in% c('Transmission','BEAT','IFNa2b treatment'),combined$study,'Other'),combined$ic50_IFNa2,chains=50,stateId=structure(1:5,.Names=c('Acute','Rebound','Chronic','Outgrowth','Post-ATI')),iter=6000))
fitB_withMix<-withAs(combined=combined[!is.na(combined$ic50_IFNb),],fitBayes2(mod4,combined$pat,combined$simpleClass,ifelse(combined$study %in% c('Transmission','BEAT','IFNa2b treatment'),combined$study,'Other'),combined$ic50_IFNb,chains=50,stateId=structure(1:5,.Names=c('Acute','Rebound','Chronic','Outgrowth','Post-ATI')),iter=6000))
betaAdjust<-exp(mean(as.matrix(fitB_withMix$fit)[,sprintf('studyMeans[%d]',fitB_withMix$studies['Transmission'])]))

pdf('test.pdf');plotSummary(fitA_withMix,addAcute=FALSE);plotSummary(fitB_withMix,addAcute=FALSE);dev.off()



ordering<-c('Pre-ATI A06','Post-ATI A06','Pre-ATI A08','Rebound A08','Post-ATI A08','Pre-ATI A09','Rebound A09','Post-ATI A09','Pre-ATI 9201','Rebound 9201','Pre-ATI 9202','Rebound 9202','Pre-ATI 9203','Rebound 9203','Pre-ATI 9207','Rebound 9207','Rebound A08','Rebound S-22','Rebound S-23','Rebound S-30','Rebound 601','Rebound BEAT-004','Rebound BEAT-030','Rebound BEAT-044','Outgrowth B106','Outgrowth B199','Outgrowth MM14','Outgrowth MM15','Outgrowth MM23','Outgrowth MM34','Outgrowth MM40','Outgrowth MM34','Acute','1 Year','Nadir','Last','Acute Recipients','Chronic Donors')
pos<-structure(1:length(unique(combined$label[!is.na(combined$label)])),.Names=unique(combined$label[!is.na(combined$label)][orderIn(combined$label[!is.na(combined$label)],ordering)]))
posStudy<-sapply(names(pos),function(xx)combined[combined$label==xx&!is.na(combined$label),'study'][1])
studySpace<-.5
pos<-pos+cumsum(c(0,posStudy[-length(posStudy)]!=posStudy[-1]))*studySpace

plotSummary<-function(fit,ylab='IFNa2 IC50 (pg/ml)',mar=c(6.9,4,.1,.1),xWidth=.4,addAcute=TRUE,xaxis=TRUE,logYAxis=TRUE,reps=2,cols=structure(rep('#00000033',nStates),.Names=stateNames),combine24=FALSE,cols2=cols,quantRange=c(.025,.975)){
  if(combine24){
    extraState<-tail(names(fit$states),1)
    fit$states<-fit$states[-length(fit$states)]
  }
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
  if(combine24){
    postProp<-mat[,'postProp']
    preProp<-mat[,'preProp']
    propRand<-runif(nrow(states)*reps,0,1)
    propRand2<-runif(nrow(states)*reps,0,1)
    states4<-states[,3]
    statesPat4<-statesPat[,3]
    statesIso4<-statesIso[,3]
    states[,3]<-ifelse(runif(nrow(states),0,1)<preProp,states[,1],states4)
    statesPat[,3]<-ifelse(propRand<preProp,statesPat[,1],statesPat4)
    statesIso[,3]<-ifelse(propRand<preProp,statesIso[,1],statesIso4)
    states<-cbind(states,states[,3])
    statesPat<-cbind(statesPat,statesPat[,3])
    statesIso<-cbind(statesPat,statesIso[,3])
    states[,4]<-ifelse(runif(nrow(states),0,1)<postProp,states[,1],states4)
    statesPat[,4]<-ifelse(propRand<postProp,statesPat[,1],statesPat4)
    statesIso[,4]<-ifelse(propRand<postProp,statesIso[,1],statesIso4)
    nStates<-4
    stateNames<-c(stateNames,extraState)
  }
  if(!addAcute)ylim<-exp(range(c(-1,1,apply(statesIso,2,quantile,c(quantRange[1],quantRange[2])))))
  makeDense<-function(xx)density(xx,from=quantile(xx,quantRange[1]),to=quantile(xx,quantRange[2]))
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
pdf('out/voaRebound_bayesSummary.pdf',height=3.5,width=9)
  layout(matrix(1:2,ncol=2),width=c(8,2.1))
  classCols<-structure(c(rep("#9EC0E1E6",3),rep("#77BCA9B3",2), rep("#9FB755B3",2), "#84C47DB3", "#B99A4BB3", "#C77C62B3", "#E581A0E6"), .Names = c("Outgrowth","Pre-ATI","Post-ATI","Acute Recipients","Acute", "Chronic Donors", "Chronic","1 Year", "Nadir", "Last", "Rebound")) 
  stCols<-sprintf('%s99',substring(classCols,1,7))
  stCols2<-sprintf('%s33',substring(classCols,1,7))
  names(stCols)<-names(stCols2)<-names(classCols)
  names(stCols2)[names(stCols2)=='Donor']<-names(stCols)[names(stCols)=='Donor']<-'Chronic'
  out<-withAs(combined=combined[!is.na(combined$ic50_IFNa2)&!is.na(combined$label),],
    plotQvoa2(combined$ic50_IFNa2,combined$label,pos,combined$displayClass,combined$study,combined$speed,ylab='IFNa2 IC50 (pg/ml)',mar=c(5.8,3.5,.1,.1),cex.axis=.9,startDown=TRUE,pats=ifelse(combined$study %in% c('Transmission','MM'),NA,combined$pat),classCols=classCols,labelXAxis=FALSE)
  )
  patPos<-tapply(out$pos,sub('BEAT-','',sub('.* ','',names(out$pos))),mean)
  patPos<-patPos[!grepl('^Acute|^Month|^Recipient|^Donor|^Nadir|^Year|^Last',names(patPos))]
  slantAxis(1,patPos,names(patPos),cex=.8)
  text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.0025,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'A',xpd=NA,adj=c(0,1),cex=2)
  #plotSummary(fitA)
  states<-plotSummary(fitA_withMix,addAcute=FALSE,ylab='Fold change from acute',mar=c(4.5,4,.1,1.6),cols=stCols,cols2=stCols2,xaxis=FALSE,combine24=TRUE)
  slantAxis(1,1:length(states),states,textOffsets=c(-.6,-.4,-.2,0),location=.7,axisArgs=list(tcl=-.4),cex=.9)
  #Add B
  text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.03,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'B',xpd=NA,adj=c(0,1),cex=2)
  #remove adjustment based on single IC50 and use bayesian estimated
  #combined$ic50_IFNb[combo$study=='Transmission']<-combo$beta[combo$study=='Transmission']/6386*2230
  #plotSummary(fitB,ylab='IFNb IC50 (pg/ml)')
  out<-withAs(combined=combined[!is.na(combined$ic50_IFNb)&!is.na(combined$label),],
    plotQvoa2(combined$ic50_IFNb/ifelse(combined$study=='Transmission',betaAdjust,1),combined$label,pos,combined$displayClass,combined$study,combined$speed,ylab='IFNb IC50 (pg/ml)',mar=c(5.8,3.5,.1,.1),cex.axis=.9,startDown=TRUE,pats=ifelse(combined$study %in% c('Transmission','MM'),NA,combined$pat),classCols=classCols,labelXAxis=FALSE)
  )
  patPos<-tapply(out$pos,sub('BEAT-','',sub('.* ','',names(out$pos))),mean)
  patPos<-patPos[!grepl('^Acute|^Month|^Recipient|^Donor|^Nadir|^Last',names(patPos))]
  slantAxis(1,patPos,names(patPos),cex=.8)
  text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.0025,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'C',xpd=NA,adj=c(0,1),cex=2)
  states<-plotSummary(fitB_withMix,addAcute=FALSE,ylab='Fold change from acute',mar=c(4.5,4,.1,1.6),cols=stCols,cols2=stCols2,xaxis=FALSE,combine24=TRUE)
  slantAxis(1,1:length(states),states,textOffsets=c(-.6,-.4,-.2,0),location=.7,axisArgs=list(tcl=-.4),cex=.9)
  text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.03,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),'D',xpd=NA,adj=c(0,1),cex=2)
dev.off()







fitBayes<-function(model,patient,states,studies,ic50,baseState='Acute',chains=50,baseStudy='Other',logFunc=log,...){
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
  fit <- sampling(model, data = dat, iter=6000, chains=chains,thin=2,control=list(adapt_delta=.99,max_treedepth=15),...)
  return(list('fit'=fit,pats=patientId,states=stateId,studies=studyId,dat=dat))
}
fitA2<-withAs(combined=combined[!is.na(combined$ic50_IFNa2),],fitBayes(mod3,combined$pat,combined$simpleClass,ifelse(combined$study %in% c('Transmission','BEAT','IFNa2b treatment'),combined$study,'Other'),combined$ic50_IFNa2,chains=50))
fitB2<-withAs(combined=combined[!is.na(combined$ic50_IFNb),],fitBayes(mod3,combined$pat,combined$simpleClass,ifelse(combined$study %in% c('Transmission','BEAT','IFNa2b treatment'),combined$study,'Other'),combined$ic50_IFNb,chains=50))
betaAdjust<-exp(mean(as.matrix(fitB2$fit)[,sprintf('studyMeans[%d]',fitB2$studies['Transmission'])]))

save(fitA2,fitB2,file='work/qvoaBayes_20200506.Rdat')


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
  withAs(combined=combined[!is.na(combined$ic50_IFNa2),],comparePreds(fitA2,combined$ic50_IFNa2,combined$virus))
  withAs(combined=combined[!is.na(combined$ic50_IFNb),],comparePreds(fitB2,combined$ic50_IFNb,combined$virus))
dev.off()



