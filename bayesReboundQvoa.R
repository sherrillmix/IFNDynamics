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
fitA_withMix<-withAs(combined=combined[!is.na(combined$ic50_IFNa2),],fitBayes2(mod4,combined$pat,combined$simpleClass,ifelse(combined$study %in% c('Transmission','BEAT','IFNa2b treatment'),combined$study,'Other'),combined$ic50_IFNa2,chains=50,stateId=structure(1:5,.Names=c('Acute','Rebound','Chronic','Outgrowth','Post-ATI')),iter=20000))
fitB_withMix<-withAs(combined=combined[!is.na(combined$ic50_IFNb),],fitBayes2(mod4,combined$pat,combined$simpleClass,ifelse(combined$study %in% c('Transmission','BEAT','IFNa2b treatment'),combined$study,'Other'),combined$ic50_IFNb,chains=50,stateId=structure(1:5,.Names=c('Acute','Rebound','Chronic','Outgrowth','Post-ATI')),iter=20000))
betaAdjust<-exp(mean(as.matrix(fitB_withMix$fit)[,sprintf('studyMeans[%d]',fitB_withMix$studies['Transmission'])]))
#save(fitA_withMix,fitB_withMix,file='work/mixFits_2020-05-14.Rdat')
load('work/mixFits_2020-05-14.Rdat')

pdf('test.pdf');plotSummary(fitA_withMix,addAcute=FALSE);plotSummary(fitB_withMix,addAcute=FALSE);dev.off()


zz<-as.matrix(fitB_withMix$fit)
meanCrI<-function(xx)c(quantile(xx,c(.025,.975)),mean(xx))[c(1,3,2)]
exp(-apply(zz[,grepl('stateMeans\\[4\\]',colnames(zz)),drop=FALSE],2,meanCrI))
exp(-apply(zz[,'postProp',drop=FALSE]*zz[,'stateMeans[2]']+(1-zz[,'postProp',drop=FALSE])*zz[,'stateMeans[4]'],2,meanCrI))
exp(-apply(zz[,'preProp',drop=FALSE]*zz[,'stateMeans[2]']+(1-zz[,'preProp',drop=FALSE])*zz[,'stateMeans[4]'],2,meanCrI))
mean(zz[,'preProp',drop=FALSE]*zz[,'stateMeans[2]']+(1-zz[,'preProp',drop=FALSE])*zz[,'stateMeans[4]']>zz[,'postProp',drop=FALSE]*zz[,'stateMeans[2]']+(1-zz[,'postProp',drop=FALSE])*zz[,'stateMeans[4]'])


#ordering<-c('Pre-ATI A06','Post-ATI A06','Pre-ATI A08','Rebound A08','Post-ATI A08','Pre-ATI A09','Rebound A09','Post-ATI A09','Pre-ATI 9201','Rebound 9201','Pre-ATI 9202','Rebound 9202','Pre-ATI 9203','Rebound 9203','Pre-ATI 9207','Rebound 9207','Rebound A08','Rebound S-22','Rebound S-23','Rebound S-30','Rebound 601','Rebound BEAT-004','Rebound BEAT-030','Rebound BEAT-044','Outgrowth B106','Outgrowth B199','Outgrowth MM14','Outgrowth MM15','Outgrowth MM23','Outgrowth MM34','Outgrowth MM40','Outgrowth MM34','Acute','1 Year','Nadir','Last','Acute Recipients','Chronic Donors')
ordering<-c('Pre-ATI A06','Post-ATI A06','Pre-ATI A08','Rebound A08','Post-ATI A08','Pre-ATI A09','Rebound A09','Post-ATI A09','Pre-ATI 9241','Rebound 9241','Pre-ATI 9242','Rebound 9242','Pre-ATI 9243','Rebound 9243','Pre-ATI 9244','Rebound 9244','Rebound A08','Rebound S22','Rebound S23','Rebound S30','Rebound 601','Rebound 004','Rebound 030','Rebound 044','Outgrowth B106','Outgrowth B199','Outgrowth MM14','Outgrowth MM15','Outgrowth MM23','Outgrowth MM34','Outgrowth MM40','Outgrowth MM34','Acute','Chronic','Acute Recipients','Chronic Donors')
pos<-structure(1:length(unique(combined$label[!is.na(combined$label)])),.Names=unique(combined$label[!is.na(combined$label)][orderIn(combined$label[!is.na(combined$label)],ordering)]))
posStudy<-sapply(names(pos),function(xx)combined[combined$label==xx&!is.na(combined$label),'study'][1])
studySpace<-.5
pos<-pos+cumsum(c(0,posStudy[-length(posStudy)]!=posStudy[-1]))*studySpace
plotSummary<-function(fit,ylab='IFNa2 IC50 (pg/ml)',mar=c(6.9,4,.1,.1),xWidth=.4,addAcute=TRUE,xaxis=TRUE,logYAxis=TRUE,reps=2,cols=structure(rep('#00000033',nStates),.Names=stateNames),combine24=FALSE,cols2=cols,quantRange=c(.025,.975),subtract12=FALSE){
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
  }else if(subtract12){
    stateNames<-names(fit$states)[-2]
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
  if(subtract12){
    states[,1]<-states[,2]-states[,1]
    statesPat[,1]<-statesPat[,2]-statesPat[,1]
    statesIso[,1]<-statesIso[,2]-statesIso[,1]
    colName<-colnames(states)[2]
    states<-states[,-2,drop=FALSE]
    statesPat<-statesPat[,-2,drop=FALSE]
    statesIso<-statesIso[,-2,drop=FALSE]
    nStates<-nStates-1
  }
  if(!addAcute)ylim<-exp(range(c(-1,1,apply(statesIso,2,quantile,c(quantRange[1],quantRange[2])))))
  makeDense<-function(xx)density(xx,from=quantile(xx,quantRange[1]),to=quantile(xx,quantRange[2]))
  denses<-apply(states,2,makeDense)
  denses2<-apply(statesPat,2,makeDense)
  denses3<-apply(statesIso,2,makeDense)
  #denses2<-lapply(1:ncol(expects),function(ii)makeDense(rnorm(n*10,expects[,ii],sds[,fit$dat$states[ii]])))
  par(mar=mar)
  plot(1,1,log='y',yaxt='n',ylab=ylab,xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim,xlim=c(.5,nStates+.5),las=1,mgp=c(2.45,1,0))
  if(logYAxis)logAxis(las=1,mgp=c(3,.6,0))
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
pdf('out/voaRebound_bayesSummary.pdf',height=7,width=7)
  layout(rbind(1:2,0,3:4),width=c(8,2.2),height=c(1,.01,1))
  #classCols<-structure(c(rep("#9EC0E1E6",3),rep("#77BCA9B3",2), rep("#9FB755B3",2), "#84C47DB3", "#B99A4BB3", "#C77C62B3", "#E581A0E6"), .Names = c("Outgrowth","Pre-ATI","Post-ATI","Acute Recipients","Acute", "Chronic Donors", "Chronic","1 Year", "Nadir", "Last", "Rebound")) 
  classCols<-structure(c(rep(classCol['qvoa'],3),rep("#aaaaaa",2), rep("#DDDDDD",2), "#aaaaaa", "#aaaaaa", "#aaaaaa", classCol['rebound']), .Names = c("Outgrowth","Pre-ATI","Post-ATI","Acute Recipients","Acute", "Chronic Donors", "Chronic","1 Year", "Nadir", "Last", "Rebound")) 
  stCols<-sprintf('%s99',substring(classCols,1,7))
  stCols2<-sprintf('%s33',substring(classCols,1,7))
  names(stCols)<-names(stCols2)<-names(classCols)
  names(stCols2)[names(stCols2)=='Donor']<-names(stCols)[names(stCols)=='Donor']<-'Chronic'
  plotFunc<-function(combined,ic50Col,fit,ylab='IFNa2 IC50 (pg/ml)',letters=LETTERS[1:2]){
    cex.axis<-1
    out<-withAs(combined=combined[!is.na(combined[,ic50Col])&!is.na(combined$label),],
      plotQvoa2(combined[,ic50Col],combined$label,pos,combined$displayClass,combined$study,combined$speed,ylab=ylab,mar=c(5.5,3.5,.1,.1),cex.axis=cex.axis,startDown=TRUE,pats=ifelse(combined$study %in% c('Transmission','MM'),NA,combined$pat),classCols=classCols,labelXAxis=FALSE)
    )
    patPos<-tapply(out$pos,sub('BEAT-','',sub('.* ','',names(out$pos))),mean)
    patPos<-patPos[!grepl('^Acute|^Month|^Recipient|^Donor|^Nadir|^Year|^Last|^Chronic',names(patPos))]
    slantAxis(1,patPos,names(patPos),cex=cex.axis,location=.8)
    text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.0025,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),letters[1],xpd=NA,adj=c(0,1),cex=2)
    #plotSummary(fitA)
    states<-plotSummary(fit,addAcute=FALSE,ylab='Fold change from acute',mar=c(4.5,4,.1,1.35),cols=stCols,cols2=stCols2,xaxis=FALSE,combine24=TRUE)
    slantAxis(1,1:length(states),sub('Outgrowth','Pre-ATI',states),textOffsets=c(-.2,-.2,-.2,-.2),location=.7,axisArgs=list(tcl=-.4),srt=-45,cex=cex.axis)
    text(grconvertX(par('fig')[1]+diff(par('fig')[1:2])*.03,from='ndc'),grconvertY(par('fig')[4]-diff(par('fig')[3:4])*.005,from='ndc'),letters[2],xpd=NA,adj=c(0,1),cex=2)
  }
  plotFunc(combined,'ic50_IFNa2',fitA_withMix)
  tmp<-combined
  tmp$betaAdjust<-tmp$ic50_IFNb/ifelse(tmp$study=='Transmission',betaAdjust,1)
  plotFunc(tmp,'betaAdjust',fitB_withMix,ylab='IFNb IC50 (pg/ml)',letters=LETTERS[3:4])
  #remove adjustment based on single IC50 and use bayesian estimated
  #combined$ic50_IFNb[combo$study=='Transmission']<-combo$beta[combo$study=='Transmission']/6386*2230
  #plotSummary(fitB,ylab='IFNb IC50 (pg/ml)')
dev.off()
file.copy('out/voaRebound_bayesSummary.pdf','out/Fig._4.pdf',overwrite=TRUE)
#system('pdfjam out/voaRebound_bayesSummary.pdf --nup 1x2 --outfile tmp.pdf;pdfcrop tmp.pdf out/Fig._4.pdf')




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
