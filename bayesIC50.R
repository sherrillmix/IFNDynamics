library('rstan')
library(dnar)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if(!exists('dat'))source('readNewData.R')

ic50Code<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nArt;
    real ic50[nVirus];
    int<lower=0> patients[nVirus];
    int<lower=0> artIds[nVirus];
    real<lower=0> days[nVirus];
    real daysBeforeArt[nVirus];
    real artStart[nArt];
    real<lower=0> hasArt[nVirus];
  }
  parameters {
    vector[nPatient] acuteRaw;
    vector[nPatient] nadirTimeRaw;
    vector[nPatient] nadirChangeRaw;
    vector[nArt] riseTimeRaw;
    vector[nArt] riseChangeRaw;
    real<lower=0> sigma;
    real nadirTimeMean;
    real<lower=0> nadirTimeSD;
    real riseTimeMean;
    real<lower=0> riseTimeSD;
    real acuteMean;
    real<lower=0> acuteSD;
    real riseChangeMean;
    real<lower=0> riseChangeSD;
    real nadirChangeMean;
    real<lower=0> nadirChangeSD;
  }
  transformed parameters{
    vector[nPatient] acute;
    real expectedIC50[nVirus];
    vector[nPatient] nadirTime;
    vector[nPatient] nadirChange;
    vector[nArt] riseChange;
    vector[nArt] riseTime;
    acute=acuteMean+acuteRaw*acuteSD;
    nadirChange=nadirChangeMean+nadirChangeRaw*nadirChangeSD;
    nadirTime=nadirTimeMean+nadirTimeRaw*nadirTimeSD;
    riseTime=riseTimeMean+riseTimeRaw*riseTimeSD;
    riseChange=riseChangeMean+riseChangeRaw*riseChangeSD;
    for(ii in 1:nVirus){
      expectedIC50[ii]=acute[patients[ii]];
      if(days[ii]<exp(nadirTime[patients[ii]]))expectedIC50[ii]=expectedIC50[ii]+nadirChange[patients[ii]]*days[ii]/exp(nadirTime[patients[ii]]);
      else expectedIC50[ii]=expectedIC50[ii]+nadirChange[patients[ii]];
      if(hasArt[ii]){ // && days[ii]>exp(nadirTime[patients[ii]])
        if(daysBeforeArt[ii]<0){
          expectedIC50[ii]=expectedIC50[ii]+riseChange[artIds[ii]];
        }else{
          if(daysBeforeArt[ii]<exp(riseTime[artIds[ii]]))expectedIC50[ii]=expectedIC50[ii]+riseChange[artIds[ii]]*(1-daysBeforeArt[ii]/exp(riseTime[artIds[ii]]));
        }
      }
    }
  }
  model {
    ic50~normal(expectedIC50,sigma);
    sigma~gamma(1,.1);
    nadirTimeSD~gamma(1,.1);
    nadirTimeRaw~normal(0,1);
    riseTimeSD~gamma(1,.1);
    riseTimeRaw~normal(0,1);
    acuteSD~gamma(1,.1);
    acuteRaw~normal(0,1);
    nadirChangeSD~gamma(1,.1);
    nadirChangeRaw~normal(0,1);
    riseChangeSD~gamma(1,.1);
    riseChangeRaw~normal(0,1);
    riseChangeMean~normal(0,30);
    nadirChangeMean~normal(0,30);
    nadirTimeMean~normal(0,30);
    riseTimeMean~normal(0,30);
  }
'
ic50Mod <- stan_model(model_code = ic50Code)

bayesIC50<-function(mod,ic50,time,timePreArt,patient,chains=50,...){
  patientId<-structure(1:length(unique(patient)),.Names=sort(unique(patient)))
  artId<-structure(1:length(unique(patient[!is.na(timePreArt)])),.Names=sort(unique(patient[!is.na(timePreArt)])))
  artStart<-tapply(time+timePreArt,patient,unique)
  dat=list(
    nVirus=length(ic50),
    nArt=max(artId),
    nPatient=max(patientId),
    ic50=log(ic50),
    patients=patientId[patient],
    artIds=ifelse(is.na(artId[patient]),9999,artId[patient]),
    days=time,
    artStart=artStart[names(artId)],
    daysBeforeArt=ifelse(is.na(timePreArt),9999,timePreArt),
    hasArt=!is.na(timePreArt)
  )
  fit <- sampling(mod, data = dat, iter=3000, chains=chains,thin=3,control=list(adapt_delta=.9,max_treedepth=15),...)
  return(list('fit'=fit,pats=patientId,arts=artId))
}



calcSims<-function(fit,dat){
  mat<-as.matrix(fit[['fit']])
  pats<-fit$pats
  arts<-fit$arts
  maxTimes<-tapply(dat$time,dat$pat,max)
  artStart<-tapply(dat$time+dat$timeBeforArt,dat$pat,unique)
  sims<-lapply(names(pats),function(pat){
    fakeTimes<-1:maxTimes[[pat]]
    acute<-mat[,sprintf('acute[%d]',pats[pat])]
    nadirChange<-mat[,sprintf('nadirChange[%d]',pats[pat])]
    nadirTime<-exp(mat[,sprintf('nadirTime[%d]',pats[pat])])
    nadirProp<-fakeTimes %*% t(1/nadirTime)
    nadirProp[nadirProp>1]<-1
    predIc50<-acute+ nadirChange * t(nadirProp)
    if(pat %in% names(artStart[!is.na(artStart)])){
      riseChange<-mat[,sprintf('riseChange[%d]',arts[pat])]
      riseTime<-exp(mat[,sprintf('riseTime[%d]',arts[pat])])
      riseTime[riseTime==0]<-1e-9
      timeBeforeArt<-artStart[pat]-fakeTimes
      riseProp<-1-matrix(timeBeforeArt,nrow=length(riseTime),ncol=length(timeBeforeArt),byrow=TRUE)/riseTime
      riseProp[riseProp<0]<-0
      riseProp[riseProp>1]<-1
      nadTimes<-matrix(nadirTime,nrow=length(riseTime),ncol=length(fakeTimes))
      timeMat<-matrix(fakeTimes,nrow=length(riseTime),ncol=length(fakeTimes),byrow=TRUE)
      riseProp[timeMat<nadTimes]<-0
      predIc50<-predIc50 + riseChange * riseProp
    }
    sigmas<-mat[,'sigma']
    if(any(is.na(predIc50)))browser()
    predInt<-rbind(apply(predIc50-sigmas,2,quantile,.025),apply(predIc50+sigmas,2,quantile,.975))
    summaries<-apply(predIc50,2,function(xx)c('mean'=mean(xx),quantile(xx,c(.025,.975))))
    summaries<-rbind(summaries,predInt)
    if(any(summaries['mean',]< -1e+08))browser()
    rownames(summaries)<-c('mean','lowCI','highCI','lowPred','highPred')
    return(rbind('time'=fakeTimes,summaries))
  })
  names(sims)<-names(pats)
  return(sims)
}

fit<-withAs(xx=dat[!is.na(dat$ic50)&!dat$qvoa,],bayesIC50(ic50Mod,xx$ic50,xx$time,xx$timeBeforArt,xx$pat))
sims<-calcSims(fit,dat)

plotFit<-function(sims,dat){
  artStart<-tapply(dat$time+dat$timeBeforArt,dat$pat,unique)
  ylim<-range(sapply(sims,function(sim)range(exp(sim[-1,]))))
  xlim<-c(1,max(sapply(sims,function(sim)max(sim[1,]))))
  print(artStart)
  for(ii in names(sims)){
    sim<-sims[[ii]]
    plot(sim['time',],exp(sim['mean',]),type='l',ylim=ylim,log='y',ylab='IC50',xlab='Days after onset of symptoms',yaxt='n',main=ii,xlim=xlim)
    dnar::logAxis(las=1)
    polygon(c(sim['time',],rev(sim['time',])),exp(c(sim['lowCI',],rev(sim['highCI',]))),border=NA,col='#00000022')
    polygon(c(sim['time',],rev(sim['time',])),exp(c(sim['lowPred',],rev(sim['highPred',]))),border=NA,col='#00000022')
    dnar::withAs(xx=dat[dat$pat==ii,],points(xx$time,xx$ic50))
    abline(v=artStart[ii],lty=2)
    print(artStart[ii])
  }
}
pdf('test.pdf')
plotFit(sims,dat)
dev.off()

mat<-as.matrix(fit$fit)
cbind(exp(apply(mat[,grepl('riseTime\\[',colnames(mat))],2,mean)),c(artDfosx,'WEAU'=391))
exp(apply(mat[,grepl('nadirTime\\[',colnames(mat))],2,mean))
