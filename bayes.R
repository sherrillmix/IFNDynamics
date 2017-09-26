if(!exists('dat'))source('readNewData.R')
library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
nThreads<-40

stanCode<-'
  data {
    int<lower=0> alphaN;
    int<lower=0> betaN;
    int<lower=0> nPat;
    int<lower=0> alphaPatientId[alphaN];
    int<lower=0> betaPatientId[betaN];
    int<lower=0> alphaDay[alphaN];
    int<lower=0> betaDay[betaN];
    real ic50Alpha[alphaN];
    real ic50Beta[betaN];
    vector[alphaN] alphaVl;
    vector[betaN] betaVl;
    vector<lower=0>[alphaN] alphaCd4;
    vector<lower=0>[betaN] betaCd4;
    int<lower=0> nTotals[nPat];
    int<lower=0> totalStarts[nPat];
    int<lower=0> alphaStartObs[nPat];
    int<lower=0> alphaNObs[nPat];
    //int<lower=0> betaStartObs[nPat];
    //int<lower=0> betaNObs[nPat];
    int<lower=0> alphaDayIndex[alphaN];
    //int<lower=0> betaDayIndex[betaN];
  }
  parameters {
    vector[sum(nTotals)] trueAlpha;
    vector[sum(nTotals)] trueBeta;
    real<lower=0> alphaAutoCor;
    real<lower=0> alphaError;
    real<lower=0> betaAutoCor;
    real<lower=0> betaError;
    //real<lower=0> ic50AutoCor;
    //real<lower=0> ic50BetaAutoCor;
    real metaAlpha;
    real metaAlphaSd;
    vector[nPat] rawAlphas;
    real metaAlpha2;
    real metaAlpha2Sd;
    vector[nPat] rawAlphas2;
  }
  transformed parameters{
    vector[nPat] betaAlphas;
    vector[nPat] betaAlphas2;
    betaAlphas = metaAlpha + rawAlphas*metaAlphaSd;
    betaAlphas2 = metaAlpha2 + rawAlphas2*metaAlpha2Sd;
  }
  model {
    rawAlphas~normal(0,1);
    rawAlphas2~normal(0,1);
    alphaAutoCor~gamma(1,1);
    alphaError~gamma(1,1);
    betaAutoCor~gamma(1,1);
    betaError~gamma(1,1);
    trueAlpha[totalStarts]~normal(0,3);
    trueBeta[totalStarts]~normal(0,3);
    for(ii in 1:nPat){
      segment(trueAlpha,totalStarts[ii]+1,nTotals[ii]-1)~normal(segment(trueAlpha,totalStarts[ii],nTotals[ii]-1),alphaAutoCor);
      segment(trueBeta,totalStarts[ii]+1,nTotals[ii]-1)~normal(segment(trueBeta,totalStarts[ii],nTotals[ii]-1),betaAutoCor);
      //should linear interpolate
      segment(ic50Alpha,alphaStartObs[ii],alphaNObs[ii])~normal(betaAlphas[ii]*segment(dayTotalIndex,alphaStartObs[ii],alphaNObs[ii])+trueAlpha[segment(dayTotalIndex,alphaStartObs[ii],alphaNObs[ii])],alphaError);
      //segment(ic50Beta,betaStartObs[ii],betaNObs[ii])~normal(trueBeta[segment(dayTotalIndex,sStartObstartObs[ii],betaNObs[ii])],betaError);
    }
  }
'

dat$patId<-as.numeric(as.factor(dat$pat))
nObs<-tapply(dat$time,dat$patId,length)
dat$lookupId<-dat$time+dayStart[dat$patId]-1
dat$week<-ceiling(dat$time/7)
nTotals<-tapply(dat$week,dat$patId,max)
totalStarts<-cumsum(c(1,nTotals))[1:max(dat$patId)]


days<-dat[!duplicated(dat$sample),c('patId','time','CD4','vl')]
days$week<-ceiling(days$time/7)
days$totalIndex<-totalStarts[days$patId]+days$week-1
nDays<-tapply(days$patId,days$patId,length)
dayStart<-cumsum(c(1,nDays))[1:max(dat$patId)]
allDays<-unlist(lapply(nTotals,function(xx)1:xx))

alpha<-dat[!is.na(data$ic50),]
beta<-dat[!is.na(data$beta),]

input<-withAs('xx'=dat[order(dat$patId,dat$time),],list(
  N=nrow(xx),
  nPat=max(xx$patId),
  patientId=xx$patId,
  day=xx$time,
  ic50=log10(xx$ic50),
  ic50Beta=log10(xx$beta),
  starts=cumsum(c(1,nObs))[1:max(dat$patId)],
  nObs=nObs,
  #lookupIds=xx$lookupId,
  nTotals=nTotals,
  nDays=nDays,
  #days=days$time,
  startDays=cumsum(c(1,nDays))[1:max(dat$patId)],
  dayTotalIndex=days$totalIndex,
  days=days$week,
  totalStarts=totalStarts,
  #dayPatientId=days$patId,
  vl=log10(days$vl),
  cd4=days$CD4
))

fit <- stan(model_code = stanCode, data = input, iter=2000, chains=nThreads,thin=10)

sims<-extract(fit)
pdf('test.pdf',width=20)
plot(apply(sims[['trueVl']],2,mean),type='l')
lines(apply(sims[['trueVl']],2,quantile,.95),col='blue')
lines(apply(sims[['trueVl']],2,quantile,.05),col='blue')
abline(v=input$totalStarts,col='red')
plot(apply(sims[['trueCd4']],2,mean),type='l')
lines(apply(sims[['trueCd4']],2,quantile,.95),col='blue')
lines(apply(sims[['trueCd4']],2,quantile,.05),col='blue')
abline(v=input$totalStarts,col='red')
dev.off()
