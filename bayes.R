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
    vector<lower=0>[alphaN] alphaDayReal;
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
    int<lower=0> betaStartObs[nPat];
    int<lower=0> betaNObs[nPat];
  }
  parameters {
    vector[sum(nTotals)] trueAlpha;
    vector[sum(nTotals)] trueBeta;
    real<lower=0> alphaAutoCor;
    real<lower=0> alphaError;
    real<lower=0> betaAutoCor;
    real<lower=0> betaError;
    real metaAlpha;
    real<lower=0> metaAlphaSd;
    vector[nPat] rawAlphas;
    real metaAlpha2;
    real<lower=0> metaAlpha2Sd;
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
    metaAlphaSd~gamma(1,.01);
    metaAlpha2Sd~gamma(1,.01);
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
      segment(ic50Alpha,alphaStartObs[ii],alphaNObs[ii])~normal(betaAlphas[ii]*segment(alphaDayReal,alphaStartObs[ii],alphaNObs[ii])+betaAlphas2[ii]*segment(alphaDayReal,alphaStartObs[ii],alphaNObs[ii]).*segment(alphaDayReal,alphaStartObs[ii],alphaNObs[ii])+trueAlpha[segment(alphaDay,alphaStartObs[ii],alphaNObs[ii])],alphaError);
      segment(ic50Beta,betaStartObs[ii],betaNObs[ii])~normal(trueBeta[segment(betaDay,betaStartObs[ii],betaNObs[ii])],betaError);
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

alpha<-dat[!is.na(dat$ic50),]
beta<-dat[!is.na(dat$beta),]

alphaNObs<-tapply(alpha$time,alpha$patId,length)
alphaStartObs<-cumsum(c(1,alphaNObs))[1:max(dat$patId)]
betaNObs<-tapply(beta$time,beta$patId,length)
betaStartObs<-cumsum(c(1,betaNObs))[1:max(dat$patId)]

input<-list(
  alphaN=nrow(alpha),
  betaN=nrow(beta),
  nPat=max(alpha$patId),
  alphaPatientId=alpha$patId,
  betaPatientId=beta$patId,
  alphaDay=alpha$week,
  alphaDayReal=alpha$week,
  betaDay=beta$week,
  ic50Alpha=log10(alpha$ic50),
  ic50Beta=log10(beta$beta),
  alphaVl=log10(alpha$vl),
  betaVl=log10(beta$vl),
  alphaCd4=alpha$CD4,
  betaCd4=beta$CD4,
  nTotals=nTotals,
  totalStarts=totalStarts,
  alphaStartObs=alphaStartObs,
  alphaNObs=alphaNObs,
  betaStartObs=betaStartObs,
  betaNObs=betaNObs
)

fit <- stan(model_code = stanCode, data = input, iter=10000, chains=nThreads,thin=20)

sims<-extract(fit)
pdf('test.pdf',width=20)
plot(apply(sims[['trueAlpha']],2,mean),type='l')
lines(apply(sims[['trueAlpha']],2,quantile,.95),col='blue')
lines(apply(sims[['trueAlpha']],2,quantile,.05),col='blue')
abline(v=input$totalStarts,col='red')
plot(apply(sims[['trueBeta']],2,mean),type='l')
lines(apply(sims[['trueBeta']],2,quantile,.95),col='blue')
lines(apply(sims[['trueBeta']],2,quantile,.05),col='blue')
abline(v=input$totalStarts,col='red')
dev.off()
