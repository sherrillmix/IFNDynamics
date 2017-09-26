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
    vector<lower=0>[alphaN] alphaDay;
    vector<lower=0>[betaN] betaDay;
    vector<lower=0>[alphaN] alphaDay2;
    vector<lower=0>[betaN] betaDay2;
    real ic50Alpha[alphaN];
    real ic50Beta[betaN];
    //vector[alphaN] alphaVl;
    //vector[betaN] betaVl;
    //vector<lower=0>[alphaN] alphaCd4;
    //vector<lower=0>[betaN] betaCd4;
  }
  parameters {
    real<lower=0> alphaError;
    real<lower=0> betaError;
    real metaAlpha;
    real metaAlpha2;
    real<lower=0> metaAlphaSd;
    real<lower=0> metaAlpha2Sd;
    vector[nPat] rawAlphas;
    vector[nPat] rawAlphas2;
    real metaBeta;
    real metaBeta2;
    real<lower=0> metaBetaSd;
    real<lower=0> metaBeta2Sd;
    vector[nPat] rawBetas;
    vector[nPat] rawBetas2;
    vector[nPat] alphaIntercepts;
    vector[nPat] betaIntercepts;
  }
  transformed parameters{
    vector[nPat] betaAlphas;
    vector[nPat] betaAlphas2;
    vector[nPat] betaBetas;
    vector[nPat] betaBetas2;
    vector[alphaN] alphaMus;
    vector[betaN] betaMus;
    betaAlphas = metaAlpha + rawAlphas*metaAlphaSd;
    betaAlphas2 = metaAlpha2 + rawAlphas2*metaAlpha2Sd;
    betaBetas = metaBeta + rawBetas*metaBetaSd;
    betaBetas2 = metaBeta2 + rawBetas2*metaBeta2Sd;
    alphaMus=alphaIntercepts[alphaPatientId];//betaAlphas2[alphaPatientId].*alphaDay2+betaAlphas[alphaPatientId].*alphaDay;
    betaMus=betaIntercepts[betaPatientId];//+betaBetas2[betaPatientId].*betaDay2+betaBetas[betaPatientId].*betaDay;
  }
  model {
    rawAlphas~normal(0,1);
    rawAlphas2~normal(0,1);
    rawBetas~normal(0,1);
    rawBetas2~normal(0,1);
    metaAlphaSd~gamma(1,1);
    metaAlpha2Sd~gamma(1,1);
    metaBetaSd~gamma(1,1);
    metaBeta2Sd~gamma(1,1);
    alphaError~gamma(1,10);
    betaError~gamma(1,10);
    ic50Alpha~normal(alphaMus,alphaError);
    ic50Beta~normal(betaMus,betaError);
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
  alphaDay=alpha$time,
  betaDay=beta$time,
  alphaDay2=alpha$time^2,
  betaDay2=beta$time^2,
  ic50Alpha=log10(alpha$ic50),
  ic50Beta=log10(beta$beta)
  #alphaVl=log10(alpha$vl),
  #betaVl=log10(beta$vl),
  #alphaCd4=alpha$CD4,
  #betaCd4=beta$CD4
)

fit <- stan(model_code = stanCode, data = input, iter=2000, chains=nThreads,thin=10)

sims<-extract(fit)
pdf('test.pdf',width=20)
print(traceplot(fit,c('metaBeta','metaAlpha','metaBeta2','metaAlpha2')))
dev.off()
