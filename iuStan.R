library('rstan')
library(dnar)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

iuCode<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nCondition;
    int<lower=0> nCount;
    int<lower=0> counts[nCount];
    real<lower=0> dilutions[nCount];
    int<lower=1,upper=nCondition> conditions[nCount];
    int<lower=1,upper=nVirus> viruses[nCount];
    //int<lower=0> nBackground;
    //int<lower=0> backgrounds[nBackground];
    int maxCount;
    real background;
  }
  parameters {
    matrix[nVirus,nCondition] baseIU;
    //real<lower=100> maxCount;
    //real<lower=0> background;
  }
  transformed parameters{
    real<lower=0> expectedCount[nCount];
    for(ii in 1:nCount){
      expectedCount[ii]=exp(baseIU[viruses[ii],conditions[ii]])/dilutions[ii]+background;
    }
  }
  model {
    for(jj in 1:nVirus)baseIU[jj,]~normal(6,3);
    for(ii in 1:nCount){
      if(expectedCount[ii]<maxCount) counts[ii]~poisson(expectedCount[ii]);
      else target += poisson_lccdf(counts[ii]|expectedCount[ii]);
    }
    //maxCount~normal(1000,2000);
    //backgrounds~poisson(background);
    //background~normal(2,2);
  }
'
#iuMod <- stan_model(model_code = iuCode)

iuCodeSimple<-'
  data {
    int<lower=0> nCount;
    int<lower=0> counts[nCount];
    vector<lower=0>[nCount] dilutions;
    real background;
  }
  parameters {
    real<lower=0> baseIU;
    vector<upper=0>[nCount] overKill;
  }
  transformed parameters{
    vector<lower=0>[nCount] expectedCount;
    for(ii in 1:nCount)expectedCount[ii]=baseIU/dilutions[ii]+background;
  }
  model {
    //for(ii in 1:nCount){
      //if(expectedCount[ii]<maxCount || ii>nCount-3) counts[ii]~poisson(expectedCount[ii]);
      //else counts[ii]~poisson(expectedCount[ii]-overKill[ii]);
      //if(expectedCount[ii]-overKill[ii]>0)counts[ii]~poisson(expectedCount[ii]-overKill[ii]+background);
      //else counts[ii]~poisson(background);
      //if(expectedCount[ii]>maxCount){
        //overKill[ii]~double_exponential(0,expectedCount[ii]-maxCount);
      //}else{
        //overKill[ii]~double_exponential(0,1);
      //}
    //}
    baseIU~normal(0,100000);
    counts~poisson(expectedCount .*exp(overKill)+background);
    overKill~double_exponential(0,.1);
  }
'
iuModSimple <- stan_model(model_code = iuCodeSimple)


countIU<-function(mod,counts,dilutions,virus,background,conditions=rep("Base",length(virus)),chains=50,...){
  virusId<-structure(1:length(unique(virus)),.Names=unique(virus))
  conditionId<-structure(1:length(unique(conditions)),.Names=unique(conditions))
  dat=list(
    nVirus=max(virusId),
    nCondition=max(conditionId),
    nCount=length(counts),
    counts=counts,
    dilutions=dilutions,
    conditions=conditionId[conditions],
    viruses=virusId[virus],
    nBackground=length(background),
    backgrounds=background,
    maxCount=500,
    background=mean(background)
  )
  #,
  fit <- sampling(mod, data = dat, iter=5000, chains=chains,thin=10,control=list(adapt_delta=.99,max_treedepth=15),...)
}

simpleCountIU<-function(mod,counts,dilutions,background,chains=20,...){
  if(length(counts)==0)stop('No counts given')
  if(length(counts)!=length(dilutions))stop('Length of counts and dilutions do not match')
  dat=list(
    nCount=length(counts),
    counts=counts,
    dilutions=dilutions,
    background=mean(background)
  )
  #,control=list(adapt_delta=.95,max_treedepth=15)
  fit <- sampling(mod, data = dat, iter=4000, chains=chains,thin=2,control=list(adapt_delta=.95,max_treedepth=15),...)
  return(fit)
}



#print(system.time(fit<-withAs(xx=tit[!tit$red&tit$dextran=='Media'&tit$virus=='MM33.14.11B1'&!is.na(tit$virus),],simpleCountIU(iuModSimple,xx$n,xx$dil,tit[!tit$red&is.na(tit$dil),'n']))))



