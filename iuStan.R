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
iuMod <- stan_model(model_code = iuCode)

iuCodeSimple<-'
  data {
    int<lower=0> nCount;
    int<lower=0> counts[nCount];
    vector<lower=0>[nCount] dilutions;
    int maxCount;
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
    maxCount=2000,
    background=mean(background)
  )
  #,control=list(adapt_delta=.95,max_treedepth=15)
  fit <- sampling(mod, data = dat, iter=4000, chains=chains,thin=2,control=list(adapt_delta=.95,max_treedepth=15),...)
  return(fit)
}



print(system.time(fit<-withAs(xx=tit[!tit$red&tit$dextran=='Media'&tit$virus=='MM33.14.11B1'&!is.na(tit$virus),],simpleCountIU(iuModSimple,xx$n,xx$dil,tit[!tit$red&is.na(tit$dil),'n']))))

fit<-withAs(xx=tit[!tit$red&tit$dextran=='Dextran'&tit$virus=='YU2'&!is.na(tit$virus),],simpleCountIU(iuModSimple,xx$n,xx$dil,tit[!tit$red&is.na(tit$dil),'n']));mean(fit)

fit<-withAs(tit=tit[!tit$red,],countIU(iuMod,tit$n,tit$dil,tit$virus,tit[grepl('Media',tit$virus)&!tit$red,'n'],tit$treatSpin,chains=40))

fit2<-withAs(tit=tit[!tit$red&tit$virus %in% c('MM33.13.1D6','Media 1','Media 2'),],countIU(iuMod,tit$n,tit$dil,tit$virus,tit[grepl('Media',tit$virus)&!tit$red,'n'],tit$treatSpin,chains=40))

fit3<-withAs(xx=tit[!tit$red&tit$virus %in% c('YU2'),],countIU(iuMod,xx$n,xx$dil,xx$virus,tit[grepl('Media',tit$virus)&!tit$red,'n'],xx$treatSpin,chains=40))
#,'maxCount','background'
print(fit,pars=c('baseIU'))

rt<-read.csv('data/RTs 01.05.2018.csv',skip=1,stringsAsFactors=FALSE)
rt$old<-suppressWarnings(as.numeric(rt$old.RT)
rt$new<-suppressWarnings(as.numeric(rt$New.RT))

rt2<-read.csv('data/20190109_virusSubset.csv',stringsAsFactors=FALSE)
rt2<-rt2[!is.na(rt2$number),]
rt2$old<-as.numeric(rt2$old.RT)
rt2$new<-as.numeric(rt2$New.RT)

lims<-range(c(rt$old,rt$new,rt2$old,rt2$new),na.rm=TRUE)
pdf('out/oldNewRT.pdf',width=4,height=4)
  par(mar=c(3.75,3.75,.2,.2))
  dnar::withAs(rt=rt[grep('MM',rt$X.2),],plot(rt$old,rt$new,xlab='Old RT',ylab='New RT',log='xy',xlim=lims,ylim=lims,yaxt='n',xaxt='n',mgp=c(2.75,1,0)))
  dnar::logAxis(las=1)
  dnar::logAxis(1)
  points(rt2$old,rt2$new,col='blue')
  abline(0,1,lty=2)
  fit<-lm(new~0+old,dat=log(rbind(rt[,c('old','new')],rt2[,c('old','new')])))
  fakeDat<-data.frame(old=seq(log(lims[1]),log(lims[2]),length.out=1000))
  pred<-predict(fit,fakeDat,interval='conf')
  pred2<-predict(fit,fakeDat,interval='pred')
  lines(exp(fakeDat[,'old']),exp(pred[,'fit']))
  polygon(exp(c(fakeDat$old,rev(fakeDat$old))),exp(c(pred[,'lwr'],rev(pred[,'upr']))),col='#00000022',border=NA)
  polygon(exp(c(fakeDat$old,rev(fakeDat$old))),exp(c(pred2[,'lwr'],rev(pred2[,'upr']))),col='#00000011',border=NA)
dev.off()
