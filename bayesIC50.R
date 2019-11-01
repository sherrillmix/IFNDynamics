library('rstan')
library(dnar)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('functions.R')

#if(!exists('dat'))source('readNewData.R')
dat<-read.csv('out/allLongitudinal.csv',stringsAsFactors=FALSE)


ic50CodeWithFast<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nArt;
    real ic50[nVirus];
    int<lower=0,upper=nPatient> patients[nVirus];
    vector<lower=0,upper=1>[nPatient] isFast;
    int<lower=0> artIds[nVirus];
    real<lower=0> days[nVirus];
    real daysBeforeArt[nVirus];
    real artStart[nArt];
    real<lower=0> hasArt[nVirus];
    int<lower=0> nSample;
    int<lower=0,upper=nSample> sample[nVirus];
  }
  parameters {
    vector[nPatient] acuteRaw;
    vector[nPatient] nadirTimeRaw;
    vector[nPatient] nadirChangeRaw;
    //DEFSSUB
    //vector[nArt] riseTimeRaw;
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
    real fastChangeMean;
  }
  transformed parameters{
    vector[nPatient] acute;
    real expectedIC50[nVirus];
    vector[nPatient] nadirTime;
    vector[nPatient] nadirChange;
    vector[nArt] riseChange;
    vector[nArt] riseTime;
    acute=acuteMean+acuteRaw*acuteSD;
    nadirChange=fastChangeMean*isFast+nadirChangeMean+nadirChangeRaw*nadirChangeSD;
    nadirTime=nadirTimeMean+nadirTimeRaw*nadirTimeSD;
    //TRANSSUB
    //riseTime=riseTimeMean+riseTimeRaw*riseTimeSD;
    riseChange=riseChangeMean+riseChangeRaw*riseChangeSD;
    for(ii in 1:nVirus){
      expectedIC50[ii]=acute[patients[ii]];
      if(days[ii]<exp(nadirTime[patients[ii]]))expectedIC50[ii]=expectedIC50[ii]+nadirChange[patients[ii]]*days[ii]/exp(nadirTime[patients[ii]]);
      else{
        expectedIC50[ii]=expectedIC50[ii]+nadirChange[patients[ii]];
      }
      if(hasArt[ii]){
        if(daysBeforeArt[ii]<0){
          expectedIC50[ii]=expectedIC50[ii]+riseChange[artIds[ii]];
        }else{
          if(daysBeforeArt[ii]<exp(riseTime[artIds[ii]]))expectedIC50[ii]=expectedIC50[ii]+riseChange[artIds[ii]]*(1-daysBeforeArt[ii]/exp(riseTime[artIds[ii]]));
        }
      }
    }
  }
  model {
    //MODELSUB
    //riseTimeRaw~normal(0,1);
    ic50~normal(expectedIC50,sigma);
    sigma~gamma(1,.1);
    nadirTimeSD~gamma(1,.1);
    nadirTimeRaw~normal(0,1);
    riseTimeSD~gamma(1,.1);
    acuteSD~gamma(1,.1);
    acuteRaw~normal(0,1);
    nadirChangeSD~gamma(1,.1);
    nadirChangeRaw~normal(0,1);
    riseChangeSD~gamma(1,.1);
    riseChangeRaw~normal(0,1);
    riseChangeMean~normal(0,10);
    nadirChangeMean~normal(0,10);
    nadirTimeMean~normal(0,10);
    riseTimeMean~normal(0,10);
    fastChangeMean~normal(0,10);
  }
'

#ic50Mod <- stan_model(model_code = ic50Code)

bayesIC50<-function(mod,ic50,time,timePreArt,patient,chains=50,fastProgressors=c(),ic50Code,...){
  patientId<-structure(1:length(unique(patient)),.Names=sort(unique(patient)))
  artId<-structure(1:length(unique(patient[!is.na(timePreArt)])),.Names=sort(unique(patient[!is.na(timePreArt)])))
  artStart<-tapply(time+timePreArt,patient,unique)
  sample<-paste(patient,time,sep='_')
  sampleId<-structure(1:length(unique(sample)),.Names=unique(sample[order(patient,time)]))
  defs<-paste(sprintf('real<upper=%f> riseTimeRaw%d;',log(artStart[names(artId)]),1:length(artId)),collapse='\n')
  mods<-paste(sprintf('riseTimeRaw%d~normal(riseTimeMean,riseTimeSD);',1:length(artId)),collapse='\n')
  trans<-paste(sprintf('riseTime[%d]=riseTimeRaw%d;',1:length(artId),1:length(artId)),collapse='\n')
  filledCode<-sub('//MODELSUB',mods,sub('//TRANSSUB',trans,sub('//DEFSSUB',defs,ic50Code)))
  ic50Mod <- stan_model(model_code = filledCode)
  dat=list(
    nVirus=length(ic50),
    nArt=max(artId),
    nPatient=max(patientId),
    ic50=log(ic50),
    patients=patientId[patient],
    artIds=ifelse(is.na(artId[patient]),9999,artId[patient]),
    days=ifelse(time<0,0,time),
    artStart=artStart[names(artId)],
    nSample=max(sampleId),
    sample=sampleId[sample],
    daysBeforeArt=ifelse(is.na(timePreArt),9999,timePreArt),
    hasArt=!is.na(timePreArt),
    isFast=names(patientId) %in% fastProgressors
  )
  fit <- sampling(ic50Mod, data = dat, iter=30000, chains=chains,thin=2,control=list(adapt_delta=.99,max_treedepth=15),...)
  return(list('fit'=fit,pats=patientId,arts=artId,code=filledCode,artStart=artStart,sample=sampleId,dat=dat))
}
#fit<-withAs(xx=dat[!is.na(dat$ic50)&!dat$qvoa,],bayesIC50(ic50Mod,xx$ic50,xx$time,xx$timeBeforeArt,xx$pat,ic50Code=ic50CodeWithFast,fastProgressors=c('MM15','WEAU')))
#sims<-calcSims(fit,dat)
#pdf('test.pdf',width=4,height=8)
#noQvoa<-dat[!dat$qvoa,]
#plotCondenseIfn(noQvoa,noQvoa$ic50,ylab='IFNa2 IC50 (pg/ml)',sims=sims)
#dev.off()


calcSims<-function(fit,dat,riseAfter=FALSE){ #,fastProgressors=c()
  mat<-as.matrix(fit[['fit']])
  pats<-fit$pats
  arts<-fit$arts
  maxTimes<-tapply(dat$time[!dat$qvoa],dat$pat[!dat$qvoa],max)
  artStart<-fit$artStart
  sims<-lapply(names(pats),function(pat){
    fakeTimes<-1:max(c(maxTimes[[pat]],artStart[[pat]]),na.rm=TRUE)
    acute<-mat[,sprintf('acute[%d]',pats[pat])]
    nadirChange<-mat[,sprintf('nadirChange[%d]',pats[pat])]
    #if(!is.null(fastProgressors)&pat %in% fastProgressors){
      #fastChange<-mat[,'fastChangeMean']
      #nadirChange<-nadirChange+fastChange
    #}
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
      if(!riseAfter)riseProp[riseProp>1]<-1
      nadTimes<-matrix(nadirTime,nrow=length(riseTime),ncol=length(fakeTimes))
      timeMat<-matrix(fakeTimes,nrow=length(riseTime),ncol=length(fakeTimes),byrow=TRUE)
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

browser()
fit<-withAs(xx=dat[!is.na(dat$ic50)&!dat$qvoa,],bayesIC50(ic50Mod,xx$ic50,xx$time,xx$timeBeforeArt,xx$pat,ic50Code=ic50CodeWithFast,fastProgressors=c('MM15','WEAU')))
fitB<-withAs(xx=dat[!is.na(dat$beta)&!dat$qvoa,],bayesIC50(ic50Mod,xx$beta,xx$time,xx$timeBeforeArt,xx$pat,ic50Code=ic50CodeWithFast,fastProgressors=c('MM15','WEAU')))

fitOld<-withAs(xx=dat[!is.na(dat$ic50)&!dat$qvoa,],bayesIC50(ic50Mod,xx$ic50,xx$time,xx$timeBeforeArt,xx$pat,ic50Code=ic50Code,fastProgressors=c('MM15','WEAU')))

fitRep<-withAs(xx=dat[!is.na(dat$replication)&!dat$qvoa,],bayesIC50(ic50Mod,xx$replication,xx$time,xx$timeBeforeArt,xx$pat,ic50Code=ic50CodeWithFast,fastProgressors=c('MM15','WEAU')))
fitInf<-withAs(xx=dat[!is.na(dat$infectivityDextran)&!dat$qvoa,],bayesIC50(ic50Mod,xx$infectivityDextran,xx$time,xx$timeBeforeArt,xx$pat,ic50Code=ic50CodeWithFast,fastProgressors=c('MM15','WEAU')))
fitInf2<-withAs(xx=dat[!is.na(dat$infectivityMedia)&!dat$qvoa,],bayesIC50(ic50Mod,xx$infectivityMedia,xx$time,xx$timeBeforeArt,xx$pat,ic50Code=ic50CodeWithFast,fastProgressors=c('MM15','WEAU')))

save(fit,fitB,fitRep,fitInf,fitInf2,file='work/bayesIC50.Rdat')
sims<-calcSims(fit,dat) #,fastProgressors=c('MM15','WEAU')
simsB<-calcSims(fitB,dat)
simsR<-calcSims(fitRep,dat)
simsInf<-calcSims(fitInf,dat)
simsInf2<-calcSims(fitInf2,dat)

artStart<-fit$artStart
plotPointsLine<-function(dat,ic50,ii,ylab,addTitle=TRUE,sims=NULL,addFit=TRUE,filterAfter=TRUE){
  plot(dat$time/7,ic50,yaxt='n',log='y',bg=patCols[dat$pat],pch=21,type='n',xlab='',ylab=ylab,xaxt='n',cex=1.4)
  if(addTitle)title(ii,line=-1)
  thisDat<-dat[dat$pat==ii,]
  thisIc50<-ic50[dat$pat==ii]
  if(sum(!is.na(thisIc50))==0)return()
  if(addFit){
    if(is.null(sims)){
      thisFit<-lm(I(log(thisIc50))~time+time2,dat=thisDat)
      fakeDays<-(min(thisDat$time)):(max(thisDat$time)+50)
      fakeDf<-data.frame('time'=fakeDays,'time2'=fakeDays^2,'logTime'=log(fakeDays),'logTime2'=log(fakeDays)^2,'logTime3'=log(fakeDays)^3,'logTime4'=log(fakeDays)^4)
      predIc50<-predict(thisFit,fakeDf,interval='confidence')
      lines(fakeDays/7,exp(predIc50[,'fit']),col=patCols[ii])
      polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols2[ii],border=NA)
      predIc50<-predict(thisFit,fakeDf,interval='prediction')
      polygon(c(fakeDays/7,rev(fakeDays)/7),c(exp(predIc50[,'lwr']),rev(exp(predIc50[,'upr']))),col=patCols3[ii],border=NA)
    }else{
      sim<-sims[[ii]][,sims[[ii]]['time',]<ifelse(is.na(artStart[ii])||!filterAfter,Inf,artStart[ii])]
      polygon(c(sim['time',],rev(sim['time',]))/7,exp(c(sim['lowCI',],rev(sim['highCI',]))),border=NA,col=patCols2[ii])
      polygon(c(sim['time',],rev(sim['time',]))/7,exp(c(sim['lowPred',],rev(sim['highPred',]))),border=NA,col=patCols3[ii])
      lines(sim['time',]/7,exp(sim['mean',]),col=patCols[ii])
      if(ii!='WEAU'||!filterAfter)abline(v=artStart[ii]/7,lty=2) #MAGIC NUMBER. SUPPRESSING WEAU VERTICAL LINE SINCE DIFFERS FROM OTHERS
    }
  }
  points(thisDat$time/7,thisIc50,pch=21+thisDat$bulk,bg=patCols[ii])
}
plotCondenseIfn<-function(dat,ic50,ylab,showLegend=TRUE,sims=NULL,addFit=TRUE,filterAfter=TRUE){
  par(mar=c(0,0,0,0))
  layout(lay2,width=c(.5,rep(1,2),.3),height=c(.01,c(1,1,1,.2,1,.2,1),1.3))
  counter<-1
  for(ii in patOrder){
    plotPointsLine(dat,ic50,ii,ylab,sims=sims,addFit=addFit,filterAfter=filterAfter)
    if(counter>8)axis(1,(0:3)*100,cex.axis=1.2,mgp=c(2.75,.7,0))
    if(counter>8)axis(1,(0:2)*100+50,rep('',3),cex.axis=1.2,mgp=c(2.75,.7,0))
    if(counter%%2==1)logAxis(2,las=1,cex.axis=1.1,mgp=c(3,.7,0))
    if(counter==5)text(par('usr')[1]-.33*diff(par('usr')[1:2]),10^mean(par('usr')[3:4]),ylab,srt=90,xpd=NA,cex=2)
    if(counter==9)text(max(par('usr')[1:2]),10^(par('usr')[3]-.27*diff(par('usr')[3:4])),'Weeks after onset of symptoms',xpd=NA,cex=2)
    if(counter==9&showLegend)legend(par('usr')[1]-diff(par('usr')[1:2])*.3,10^(par('usr')[3]-diff(par('usr')[3:4])*.45),c(ifelse(is.null(sims),'Quadratic regression','Bayesian model'),ifelse(is.null(sims),'95% confidence interval','95% credible interval'),'95% prediction interval','Limiting dilution isolate','Bulk isolate'),col=c(patCols[1],NA,NA,'black','black'),pt.bg=c(NA,patCols2[1],patCols3[1],patCols[1],patCols[1]),lty=c(1,NA,NA,NA,NA),pch=c(NA,22,22,21,22),border=NA,pt.cex=c(3.2,3.2,3.2,1.4,1.4),cex=1.1,xjust=0,yjust=1,xpd=NA)
    counter<-counter+1
  }
}
plotFit<-function(sims,dat,var='ic50',ylab='IC50',filterAfter=TRUE){
  ylim<-range(dat[,var],na.rm=TRUE)
  xlim<-c(1,max(sapply(sims,function(sim)max(sim[1,]))))
  for(ii in names(sims)){
    if(filterAfter)sim<-sims[[ii]][,sims[[ii]]['time',]<ifelse(is.na(artStart[ii]),Inf,artStart[ii])]
    plot(sim['time',]/7,exp(sim['mean',]),type='l',ylim=ylim,log='y',ylab=ylab,xlab='Days after onset of symptoms',yaxt='n',main=ii,xlim=xlim/7,mgp=c(2.25,1,0))
    dnar::logAxis(las=1,mgp=c(1,.6,0))
    polygon(c(sim['time',],rev(sim['time',]))/7,exp(c(sim['lowCI',],rev(sim['highCI',]))),border=NA,col='#00000022')
    polygon(c(sim['time',],rev(sim['time',]))/7,exp(c(sim['lowPred',],rev(sim['highPred',]))),border=NA,col='#00000022')
    dnar::withAs(xx=dat[dat$pat==ii,],points(xx$time/7,xx[,var]))
    abline(v=artStart[ii]/7,lty=2)
  }
}
plotSummaries<-function(fit,var='IFN IC50'){ #,fastProgressors=c()
  mat<-as.matrix(fit$fit)
  artStart<-fit$artStart
  #,'fastChange'=sprintf('Fast progressor nadir fold change in\n %s from other progression',var)
  vars<-c('nadirTime'='Nadir time\n(weeks after onset of infection)','nadirChange'=sprintf('Nadir fold change in\n %s from acute',var),'riseTime'='Rise time\n(weeks before ART initiation)','riseChange'=sprintf('Rise fold change in\n %s from nadir',var),'acute'=sprintf('Base level in %s\nat acute infection',var)) 
  for(ii in names(vars)){
    nadirTimes<-c(colnames(mat)[grep(sprintf('%s\\[',ii),colnames(mat))],sprintf('%sMean',ii))
    nads<-mat[,nadirTimes,drop=FALSE]
    isFast<-any(grepl('fastChangeMean',colnames(mat)))&&ii=='nadirChange'
    if(isFast){
      nadirTimes<-c(nadirTimes,'fastProgressor')
      nads<-cbind(nads,mat[,'nadirChangeMean']+mat[,'fastChangeMean'])
      colnames(nads)<-nadirTimes
    }
    nadMeans<-apply(nads,2,mean)
    nadRanges<-apply(nads,2,quantile,c(.025,.975))
    if(max(nadRanges)>10000)browser()
    if(grepl('rise',ii)) pats<-sort(fit$arts)
    else pats<-sort(fit$pats)
    if(grepl('Change|acute',ii)){
      xlim<-exp(max(abs(nadRanges))*c(-1,1))
      scale<-1
    } else{
      scale<-7
      xlim<-exp(range(nadRanges))/scale
      xlim[2]<-min(1000,xlim[2])
    }
    plot(1,1,type='n',xlim=xlim,ylim=c(1,length(nadirTimes)),xlab=vars[ii],yaxt='n',ylab='',log=ifelse(grepl('Change',ii),'x',''),xaxt=ifelse(grepl('Change',ii),'n','s'))
    if(grepl('Change',ii))logAxis(1)
    axis(2,c(pats,length(pats)+1:(length(nadirTimes)-length(pats))),c(names(pats),if(isFast)c('Standard','Fast')else 'Overall'),las=1,adj=.5)
    segments(exp(nadRanges[1,])/scale,1:length(nadirTimes),exp(nadRanges[2,])/scale,1:length(nadirTimes),col=c('Overall'='black',patCols)[c(names(pats),rep('Overall',length(nadirTimes)-length(pats)))])
    points(exp(nadMeans)/scale,1:length(nadirTimes),pch=21,bg=c('Overall'='black',patCols)[c(names(pats),rep('Overall',length(nadirTimes)-length(pats)))])
    if(grepl('Change',ii))abline(v=1,lty=2)
  }
}
pdf('out/bayesSummary.pdf',width=4,height=8)
  par(mar=c(4.1,4.5,.1,.1))
  plotSummaries(fit,var='IFNa2 IC50') #,fastProgressors=c('MM15','WEAU')
  plotSummaries(fitB,var='IFNb IC50')
  plotSummaries(fitRep,var='replicative capacity')
  plotSummaries(fitInf,var='infectivity (with dextran)')
  #plotSummaries(fitInf2,var='infectivity (without dextran)')
dev.off()

fit350<-withAs(xx=dat[!is.na(dat$ic50)&!dat$qvoa,],bayesIC50(ic50CodeRiseAfter,xx$ic50,xx$time,xx$timeBefore350,xx$pat))
sims350<-calcSims(fit350,dat,riseAfter=TRUE)
artStart<-fit350$artStart
pdf('out/bayes350.pdf',width=4,height=8)
  plotSummaries(fit350,var='IFNa2 IC50')
  noQvoa<-dat[!dat$qvoa,]
  plotCondenseIfn(noQvoa,noQvoa$ic50,ylab='IFNa2 IC50 (pg/ml)',sims=sims350,filterAfter=FALSE)
dev.off()

pdf('out/bayesFit.pdf',width=4,height=8)
noQvoa<-dat#dat[!dat$qvoa,]
plotCondenseIfn(noQvoa,noQvoa$ic50,ylab='IFNa2 IC50 (pg/ml)',sims=sims,filterAfter=TRUE)
  text(grconvertX(.01,from='ndc'),grconvertY(.99,from='ndc'),'A',xpd=NA,adj=c(0,1),cex=2.5)
plotCondenseIfn(noQvoa,noQvoa$beta,ylab='IFNb IC50 (pg/ml)',sims=simsB,filterAfter=TRUE)
  text(grconvertX(.01,from='ndc'),grconvertY(.99,from='ndc'),'B',xpd=NA,adj=c(0,1),cex=2.5)
#plotCondenseIfn(noQvoa,noQvoa$replication,ylab='Replicative capacity (day 7 p24 ng/ml)',sims=simsR)
#plotCondenseIfn(noQvoa,noQvoa$infectivityDextran,ylab='Infectivity (IU/ng RT with dextran)',sims=simsInf)
#plotCondenseIfn(noQvoa,noQvoa$infectivityMedia,ylab='Infectivity (IU/ng RT without dextran)',sims=simsInf2)
dev.off()
pdf('out/bayesFitNoPred.pdf',width=4,height=8)
noQvoa<-dat[!dat$qvoa,]
plotCondenseIfn(noQvoa,noQvoa$ic50,ylab='IFNa2 IC50',sims=sims,addFit=FALSE)
plotCondenseIfn(noQvoa,noQvoa$beta,ylab='IFNb IC50',sims=simsB,addFit=FALSE)
#plotCondenseIfn(noQvoa,noQvoa$replication,ylab='Replicative capacity (day 7 p24 ng/ml)',sims=simsR,addFit=FALSE)
#plotCondenseIfn(noQvoa,noQvoa$infectivityDextran,ylab='Infectivity (IU/ng RT with dextran)',sims=simsInf,addFit=FALSE)
#plotCondenseIfn(noQvoa,noQvoa$infectivityMedia,ylab='Infectivity (IU/ng RT without dextran)',sims=simsInf2,addFit=FALSE)
dev.off()


exampleCurve<-function(fit,reps=10,isFast=FALSE,timeToArt=5,times=1:(timeToArt*365),confInt=.05){
  mat<-as.matrix(fit)
  generateMat<-function(mat,base,reps,addFast=FALSE){
    thisMean<-sprintf('%sMean',base)
    thisSd<-sprintf('%sSD',base)
    fastMean<-'fastChangeMean'
    isoSd<-'sigma'
    out<-data.frame('mean'=rep(mat[,thisMean],each=reps))
    if(addFast)out$mean<-out$mean+rep(mat[,fastMean],each=reps)
    out$pat<-rnorm(nrow(out),out$mean,rep(mat[,thisSd],each=reps))
    #out$iso<-rnorm(out$pat,rep(mat[,isoSd],each=reps))
    return(out)
  }
  acute<-generateMat(mat,'acute')
  nadirChange<-generateMat(mat,'nadirChange',reps,isFast)
  riseChange<-generateMat(mat,'riseChange',reps)
  nadirTime<-exp(generateMat(mat,'nadirTime',reps))
  riseTime<-exp(generateMat(mat,'riseTime',reps))
  sigma<-rep(mat[,'sigma'],each=reps)
  #acute<-data.frame(''=rep(mat[,'acuteMean'],'pat'=,'iso'=rnorm(nrow(mat)*reps,mat[,'acuteMean'],mat[,'acuteSD'])
  #nadirChange<-rnorm(nrow(mat)*reps,mat[,'nadirChangeMean']+if(isFast)mat[,'fastChangeMean'] else 0,mat[,'nadirChangeSD'])
  #riseChange<-rnorm(nrow(mat)*reps,mat[,'riseChangeMean'],mat[,'riseChangeSD'])
  #nadirTime<-exp(rnorm(nrow(mat)*reps,mat[,'nadirTimeMean'],mat[,'nadirTimeSD']))
  #riseTime<-exp(rnorm(nrow(mat)*reps,mat[,'riseTimeMean'],mat[,'riseTimeSD']))
  out<-lapply(structure(c('mean','pat','iso'),.Names=c('mean','pat','iso')),function(ii){
    if(ii=='iso')col<-'pat'
    else col<-ii
    nadirProp<-times %*% t(1/nadirTime[,col])
    nadirProp[nadirProp>1]<-1
    predIc50<-acute[,col] + nadirChange[,col] * t(nadirProp)
    if(!is.na(timeToArt)){
      #riseTime[riseTime==0]<-1e-9
      timeBeforeArt<-timeToArt*365-times
      riseProp<-1-matrix(timeBeforeArt,nrow=length(riseTime[,col]),ncol=length(timeBeforeArt),byrow=TRUE)/riseTime[,col]
      riseProp[riseProp<0]<-0
      riseProp[riseProp>1]<-1
      predIc50<-predIc50 + riseChange[,col] * riseProp
    }
    if(ii=='iso')predIc50[,]<-rnorm(prod(dim(predIc50)),unlist(predIc50),sigma)
    preds<-data.frame('time'=times,'mean'=exp(apply(predIc50,2,mean)),'lower'=exp(apply(predIc50,2,quantile,confInt/2)),'upper'=exp(apply(predIc50,2,quantile,1-confInt/2)))
    return(preds)
  })
  #return(list('times'=out,'params'=list('acute'=density(acute$mean),'nadirChange'=density(nadirChange$mean),'riseChange'=density(riseChange$mean),'nadirTime'=density(nadirTime$mean),'riseTime'=density(riseTime$mean))))
  return(out)
}
predIc50<-exampleCurve(fit$fit,reps=2)
predIc50Fast<-exampleCurve(fit$fit,isFast=TRUE,timeToArt=1.5,reps=2)
predIc50Slow<-exampleCurve(fit$fit,timeToArt=NA,times=seq(1,5*365),reps=2)

predIc50Beta<-exampleCurve(fitB$fit,reps=2)
predIc50FastBeta<-exampleCurve(fitB$fit,isFast=TRUE,timeToArt=1.5,reps=2)
predIc50SlowBeta<-exampleCurve(fitB$fit,timeToArt=NA,times=seq(1,5*365),reps=2)

plotExamples<-function(predIc50,predIc50Fast=NULL,predIc50Slow=NULL,ylab='IFNa2 IC50'){
  plot(1,1,type='n',xlim=c(1,5*365)/7,ylim=range(predIc50$pat[,c('mean','lower','upper')]),las=1,log='y',yaxt='n',ylab=ylab,xlab='',mgp=c(2.6,.7,0))
  title(xlab='Time after infection (weeks)',mgp=c(1.9,.7,0))
  logAxis(las=1)
  lines(predIc50$mean$time/7,predIc50$mean$mean,col='#FF8000')
  polygon(c(predIc50$mean$time,rev(predIc50$mean$time))/7,c(predIc50$mean$lower,rev(predIc50$mean$upper)),border='#FF800011',col='#FF800033')
  #polygon(c(predIc50$pat$time,rev(predIc50$pat$time))/7,c(predIc50$pat$lower,rev(predIc50$pat$upper)),border='#FF800011',col='#FF800011')
  polygon(c(predIc50$iso$time,rev(predIc50$iso$time))/7,c(predIc50$iso$lower,rev(predIc50$iso$upper)),border='#FF800011',col='#FF800011')
  if(!is.null(predIc50Fast)){
    lines(predIc50Fast$mean$time/7,predIc50Fast$mean$mean,col='#FF0000')
    polygon(c(predIc50Fast$mean$time,rev(predIc50Fast$mean$time))/7,c(predIc50Fast$mean$lower,rev(predIc50Fast$mean$upper)),border='#FF000011',col='#FF000033')
    #polygon(c(predIc50Fast$pat$time,rev(predIc50Fast$pat$time))/7,c(predIc50Fast$pat$lower,rev(predIc50Fast$pat$upper)),border='#FF000011',col='#FF000011')
    polygon(c(predIc50Fast$iso$time,rev(predIc50Fast$iso$time))/7,c(predIc50Fast$iso$lower,rev(predIc50Fast$iso$upper)),border='#FF000011',col='#FF000011')
  }
  if(!is.null(predIc50Slow)){
    lines(predIc50Slow$mean$time/7,predIc50Slow$mean$mean,col='#0000FF')
    polygon(c(predIc50Slow$mean$time,rev(predIc50Slow$mean$time))/7,c(predIc50Slow$mean$lower,rev(predIc50Slow$mean$upper)),border='#0000FF11',col='#0000FF11')
    #polygon(c(predIc50Slow$pat$time,rev(predIc50Slow$pat$time))/7,c(predIc50Slow$pat$lower,rev(predIc50Slow$pat$upper)),border='#0000FF11',col='#0000FF11')
    polygon(c(predIc50Slow$iso$time,rev(predIc50Slow$iso$time))/7,c(predIc50Slow$iso$lower,rev(predIc50Slow$iso$upper)),border='#0000FF11',col='#0000FF11')
  }
}
pdf('out/bayesExample.pdf',width=3.5,height=3.5)
  par(mar=c(3,3.5,.1,.1))
  plotExamples(predIc50,predIc50Fast,predIc50Slow)
  plotExamples(predIc50Beta,predIc50FastBeta,predIc50SlowBeta,ylab='IFNb IC50')
dev.off()

plotParams<-function(fit,totalTime=6,timeAfterInf=1,ylab='IFNa2 IC50'){
  mat<-as.matrix(fit)
  acute<-mat[,'acuteMean']
  fast<-mat[,'fastChangeMean']
  nadirChange<-mat[,'nadirChangeMean']
  riseChange<-mat[,'riseChangeMean']
  nadirTime<-exp(mat[,'nadirTimeMean'])/365
  riseTime<-exp(mat[,'riseTimeMean'])/365
  plot(1,1,type='n',xlim=c(-.3,totalTime),ylim=exp(quantile(c(acute,acute+nadirChange),c(.01,.99))),bty='n',xlab='',ylab=ylab,yaxt='n',log='y',xaxt='n')
  logAxis(las=1)
  mtext("Time after onset\nof symptoms",1,at=mean(c(0,timeAfterInf)),3.75)
  axis(1,0:timeAfterInf,mgp=c(2,1.75,1))
  axis(1,(timeAfterInf+1):totalTime,rev(totalTime-totalTime:(timeAfterInf+1)),mgp=c(2,1.75,1))
  mtext("Time prior to\n ART initiation",1,at=mean(c(timeAfterInf+1,totalTime)),3.75)
  makeDense<-function(xx)density(xx,from=quantile(xx,.025),to=quantile(xx,.975))
  acuteDen<-makeDense(exp(acute))
  nadirDen<-makeDense(exp(acute+nadirChange))
  nadirFastDen<-makeDense(exp(acute+nadirChange+fast))
  nadirTimeDen<-makeDense(nadirTime)
  riseTimeDen<-makeDense(riseTime)
  riseDen<-makeDense(exp(acute+nadirChange+riseChange))
  riseFastDen<-makeDense(exp(acute+nadirChange+fast+riseChange))
  xWidth<-diff(par('usr')[1:2])*.02
  yWidth<-diff(par('usr')[3:4])*.02
  riseMean<-exp(mean(log(riseTime)))
  nadirMean<-exp(mean(log(nadirTime)))
  nadChangeMean<-exp(mean(acute+nadirChange))
  nadChangeFastMean<-exp(mean(acute+nadirChange+fast))
  acuteMean<-exp(mean(acute))
  scaleFunc<-max
  riseChangeMean<-exp(mean(acute+nadirChange+riseChange))
  riseChangeFastMean<-exp(mean(acute+nadirChange+fast+riseChange))
  normCol='#FF800077'
  fastCol='#FF000077'
  bothCol='#00000033'
  #acute
  polygon(c(acuteDen$y,-rev(acuteDen$y))/scaleFunc(acuteDen$y)*xWidth,c(acuteDen$x,rev(acuteDen$x)),col=bothCol)
  segments(-xWidth,acuteMean,xWidth,acuteMean)
  #acute to nadir
  segments(0,acuteMean,nadirMean,nadChangeMean)
  #nadir
  polygon(nadirMean+c(nadirDen$y,-rev(nadirDen$y))/scaleFunc(nadirDen$y)*xWidth,c(nadirDen$x,rev(nadirDen$x)),col=normCol)
  segments(0,acuteMean,nadirMean,nadChangeFastMean)
  segments(-xWidth+nadirMean,nadChangeMean,xWidth+nadirMean,nadChangeMean)
  #nadirFast
  polygon(nadirMean+c(nadirFastDen$y,-rev(nadirFastDen$y))/scaleFunc(nadirFastDen$y)*xWidth,c(nadirFastDen$x,rev(nadirFastDen$x)),col=fastCol)
  segments(0,acuteMean,nadirMean,nadChangeFastMean)
  segments(-xWidth+nadirMean,nadChangeFastMean,xWidth+nadirMean,nadChangeFastMean)
  #nadirTime
  polygon(c(nadirTimeDen$x,rev(nadirTimeDen$x)),10^(par('usr')[3]-diff(par('usr')[3:4])*.03+ yWidth*c(nadirTimeDen$y,-rev(nadirTimeDen$y))/scaleFunc(nadirTimeDen$y)),xpd=NA,col=bothCol)
  segments(nadirMean,10^(par('usr')[3]-diff(par('usr')[3:4])*.03- yWidth),nadirMean,10^(par('usr')[3]-diff(par('usr')[3:4])*.03+ yWidth),xpd=NA)
  #connecting
  segments(nadirMean,nadChangeFastMean,totalTime-riseMean,nadChangeFastMean,lty=2)
  segments(nadirMean,nadChangeMean,totalTime-riseMean,nadChangeMean,lty=2)
  #riseTime
  polygon(totalTime-c(riseTimeDen$x,rev(riseTimeDen$x)),10^(par('usr')[3]-diff(par('usr')[3:4])*.03+ yWidth*c(riseTimeDen$y,-rev(riseTimeDen$y))/scaleFunc(riseTimeDen$y)),xpd=NA,col=bothCol)
  segments(totalTime-riseMean,10^(par('usr')[3]-diff(par('usr')[3:4])*.03- yWidth),totalTime-riseMean,10^(par('usr')[3]-diff(par('usr')[3:4])*.03+ yWidth),xpd=NA)
  #rise
  polygon(totalTime+c(riseDen$y,-rev(riseDen$y))/scaleFunc(riseDen$y)*xWidth,c(riseDen$x,rev(riseDen$x)),col=normCol)
  segments(totalTime-xWidth,riseChangeMean,totalTime+xWidth,riseChangeMean)
  segments(totalTime-riseMean,nadChangeMean,totalTime,riseChangeMean)
  #riseFast
  polygon(totalTime+c(riseFastDen$y,-rev(riseFastDen$y))/scaleFunc(riseFastDen$y)*xWidth,c(riseFastDen$x,rev(riseFastDen$x)),xpd=NA,col=fastCol)
  segments(totalTime-xWidth,riseChangeFastMean,totalTime+xWidth,riseChangeFastMean)
  segments(totalTime-riseMean,nadChangeFastMean,totalTime,riseChangeFastMean)
}
pdf('out/bayesParams.pdf',width=3.5,height=3.5)
  par(mar=c(4.6,4,1,.5))
  plotParams(fit$fit)
  plotParams(fitB$fit,timeAfterInf=2,totalTime=5,ylab='IFNb IC50')
dev.off()

"
  transformed parameters{
    vector[nPatient] acute;
    real expectedIC50[nVirus];
    vector[nPatient] nadirTime;
    vector[nPatient] nadirChange;
    vector[nArt] riseChange;
    vector[nArt] riseTime;
    acute=acuteMean+acuteRaw*acuteSD;
    nadirChange=fastChangeMean*isFast+nadirChangeMean+nadirChangeRaw*nadirChangeSD;
    nadirTime=nadirTimeMean+nadirTimeRaw*nadirTimeSD;
    //TRANSSUB
    //riseTime=riseTimeMean+riseTimeRaw*riseTimeSD;
    riseChange=riseChangeMean+riseChangeRaw*riseChangeSD;
    for(ii in 1:nVirus){
      expectedIC50[ii]=acute[patients[ii]];
      if(days[ii]<exp(nadirTime[patients[ii]]))expectedIC50[ii]=expectedIC50[ii]+nadirChange[patients[ii]]*days[ii]/exp(nadirTime[patients[ii]]);
      else{
        expectedIC50[ii]=expectedIC50[ii]+nadirChange[patients[ii]];
      }
      if(hasArt[ii]){
        if(daysBeforeArt[ii]<0){
          expectedIC50[ii]=expectedIC50[ii]+riseChange[artIds[ii]];
        }else{
          if(daysBeforeArt[ii]<exp(riseTime[artIds[ii]]))expectedIC50[ii]=expectedIC50[ii]+riseChange[artIds[ii]]*(1-daysBeforeArt[ii]/exp(riseTime[artIds[ii]]));
        }
      }
    }
  }
  model {
    //MODELSUB
    //riseTimeRaw~normal(0,1);
    ic50~normal(expectedIC50,sigma);
    sigma~gamma(1,.1);
    nadirTimeSD~gamma(1,.1);
    nadirTimeRaw~normal(0,1);
    riseTimeSD~gamma(1,.1);
    acuteSD~gamma(1,.1);
    acuteRaw~normal(0,1);
    nadirChangeSD~gamma(1,.1);
    nadirChangeRaw~normal(0,1);
    riseChangeSD~gamma(1,.1);
    riseChangeRaw~normal(0,1);
    riseChangeMean~normal(0,10);
    nadirChangeMean~normal(0,10);
    nadirTimeMean~normal(0,10);
    riseTimeMean~normal(0,10);
    fastChangeMean~normal(0,10);
  }
"

mat<-as.matrix(fit$fit)
cbind(exp(apply(mat[,grepl('riseTime\\[',colnames(mat))],2,mean)),c(artDfosx,'WEAU'=391))
exp(apply(mat[,grepl('nadirTime\\[',colnames(mat))],2,mean))


mat<-as.matrix(fit$fit)
matB<-as.matrix(fitB$fit)
message('Mean IFNa2 acute:',format(exp(mean(mat[,'acuteMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(mat[,'acuteMean'],c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNb acute:',format(exp(mean(matB[,'acuteMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(matB[,'acuteMean'],c(.025,.975))),digits=4),collapse=' - '))
randomAcute<-rnorm(nrow(mat)*3,mat[,'acuteMean'],sqrt(mat[,'acuteSD']^2+mat[,'sigma']^2))
randomAcuteB<-rnorm(nrow(matB)*3,mat[,'acuteMean'],sqrt(matB[,'acuteSD']^2+matB[,'sigma']^2))
message('95% prediction interval IFNa2 acute: ',paste(format(exp(quantile(randomAcute,c(.025,.975))),digits=4),collapse=' - '))
message('99% prediction interval IFNa2 acute: ',paste(format(exp(quantile(randomAcute,c(.005,.995))),digits=4),collapse=' - '))
message('95% prediction interval IFNb acute: ',paste(format(exp(quantile(randomAcuteB,c(.025,.975))),digits=4),collapse=' - '))
message('99% prediction interval IFNb acute: ',paste(format(exp(quantile(randomAcuteB,c(.005,.995))),digits=4),collapse=' - '))

message('Mean IFNa2 nadir drop:',format(exp(-mean(mat[,'nadirChangeMean'])),digits=4),' 95% credible interval: ',paste(format(exp(-quantile(mat[,'nadirChangeMean'],c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNb nadir drop:',format(exp(-mean(matB[,'nadirChangeMean'])),digits=4),' 95% credible interval: ',paste(format(exp(-quantile(matB[,'nadirChangeMean'],c(.025,.975))),digits=4),collapse=' - '))

message('Mean IFNa2 nadir time:',format(exp(mean(mat[,'nadirTimeMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(mat[,'nadirTimeMean'],c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNb nadir time:',format(exp(mean(matB[,'nadirTimeMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(matB[,'nadirTimeMean'],c(.025,.975))),digits=4),collapse=' - '))


fastAN<-mat[,'nadirChangeMean']+mat[,'fastChangeMean']
fastBN<-matB[,'nadirChangeMean']+matB[,'fastChangeMean']
fastA<-mat[,'fastChangeMean']
fastB<-matB[,'fastChangeMean']
message('Mean IFNa2 nadir level fast:',format(exp(-mean(fastAN)),digits=4),' 95% credible interval: ',paste(format(exp(-quantile(fastAN,c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNb nadir level fast:',format(exp(-mean(fastBN)),digits=4),' 95% credible interval: ',paste(format(exp(-quantile(fastBN,c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNa2 fast increase:',format(exp(mean(fastA)),digits=4),' 95% credible interval: ',paste(format(exp(quantile(fastA,c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNb fast increase:',format(exp(mean(fastB)),digits=4),' 95% credible interval: ',paste(format(exp(quantile(fastB,c(.025,.975))),digits=4),collapse=' - '))

message('Mean IFNa2 rise time:',format(exp(mean(mat[,'riseTimeMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(mat[,'riseTimeMean'],c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNb rise time:',format(exp(mean(matB[,'riseTimeMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(matB[,'riseTimeMean'],c(.025,.975))),digits=4),collapse=' - '))


message('Mean IFNa2 rise mean:',format(exp(mean(mat[,'riseChangeMean'])),deigits=4),' 95% credible interval: ',paste(format(exp(quantile(mat[,'riseChangeMean'],c(.025,.975))),digits=4),collapse=' - '))
message('Mean IFNb rise mean:',format(exp(mean(matB[,'riseChangeMean'])),digits=4),' 95% credible interval: ',paste(format(exp(quantile(matB[,'riseChangeMean'],c(.025,.975))),digits=4),collapse=' - '))
