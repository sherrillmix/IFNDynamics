library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat<-read.csv('out/allLongitudinal.csv',stringsAsFactors=FALSE)
meta<-read.csv('out/allLongitudinalMeta.csv')
founders<-read.csv('out/founders.csv',stringsAsFactors=FALSE,row.names=1)
superTimes<-structure(founders$superTime,.Names=rownames(founders))

ic50_bp_mar<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nArt;
    real cd4[nVirus];
    real vl[nVirus];
    real ic50[nVirus];
    int<lower=0,upper=nPatient> patients[nVirus];
    vector<lower=0,upper=1>[nPatient] isFast;
    vector<lower=0,upper=1>[nPatient] isNon;
    int<lower=0> artIds[nVirus];
    int<lower=0> days[nVirus];
    real weeks[nVirus];
    //real daysBeforeArt[nVirus];
    real artStart[nArt];
    real<lower=0> hasArt[nVirus];
    int<lower=0> nSample;
    int<lower=0,upper=nSample> sample[nVirus];
    int<lower=0> tMax;
    matrix[nPatient,tMax] vlMat;
    matrix[nPatient,tMax] cd4Mat;
  }
  parameters {
    vector[nPatient] acuteRaw;
    real acuteMean;
    real<lower=0> acuteSD;
    //vector[nPatient] nadirTimeRaw;
    real<lower=0,upper=tMax> nadirTimeMean;
    real<lower=0> nadirTimeSD;
    vector[nPatient] nadirChangeRaw;
    real nadirChangeMean;
    real<lower=0> nadirChangeSD;
    real<lower=0> sigma;
    real fastChangeMean;
    real nonChange;
    vector[nPatient] cd4BetaRaw;
    real cd4BetaMean[2];
    real<lower=0> cd4BetaSD;
    //vector[nPatient] vlBetaRaw;
    //real vlBetaMean[2];
    //real<lower=0> vlBetaSD;
  }
  transformed parameters{
    vector[nPatient] acute;
    //real expectedIC50[nVirus];
    //vector<lower=0,upper=tMax>[nPatient] nadirTime;
    vector[nPatient] nadirChange;
    vector[nPatient] cd4Beta;
    //vector[nPatient] vlBeta;
    vector[tMax] lp;
    acute=acuteMean*(1-isNon)+acuteRaw*acuteSD+isNon*nonChange;
    nadirChange=nadirChangeMean*(1-isFast)+nadirChangeRaw*nadirChangeSD+fastChangeMean*isFast;
    cd4Beta=(cd4BetaMean[1]*(1-isFast)+cd4BetaMean[2]*isFast)+cd4BetaRaw*cd4BetaSD;
    //vlBeta=(vlBetaMean[1]*(1-isFast)+vlBetaMean[2]*isFast)+vlBetaRaw*vlBetaSD;
    //nadirChange=nadirChangeMean+nadirChangeRaw*nadirChangeSD;
    //nadirTime=exp(nadirTimeMean+nadirTimeRaw*nadirTimeSD);
    for (s in 1:tMax){
      lp[s]=neg_binomial_2_lpmf(s|nadirTimeMean,nadirTimeMean^2/exp(nadirTimeSD));
      for (ii in 1:nVirus){
        lp[s] = lp[s] + normal_lpdf(ic50[ii] | weeks[ii] < s ? 
          acute[patients[ii]]+nadirChange[patients[ii]]*weeks[ii]/s  : 
          acute[patients[ii]]+nadirChange[patients[ii]] +(cd4Beta[patients[ii]])*(cd4[ii]-cd4Mat[patients[ii],s])
          ,sigma);
      }
    }
  }
  model {
    target += log_sum_exp(lp);
    //ic50~normal(expectedIC50,sigma);
    sigma~gamma(1,.1);
    nadirTimeSD~normal(0,10); //gamma(1,.1);
    //nadirTimeRaw~normal(0,1);
    nadirChangeSD~gamma(1,.1);
    nadirChangeRaw~normal(0,1);
    acuteSD~gamma(1,.1);
    acuteRaw~normal(0,1);
    nadirChangeMean~normal(0,10);
    fastChangeMean~normal(0,10);
    nonChange~normal(0,10);
    cd4BetaMean~normal(0,10);
    cd4BetaSD~gamma(1,.1);
    cd4BetaRaw~normal(0,1);
    //vlBetaMean~normal(0,10);
    //vlBetaSD~gamma(1,.1);
    //vlBetaRaw~normal(0,1);
  }\n
'
ic50Mod <- stan_model(model_code = ic50_bp_mar)

bayesIC50_3<-function(mod,ic50,time,timePreArt,patient,chains=50,nIter=30000,fastProgressors=c(),cd4=c(),vl=c(),baseVl=c(),meta=NULL,superTimes=NULL,nonProgressors=FALSE,...){
  #treating all postART   as one group (currently only WEAU)
  patientId<-structure(1:length(unique(patient)),.Names=sort(unique(patient)))
  artId<-structure(1:length(unique(patient[!is.na(timePreArt)])),.Names=sort(unique(patient[!is.na(timePreArt)])))
  artStart<-tapply(time+timePreArt,patient,unique)
  sample<-paste(patient,time,sep='_')
  sampleId<-structure(1:length(unique(sample)),.Names=unique(sample[order(patient,time)]))
  ####patient-week matrix of inferred CD4 and VL
  weekMax<-tapply(round(time/7),patient,max)
  newDat<-do.call(rbind,mapply(function(xx,yy){data.frame('pat'=yy,'week'=min(round(time/7)):xx,stringsAsFactors=FALSE)},weekMax,names(weekMax),SIMPLIFY=FALSE))
  newDat$day<-newDat$week*7-3.5
  newDat$vl<-newDat$cd4<-NA
  for(ii in unique(newDat$pat)){
    thisDat<-meta[meta$mm==ii&meta$time<=weekMax[ii]*7,]
    newDat[newDat$pat==ii,'vl']<-approx(thisDat$time,thisDat$vl,newDat[newDat$pat==ii,'day'],rule=2)$y
    newDat[newDat$pat==ii,'cd4']<-approx(thisDat$time,thisDat$cd4,newDat[newDat$pat==ii,'day'],rule=2)$y
  }
  cd4Mat<-tapply(newDat$cd4,list(newDat$pat,newDat$week),c)
  vlMat<-tapply(newDat$vl,list(newDat$pat,newDat$week),c)
  weekOffset<-1-min(round(time/7))
  tMax<-100 #MM15 ends at 103
  superIds<-structure(1:(sum(!is.na(superTimes))+1),.Names=c('NA',names(superTimes)[!is.na(superTimes)]))
  superId<-ifelse(time<superTimes[patient],NA,superIds[patient])
  superId[is.na(superId)]<-1
  newDat$isSuper<-newDat$day>=ifelse(is.na(superTimes),Inf,superTimes)[newDat$pat]
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
    isFast=names(patientId) %in% fastProgressors,
    isNon=names(patientId) %in% nonProgressors,
    cd4=as.numeric(cd4),
    vl=as.numeric(vl),
    baseVl=baseVl[names(patientId)],
    weeks=round(time/7),
    vlMat=log(vlMat[names(patientId),weekOffset+1:tMax]),
    cd4Mat=cd4Mat[names(patientId),weekOffset+1:tMax]/100,
    tMax=tMax,
    superId=superId,
    isSuper=superId>1,
    nSuper=sum(!is.na(superTimes))
  )
  fit <- sampling(mod, data = dat, iter=nIter, chains=chains,thin=2,control=list(adapt_delta=.98,max_treedepth=12),...)
  return(list('fit'=fit,pats=patientId,arts=artId,artStart=artStart,sample=sampleId,dat=dat,superId=superIds,simDat=newDat))
}
baseVl<-dnar::withAs(xx=unique(meta[meta$time<365*2&meta$time>180,c('time','vl','mm')]),tapply(xx$vl,xx$mm,function(xx)mean(log(xx),na.rm=TRUE)))
dat$week<-as.integer(round(dat$time/7))

fit<-dnar::withAs(xx=dat[!is.na(dat$ic50)&!dat$qvoa,],bayesIC50_3(ic50Mod,xx$ic50,xx$time,xx$timeBeforeArt,xx$pat,fastProgressors=c('MM15','WEAU'),nonProgressors=c('MM55','MM62'),cd4=(xx$fillCD4)/100,vl=log(xx$fillVl),baseVl=baseVl,nIter=1000,chains=50,meta=meta,superTimes=superTimes<50))

fitB<-dnar::withAs(xx=dat[!is.na(dat$beta)&!dat$qvoa,],bayesIC50_3(ic50Mod,xx$beta,xx$time,xx$timeBeforeArt,xx$pat,fastProgressors=c('MM15','WEAU'),nonProgressors=c('MM55','MM62'),cd4=(xx$fillCD4)/100,vl=log(xx$fillVl),baseVl=baseVl,nIter=1000,chains=50,meta=meta,superTimes=superTimes))
#save(fit,file='out/bayesFit_20200602.Rdat') #CD4 marginalized
#save(fit,file='out/bayesFit_20200602.Rdat') #CD4 marginalized
#save(fit,file='out/bayesFit_20200602c.Rdat')# CD4+VL marginalized 2000 rep, 125 weeks
#save(fit,file='out/bayesFit_20200604.Rdat')# CD4+VL marginalized 2000 rep, 100 weeks, setpoint
#save(fit,file='out/bayesFit_20200605.Rdat')# CD4+VL (split out fast), marginalized 1000 rep, super, 100 weeks, setpoint
#load(file='out/bayesFit_20200604.Rdat')
#
print(fit$fit,pars=c('nadirTimeRaw','nadirChangeRaw','expectedIC50','acuteRaw','lp','vlBetaRaw','cd4BetaRaw','superBetaRaw'),include=FALSE)
print(fitB$fit,pars=c('nadirTimeRaw','nadirChangeRaw','expectedIC50','acuteRaw','lp','vlBetaRaw','cd4BetaRaw','superBetaRaw'),include=FALSE)

softMax<-function(xx)exp(xx)/sum(exp(xx))
logsumexp <- function (x) {
    y = max(x)
  y + log(sum(exp(x - y)))
}
softmax <- function (x) {
    exp(x - logsumexp(x))
}

mat<-as.matrix(fit$fit)
matB<-as.matrix(fitB$fit)
calcPreds<-function(mat,dat,newDat){
  newDat$scaleCd4<-newDat$cd4/100
  newDat$scaleVl<-log(newDat$vl)
  weekMax<-tapply(round(dat$days/7),names(dat$patients),max)
  vlCd4<-unique(data.frame('pat'=names(dat$patients),'day'=dat$days,'vl'=dat$vl,'cd4'=dat$cd4))
  newDat$patId<-sapply(newDat$pat,function(xx)dat$patients[names(dat$patients)==xx][1])
  weeks<-outer(newDat$pat,1:dat$tMax,function(xx,yy)yy)
  isNadir<-newDat$week<weeks
  propNadir<-newDat$week/weeks
  vlMat<-dat$vlMat[newDat$pat,]
  cd4Mat<-dat$cd4Mat[newDat$pat,]
  preds<-do.call(cbind,parallel::mclapply(split(mat,sort(rep(1:30,length.out=nrow(mat)))),function(xx,...){
    thisMat<-matrix(xx,ncol=ncol(mat))
    colnames(thisMat)<-colnames(mat)
    apply(thisMat,1,function(xx){
      pS<-softmax(xx[grep('lp\\[',names(xx))])    
      #cd4Bs<-xx[sprintf('cd4Beta[%d]',dat$patients)]
      #vlBs<-xx[sprintf('vlBeta[%d]',dat$patients)]
      #acutes<-xx[sprintf('acute[%d]',dat$patients)]
      #nadirs<-xx[sprintf('nadirChange[%d]',dat$patients)]
      cd4Bs<-xx[sprintf('cd4Beta[%d]',newDat$patId)]
      acutes<-xx[sprintf('acute[%d]',newDat$patId)]
      nadirs<-xx[sprintf('nadirChange[%d]',newDat$patId)]
      pred<-acutes+ifelse(isNadir,propNadir*nadirs,nadirs+(newDat$scaleCd4-cd4Mat)*cd4Bs)
      out<-apply(pred,1,function(xx)sum(xx*pS))
    })
  },mc.cores=40))
  sigmas<-mat[,'sigma']
  isos<-apply(rbind(mat[,'sigma'],preds),2,function(xx)rnorm(length(xx[-1]),xx[-1],xx[1]))
  list('preds'=preds,'isos'=isos,'simData'=newDat)
}
predIc50<-calcPreds(mat,fit$dat,fit$simDat)
predIc50B<-calcPreds(matB,fitB$dat,fitB$simDat)

plotIfn<-function(fit,predIc50,ylab='IFNa2 IC50'){
  for(ii in sort(unique(names(fit$dat$patients)))){
    plot(1,1,type='n',xlim=range(fit$dat$days),ylim=range(exp(fit$dat$ic50))*c(.8,1),main=ii,xlab='DFOSx',ylab='',log='y',yaxt='n',mgp=c(1.8,.5,0),tcl=-.3)
    mtext(ylab,2,line=2.9,cex=.83)
    selector<-names(fit$dat$patients)==ii
    points(fit$dat$days[selector],exp(fit$dat$ic50[selector]))
    dnar::logAxis(las=1)
    selector2<-predIc50$simData$pat==ii
    lines(predIc50$simData[selector2,'day'],exp(apply(predIc50$preds[selector2,],1,mean)))
    quants<-(apply(predIc50$preds[selector2,],1,quantile,c(.025,.975),na.rm=TRUE))
    quants2<-(apply(predIc50$isos[selector2,],1,quantile,c(.025,.975),na.rm=TRUE))
    polygon(c(predIc50$simData[selector2,'day'],rev(predIc50$simData[selector2,'day'])),exp(c(quants[1,],rev(quants[2,]))),col='#00000033',border=NA)
    polygon(c(predIc50$simData[selector2,'day'],rev(predIc50$simData[selector2,'day'])),exp(c(quants2[1,],rev(quants2[2,]))),col='#00000033',border=NA)
    abline(v=superTimes[ii],col='red')
  }
}
pdf('Rplots.pdf',width=17,height=7)
  par(mfrow=c(2,5),mar=c(3.5,4,1,.1))
  plotIfn(fit,predIc50)
  plotIfn(fitB,predIc50B,'IFNb IC50')
  for(ii in sort(unique(names(fit$dat$patients)))){
    selector2<-predIc50$simData$pat==ii
    plot(predIc50$simData[selector2,'day'],(predIc50$simData[selector2,'cd4']/100),xlim=range(fit$dat$days),ylim=range((fit$dat$cd4)),main=ii,xlab='DFOSx',ylab='CD4',mgp=c(1.8,.5,0),tcl=-.3)
    abline(v=superTimes[ii],col='red')
  }
  for(ii in sort(unique(names(fit$dat$patients)))){
    selector2<-predIc50$simData$pat==ii
    plot(predIc50$simData[selector2,'day'],log(predIc50$simData[selector2,'vl']),xlim=range(fit$dat$days),ylim=range((fit$dat$vl)),main=ii,xlab='DFOSx',ylab='VL',mgp=c(1.8,.5,0),tcl=-.3)
    abline(v=superTimes[ii],col='red')
  }
dev.off()

apply(mat[,grep('lp\\[',colnames(mat))],2,mean)
pdf('Rplots.pdf',height=20,width=20);plot(apply(mat[,grep('lp\\[',colnames(mat))],2,mean));dev.off()
pdf('Rplots.pdf',height=20,width=20);traceplot(fit$fit,'lp');dev.off()
