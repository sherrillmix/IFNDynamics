stanCode<-'
  data {
    int<lower=0> N;
    int<lower=0> nPat;
    int<lower=0> patientId[N];
    int<lower=0> day[N];
//    real ic50[N];
 //   real ic50Beta[N];
    int<lower=0> nObs[nPat];
    int<lower=0> starts[nPat];
    int<lower=0> startDays[nPat];
    //int<lower=0> lookupIds[N];
    int<lower=0> nDays[nPat];
    int<lower=0> nTotals[nPat];
    int<lower=0> totalStarts[nPat];
    vector<lower=0>[sum(nDays)] days;
    //int<lower=0> dayPatientId[sum(nDays)];
    int<lower=0> dayTotalIndex[sum(nDays)];
    vector<lower=0>[sum(nDays)] vl;
    vector<lower=0>[sum(nDays)] cd4;
  }
  parameters {
    vector[sum(nTotals)] trueVl;
    vector<lower=0>[sum(nTotals)] trueCd4;
    //vector[sum(nDays)]<lower=0> ic50;
    //vector[sum(nDays)]<lower=0> ic50Beta;
    real<lower=0> vlAutoCor;
    real<lower=0> vlError;
    real<lower=0> cd4AutoCor;
    real<lower=0> cd4Error;
    //real<lower=0> ic50AutoCor;
    //real<lower=0> ic50BetaAutoCor;
    real metaCd4;
    real metaCd4Sd;
    vector[nPat] rawCd4s;
  }
  transformed parameters{
    vector[nPat] betaCd4s;
    betaCd4s = metaCd4 + rawCd4s*metaCd4Sd;
  }
  model {
    rawCd4s~normal(0,1);
    vlAutoCor~gamma(1,10);
    //vlError~gamma(1,10);
    vlError~normal(.05,.05);
    cd4AutoCor~gamma(1,.5);
    //cd4Error~gamma(1,.5);
    cd4Error~normal(5,5);
    metaCd4Sd~gamma(1,1);
    trueVl[totalStarts]~normal(5,10);
    trueCd4[totalStarts]~normal(750,1500);
    for(ii in 1:nPat){
      segment(trueVl,totalStarts[ii]+1,nTotals[ii]-1)~normal(segment(trueVl,totalStarts[ii],nTotals[ii]-1),vlAutoCor);
      segment(trueCd4,totalStarts[ii]+1,nTotals[ii]-1)~normal(segment(trueCd4,totalStarts[ii],nTotals[ii]-1),cd4AutoCor);
      //segment(trueCd4,totalStarts[ii]+1,nTotals[ii]-1)~normal(betaCd4s[ii]*segment(days,totalStarts[ii]+1,nTotals[ii]-1)+segment(trueCd4,totalStarts[ii],nTotals[ii]-1),cd4AutoCor);
      //should linear interpolate
      segment(vl,startDays[ii],nDays[ii])~normal(trueVl[segment(dayTotalIndex,startDays[ii],nDays[ii])],vlError);
      segment(cd4,startDays[ii],nDays[ii])~normal(betaCd4s[ii]*segment(days,startDays[ii],nDays[ii])+trueCd4[segment(dayTotalIndex,startDays[ii],nDays[ii])],cd4Error);
    }
  }
'

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


ic50CodeSample<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nArt;
    real ic50[nVirus];
    int<lower=0,upper=nPatient> patients[nVirus];
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
    real sampleNoiseRaw[nSample];
    real<lower=0> sampleSD;
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
    //TRANSSUB
    //riseTime=riseTimeMean+riseTimeRaw*riseTimeSD;
    riseChange=riseChangeMean+riseChangeRaw*riseChangeSD;
    for(ii in 1:nVirus){
      expectedIC50[ii]=acute[patients[ii]]+sampleNoiseRaw[sample[ii]]*sampleSD;
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
    sampleNoiseRaw~normal(0,1);
    sampleSD~gamma(1,.1);
  }
'

ic50Code<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nArt;
    real ic50[nVirus];
    int<lower=0,upper=nPatient> patients[nVirus];
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
  }
'


ic50CodeRiseAfter<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nPatient;
    int<lower=0> nArt;
    real ic50[nVirus];
    int<lower=0,upper=nPatient> patients[nVirus];
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
        if(daysBeforeArt[ii]<exp(riseTime[artIds[ii]]))expectedIC50[ii]=expectedIC50[ii]+riseChange[artIds[ii]]*(1-daysBeforeArt[ii]/exp(riseTime[artIds[ii]]));
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
  }
'


