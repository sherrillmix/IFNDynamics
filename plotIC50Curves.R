library(xlsx)
library(dnar)

concAlpha<-c(0, 0.001, 0.003, 0.008, 0.023, 0.068, 0.204, 0.611, 1.833, 5.5)
concBeta<-c(0, 4.4E-05, 0.00044, 0.0044, 0.044, 0.44, 4.4, 44, 440, 4400)
lowerLimit<-50

#https://www.myassays.com/four-parameter-logistic-regression.html
#http://stats.stackexchange.com/questions/61144/how-to-do-4-parametric-regression-for-elisa-data-in-r
#a = the minimum value that can be obtained (i.e. what happens at 0 dose)
#b = Hill’s slope of the curve (i.e. this is related to the steepness of the curve at point c).
#c = the point of inflection (i.e. the point on the S shaped curve halfway between a and d)
#d = the maximum value that can be obtained (i.e. what happens at infinite dose)
#f<-function(B,x,transformX=function(x)x){if(any(x<=0))browser();exp(log(B[1]-B[4])-log(1+(transformX(x)/B[3])^B[2]))+B[4]}
f5<-function(B,x,transformX=function(x)x){
  #if(x==-Inf)return(B[1])
  exp(log(B[1]-B[4])-B[5]*log(1+(transformX(x)/B[3])^B[2]))+B[4]
}
#f5<-function(B,x,transformX=log){exp(log(B[1]-B[4])-B[5]*log(1+(transformX(x)/B[3])^B[2]))+B[4]}
f4<-function(B,x,transformX=function(x)x){exp(log(B[1]-B[4])-log(1+(transformX(x)/B[3])^B[2]))+B[4]}
LS<-function(B,y,x){sum((log(y)-log(fitFunc(B,x)))^2)}
optIc50<-function(B)optim(1,function(x,B){(B[1]*.5-fitFunc(B,x))^2},B=B,method='Brent',lower=10^-10,upper=10^2)$par
fit4par<-function(concs,p24s){
  #p24s<-p24s[concs>0]
  #concs<-concs[concs>0]
  fits<-lapply(list(c(max(p24s),1,1,min(p24s)),c(max(p24s)*2,1,1,min(p24s)/10),c(mean(p24s),1,1,mean(p24s)),c(max(p24s),.01,1,1)),function(starts)suppressWarnings(nlminb(starts,LS,x=concs,y=p24s,lower=c(0,-Inf,-Inf,0))))
  #fit<-suppressWarnings(nlminb(c(max(p24s),1,1,min(p24s)),LS,x=concs,y=p24s,lower=c(0,-Inf,-Inf,0))$par)
  obj<-sapply(fits,'[[','objective')
  return(fits[[which.min(obj)]]$par)
}
fit5par<-function(concs,p24s){
  #p24s<-p24s[concs>0]
  #concs<-concs[concs>0]
  fits<-lapply(list(c(max(p24s),1,1,min(p24s),1),c(max(p24s)*2,1,1,min(p24s)/10),c(mean(p24s),1,1,mean(p24s),1),c(max(p24s),.01,1,1),1),function(starts)suppressWarnings(nlminb(starts,LS5,x=concs,y=p24s,lower=c(0,-Inf,-Inf,0,-Inf))))
  #fit<-suppressWarnings(nlminb(c(max(p24s),1,1,min(p24s)),LS,x=concs,y=p24s,lower=c(0,-Inf,-Inf,0))$par)
  obj<-sapply(fits,'[[','objective')
  return(fits[[which.min(obj)]]$par)
}
fitFunc<-f5
fitter<-fit5par

findVresIc50<-function(concAlpha,p24s){
  origConcs<-rep(concAlpha,each=2*nrow(p24s))
  p24s<-unlist(p24s)
  fit<-fitter(origConcs,p24s)
  vres<-fitFunc(fit,max(origConcs))
  ic50<-optIc50(fit)
  percVres<-vres/fit[1]*100
  return(c('ic50'=ic50,'vres'=vres,'percVres'=percVres,'max'=fit[1]))
}
calcIc50s<-function(dilutes,concAlpha){
  out<-do.call(rbind,lapply(unique(dilutes$sample),function(thisSample){
    xx<-dilutes[dilutes$sample==thisSample,1:20]
    findVresIc50(concAlpha,xx)
  }))
  rownames(out)<-unique(dilutes$sample)
  return(out)
}

if(!exists('dilutes')){
  dilutions<-read.csv('VOA MM-cohort/dilutions.csv',stringsAsFactors=FALSE)
  readIfns<-function(file){
    wb <- loadWorkbook(file)
    counter<-0
    ic50Curves<-lapply(getSheets(wb),function(sheet){
      counter<<-counter+1
      if(counter==1)return(NULL)
      rows<-getCells(getRows(sheet)[1:40],simplify=FALSE)
      vals<-lapply(rows,function(row){
        tmp<-sapply(row,function(xx)ifelse(is.null(xx),NA,getCellValue(xx)))
        #100 is arbitrary number to make all same width
        out<-rep(NA,100)
        names(out)<-1:100
        out[names(tmp)[names(tmp)!='']]<-tmp[names(tmp)!='']
        return(out)
      })
      lastRow<-2+max(which(!is.na(sapply(vals[3:12],'[[',4))))
      dat<-as.data.frame(apply(do.call(rbind,lapply(vals[3:lastRow],function(zz)zz[3:22])),2,function(xx){ifelse(xx==''|is.na(xx),NA,as.numeric(sub('[><]','',xx)))}),stringsAsFactors=FALSE)
      dat$sample<-fillDown(sapply(vals[3:lastRow],'[[',2))
      return(dat)
    })
    names(ic50Curves)<-names(getSheets(wb))
    dilutes<-do.call(rbind,ic50Curves)
    return(dilutes)
  }
  dilutes<-readIfns("VOA MM-cohort/IC50 Alpha for New patients BULK isol. .xls")
  dilutes[,1:20][dilutes[,1:20]<lowerLimit]<-lowerLimit
  dilutes$dilution<-sapply(dilutes$sample,function(xx)if(xx %in% dilutions$sample[dilutions$ifn=='alpha']) dilutions[dilutions$ifn=='alpha'&dilutions$sample==xx,'dilution'] else 50)
  dilutes[,1:20]<-dilutes[,1:20]*dilutes$dilution/1000
  dilutesBeta<-readIfns("VOA MM-cohort/IC50 Beta for all BULK isolates .xlsx")
  dilutesBeta[,1:20][dilutesBeta[,1:20]<lowerLimit]<-lowerLimit
  dilutesBeta$dilution<-sapply(dilutesBeta$sample,function(xx)if(xx %in% dilutions$sample[dilutions$ifn=='beta']) dilutions[dilutions$ifn=='beta'&dilutions$sample==xx,'dilution'] else 50)
  dilutesBeta[,1:20]<-dilutesBeta[,1:20]*dilutesBeta$dilution/1000
}

plotIfn<-function(concAlpha,p24s,main='',xlab='',ylims=range(p24s),log='xy',scaleMax=FALSE){
  origConcs<-concs<-rep(concAlpha,each=2*nrow(p24s))
  if(scaleMax)maxs<-unlist(p24s[,1:2])
  else maxs<-1
  vresIc50<-findVresIc50(concAlpha,p24s)
  p24s<-unlist(p24s)
  zeroOffset<-.01
  concs[concs==0]<-min(concs[concs>0])*zeroOffset
  fit<-fitter(origConcs,p24s)
  cols<-rainbow.lab(12)[c(1:2,11:12)]
  names(cols)<-c('Bio replicate 1\nTech replicate 1','Bio replicate 1\nTech replicate 2','Bio replicate 2\nTech replicate 1','Bio replicate 2\nTech replicate 2')
  plot(concs,p24s/maxs,xlab=xlab,ylab='',log=log,las=1,xaxt='n',main=main,bg=cols,pch=21,mgp=c(2,1,0),ylim=ylims/mean(maxs),yaxt='n')
  #fit2<-suppressWarnings(nlminb(c(max(p24s),1,1,min(p24s)),LS,x=concs[origConcs!=0],y=p24s[origConcs!=0],lower=c(0,-Inf,-Inf,0))$par)
  fakeConc<-10^seq(-10,10,.001)
  fitLine<-fitFunc(fit,fakeConc)
  lines(fakeConc,fitLine/mean(maxs),col='#00000066',lwd=3)
  if(scaleMax)title(ylab='Proportion maximum p24',mgp=c(2,1,0))
  else title(ylab='p24 concentration (ng/ml)',mgp=c(2.25,1,0))
  if(grepl('x',log))logAxis(1,axisMin=min(origConcs[origConcs>0]),mgp=c(2.5,.7,0))
  else axis(1,pretty(par('usr')[1:2]),mgp=c(2.5,.7,0))
  if(grepl('y',log))logAxis(2,las=1,mgp=c(2.5,.7,0))
  else axis(2,pretty(par('usr')[3:4]),las=1,mgp=c(2.5,.7,0))
  axis(1,min(origConcs[origConcs>0])*zeroOffset,0,mgp=c(2.5,.7,0))
  #abline(h=c(fit[c(1,4)]),lty=3,col='#00000055')
  #abline(v=fit[3],lty=3,col='#00000055')
  isAbove<-vresIc50['ic50']>10^(par('usr')[3]*.75+par('usr')[4]*.25)
  abline(v=vresIc50['ic50'],lty=2,col='#FF000055')
  #if(fit[4]<par('usr')[4])
  vres<-fitFunc(fit,max(origConcs))
  abline(h=vresIc50['vres']/mean(maxs),lty=2,col='#0000FF55')
  yCoord<-vresIc50['vres']/2^ifelse(isAbove,1,-1)
  xCoord<-10^(par('usr')[1]*.72+par('usr')[2]*.28)
  text(xCoord,yCoord,sprintf('Vres=%s',format(vresIc50['percVres'],digits=2,width=3)),adj=c(1,0.5),col='blue')
  segments(xCoord*1.1,yCoord,10^(par('usr')[1]*.65+par('usr')[2]*.35),vresIc50['vres']/mean(maxs),col='blue')
  #ic50 label
  isAbove<-vresIc50['vres']>10^mean(par('usr')[3:4])
  yCoord<-ifelse(isAbove,10^(par('usr')[3]*.2+par('usr')[4]*.8),10^(par('usr')[3]*.55+par('usr')[4]*.45)*1.2)
  text(vresIc50['ic50']/1.85,yCoord,sprintf('IC50=%s',format(vresIc50['ic50'],digits=2)),adj=c(1,0.5),col='red')
  segments(vresIc50['ic50']/1.75,yCoord,vresIc50['ic50'],yCoord*ifelse(isAbove,1.2,.8),col='red')
  axis(4,vresIc50['max']/c(1,2,10,100,1000),as.character(c(1,.5,.1,.01,.001)*100),las=1,tcl=-.2,mgp=c(0,.6,0))
  abline(h=vresIc50['max']/2,lty=3,col='#00000055')
  text(convertLineToUser(2.6,4),10^mean(par('usr')[3:4]),'Percent of maximum',xpd=NA,srt=-90)
  par('lheight'=.72)
  legend('topright',names(cols),pch=21,pt.bg=cols,inset=.01,cex=.5,pt.cex=1,y.intersp=1.5)
  #if(thisSample=='MM55.12.2B1 bulk')browser()
  return(fit)
}
plotIfns<-function(dilutes,concAlpha,xlab='',...){
  ylims<-range(dilutes[,1:20],na.rm=TRUE)
  par(mar=c(3,3.25,.25,3))
  for(thisSample in unique(dilutes$sample)){
    xx<-dilutes[dilutes$sample==thisSample,1:20]
    plotIfn(concAlpha,xx,'',xlab=xlab,ylims=ylims,...)
  }
}
plotDualIfns<-function(dilutes,dilutesBeta,concAlpha,concBeta){
  ylims<-range(dilutes[,1:20],na.rm=TRUE)
  ylims2<-range(dilutesBeta[,1:20],na.rm=TRUE)
  par(mar=c(3.5,4.5,1.5,1),mfrow=c(1,2))
  uniqs<-unique(c(dilutes$sample,dilutesBeta$sample))
  uniqs<-uniqs[uniqs %in% dilutes$sample & uniqs %in% dilutesBeta$sample]
  for(thisSample in uniqs){
    xx<-dilutes[dilutes$sample==thisSample,1:20]
    yy<-dilutesBeta[dilutesBeta$sample==thisSample,1:20]
    plotIfn(concAlpha,xx,thisSample,xlab='IFNa2 concentration (pg/ml)',ylims=ylims)
    plotIfn(concBeta,yy,thisSample,xlab='IFNb concentration (pg/ml)',ylims=ylims2)
  }
}
pdf('out/allCurves.pdf',width=6,height=4)
  plotIfns(dilutes,concAlpha,'IFNa2 concentration (pg/ml)')
dev.off()
pdf('out/allCurvesBeta.pdf',width=6,height=4)
  plotIfns(dilutesBeta,concBeta,'IFNb concentration (pg/ml)')
dev.off()

pdf('out/allCurvesCombo.pdf',width=10,height=5)
  plotDualIfns(dilutes,dilutesBeta,concAlpha,concBeta)
dev.off()

ic50Fits<-calcIc50s(dilutes,concAlpha)
ic50FitsBeta<-calcIc50s(dilutesBeta,concBeta)

