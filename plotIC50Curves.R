library(xlsx)
library(dnar)

concAlpha<-c(0, 0.001, 0.003, 0.008, 0.023, 0.068, 0.204, 0.611, 1.833, 5.5)
concBeta<-c(0, 4.4E-05, 0.00044, 0.0044, 0.044, 0.44, 4.4, 44, 440, 4400)

#https://www.myassays.com/four-parameter-logistic-regression.html
#http://stats.stackexchange.com/questions/61144/how-to-do-4-parametric-regression-for-elisa-data-in-r
#a = the minimum value that can be obtained (i.e. what happens at 0 dose)
#b = Hill’s slope of the curve (i.e. this is related to the steepness of the curve at point c).
#c = the point of inflection (i.e. the point on the S shaped curve halfway between a and d)
#d = the maximum value that can be obtained (i.e. what happens at infinite dose)
f<-function(B,x)exp(log(B[1]-B[4])-log(1+(x/B[3])^B[2]))+B[4]
LS<-function(B,y,x){sum((log(y)-log(f(B,x)))^2)}
fit4par<-function(concs,p24s){
  fits<-lapply(list(c(max(p24s),1,1,min(p24s)),c(max(p24s)*2,1,1,min(p24s)/10),c(mean(p24s),1,1,mean(p24s)),c(max(p24s),.01,1,1)),function(starts)suppressWarnings(nlminb(starts,LS,x=concs,y=p24s,lower=c(0,-Inf,-Inf,0))))
  #fit<-suppressWarnings(nlminb(c(max(p24s),1,1,min(p24s)),LS,x=concs,y=p24s,lower=c(0,-Inf,-Inf,0))$par)
  obj<-sapply(fits,'[[','objective')
  return(fits[[which.min(obj)]]$par)
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
  dilutes$dilution<-sapply(dilutes$sample,function(xx)if(xx %in% dilutions$sample[dilutions$ifn=='alpha']) dilutions[dilutions$ifn=='alpha'&dilutions$sample==xx,'dilution'] else 50)
  dilutes[,1:20]<-dilutes[,1:20]*dilutes$dilution/1000
  dilutesBeta<-readIfns("VOA MM-cohort/IC50 Beta for all BULK isolates .xlsx")
  dilutesBeta$dilution<-sapply(dilutesBeta$sample,function(xx)if(xx %in% dilutions$sample[dilutions$ifn=='beta']) dilutions[dilutions$ifn=='beta'&dilutions$sample==xx,'dilution'] else 50)
  dilutesBeta[,1:20]<-dilutesBeta[,1:20]*dilutesBeta$dilution/1000
}

plotIfn<-function(concAlpha,p24s,main='',xlab='',ylims=range(p24s)){
  origConcs<-concs<-rep(concAlpha,each=2*nrow(xx))
  concs[concs==0]<-min(concs[concs>0])*.1
  plot(concs,p24s,xlab=xlab,ylab='',log='xy',las=1,xaxt='n',main=main,bg=c('red','green','blue','orange'),pch=21,mgp=c(2.5,1,0),ylim=ylims,yaxt='n')
  fit<-fit4par(origConcs,p24s)
  fit2<-suppressWarnings(nlminb(c(max(p24s),1,1,min(p24s)),LS,x=concs[origConcs!=0],y=p24s[origConcs!=0],lower=c(0,-Inf,-Inf,0))$par)
  fakeConc<-10^seq(-10,10,.001)
  fitLine<-f(fit,fakeConc)
  lines(fakeConc,fitLine,col='#FF000066',lwd=3)
  title(ylab='p24 concentration (ng/ml)',mgp=c(3.5,1,0))
  logAxis(1,axisMin=min(origConcs[origConcs>0]))
  logAxis(2,las=1)
  axis(1,min(origConcs[origConcs>0])*.1,0)
  abline(h=c(fit[c(1,4)]),lty=3,col='#00000055')
  abline(v=fit[3],lty=3,col='#00000055')
  #if(fit[4]<par('usr')[4])
  text(10^par('usr')[1]*2,10^par('usr')[3]*2,sprintf('Vres=%s IC50=%s',format(fit[4],digits=2),format(fit[3],digits=2)),adj=0)
  #if(thisSample=='MM55.12.2B1 bulk')browser()
}
plotIfns<-function(dilutes,concAlpha,xlab=''){
  ylims<-range(dilutes[,1:20],na.rm=TRUE)
  par(mar=c(3.5,4.5,1.5,1))
  for(thisSample in unique(dilutes$sample)){
    xx<-dilutes[dilutes$sample==thisSample,1:20]
    plotIfn(concAlpha,unlist(xx),thisSample,xlab=xlab,ylims=ylims)
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
    plotIfn(concAlpha,unlist(xx),thisSample,xlab='IFNa2 concentration (pg/ml)',ylims=ylims)
    plotIfn(concBeta,unlist(yy),thisSample,xlab='IFNb concentration (pg/ml)',ylims=ylims2)
  }
}

pdf('out/allCurves.pdf',width=5,height=5)
  plotIfns(dilutes,concAlpha,'IFNa2 concentration (pg/ml)')
dev.off()

pdf('out/allCurvesBeta.pdf',width=5,height=5)
  plotIfns(dilutesBeta,concBeta,'IFNb concentration (pg/ml)')
dev.off()

pdf('out/allCurvesCombo.pdf',width=10,height=5)
  plotDualIfns(dilutes,dilutesBeta,concAlpha,concBeta)
dev.off()
