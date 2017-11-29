
if(!exists('dat'))source('readNewData.R')
if(!exists('ic50Fits'))source('plotIC50Curves.R')

alphaIds<-lapply(rownames(ic50Fits),function(xx)which(dat$ID.for.Publications==xx))
if(any(sapply(alphaIds,length)!=1))stop('Problem matching alphas')
alphaIds<-unlist(alphaIds)
betaIds<-lapply(rownames(ic50FitsBeta),function(xx)which(dat$ID.for.Publications==xx))
if(any(sapply(betaIds,length)!=1))stop('Problem matching betas')
betaIds<-unlist(betaIds)

marvinAlpha<-dat[alphaIds,c('IFNa2.....Sept.2017.....IFNa2.Pooled.Donor.cells.IC50..pg.ml.','IFNa2.....Sept.2017.IFNa2.Pooled.Donor.Vres.at.5.5pg.ml...UT.')]
colnames(marvinAlpha)<-c('ic50','vres')
rownames(marvinAlpha)<-rownames(ic50Fits)

marvinBeta<-dat[betaIds,c('IFNbeta..Pooled.Donor.cells.IC50..pg.ml.','IFNbeta..Pooled.Donor.Vres.at.4.400.pg.ml...UT.')]
colnames(marvinBeta)<-c('ic50','vres')
rownames(marvinBeta)<-rownames(ic50FitsBeta)

pdf('out/marvin_vs_regression.pdf')
  plot(marvinAlpha$ic50,ic50Fits[,'ic50'],xlab='Marvin IC50',ylab='Regression IC50',main='IFNa2',log='xy',xaxt='n',yaxt='n')
  text(marvinAlpha$ic50,ic50Fits[,'ic50'],rownames(marvinAlpha),cex=.25,col='#FF000099')
  logAxis(1);logAxis(2,las=1)
  abline(0,1)
  plot(marvinAlpha$vres,ic50Fits[,'percVres'],xlab='Marvin Vres',ylab='Regression Vres',main='IFNa2',log='xy',xaxt='n',yaxt='n')
  text(marvinAlpha$vres,ic50Fits[,'percVres'],rownames(marvinAlpha),cex=.25,col='#FF000099')
  logAxis(1);logAxis(2,las=1)
  abline(0,1)
  plot(marvinBeta$ic50,ic50FitsBeta[,'ic50'],xlab='Marvin IC50',ylab='Regression IC50',main='IFNb',log='xy',xaxt='n',yaxt='n')
  text(marvinBeta$ic50,ic50FitsBeta[,'ic50'],rownames(marvinBeta),cex=.25,col='#FF000099')
  logAxis(1);logAxis(2,las=1)
  abline(0,1)
  plot(marvinBeta$vres,ic50FitsBeta[,'percVres'],xlab='Marvin Vres',ylab='Regression Vres',main='IFNb',log='xy',xaxt='n',yaxt='n')
  text(marvinBeta$vres,ic50FitsBeta[,'percVres'],rownames(marvinBeta),cex=.25,col='#FF000099')
  logAxis(1);logAxis(2,las=1)
  abline(0,1)
dev.off()

plotIfnsMarvin<-function(dilutes,concAlpha,xlab='',marvins){
  ylims<-range(dilutes[,1:20],na.rm=TRUE)
  par(mar=c(3.5,4.5,1.5,1),mfrow=c(1,2))
  for(thisSample in unique(dilutes$sample)){
    xx<-dilutes[dilutes$sample==thisSample,1:20]
    thisMax<-mean(unlist(xx[,1:2]))
    plotIfn(concAlpha,xx,thisSample,xlab=xlab,ylims=ylims)
    abline(v=marvins[thisSample,'ic50'],col='red',lty=2)
    abline(h=marvins[thisSample,'vres']*thisMax/100,col='red',lty=2)
    text(10^par('usr')[1]*2,10^par('usr')[3]*1.3,sprintf('Vres=%s IC50=%s (Marvin)',format(marvins[thisSample,'vres'],digits=2,width=3),format(marvins[thisSample,'ic50'],digits=2)),adj=0)
    plotIfn(concAlpha,xx,thisSample,xlab=xlab,log='x',scaleMax=TRUE,ylims=c(0,thisMax))
    abline(v=marvins[thisSample,'ic50'],col='red',lty=2)
    abline(h=marvins[thisSample,'vres']/100,col='red',lty=2)
    abline(h=0.5,lty=2)
  }
}
pdf('out/allCurvesMarvin.pdf',width=10,height=5)
  plotIfnsMarvin(dilutes,concAlpha,'IFNa2 concentration (pg/ml)',marvinAlpha)
dev.off()
pdf('out/allCurvesBetaMarvin.pdf',width=10,height=5)
  plotIfnsMarvin(dilutesBeta,concBeta,'IFNb concentration (pg/ml)',marvinBeta)
dev.off()


ic50Diffs<-list('IFNa2 IC50'=log2(marvinAlpha$ic50/ic50Fits[,'ic50']),'IFNb IC50'=log2(marvinBeta$ic50/ic50FitsBeta[,'ic50']))
vresDiffs<-list('IFNa2 Vres'=marvinAlpha$vres-ic50Fits[,'percVres'],'IFNb Vres'=marvinBeta$vres-ic50FitsBeta[,'percVres'])
lapply(ic50Diffs,function(xx)names(which(abs(xx)>2)))
