#!/usr/bin/env Rscript
library(dnar)
#args = commandArgs(trailingOnly=TRUE)
#if(length(args)!=2)stop("Incorrect arguments. Please specify:\n 1) a no header csv file with first column sample names and all additional columns giving sample1Conc1 sample2Conc1 sample1Conc2 ... \n 2) a no header csv file giving the concentration for each column of 1) (the first concentration should be 0)")

#https://www.myassays.com/four-parameter-logistic-regression.html
#http://stats.stackexchange.com/questions/61144/how-to-do-4-parametric-regression-for-elisa-data-in-r
#a = the minimum value that can be obtained (i.e. what happens at 0 dose)
#b = Hill’s slope of the curve (i.e. this is related to the steepness of the curve at point c).
#c = the point of inflection (i.e. the point on the S shaped curve halfway between a and d)
#d = the maximum value that can be obtained (i.e. what happens at infinite dose)
f<-function(B,x)(B[1]-B[4])/(1+(x/B[3])^B[2])+B[4]
fFixMax<-function(B,x,max=1)f(c(B,max),x)
LS<-function(B,y,x)sum((y-f(B,x))^2)

args<-c('data/ic50.csv','data/conc.csv')
ic50s<-read.csv(args[[1]],stringsAsFactors=FALSE,header=FALSE)
conc<-unlist(read.csv(args[[2]],header=FALSE)[1,])
samples<-ic50s[,1]
ic50s<-ic50s[,-1]
splits<-list(as.matrix(ic50s[,seq(1,ncol(ic50s),2)]),as.matrix(ic50s[,seq(2,ncol(ic50s),2)]))
props<-lapply(splits,function(xx)sweep(xx,1,xx[,1],'/'))
if(any(dim(props[[1]])!=dim(props[[2]])))stop("Mismatch in IC50 dimensions between replicates")
fakeConc<-10^seq(-10,10,.001)
out<-data.frame('sample'=samples,'IC50_1'=NA,'IC50_2'=NA,'Vres_1'=NA,'Vres_2'=NA,'analVres_1'=NA,'analVres_2'=NA,row.names=samples)
pdf('out/check.pdf',height=5,width=10)
for(ii in 1:nrow(props[[1]])){
  message(samples[ii])
  fit1<-suppressWarnings(nlminb(c(1,1,1,0),LS,x=conc[-1],y=props[[1]][ii,-1],lower=c(0,-Inf,-Inf,0),upper=c(1,Inf,Inf,1))$par)
  fit2<-suppressWarnings(nlminb(c(1,1,1,0),LS,x=conc[-1],y=props[[2]][ii,-1],lower=c(0,-Inf,-Inf,0),upper=c(1,Inf,Inf,1))$par)
  fitLine1<-f(fit1,fakeConc)
  fitLine2<-f(fit2,fakeConc)
  ic1<-approx(fitLine1,fakeConc,.5)$y
  ic2<-approx(fitLine2,fakeConc,.5)$y
  plot(1,1,type='n',log='x',xlim=range(conc[-1])*c(1e-3,1e3),ylim=c(0,1),xlab='Concentration',ylab='Proportion',xaxt='n',main=samples[ii],las=1)
  logAxis(1)
  points(conc[conc>0],props[[1]][ii,-1],pch=21,bg='#FF000099',cex=1.5,col=NA)
  points(conc[conc>0],props[[2]][ii,-1],pch=21,bg='#0000FF99',cex=1.5,col=NA)
  lines(fakeConc,fitLine1,col='#FF000066')
  lines(fakeConc,fitLine2,col='#0000FF66')
  abline(h=.5,lty=2)
  segments(ic1,par('usr')[3],ic1,.5,col='#FF000066',lty=2)
  segments(ic2,par('usr')[3],ic2,.5,col='#0000FF66',lty=2)
  out[samples[ii],c('IC50_1','IC50_2','Vres_1','Vres_2','analVres_1','analVres_2')]<-c(ic1,ic2,props[[1]][ii,length(conc)],props[[2]][ii,length(conc)],fit1[4],fit2[4])
}
dev.off()

out2<-data.frame('sample'=samples,'IC50_1'=NA,'IC50_2'=NA,'Vres_1'=NA,'Vres_2'=NA,'analVres_1'=NA,'analVres_2'=NA,row.names=samples)
pdf('out/check2.pdf',height=5,width=10)
for(ii in 1:nrow(splits[[1]])){
  message(samples[ii])
  fit1<-suppressWarnings(nlminb(c(max(splits[[1]][ii,]),1,1,1),LS,x=conc,y=splits[[1]][ii,],lower=c(0,-Inf,-Inf,0))$par)
  fit2<-suppressWarnings(nlminb(c(max(splits[[1]][ii,]),1,1,1),LS,x=conc,y=splits[[2]][ii,],lower=c(0,-Inf,-Inf,0))$par)
  fitLine1<-f(fit1,fakeConc)
  fitLine2<-f(fit2,fakeConc)
  ic2<-approx(fitLine2,fakeConc,.5)$y
  ic1<-fit1[3]
  fifty1<-f(fit1,ic1)
  ic2<-fit2[3]
  fifty2<-f(fit2,ic2)
  plot(1,1,type='n',log='x',xlim=range(conc[-1])*c(1e-3,1e3),ylim=c(0,max(c(splits[[1]][ii,],splits[[2]][ii,]))),xlab='Concentration',ylab='p24',xaxt='n',main=samples[ii],las=1)
  logAxis(1)
  points(conc[conc>0],splits[[1]][ii,-1],pch=21,bg='#FF000099',cex=1.5,col=NA)
  points(conc[conc>0],splits[[2]][ii,-1],pch=21,bg='#0000FF99',cex=1.5,col=NA)
  lines(fakeConc,fitLine1,col='#FF000066')
  lines(fakeConc,fitLine2,col='#0000FF66')
  segments(10^par('usr')[1],fifty1,ic1,fifty1,col='#FF000066',lty=2)
  segments(10^par('usr')[1],fifty2,ic2,fifty2,col='#0000FF66',lty=2)
  segments(ic1,par('usr')[3],ic1,fifty1,col='#FF000066',lty=2)
  segments(ic2,par('usr')[3],ic2,fifty2,col='#0000FF66',lty=2)
  out2[samples[ii],c('IC50_1','IC50_2','Vres_1','Vres_2','analVres_1','analVres_2')]<-c(ic1,ic2,props[[1]][ii,length(conc)],props[[2]][ii,length(conc)],fit1[4],fit2[4])
}
dev.off()


