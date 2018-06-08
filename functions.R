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
LS<-function(B,y,x){sum(abs((log(y)-log(fitFunc(B,x)))))}
optIc50<-function(B)exp(optim(1,function(x,B){abs(B[1]*.5-fitFunc(B,exp(x)))},B=B,method='Brent',lower=-20,upper=20)$par)
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
  fits<-lapply(list(c(max(p24s),1,1,min(p24s),1),c(max(p24s)*2,1,1,min(p24s)/10),c(mean(p24s),1,1,mean(p24s),1),c(max(p24s),.01,1,1,1)),function(starts)suppressWarnings(nlminb(starts,LS,x=concs,y=p24s,lower=c(0,-Inf,-Inf,0,-Inf))))
  #fit<-suppressWarnings(nlminb(c(max(p24s),1,1,min(p24s)),LS,x=concs,y=p24s,lower=c(0,-Inf,-Inf,0))$par)
  obj<-sapply(fits,'[[','objective')
  return(fits[[which.min(obj)]]$par)
}
fitFunc<-f5
fitter<-fit5par

findVresIc50<-function(concAlpha,p24s){
  origConcs<-rep(concAlpha,each=ncol(p24s)/length(concAlpha)*nrow(p24s))
  p24s<-unlist(p24s)
  fit<-fitter(origConcs[!is.na(p24s)],p24s[!is.na(p24s)])
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

readIfns<-function(file,exclude=1,firstCol=3,nameCol=firstCol-1,ifnCol=0){
  wb <- loadWorkbook(file)
  counter<-0
  ic50Curves<-lapply(getSheets(wb),function(sheet){
    counter<<-counter+1
    if(counter %in% exclude)return(NULL)
    message(' Sheet ',counter)
    #can get weird error otherwise
    #rows<-lapply(1:40,function(xx)tryCatch(getCells(getRows(sheet)[xx])[[1]],error=function(e)return(NULL)))
    if(length(getRows(sheet))==0)return(NULL)
    rows<-getCells(getRows(sheet,1:40),1:50,simplify=FALSE)
    vals<-lapply(rows,function(row){
      tmp<-sapply(row,function(xx)ifelse(is.null(xx),NA,getCellValue(xx)))
      out<-rep(NA,50)
      names(out)<-1:50
      out[names(tmp)[names(tmp)!='']]<-tmp[names(tmp)!='']
      return(out)
    })
    if(any(sapply(vals[1:10],function(xx)is.null(xx[[3]])))){message('  Not found');return(NULL)}
    #assuming col 3, first row contains 1
    firstRow<-min(c(Inf,which(!is.na(sapply(vals[1:10],'[[',firstCol))&sapply(vals[1:10],'[[',firstCol)=='1')))+1
    if(firstRow==Inf){message('  Not found');return(NULL)}
    lastRow<-firstRow+min(c(Inf,which(sapply(vals[firstRow:40],function(xx)is.null(xx)||all(is.na(xx[firstCol+0:5])|is.null(xx[firstCol+0:5]))))))-1-1
    if(lastRow==Inf)return(NULL)
    #if((lastRow-firstRow+1) %%2 !=0)warning('Number of rows not a multiple of 2 on sheet ',counter,' of ',file)
    dat<-as.data.frame(apply(do.call(rbind,lapply(vals[firstRow:lastRow],function(zz)zz[firstCol+0:19])),2,function(xx){as.numeric(ifelse(xx==''|is.na(xx)|grepl('\\?\\?',xx),NA,sub('[><]','',xx)))}),stringsAsFactors=FALSE)
    dat$sample<-fillDown(sapply(vals[firstRow:lastRow],'[[',nameCol))
    dat$ifn<-fillDown(sapply(vals[firstRow:lastRow],'[[',ifnCol),errorIfFirstEmpty=FALSE)
    return(dat)
  })
  names(ic50Curves)<-names(getSheets(wb))
  dilutes<-do.call(rbind,ic50Curves)
  dilutes$sheet<-rep(names(ic50Curves),sapply(ic50Curves,function(xx)ifelse(is.null(xx),0,nrow(xx))))
  dilutes$sample<-sprintf('%s (%s)',dilutes$sample,dilutes$sheet)
  return(dilutes)
}

#assuming 
#bio#1-tech#1 bio#1-tech#2
#bio#2-tech#1 bio#2-tech#2
plotIfn<-function(concAlpha,p24s,main='',xlab='',ylims=range(p24s),log='xy',scaleMax=FALSE,findVresIc50=findVresIc50){
  origConcs<-concs<-rep(concAlpha,each=ncol(p24s)/length(concAlpha)*nrow(p24s))
  if(scaleMax)maxs<-unlist(p24s[,1:2])
  else maxs<-1
  vresIc50<-findVresIc50(concAlpha,p24s)
  zeroOffset<-.01
  concs[concs==0]<-min(concs[concs>0])*zeroOffset
  if(nrow(p24s)*ncol(p24s)/length(concAlpha)==8){
    cols<-rainbow.lab(17)[c(1:2,6:7,11:12,16:17)]
    names(cols)<-c('Bio replicate 1\nTech replicate 1','Bio replicate 1\nTech replicate 2','Bio replicate 2\nTech replicate 1','Bio replicate 2\nTech replicate 2','Bio replicate 3\nTech replicate 1','Bio replicate 3\nTech replicate 2','Bio replicate 4\nTech replicate 1','Bio replicate 4\nTech replicate 2')
    cols<-cols[c(seq(1,8,2),seq(2,8,2))]
  }else if(nrow(p24s)*ncol(p24s)/length(concAlpha)==4){
    cols<-rainbow.lab(12)[c(1:2,11:12)]
    names(cols)<-c('Bio replicate 1\nTech replicate 1','Bio replicate 1\nTech replicate 2','Bio replicate 2\nTech replicate 1','Bio replicate 2\nTech replicate 2')
    cols<-cols[c(seq(1,4,2),seq(2,4,2))]
  }else if(nrow(p24s)*ncol(p24s)/length(concAlpha)==2){
    cols<-rainbow.lab(2)
    names(cols)<-c('Replicate 1','Replicate 2')
  }else if(nrow(p24s)*ncol(p24s)/length(concAlpha)==1){
    cols<-rainbow.lab(2)[1]
    names(cols)<-c('Replicate 1')
  }else{
    browser()
    stop('p24s not 2 or 4 rows')
  }
  p24s<-unlist(p24s)
  fit<-fitter(origConcs[!is.na(p24s)],p24s[!is.na(p24s)])
  plot(concs,p24s/maxs,xlab=xlab,ylab='',log=log,las=1,xaxt='n',main=main,bg=cols,pch=21,mgp=c(2,1,0),ylim=ylims/mean(maxs),yaxt='n')
  #fit2<-suppressWarnings(nlminb(c(max(p24s),1,1,min(p24s)),LS,x=concs[origConcs!=0],y=p24s[origConcs!=0],lower=c(0,-Inf,-Inf,0))$par)
  fakeConc<-10^seq(-10,10,.001)
  fitLine<-fitFunc(fit,fakeConc)
  lines(fakeConc,fitLine/mean(maxs),col='#00000066',lwd=3)
  if(scaleMax)title(ylab='Proportion maximum p24',mgp=c(2,1,0))
  else title(ylab='p24 concentration (ng/ml)',mgp=c(ifelse(grepl('y',log),2.25,2.9),1,0))
  if(grepl('x',log))logAxis(1,axisMin=min(origConcs[origConcs>0]),mgp=c(2.5,.7,0))
  else axis(1,pretty(par('usr')[1:2]),mgp=c(2.5,.7,0))
  if(grepl('y',log))logAxis(2,las=1,mgp=c(2.5,.7,0))
  else axis(2,pretty(par('usr')[3:4]),las=1,mgp=c(2.5,.7,0))
  axis(1,min(origConcs[origConcs>0])*zeroOffset,0,mgp=c(2.5,.7,0))
  #abline(h=c(fit[c(1,4)]),lty=3,col='#00000055')
  #abline(v=fit[3],lty=3,col='#00000055')
  #isAbove<-vresIc50['ic50']>10^(par('usr')[3]*.75+par('usr')[4]*.25)
  isAbove<-vresIc50['vres']>10^mean(par('usr')[3:4])
  abline(v=vresIc50['ic50'],lty=2,col='#FF000055')
  abline(h=vresIc50['vres']/mean(maxs),lty=2,col='#0000FF55')
  yCoord<-ifelse(grepl('y',log),vresIc50['vres']/2^ifelse(isAbove,1,-1),vresIc50['vres']+diff(par('usr')[3:4])*ifelse(isAbove,-.1,.1))
  xCoord<-10^(par('usr')[1]*.72+par('usr')[2]*.28)
  text(xCoord,yCoord,sprintf('Vres=%s%%',format(vresIc50['percVres'],digits=2,width=3)),adj=c(1,0.5),col='blue')
  segments(xCoord*1.1,yCoord,10^(par('usr')[1]*.65+par('usr')[2]*.35),vresIc50['vres']/mean(maxs),col='blue')
  #ic50 label
  yCoord<-ifelse(grepl('y',log),ifelse(isAbove,10^(par('usr')[3]*.8+par('usr')[4]*.2),10^(par('usr')[3]*.7+par('usr')[4]*.3)*1.2),par('usr')[3]*.7+par('usr')[4]*.3)
  text(vresIc50['ic50']/1.85,yCoord,sprintf('IC50=%s',format(vresIc50['ic50'],digits=2)),adj=c(1,0.5),col='red')
  segments(vresIc50['ic50']/1.75,yCoord,vresIc50['ic50'],yCoord*ifelse(isAbove,1.2,.8),col='red')
  if(is.na(vresIc50['max']))browser()
  if(grepl('y',log))axis(4,vresIc50['max']/c(1,2,10,100,1000),as.character(c(1,.5,.1,.01,.001)*100),las=1,tcl=-.2,mgp=c(0,.6,0))
  abline(h=vresIc50['max']/2,lty=3,col='#00000055')
  text(convertLineToUser(2.6,4),10^mean(par('usr')[3:4]),'Percent of maximum',xpd=NA,srt=-90)
  par('lheight'=.72)
  legend('topright',names(cols),pch=21,pt.bg=cols,inset=.01,cex=ifelse(length(cols)==2,1,.6),pt.cex=1,y.intersp=ifelse(length(cols)<=2,1,1.5))
  return(fit)
}
condenseReps<-function(xx){
  do.call(cbind,lapply(seq(1,ncol(xx),2),function(ii)apply(xx[,ii+0:1,drop=FALSE],1,mean,na.rm=TRUE)))
}
plotIfns<-function(dilutes,concAlpha,xlab='',condenseTechs=TRUE,log='xy',multiple=1,...){
  ylims<-range(dilutes[,1:20],na.rm=TRUE)
  dilutes[,1:20]<-dilutes[,1:20]*multiple
  par(mar=c(3,3.8,1,3))
  fits<-list()
  for(thisSample in sort(unique(dilutes$sample))){
    xx<-dilutes[dilutes$sample==thisSample,1:20]
    if(!grepl('y',log))ylims<-c(0,max(xx,na.rm=TRUE))
    fits[[thisSample]]<-plotIfn(concAlpha,xx,main=thisSample,xlab=xlab,ylims=ylims,log=log,...)
    if(!grepl('y',log))ylims<-c(0,max(condenseReps(xx),na.rm=TRUE))
    if(condenseTechs)plotIfn(concAlpha,condenseReps(xx),main=thisSample,xlab=xlab,ylims=ylims,log=log,...)
  }
  invisible(fits)
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

GROWTH<-function(x,y,newX){
  mod<-lm(log(y)~x)
  exp(predict(mod,data.frame(x=newX)))
}

calculateBasicIc50<-function(concs,p24s,vocal=FALSE){
  means<-condenseReps(p24s)
  props<-t(apply(means,1,function(xx)xx/xx[1]))
  #ic50<-mean(apply(props,1,function(xx)approx(xx,concs,.5)$y))
  if(vocal)print(t(props))
  ic50s<-apply(props,1,function(xx){
    id<-min(c(Inf,which(xx<.5)))
    if(id==Inf)return(NA)
    if(any(concs[id+-1:0]==0))return(NA)
    GROWTH(xx[id+-1:0],concs[id+-1:0],.5)
  })
  if(vocal)print(ic50s)
  ic50<-mean(ic50s)
  percVres<-mean(props[,ncol(props)])*100
  vres<-mean(means[,ncol(props)])
  max<-mean(means[,1])
  return(c('ic50'=ic50,'vres'=vres,'percVres'=percVres,'max'=max))
}
calcBasicIc50<-function(...){out<-calculateBasicIc50(...);out<-out[c('ic50','percVres')];names(out)[2]<-'vres';out}


