library(xlsx)
library(dnar)
#concAlpha<-c(0, 0.001, 0.003, 0.008, 0.023, 0.068, 0.204, 0.611, 1.833, 5.5)
concAlpha<-c(0,5.55*3^(-8:0))
#concBeta<-c(0, 4.4E-05, 0.00044, 0.0044, 0.044, 0.44, 4.4, 44, 440, 4400)
concBeta<-c(0,4400*10^(-8:0))
lowerLimit<-50

concAlpha6<-concAlpha[c(1,seq(2,length(concAlpha),2))]
concBeta6<-concBeta[c(1,seq(2,length(concBeta),2))]


#https://www.myassays.com/four-parameter-logistic-regression.html
#http://stats.stackexchange.com/questions/61144/how-to-do-4-parametric-regression-for-elisa-data-in-r
#a = the minimum value that can be obtained (i.e. what happens at 0 dose)
#b = Hill’s slope of the curve (i.e. this is related to the steepness of the curve at point c).
#c = the point of inflection (i.e. the point on the S shaped curve halfway between a and d)
#d = the maximum value that can be obtained (i.e. what happens at infinite dose)
#e assymetry
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

readIfns<-function(file,exclude=1,firstCol=3,nameCol=firstCol-1,ifnCol=0,dilCol=0,minRows=10){
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
    if(any(sapply(vals[1:minRows],function(xx)is.null(xx[[3]])))){message('  Not found');return(NULL)}
    #assuming col 3, first row contains 1
    firstRow<-min(c(Inf,which(!is.na(sapply(vals[1:10],'[[',firstCol))&sapply(vals[1:10],'[[',firstCol)=='1')))+1
    if(firstRow==Inf){message('  Not found');return(NULL)}
    lastRow<-firstRow+min(c(Inf,which(sapply(vals[firstRow:40],function(xx)is.null(xx)||all(is.na(xx[firstCol+0:5])|is.null(xx[firstCol+0:5]))))))-1-1
    if(lastRow==Inf)return(NULL)
    #if((lastRow-firstRow+1) %%2 !=0)warning('Number of rows not a multiple of 2 on sheet ',counter,' of ',file)
    dat<-as.data.frame(apply(do.call(rbind,lapply(vals[firstRow:lastRow],function(zz)zz[firstCol+0:19])),2,function(xx){as.numeric(ifelse(xx==''|is.na(xx)|grepl('\\?\\?',xx),NA,sub('[><]','',xx)))}),stringsAsFactors=FALSE)
    dat$sample<-fillDown(sapply(vals[firstRow:lastRow],'[[',nameCol))
    dat$ifn<-if(ifnCol==0)NA else fillDown(sapply(vals[firstRow:lastRow],'[[',ifnCol),errorIfFirstEmpty=FALSE)
    dat$dil<-if(dilCol==0)NA else fillDown(sapply(vals[firstRow:lastRow],'[[',dilCol),errorIfFirstEmpty=FALSE)
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
plotIfn<-function(concAlpha,p24s,main='',xlab='',ylims=range(p24s),log='xy',scaleMax=FALSE,findVresIc50Func=findVresIc50,showVres=TRUE,showMax=FALSE,showLegend=FALSE,showPercent=grepl('y',log),ylab='p24 concentration (pg/ml)'){
  origConcs<-concs<-rep(concAlpha,each=ncol(p24s)/length(concAlpha)*nrow(p24s))
  if(scaleMax){
    p24s<-p24s/unlist(p24s[,1:2])
    maxs<-unlist(p24s[,1:2])
  } else maxs<-1
  vresIc50<-findVresIc50Func(concAlpha,p24s)
  zeroOffset<-.01
  fakeZero<-min(concs[concs>0])*zeroOffset
  concs[concs==0]<-fakeZero
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
  repMeans<-apply(p24s,1,function(xx)xx[seq(1,length(xx),2)]+xx[seq(2,length(xx),2)])/2
  p24s<-unlist(p24s)
  fit<-fitter(origConcs[!is.na(p24s)],p24s[!is.na(p24s)])
  plot(concs,p24s,xlab=xlab,ylab='',log=log,las=1,xaxt='n',main=main,bg=cols,pch=21,mgp=c(2,1,0),ylim=ylims,yaxt='n')
  #fit2<-suppressWarnings(nlminb(c(max(p24s),1,1,min(p24s)),LS,x=concs[origConcs!=0],y=p24s[origConcs!=0],lower=c(0,-Inf,-Inf,0))$par)
  if(scaleMax){
    apply(repMeans,2,function(xx)lines(ifelse(concAlpha==0,fakeZero,concAlpha),xx,col='#00000066'))
  }else{
    fakeConc<-10^seq(-10,10,.001)
    fitLine<-fitFunc(fit,fakeConc)
    lines(fakeConc,fitLine/mean(maxs),col='#00000066',lwd=3)
  }
  if(scaleMax)title(ylab='Proportion maximum p24',mgp=c(2,1,0))
  else title(ylab=ylab,mgp=c(ifelse(grepl('y',log),2.25,2.9),1,0))
  if(grepl('x',log))dnar::logAxis(1,axisMin=min(origConcs[origConcs>0]),mgp=c(2.5,.7,0))
  else axis(1,pretty(par('usr')[1:2]),mgp=c(2.5,.7,0))
  if(grepl('y',log))dnar::logAxis(2,las=1,mgp=c(2.5,.7,0))
  else axis(2,pretty(par('usr')[3:4]),las=1,mgp=c(2.5,.7,0))
  axis(1,min(origConcs[origConcs>0])*zeroOffset,0,mgp=c(2.5,.7,0))
  #abline(h=c(fit[c(1,4)]),lty=3,col='#00000055')
  #abline(v=fit[3],lty=3,col='#00000055')
  #isAbove<-vresIc50['ic50']>10^(par('usr')[3]*.75+par('usr')[4]*.25)
  isAbove<-vresIc50['vres']>10^mean(par('usr')[3:4])
  abline(v=vresIc50['ic50'],lty=2,col='#FF000055',lwd=1.5)
  xCoord<-10^(par('usr')[1]*.72+par('usr')[2]*.28)
  if(showVres){
    abline(h=vresIc50['vres']/mean(maxs),lty=2,col='#0000FF55',lwd=1.5)
    yCoord<-ifelse(grepl('y',log),vresIc50['vres']/2^ifelse(isAbove,1,-1),vresIc50['vres']+diff(par('usr')[3:4])*ifelse(isAbove,-.1,.1))
    text(xCoord,yCoord,sprintf('Vres=%s%%',format(vresIc50['percVres'],digits=2,width=3)),adj=c(1,0.5),col='blue')
    segments(xCoord*1.1,yCoord,10^(par('usr')[1]*.65+par('usr')[2]*.35),vresIc50['vres']/mean(maxs),col='blue')
  }
  if(showMax){
    abline(h=vresIc50['max'],lty=2,lwd=1.5,col='#00000055')
  }
  #ic50 label
  yCoord<-ifelse(grepl('y',log),ifelse(isAbove,10^(par('usr')[3]*.8+par('usr')[4]*.2),10^(par('usr')[3]*.7+par('usr')[4]*.3)*1.2),par('usr')[3]*.7+par('usr')[4]*.3)
  text(vresIc50['ic50']/1.85,yCoord,sprintf('IC50=%s',format(vresIc50['ic50'],digits=2)),adj=c(1,0.5),col='red')
  segments(vresIc50['ic50']/1.75,yCoord,vresIc50['ic50'],yCoord*ifelse(isAbove,1.2,.8),col='red')
  if(is.na(vresIc50['max']))browser()
  if(showPercent){
    if(grepl('y',log)){
      labs<-c(1,2,10,100,1000)
      yPos<-10^mean(par('usr')[3:4])
    }else{
      labs<-c(1,2,Inf)
      yPos<-mean(par('usr')[3:4])
    }
    axis(4,vresIc50['max']/labs,as.character(1/labs*100),las=1,tcl=-.2,mgp=c(0,.6,0))
    text(dnar::convertLineToUser(2.4,4),yPos,'Percent of untreated',xpd=NA,srt=-90)
  }
  abline(h=vresIc50['max']/2,lty=2,lwd=1.5,col='#00000055')
  par('lheight'=.72)
  if(showLegend)legend('topright',names(cols),pch=21,pt.bg=cols,inset=.01,cex=ifelse(length(cols)==2,1,.6),pt.cex=1,y.intersp=ifelse(length(cols)<=2,1,1.5))
  return(fit)
}
condenseReps<-function(xx){
  do.call(cbind,lapply(seq(1,ncol(xx),2),function(ii)apply(xx[,ii+0:1,drop=FALSE],1,mean,na.rm=TRUE)))
}
plotIfns<-function(dilutes,concAlpha,xlab='',condenseTechs=TRUE,log='xy',multiple=1,showMain=TRUE,ylims=range(dilutes[,1:nCol],na.rm=TRUE),...){
  nCol<-length(concAlpha)*2
  dilutes[,1:nCol]<-dilutes[,1:nCol]*multiple
  par(mar=c(3,3.8,1,3))
  fits<-list()
  for(thisSample in sort(unique(dilutes$sample))){
    xx<-dilutes[dilutes$sample==thisSample,1:nCol]
    if(!grepl('y',log))ylims<-c(0,ylims[2])
    fits[[thisSample]]<-plotIfn(concAlpha,xx,main=ifelse(showMain,thisSample,''),xlab=xlab,ylims=ylims,log=log,...)
    if(!grepl('y',log))ylims<-c(0,ylims[2])
    if(condenseTechs)plotIfn(concAlpha,condenseReps(xx),main=ifelse(showMain,thisSample,''),xlab=xlab,ylims=ylims,log=log,...)
  }
  invisible(fits)
}

plotDualIfns<-function(dilutes,dilutesBeta,concAlpha,concBeta){
  nCol<-length(concAlpha)*2
  ylims<-range(dilutes[,1:nCol],na.rm=TRUE)
  ylims2<-range(dilutesBeta[,1:nCol],na.rm=TRUE)
  par(mar=c(3.5,4.5,1.5,1),mfrow=c(1,2))
  uniqs<-unique(c(dilutes$sample,dilutesBeta$sample))
  uniqs<-uniqs[uniqs %in% dilutes$sample & uniqs %in% dilutesBeta$sample]
  for(thisSample in uniqs){
    xx<-dilutes[dilutes$sample==thisSample,1:nCol]
    yy<-dilutesBeta[dilutesBeta$sample==thisSample,1:nCol]
    plotIfn(concAlpha,xx,thisSample,xlab='IFNa2 concentration (pg/ml)',ylims=ylims)
    plotIfn(concBeta,yy,thisSample,xlab='IFNb concentration (pg/ml)',ylims=ylims2)
  }
}

GROWTH<-function(x,y,newX){
  mod<-lm(log(y)~x)
  exp(predict(mod,data.frame(x=newX)))
}

calculateBasicIc50<-function(concs,p24s,vocal=FALSE,means=condenseReps(p24s)){
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
calcBasicIc50<-function(...,dil=NULL){
  if(any(dil!=dil[1]))stop("Multiple dilutions")
  dil<-dil[1]
  out<-calculateBasicIc50(...)
  out<-out[c('ic50','percVres','max')]
  names(out)[2]<-'vres'
  if(!is.null(dil))out<-c(out,'repCap'=unname(out['max']*dil/1000))
  out
}

plotIfnStacked<-function(conc,p24,bioRep=rep(1,length(p24)),techRep=rep(1,length(p24)),main='',xlab='',ylims=range(p24),log='xy',scaleMax=FALSE,findVresIc50=findVresIc50){
  origConcs<-conc
  means<-tapply(p24,list(bioRep,conc),mean)
  vresIc50<-calculateBasicIc50(as.numeric(colnames(means)),means=means)
  zeroOffset<-.01
  conc[conc==0]<-min(conc[conc>0])*zeroOffset
  reps<-sprintf('Bio rep %d Tech rep %d',bioRep,techRep)
  repCols<-dnar::rainbow.lab(length(unique(reps)))
  names(repCols)<-sort(unique(reps))
  plot(conc,p24,xlab=xlab,ylab='',log=log,las=1,xaxt='n',main=main,bg=repCols[reps],pch=21,mgp=c(2,1,0),ylim=ylims,yaxt='n')
  title(ylab='p24 concentration (pg/ml)',mgp=c(ifelse(grepl('y',log),2.25,2.9),1,0))
  if(grepl('x',log))dnar::logAxis(1,axisMin=min(origConcs[origConcs>0]),mgp=c(2.5,.7,0))
  else axis(1,pretty(par('usr')[1:2]),mgp=c(2.5,.7,0))
  if(grepl('y',log))dnar::logAxis(2,las=1,mgp=c(2.5,.7,0))
  else axis(2,pretty(par('usr')[3:4]),las=1,mgp=c(2.5,.7,0))
  axis(1,min(origConcs[origConcs>0])*zeroOffset,0,mgp=c(2.5,.7,0))
  isAbove<-vresIc50['vres']>10^mean(par('usr')[3:4])
  abline(v=vresIc50['ic50'],lty=2,col='#FF000055')
  abline(h=vresIc50['vres'],lty=2,col='#0000FF55')
  yCoord<-ifelse(grepl('y',log),vresIc50['vres']/2^ifelse(isAbove,1,-1),vresIc50['vres']+diff(par('usr')[3:4])*ifelse(isAbove,-.1,.1))
  xCoord<-10^(par('usr')[1]*.72+par('usr')[2]*.28)
  text(xCoord,yCoord,sprintf('Vres=%s%%',format(vresIc50['percVres'],digits=2,width=3)),adj=c(1,0.5),col='blue')
  segments(xCoord*1.1,yCoord,10^(par('usr')[1]*.65+par('usr')[2]*.35),vresIc50['vres'],col='blue')
  #ic50 label
  yCoord<-ifelse(grepl('y',log),ifelse(isAbove,10^(par('usr')[3]*.8+par('usr')[4]*.2),10^(par('usr')[3]*.7+par('usr')[4]*.3)*1.2),par('usr')[3]*.7+par('usr')[4]*.3)
  text(vresIc50['ic50']/1.85,yCoord,sprintf('IC50=%s',format(vresIc50['ic50'],digits=2)),adj=c(1,0.5),col='red')
  segments(vresIc50['ic50']/1.75,yCoord,vresIc50['ic50'],yCoord*ifelse(isAbove,1.2,.8),col='red')
  if(is.na(vresIc50['max']))browser()
  if(grepl('y',log))axis(4,vresIc50['max']/c(1,2,4,10,100,1000),as.character(c(1,.5,.25,.1,.01,.001)*100),las=1,tcl=-.2,mgp=c(0,.6,0))
  abline(h=vresIc50['max']/2,lty=3,col='#00000055')
  text(dnar::convertLineToUser(2.6,4),10^mean(par('usr')[3:4]),'Percent of maximum',xpd=NA,srt=-90)
  par('lheight'=.72)
  legend('topright',names(repCols),pch=21,pt.bg=repCols,inset=.01,cex=ifelse(length(repCols)==2,1,.6),pt.cex=1,y.intersp=ifelse(length(repCols)<=2,1,1.5))
  fit<-fitter(origConcs,p24)
  fakeConc<-10^seq(-10,10,.001)
  fitLine<-fitFunc(fit,fakeConc)
  lines(fakeConc,fitLine,col='#00000066',lwd=3)
  #text(fit[3],fit[1]/2,sprintf('Line IC50: %0.3f',fit[3]))
  #abline(v=fit[3],h=fit[1]/2+fit[4]/2,lty=3)
  fifty<-fakeConc[which.min(abs(fitLine-fit[1]/2))]
  #abline(v=fifty,lty=3)
  return(c(fit[3],fifty))
}

read6ConcIfns<-function(file,exclude=1,firstCol=3,nameCol=firstCol-1,ifnCol=0,dilCol=0,minRows=8){
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
    if(any(sapply(vals[1:minRows],function(xx)is.null(xx[[3]])))){message('  Not found');return(NULL)}
    #assuming col 3, first row contains 1
    firstRow<-min(c(Inf,which(!is.na(sapply(vals[1:10],'[[',firstCol))&sapply(vals[1:10],'[[',firstCol)=='1')))+1
    if(firstRow==Inf){message('  Not found');return(NULL)}
    lastRow<-firstRow+min(c(Inf,which(sapply(vals[firstRow:40],function(xx)is.null(xx)||all(is.na(xx[firstCol+0:5])|is.null(xx[firstCol+0:5]))))))-1-1
    if(lastRow==Inf)return(NULL)
    #if((lastRow-firstRow+1) %%2 !=0)warning('Number of rows not a multiple of 2 on sheet ',counter,' of ',file)
    dat<-as.data.frame(apply(do.call(rbind,lapply(vals[firstRow:lastRow],function(zz)zz[firstCol+0:23])),2,function(xx){as.numeric(ifelse(xx==''|is.na(xx)|grepl('\\?\\?',xx),NA,sub('[><]','',xx)))}),stringsAsFactors=FALSE)
    dat$sample<-fillDown(sapply(vals[firstRow:lastRow],'[[',nameCol))
    if(ifnCol!=0)dat$ifn<-fillDown(sapply(vals[firstRow:lastRow],'[[',ifnCol),errorIfFirstEmpty=FALSE)
    if(dilCol!=0)dat$dilution<-fillDown(sapply(vals[firstRow:lastRow],'[[',dilCol),errorIfFirstEmpty=FALSE)
    return(dat)
  })
  names(ic50Curves)<-names(getSheets(wb))
  dilutes<-do.call(rbind,ic50Curves)
  dilutes$sheet<-rep(names(ic50Curves),sapply(ic50Curves,function(xx)ifelse(is.null(xx),0,nrow(xx))))
  dilutes$sample<-sprintf('%s (%s)',dilutes$sample,dilutes$sheet)
  return(dilutes)
}

convertIfn6<-function(p24s){
  nConc<-ncol(p24s)/2
  bio1<-p24s[,1:nConc,drop=FALSE]
  bio2<-p24s[,nConc+1:nConc,drop=FALSE]
  colnames(bio1)<-colnames(bio2)<-1:nConc
  do.call(rbind,lapply(1:nrow(bio1),function(xx)rbind(bio1[xx,,drop=FALSE],bio2[xx,,drop=FALSE])))
}

readCounts<-function(countFile){
  counts<-read.csv(countFile,stringsAsFactors=FALSE)
  counts$plate<-basename(counts$dir)
  counts$well<-sub('\\.CTL$','',counts$file)
  counts$col<-as.numeric(sub('^[A-Z]','',counts$well))
  counts$row<-sub('[0-9]+$','',counts$well)
  letterLookup<-structure(1:26,.Names=LETTERS)
  counts$rowNum<-letterLookup[counts$row]
  return(counts)
}
readPlateViruses<-function(plateFile,virusFile){
  plate<-read.csv(plateFile,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
  viruses<-read.csv(virusFile,stringsAsFactors=FALSE)
  plate<-apply(plate,2,trimws)
  viruses<-rbind(viruses,c('Empty','Empty'))
  rownames(viruses)<-viruses$id
  if(any(!unlist(plate) %in% viruses$id))warning('Unknown virus ',unique(unlist(plate)[!unlist(plate) %in% viruses$id]))
  plateIds<-data.frame('row'=rep(rownames(plate),ncol(plate)),'col'=rep(colnames(plate),each=nrow(plate)),'vId'=as.vector(unlist(plate)),stringsAsFactors=FALSE)
  plateIds$well<-sprintf('%s%s',plateIds$row,plateIds$col)
  rownames(plateIds)<-plateIds$well
  plateIds$virus<-viruses[plateIds$vId,'sample']
  return(plateIds)
}

#source('iuStan.R')
runIuStan<-function(virus,tit){
  times<-sort(unique(tit$time))
  do.call(rbind,lapply(structure(times,.Names=times),function(time){
    thisDat<-tit[!is.na(tit$virus)&tit$virus==virus&tit$time==time,c('n','dilution')]
    fit<-simpleCountIU(iuModSimple,thisDat$n,thisDat$dilution,tit[grepl('MEDIA|Media',tit$virus)&!is.na(tit$virus),'n'])
    means<-mean(as.matrix(fit)[,'baseIU'])/100
    sds<-sd(as.matrix(fit)[,'baseIU'])/100
    return(c(mean=means,sd=sds))
  }))
}
runAllIu<-function(viruses,tit,mc.cores=20){
  out<-parallel::mclapply(viruses,runIuStan,tit,mc.cores=20)
  means<-do.call(rbind,lapply(out,function(xx)xx[,'mean']))
  sds<-do.call(rbind,lapply(out,function(xx)xx[,'sd']))
  return(list(mean=means,sd=sds))
}
plotRawIce<-function(tit,viruses,ius,sds=NULL){
  timeCols<-structure(dnar::rainbow.lab(length(unique(tit$time))),.Names=sort(unique(tit$time)))
  fakeDils=2^seq(0,16,length.out=1000)
  iuRange<-range(ius)
  par(mfrow=c(1,2),mar=c(4,4,1,3))
  for(virus in unique(viruses)){
    thisDat<-tit[tit$virus==virus&!is.na(tit$virus),]
    dnar::withAs(zz=thisDat,plot(zz$dilution,zz$n+1,xlab='Dilution',ylab='TZMBL count',las=1,log='yx',main=sprintf('%s',virus),ylim=range(tit$n+1),xlim=c(1,max(tit$dilution,na.rm=TRUE)),pch=21,bg=timeCols[as.character(zz$time)],cex=2,yaxt='n'))
    dnar::logAxis(2,offset=1,axisMin=1,las=1)
    axis(2,1,0,las=1)
    legend('topright',pch=21,names(timeCols),pt.bg=timeCols,bty='n')
    for(ii in names(timeCols)){
      preds<-ius[virus,ii]/fakeDils*100
      lines(fakeDils,preds+1,col=timeCols[ii])
    }
    plot(as.numeric(colnames(ius)),ius[virus,],xlab='Time on ice',ylab='Estimated IU/ul',las=1,log='y',ylim=iuRange,yaxt='n',type='b')
    dnar::logAxis(2,las=1)
    prettyY<-c(1,.5,.1,.01,.001)
    axis(4,prettyY*ius[virus,'0'],prettyY,las=1)
    if(any(ius[virus,]>1)){
      abline(h=ius[virus,'0']*.5,lty=2)
      times<-as.numeric(colnames(ius))
      half<-approx(ius[virus,]/ius[virus,'0'],times,.5)$y
      if(!is.null(sds)){
        thisSd<-sds[virus,colnames(ius)]
        segments(times,pmax(10^(par('usr')[3]),ius[virus,]-2*thisSd),times,ius[virus,]+2*thisSd)
      }
      #n<-log10(ius[virus,])
      #mod<-lm(n~times)
      #fakeTime<-seq(0,60,.01) 
      #preds<-predict(mod,data.frame(times=fakeTime),interval='confidence')
      #lines(fakeTime,10^(preds[,1]))
      #polygon(c(fakeTime,rev(fakeTime)),10^(c(preds[,2],rev(preds[,3]))),border=NA,col='#00000033')
      #half<-log10(.5)/coefficients(mod)[2]
      abline(v=half,lty=3)
      title(main=sprintf('Ice half life: %s hr',ifelse(is.na(half),'>48',format(half,digits=3))))
    }
  }
}

plotStackedIfns<-function(dilutes,concAlpha,xlab='',cols=NULL,pch=NULL,...){
  uniqSamples<-unique(dilutes$sample)
  cols<-c(cols,structure(rep('black',sum(!uniqSamples %in% names(cols))),.Names=uniqSamples[!uniqSamples %in% names(cols)]))
  pch<-c(pch,structure(rep(21,sum(!uniqSamples %in% names(pch))),.Names=uniqSamples[!uniqSamples %in% names(pch)]))
  nCol<-length(concAlpha)*2 #assume two tech replicates
  par(mar=c(3,3.8,.1,.1))
  zeroOffset<-.1
  concs<-concAlpha
  fakeZero<-min(concs[concs>0])*zeroOffset
  concs[concs==0]<-fakeZero
  props<-do.call(rbind,by(dilutes[,1:nCol],dilutes$sample,function(xx){
    bioRepMeans<-(xx[seq(1,length(xx),2)]+xx[seq(2,length(xx),2)])/2
    bioRepProp<-bioRepMeans/bioRepMeans[,1]
    apply(bioRepProp,2,mean)
  }))
  plot(rep(concs,each=nrow(props)),unlist(props),log='x',xaxt='n',las=1,mgp=c(2,1,0),xlab=xlab,ylab='',type='n')
  title(ylab='Proportion of maximum p24',mgp=c(2.6,1,0))
  for(ii in rownames(props)){
    lines(concs,props[ii,],col=cols[ii],lwd=2)
    points(concs,props[ii,],bg=cols[ii],pch=pch[ii])
  }
  dnar::logAxis(1)
  abline(h=.5,lty=2)
  legend('bottomleft',names(cols),inset=.01,lty=1,col=cols,lwd=2,pch=pch[names(cols)],pt.bg=cols,...)
}
