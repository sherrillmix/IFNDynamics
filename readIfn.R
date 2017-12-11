library(xlsx)


wb <- loadWorkbook("data/EJ MM plasma cytokine data Oct 2017 CORRECTED.xlsx")
ifns<-lapply(getSheets(wb),function(sheet){
  rows<-getCells(getRows(sheet),simplify=FALSE)
  vals<-lapply(rows,function(row){
    tmp<-sapply(row,function(xx)ifelse(is.null(xx),NA,getCellValue(xx)))
    #30 is arbitrary number to make all same width
    out<-rep(NA,50)
    names(out)<-1:50
    out[names(tmp)[names(tmp)!='']]<-tmp[names(tmp)!='']
    return(out)
  })
  dates<-sapply(vals,'[',3)
  dat<-do.call(rbind,vals[1:length(vals)>2&!is.na(dates)])
  cols<-c('sample','date','dfosx','vl','cd4','diluted','XXX','ifna1','ifna2','XXX','ifnb1','ifnb2','XXX','ifno1','ifno2','XXX','ifng1','ifng2')
  ifnCols<-grep('IFN',vals[[1]])
  if(any(ifnCols!=grep('^ifn.*1$',cols)))stop('Mismatch in ifn cols')
  colnames(dat)[1:length(cols)]<-cols
  dat<-dat[,cols[cols!='XXX']]
  dat[dat=='BDL']<- -Inf
  dat<-as.data.frame(dat,stringsAsFactors=FALSE)
  dat[,grepl('^ifn',colnames(dat))]<- apply(dat[,grepl('^ifn',colnames(dat))],2,as.numeric)
  dat$dfosx<-as.numeric(dat$dfosx)
  return(dat)
})
names(ifns)<-names(getSheets(wb))


lowers<-c('ifna'=12.5,'ifnb'=50,'ifno'=5,'ifng'=8)
ifnCols<-lapply(names(lowers),function(xx)sprintf('%s%d',xx,1:2))
names(ifnCols)<-sub('ifn','IFN',names(lowers))

ifnRanges<-lapply(ifnCols,function(xx){
  allVals<-unlist(lapply(ifns,'[',xx))
  allVals[allVals==-Inf]<-0
  range(allVals,na.rm=TRUE)
})

pdf('ifns.pdf',width=8)
for(ii in names(ifns)){
  par(mfrow=c(2,2),mar=c(4,4,2,.4))
  thisDat<-ifns[[ii]]
  thisDat[thisDat==-Inf]<-0
  for(jj in names(ifnCols)){
    plot(1,1,type='n',xlim=range(thisDat$dfosx),ylim=ifnRanges[[jj]],main=jj,las=1,ylab=jj,xlab='')
    title(xlab='DFOSx',mgp=c(2.2,1,0))
    for(col in ifnCols[[jj]])points(thisDat$dfosx,thisDat[,col])
    segments(thisDat$dfosx,thisDat[,ifnCols[[jj]][1]],thisDat$dfosx,thisDat[,ifnCols[[jj]][2]],col='#00000033')
    avg<-apply(thisDat[,ifnCols[[jj]]],1,mean)
    lines(thisDat$dfosx[!is.na(avg)],avg[!is.na(avg)])
  }
  title(main=ii,outer=TRUE,line=-1)
}
dev.off()

