library(xlsx)
library(lubridate)
wb <- loadWorkbook("meta/EJ MM plasma cytokine data CORRECTED updated VL CD4 Jan2018.xlsx")

metas<-lapply(getSheets(wb),function(sheet){
  rows<-getCells(getRows(sheet),simplify=FALSE)
  vals<-lapply(rows,function(row){
    tmp<-sapply(row,function(xx)ifelse(is.null(xx),NA,getCellValue(xx)))
    #30 is arbitrary number to make all same width
    out<-rep(NA,50)
    names(out)<-1:50
    out[names(tmp)[names(tmp)!='']]<-tmp[names(tmp)!='']
    return(out)
  })
  dates<-sapply(vals,'[',2)
  goodDates<-!is.na(suppressWarnings(as.numeric(dates)))
  isStringDate<-all(!goodDates)
  if(isStringDate)goodDates<-grepl('[0-9]{2}[./][0-9]{2}[./][12]?[90]?[0-9]{2}',dates)
  dat<-do.call(rbind,vals[1:length(vals)>2&goodDates])
  if(is.null(dat))browser()
  cols<-c('sample','date','dfosx','oldViralLoad','viralLoad','cd4','diluted','XXX','ifna1','ifna2','XXX','ifnb1','ifnb2','XXX','ifno1','ifno2','XXX','ifng1','ifng2')
  if(!vals[[1]][[5]] %in% c('Corrected','corrected','Viral load (copies/ml)'))cols<-cols[-4]
  if(is.na(vals[[1]][[3]]))cols[3:4]<-c('newDate','dfosx')
  ifnCols<-grep('IFN',vals[[1]])
  colnames(dat)[1:length(cols)]<-cols
  dat<-dat[,cols[cols!='XXX']]
  dat[dat=='BDL']<- -Inf
  dat<-as.data.frame(dat,stringsAsFactors=FALSE)
  dat[,grepl('^ifn',colnames(dat))]<- apply(dat[,grepl('^ifn',colnames(dat))],2,as.numeric)
  dat$dfosx<-as.numeric(dat$dfosx)
  if(isStringDate)dat$rDate<-dmy(dat$date)
  else dat$rDate<-as.Date(as.numeric(dat$date),origin='1899-12-30')
  #print(dat$cd4)
  dat$cd4<-as.numeric(ifelse(dat$cd4 %in% c('not done','no data'),NA,dat$cd4))
  dat$vl<-as.numeric(ifelse(dat$viralLoad %in% c('not done','no data'),NA,gsub('[><,]','',dat$viralLoad)))
  if(any(colnames(dat)=='oldViralLoad')){
    oldVl<-as.numeric(gsub('[<>,]','',dat$oldViralLoad))
    if(any(!is.na(dat$vl)&(abs(log2(dat$vl/oldVl))>1)))stop('Big difference in old and new vl')
    dat[is.na(dat$vl),'vl']<-oldVl[is.na(dat$vl)]
  }
  return(dat)
})
names(metas)<-names(getSheets(wb))
lapply(metas,'[','cd4')

