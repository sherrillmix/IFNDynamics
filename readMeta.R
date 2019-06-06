library(dnar)
library(xlsx)
library(lubridate)

convertUKMM<-function(xx,mmLookup,errorOnNA=TRUE){
  splits<-strsplit(xx,'\\.')
  splits<-lapply(splits,function(xx){
    orig<-xx
    if(length(xx)==0)return('')
    xx[1]<-sub('UK','EJ',gsub(' ','',xx[1]))
    if(grepl('EJ',xx[1]))xx[1]<-mmLookup[xx[1]]
    if(is.na(xx[1])&errorOnNA)stop('Problem converting ID "',paste(orig,collapse='.'),'"')
    return(xx)
  })
  return(sapply(splits,paste,collapse='.'))
}
fixDecimals<-function(xx,maxVisit=29,checkRank=TRUE){
  splits<-strsplit(xx,'\\.')
  visits<-as.numeric(sub('[oO]','0',sapply(splits,'[',2)))
  pats<-sapply(splits,'[',1)
  fixVisit<-ave(visits,pats,FUN=function(xx){
    out<-xx
    #catch crazy visits e.g. 70
    out[out>maxVisit&out%%10]<-out[out>maxVisit&out%%10]/10
    #catch too high values
    potentialProbs<-which(xx%%10==0)
    for(ii in rev(potentialProbs)){
      if(ii==length(xx))break
      #fix too high e.g. 20 instead of 2
      if(any(out[(ii+1):length(out)]<out[ii]))out[ii]<-out[ii]/10
    }
    #catch too low value
    probs<-out<cummax(out)
    if(any(probs))out[probs]<-out[probs]*10
    if(checkRank&&any(rank(out)!=1:length(out)))stop('Problem fixing visits: ',paste(xx,collapse=', '))
    out
  })
  return(sprintf('%s.%s%s',pats,sprintf('%02d',fixVisit),sapply(splits,function(xx)ifelse(length(xx)>2,sprintf('.%s',paste(xx[-2:-1],collapse='.')),''))))
}

##MASTER META LIST##
meta1<-read.csv('data/For Scot, Complete Master table AUG.2017_meta.csv',stringsAsFactors=FALSE)[,-1:-2]
meta1<-meta1[,1:6]
#meta1$id<-fillDown(meta1$ID)
meta1$id<-sapply(strsplit(meta1$Time.Points,'\\.'),'[',1)
meta1<-meta1[meta1$Time.Points!='Total number of sequences',]
rownames(meta1)<-meta1$Time.Points

meta2<-read.csv('data/New MM cohort patients.csv',stringsAsFactors=FALSE)
#meta2<-meta2[meta2$Date!=''&!is.na(meta2$Date)&meta2$Date!='Date',]
meta2<-meta2[meta2[,2]!='',]
colnames(meta2)[1:2]<-c('ID','Time.Points')
colnames(meta2)[colnames(meta2)=='Viral.load']<-'VL'
colnames(meta2)[colnames(meta2)=='CD4.count']<-'CD4'
meta2<-meta2[,1:6]
meta2$id<-fillDown(meta2$ID)
meta2$Time.Points<-sprintf('MM%s',meta2$Time.Points)
tmp<-sub('\\.([0-9])$','.0\\1',mapply(function(xx,yy)sub('^MM[0-9]+',xx,yy),meta2$id,meta2$Time.Points))
meta2$Time.Points<-tmp
rownames(meta2)<-meta2$Time.Points

meta<-rbind(meta1,meta2)
meta$mm<-meta$id
meta$Date[meta$Date==' 12/07/2001'&meta$id=='MM14']<-'7/12/2001'
meta$Date[meta$Date=='07/0806']<-'07/08/2006'
meta$Date[meta$Date==' 08/21/08']<-'08/21/2008'
meta$rDate<-as_date(sapply(meta$Date,function(xx)ifelse(grepl('/[0-9]{4}',xx),mdy(xx),dmy(xx))))
meta$vl<-as.numeric(gsub('[><,]','',meta$VL))
meta$cd4<-as.numeric(meta$CD4)
meta$visit<-sapply(strsplit(rownames(meta),'\\.'),'[',2)


##CORRECTIONS and IFN###

wb <- loadWorkbook("meta/EJ MM plasma cytokine data CORRECTED updated VL CD4 Jan2018.xlsx")

rawMetas<-lapply(getSheets(wb),function(sheet){
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
  dat[dat=='BDL']<- 1
  dat<-as.data.frame(dat,stringsAsFactors=FALSE)
  dat[,grepl('^ifn',colnames(dat))]<- apply(dat[,grepl('^ifn',colnames(dat))],2,as.numeric)
  dat$dfosx<-as.numeric(dat$dfosx)
  if(any(colnames(dat)=='newDate')){
    if(any(!is.na(dat$newDate)&(abs(as.numeric(dat$newDate)-as.numeric(dat$oldDate))>10)))stop('Big difference in old and new date')
    dat[!is.na(dat$newDate),'date']<-dat[!is.na(dat$newDate),'newDate']
  }
  if(isStringDate)dat$rDate<-dmy(dat$date)
  else dat$rDate<-as.Date(as.numeric(dat$date),origin='1899-12-30')
  dat$cd4<-as.numeric(ifelse(dat$cd4 %in% c('not done','no data','NA'),NA,dat$cd4))
  dat$vl<-as.numeric(ifelse(dat$viralLoad %in% c('not done','no data'),NA,gsub('[><,]','',dat$viralLoad)))
  if(any(colnames(dat)=='oldViralLoad')){
    oldVl<-as.numeric(gsub('[<>,]','',dat$oldViralLoad))
    if(any(!is.na(dat$vl)&(abs(log2(dat$vl/oldVl))>1)))stop('Big difference in old and new vl')
    dat[is.na(dat$vl),'vl']<-oldVl[is.na(dat$vl)]
  }
  return(dat)
})
names(rawMetas)<-names(getSheets(wb))

colCounts<-table(unlist(sapply(rawMetas,colnames)))
targetCols<-names(colCounts)[colCounts==length(rawMetas)]
ifnMeta<-do.call(rbind,lapply(rawMetas,'[',targetCols[orderIn(targetCols,colnames(rawMetas[[1]]))]))
ifnMeta$ej<-sapply(strsplit(rownames(ifnMeta),' '),'[',1)
ifnMeta$mm<-sapply(strsplit(rownames(ifnMeta),'[ .]'),'[',2)
ifnMeta$sample<-sub('[oO]','0',ifnMeta$sample)
isVisitSample<-!is.na(ifnMeta$sample)&grepl('\\.',ifnMeta$sample)&!grepl('[a-zA-Z]$',ifnMeta$sample)
ifnMeta$oldSample<-ifnMeta$sample
ifnMeta$sample[isVisitSample]<-fixDecimals(ifnMeta$sample[isVisitSample])
ifnMeta$visit<-sapply(strsplit(ifnMeta$sample,'[.]'),'[',2)
ifnMeta<-ifnMeta[ifnMeta$mm %in% meta$mm,]
#fix inconsistent date formatting
ifnMeta[ifnMeta$ej=='EJ52'&ifnMeta$rDate=='2001-06-09','rDate']<-ymd('2001-09-06')
tmp<-unique(ifnMeta[,c('mm','ej')])
ejLookup<-structure(tmp$ej,.Names=tmp$mm)
mmLookup<-structure(tmp$mm,.Names=tmp$e)



## cell sorting ##
sorts<-read.csv('meta/AFM MM data summary Jan2018.csv',stringsAsFactors=FALSE)
colnames(sorts)[1]<-'patient'
sorts$patient<-fillDown(sorts$patient)
newColnames<-sub('\\.\\.3\\.replicate\\.values\\.|\\.\\.\\..+cells\\.','',fillDown(ifelse(grepl('^X\\.[0-9]+$',colnames(sorts)),NA,colnames(sorts))))
newColnames<-sprintf('%s%s',newColnames,ave(newColnames,newColnames,FUN=function(xx)if(length(xx)==1)'' else sprintf('__%d',1:length(xx))))
colnames(sorts)<-newColnames
sorts$ej<-sub(' ','',sapply(strsplit(sorts$Donor,'[.]'),'[',1))
paste(sorts$ej,sorts$DFOSx) %in% paste(meta$ej,meta$dfosx)
sorts$rDate<-as_date(sapply(sorts$Visit.Date,function(xx)ifelse(grepl('/',xx),mdy(xx),dmy(xx))))

## Inflammation markers ##
trans<-read.csv('meta/London Cohort inflammation markers ELISA data for Scott 07172017.csv',stringsAsFactors=FALSE)
trans$sample<-fixDecimals(as.character(trans$Sample))

## messy all combined list
if(!exists('ejs'))source('readAllPats.R')

## Post ART data
pbmc<-read.csv('meta/EJ post ART PBMC available.csv',header=FALSE,stringsAsFactors=FALSE)
pbmc[,1]<-trimws(pbmc[,1])
pbmc<-pbmc[grepl('^E?J?[0-9]+\\.[0-9]+$',pbmc[,1]),]
colnames(pbmc)<-c('sample','date','DFOSx','viralLoad','CD4','vials','postArt')
pbmc$ej<-sprintf('EJ%s',sub('EJ','',sub('\\.[0-9]+$','',pbmc$sample)))
pbmc$visit<-sub('.*\\.([0-9]+)$','\\1',pbmc$sample)
pbmc$vl<-as.numeric(sub('<','',pbmc$viralLoad))
pbmc$cd4<-as.numeric(sub('N/A','',pbmc$CD4))

## Additional data
more<-read.csv('meta/moreMetaData.csv',stringsAsFactors=FALSE)
more$vl<-as.numeric(sub('<','',more$VL))
more$cd4<-as.numeric(more$CD4)

art<-read.csv('data/artDates.csv',stringsAsFactors=FALSE)
artDates<-withAs(xx=art[!is.na(art$date)&art$mm %in% meta$mm,],structure(dmy(xx$date),.Names=xx$mm))
#art$lastDate<-ymd(apply(art[,c('lastClinic','lastSample')],1,function(xx)if(all(is.na(xx)))return(NA)else as.character(max(dmy(xx),na.rm=TRUE))))
art$lastDate<-dmy(art$lastSample)
lastDates<-withAs(xx=art[!is.na(art$lastDate)&art$mm %in% meta$mm,],structure(xx$lastDate,.Names=xx$mm))

## Joining ##

##combine ifnMeta and meta
#newMeta<-!paste(ifnMeta$rDate,ifnMeta$mm) %in% paste(meta$rDate,meta$mm)
#minDiff<-apply(ifnMeta[newMeta,c('mm','rDate')],1,function(xx)min(c(Inf,abs(meta[meta$mm==xx[1],'rDate']-ymd(xx[2])))))
#checked manually look all distinct
#ifnMeta[newMeta,][minDiff<10,]
newMeta<-!paste(meta$rDate,meta$mm) %in% paste(ifnMeta$rDate,ifnMeta$mm)
minDiff<-apply(meta[newMeta,c('mm','rDate')],1,function(xx)min(c(Inf,abs(ifnMeta[ifnMeta$mm==xx[1],'rDate']-ymd(xx[2])))))
#checked manually look distinct
meta[newMeta,][minDiff<10,]
meta$ej<-ejLookup[meta$mm]
metaMerge<-meta
metaMerge[,colnames(ifnMeta)[!colnames(ifnMeta) %in% colnames(meta)]]<-NA
metaMerge[,c('sample','date','dfosx','viralLoad')]<-metaMerge[,c('Time.Points','Date','DFOSx','VL')]
metaMerge$source<-'meta'

ifnMeta$source<-'ifn'
comboMeta<-rbind(ifnMeta,metaMerge[newMeta,colnames(ifnMeta)])

#combine ejs
thisEjs<-ejs[ejs$ej %in% comboMeta$ej&!grepl('[a-z]',ejs$date),]
thisEjs$rDate<-as_date(sapply(thisEjs$date,function(xx)ifelse(grepl('/[0-9]{4}',xx),mdy(xx),dmy(xx))))
newEj<-!paste(thisEjs$ej,thisEjs$rDate) %in% paste(comboMeta$ej,comboMeta$rDate)&!grepl('HAART',thisEjs$notes)
thisEjs<-thisEjs[newEj,]
thisEjs[,c('viralLoad','sample')]<-thisEjs[,c('vl','id')]
thisEjs$visit<-trimws(sapply(strsplit(thisEjs$sample,'\\.'),'[',2))
thisEjs$mm<-sapply(thisEjs$ej,function(xx)names(ejLookup)[ejLookup==xx])
thisEjs<-thisEjs[!thisEjs$id %in% c('EJ 85.14','EJ85.11'),]
minDiff<-apply(thisEjs[,c('mm','rDate')],1,function(xx)min(c(Inf,abs(comboMeta[comboMeta$mm==xx[1],'rDate']-ymd(xx[2])))))
if(any(minDiff<10))stop('Close date in ejs')
thisEjs[,colnames(comboMeta)[!colnames(comboMeta) %in% colnames(thisEjs)]]<-NA
thisEjs$source<-'ej'

comboMeta<-rbind(comboMeta,thisEjs[,colnames(comboMeta)])

#combine pbmc
thisPbmc<-pbmc[pbmc$ej %in% comboMeta$ej & !paste(pbmc$ej,pbmc$visit) %in% paste(comboMeta$ej,comboMeta$visit),]
thisPbmc$rDate<-dmy(thisPbmc$date)
thisPbmc$mm<-sapply(thisPbmc$ej,function(xx)names(ejLookup)[ejLookup==xx])
thisPbmc[,colnames(comboMeta)[!colnames(comboMeta) %in% colnames(thisPbmc)]]<-NA
thisPbmc$source<-'pbmc'

comboMeta<-rbind(comboMeta,thisPbmc[,colnames(comboMeta)])

##combine additional data
thisMore<-more
thisMore$ej<-ejLookup[thisMore$mm]
thisMore$rDate<-dmy(thisMore$date)
thisMore[,colnames(comboMeta)[!colnames(comboMeta) %in% colnames(thisMore)]]<-NA
thisMore$source<-'additional'

comboMeta<-rbind(comboMeta,thisMore[,colnames(comboMeta)])

##combine trans
if(any(!trans$sample %in% sprintf('%s.%s',sub('EJ','',comboMeta$ej),comboMeta$visit)))stop('Found unknown sample in trans data')
rownames(trans)<-trans$sample
transCols<-colnames(trans)[!colnames(trans) %in% c('sample','Sample')]
if(any(transCols %in% colnames(comboMeta)))stop('Duplicate column in trans')
comboMeta[,transCols]<-trans[sprintf('%s.%s',sub('EJ','',comboMeta$ej),comboMeta$visit),transCols]

##combine sort
#first two samples are controls
if(any(!paste(sorts$ej,sorts$rDate)[-1:-2] %in% paste(comboMeta$ej,comboMeta$rDate)))stop('Unknown sample in sorts')
sortCols<-colnames(sorts)[grepl('BST|HLA|CD38|__',colnames(sorts))]
if(any(sortCols %in% colnames(comboMeta)))stop('Duplicate column in trans')
rownames(sorts)<-paste(sorts$ej,sorts$rDate)
comboMeta[,sortCols]<-sorts[paste(comboMeta$ej,comboMeta$rDate),sortCols]
comboMeta<-comboMeta[order(comboMeta$mm,comboMeta$rDate),]
comboMeta$dfosx<-as.numeric(comboMeta$dfosx)
comboMeta$qvoa<-comboMeta$rDate>as_date(ifelse(comboMeta$mm %in% names(artDates),artDates[comboMeta$mm],Inf))


sapply(by(comboMeta[,c('dfosx','rDate')],comboMeta$mm,function(xx){zz<-table(xx$rDate-xx$dfosx)}),function(xx)diff(range(ymd(names(xx)))))
baseDate<-by(comboMeta[,c('dfosx','rDate')],comboMeta$mm,function(xx){zz<-table(xx$rDate-xx$dfosx);names(zz)[which.max(zz)]})
comboMeta$time<-comboMeta$rDate-ymd(baseDate[comboMeta$mm])
comboMeta[comboMeta$visit=='12 MW'&comboMeta$mm=='MM39','visit']<-'13'

comboMeta[comboMeta$vl==37611600&!is.na(comboMeta$vl),'vl']<-NA

if(any(apply(table(comboMeta$visit,comboMeta$mm)>1,2,any)))stop('Duplicate visit found')
write.csv(comboMeta,'out/combinedMeta.csv')
tmp<-comboMeta[,c('mm','ej','date','rDate','vl','cd4','source')]
tmp$dfosx<-comboMeta$time
write.csv(tmp,'out/combinedMetadata.csv',row.names=FALSE)

artDfosx<-sapply(names(artDates),function(xx)artDates[xx]-ymd(baseDate[xx]))
names(artDfosx)<-names(artDates)
lastDfosx<-sapply(names(lastDates),function(xx)lastDates[xx]-ymd(baseDate[xx]))
names(lastDfosx)<-names(lastDates)
for(ii in names(lastDfosx))lastDfosx[ii]<-max(as.numeric(comboMeta[comboMeta$mm==ii,'time']),lastDfosx[ii])

customCols<-read.csv('data/Hex color no. for MM cohort colorcode.csv',stringsAsFactors=FALSE,header=FALSE)[,1:2]
customCols<-customCols[customCols[,1]!='',]
colnames(customCols)<-c('sample','color')
customCols$name<-fixDecimals(sub(' ?\\(.*$','',customCols$sample))
rownames(customCols)<-customCols$name



wb <- loadWorkbook("meta/EJ MM CD4 VL pre and post ART 08June2018_sasm.xlsx")
vals<-lapply(getSheets(wb),function(sheet){
  rows<-getCells(getRows(sheet),simplify=FALSE)
  vals<-lapply(rows,function(row){
    tmp<-lapply(as.character(1:8),function(xx)ifelse(any(names(row)==xx),getCellValue(row[[xx]]),NA))
    if(is.na(tmp[[2]])&is.na(tmp[[3]]))return(NULL)
    if((grepl('Date',tmp[[2]])|grepl('Date',tmp[[3]])))return(NULL)
    out<-data.frame('id'='999.99','origDate'='99.99.99','date'=99999,'DFOSx'=99999,'VL'=999999999,'CD4'=9999999,'ART'='','Notes'='',stringsAsFactors=FALSE)[0,]
    out[1,]<-rep(NA,8)
    for(ii in 1:8)if(length(tmp)>=ii)out[1,ii]<-tmp[[ii]] else out[1,ii]<-NA
    return(out)
  })
  return(do.call(rbind,vals))
})
compiledMeta<-do.call(rbind,mapply(function(xx,yy){xx$pat<-yy;xx},vals,names(vals),SIMPLIFY=FALSE))
compiledMeta[compiledMeta$origDate=='05.01.12'&compiledMeta$id=='108.1','id']<-'108.10'
compiledMeta[compiledMeta$id=='85.12MW'&!is.na(compiledMeta$id),'id']<-'85.13'
compiledMeta$mm<-sub('.* ','',sub('MM ','MM',compiledMeta$pat))
compiledMeta$ej<-sub(' .*','',sub('EJ ','EJ',compiledMeta$pat))
compiledMeta<-compiledMeta[compiledMeta$mm %in% mmLookup,]
compiledMeta$rDate<-as.Date(as.numeric(compiledMeta$date),origin='1899-12-30')
compiledMeta$vl<-as.numeric(gsub(' ','',sub('<','',compiledMeta$VL)))
compiledMeta$cd4<-as.numeric(compiledMeta$CD4)
rownames(compiledMeta)<-sapply(strsplit(sub('^[^ ]+ ','',rownames(compiledMeta)),'\\.'),function(xx)sprintf('%s.%02d',xx[1],as.numeric(xx[2])))
if(any(is.na(compiledMeta$rDate)))stop('Problem interpreting date')
if(year(min(compiledMeta$rDate))<2000)stop('Year <2000 detected')
if(year(min(compiledMeta$rDate))>2015)stop('Year >2015 detected')
startDates<-tapply(compiledMeta$rDate-compiledMeta$DFOSx,compiledMeta$mm,mostAbundant)
compiledMeta$time<-compiledMeta$rDate-as.Date(startDates[compiledMeta$mm])
if(any(abs(compiledMeta$time-compiledMeta$DFOSx)>1))warning('Disagreement in dfosx')


#2nd column likely gives day from exposure
weauSymptomDate<-ymd('1990-06-04')
weauMeta<-read.csv('meta/weau.csv',stringsAsFactors=FALSE)
weauMeta$origDate<-weauMeta$date<-weauMeta$Date
weauMeta$rDate<-dmy(weauMeta$Date)
weauMeta$ID<-weauMeta$id<-weauMeta$visit<-1:nrow(weauMeta)
weauMeta$time<-weauMeta$rDate-weauSymptomDate
weauAdditional<-read.csv('meta/additionalWEAUMeta.csv')
weauAdditional$origDate<-weauAdditional$date<-weauAdditional$Date
weauAdditional$rDate<-mdy(weauAdditional$Date)
weauAdditional$ID<-weauAdditional$id<-weauAdditional$visit<-nrow(weauMeta)+1:nrow(weauAdditional)
weauAdditional$VL<-NA
weauAdditional$Time<-NA
weauAdditional$Available<-NA
weauAdditional$time<-weauAdditional$rDate-weauSymptomDate
weauAdditional<-weauAdditional[weauAdditional$time>100,]
weauMeta<-rbind(weauMeta,weauAdditional[,colnames(weauMeta)])
weauMeta$cd4<-as.numeric(ifelse(weauMeta$CD4=='nd',NA,weauMeta$CD4))
weauMeta$vl<-as.numeric(ifelse(weauMeta$VL=='nd',NA,weauMeta$VL))
weauMeta$DFOSx<-weauMeta$time
weauMeta$ART<-weauMeta$Notes<-NA
weauMeta$pat<-weauMeta$mm<-weauMeta$ej<-'WEAU'
rownames(weauMeta)<-weauMeta$Time.Points<-sprintf('WEAU.%02d',weauMeta$id)
weauMeta<-weauMeta[order(weauMeta$time),]


compiledMeta<-rbind(compiledMeta,weauMeta[,colnames(compiledMeta)])
meta<-rbind(meta,weauMeta[,colnames(meta)])
compiledMeta$visit<-sub('[^.]+\\.','',compiledMeta$id)
compiledMeta$visit<-ifelse(grepl('^[0-9]+$',compiledMeta$visit),sprintf('%02d',suppressWarnings(as.integer(compiledMeta$visit))),compiledMeta$visit)
compiledMeta$sample<-ifelse(compiledMeta$visit==''|is.na(compiledMeta$id),sprintf('XX%s',1:nrow(compiledMeta)),paste(compiledMeta$mm,compiledMeta$visit,sep='.'))
rownames(compiledMeta)<-compiledMeta$sample

#WEAU no ART but calling first record of low CD4 as when would have initiated
meta$artDay<-c(artDfosx,'WEAU'=391)[meta$mm]
meta$daysBeforeArt<-meta$artDay-as.numeric(meta$DFOSx)
compiledMeta$artDay<-c(artDfosx,'WEAU'=391)[compiledMeta$mm]
compiledMeta$daysBeforeArt<-compiledMeta$artDay-as.numeric(compiledMeta$DFOSx)

founders<-read.csv('founder.csv',stringsAsFactors=FALSE,row.names=1)
superDate<-ymd(founders$superDate)
founders$superTime<-superDate-ymd(startDates[rownames(founders)])

less350Time<-by(compiledMeta[!is.na(compiledMeta$cd4),],compiledMeta[!is.na(compiledMeta$cd4),'mm'],function(xx){
  lastInfect<-ifelse(is.na(founders[xx$mm[1],'superTime']),1,founders[xx$mm[1],'superTime'])
  xx$previousLess<-c(Inf,xx$cd4[-nrow(xx)])<350
  out<-min(c(xx[xx$time>lastInfect+180&xx$cd4<350&xx$previousLess,'time'],Inf))
  if(out==Inf)out<-NA
  if(is.na(out)&!is.na(artDfosx[xx$mm[1]]))out<-artDfosx[xx$mm[1]]
  return(out)
})
compiledMeta$day350<-less350Time[compiledMeta$mm]
compiledMeta$daysBefore350<-compiledMeta$day350-as.numeric(compiledMeta$DFOSx)

if(FALSE){
comboMeta[which(!paste(comboMeta$mm,comboMeta$rDate) %in% paste(compiledMeta$mm,compiledMeta$rDate) & !is.na(comboMeta$mm)&(!is.na(comboMeta$vl)|!is.na(comboMeta$cd4))),c('mm','date','rDate','time','vl','cd4','source')]
tmp<-comboMeta$vl
names(tmp)<-paste(comboMeta$mm,comboMeta$rDate)
tmp<-tmp[paste(compiledMeta$mm,compiledMeta$rDate)]
probs<-tmp!=sub('<','',compiledMeta$VL)&!is.na(tmp)
cbind(compiledMeta[probs,],tmp[probs])
tmp<-comboMeta$cd4
names(tmp)<-paste(comboMeta$mm,comboMeta$rDate)
tmp<-tmp[paste(compiledMeta$mm,compiledMeta$rDate)]
probs<-tmp!=compiledMeta$CD4&!is.na(tmp)
cbind(compiledMeta[probs,],tmp[probs])
}
