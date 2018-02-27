if(!exists('dat'))source('readNewData.R')
set.seed(12345)
message('Marvin worst')
mins<-dat[dat$replication<=sort(dat[!is.na(dat$replication)&!dat$qvoa,'replication'])[5]&!is.na(dat$replication),c('id','replication')]
print(mins)

message('Marvin random')
pairs<-withAs(xx=dat[!is.na(dat$replication)&!dat$qvoa,],tapply(xx$id,xx$pat,sample,2))
dummy<-lapply(pairs,sapply,message)

message('Marvin QVOA')
dummy<-sapply(dat[dat$qvoa,'id'],message)

if(!exists('hiv'))source('~/projects/hivPair/readData.R',chdir=TRUE)

message('Shilpa worst')
mins<-hiv[hiv$Replicative.capacity.Pooled.Donor.cells.p24.d7<=sort(hiv[!is.na(hiv$Replicative.capacity.Pooled.Donor.cells.p24.d7),'Replicative.capacity.Pooled.Donor.cells.p24.d7'])[5],c('Renamed','Replicative.capacity.Pooled.Donor.cells.p24.d7','isGenital','select')]
print(mins)

message('Shilpa nongenital random')
dummy<-sapply(withAs(xx=hiv[!hiv$isGenital&hiv$select=='UT',],tapply(xx$Renamed,xx$baseName,sample,1)),message)

message('Shilpa genital random')
dummy<-sapply(withAs(xx=hiv[hiv$isGenital&hiv$select=='UT',],tapply(xx$Renamed,xx$baseName,sample,2)),sapply,message)


set.seed(123456)
message('Marvin ice')
excludeTimes<-list('MM40'=c(),'MM23'=c(64),'MM34'=c(74))
num<-c('MM40'=2,'MM23'=1,'MM34'=1)
allPicks<-list()
for(ii in names(excludeTimes)){
  thisDat<-dat[dat$pat==ii&!is.na(dat$rt)&!dat$qvoa,]
  uniqTime<-unique(thisDat$time)
  uniqTime<-uniqTime[!uniqTime %in% excludeTimes[[ii]]]
  rtThresh<-quantile(thisDat$rt,.25)
  message('RT threshold for ',ii,': ',rtThresh)
  picks<-lapply(uniqTime,function(jj){
    thisTime<-thisDat[thisDat$time==jj,]
    gt005<-thisTime$rt>rtThresh
    message(sum(gt005),' viruses with rt>',rtThresh,' out of ',nrow(thisTime),' for time ',jj)
    if(any(gt005))sample(thisTime$ID.for.Publications[gt005],num[ii])
    else(tail(thisTime$ID.for.Publications[order(thisTime$rt)],num[ii]))
  })
  message(paste(unlist(picks),collapse=', '))
  allPicks<-c(allPicks,list(picks))
}
names(allPicks)<-names(excludeTimes)
message(paste(unlist(mapply(head,allPicks[['MM40']],c(2,1,2,1,2))),collapse=', '))
