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
