
plateLookup<-read.csv('ice/2020-02-01_plateIds.csv',stringsAsFactors=FALSE)
plateLookup$type<-ifelse(grepl('[A-Z]',plateLookup$plate),'plasma','PBMC')

spots<-read.csv('ice/2020-01-31_positiveSpots.csv',stringsAsFactors=FALSE)
convert<-convert96To48(spots$Well)
out<-cbind(spots,convert)
out$plate48<-mapply('[',strsplit(out$Plate,'/'),out$plate)
out$well48<-paste(out$row,out$col,sep='')
colnames(out)[colnames(out)=='Well']<-'well96'
colnames(out)[colnames(out)=='Plate']<-'plate96'
write.csv(out[,c('well48','plate48','Note','well96','plate96')],'out/2020-01-31_positiveSpotsConverted.csv')
out$date<-'2020/01/31'
out$animal<-apply(out[,c('plate48','row')],1,function(xx){out<-plateLookup[plateLookup$plate==xx[1]&grepl(xx[2],plateLookup$row),'animal'];ifelse(length(out)==1,out,NA)})
out$day<-apply(out[,c('plate48','row')],1,function(xx){out<-plateLookup[plateLookup$plate==xx[1]&grepl(xx[2],plateLookup$row),'day'];ifelse(length(out)==1,out,NA)})
out$type<-apply(out[,c('plate48','row')],1,function(xx){out<-plateLookup[plateLookup$plate==xx[1]&grepl(xx[2],plateLookup$row),'type'];ifelse(length(out)==1,out,NA)})
write.csv(out,'out/2020-01-31_positiveSpotsLabeled.csv')
bak<-comboSpots

spots2<-read.csv('ice/2020-02-03_positiveSpots.csv',stringsAsFactors=FALSE)
convert2<-convert96To48(spots2$Well)
out2<-cbind(spots2,convert2)
out2$plate48<-mapply('[',strsplit(out2$Plate,'/'),out2$plate)
out2$well48<-paste(out2$row,out2$col,sep='')
colnames(out2)[colnames(out2)=='Well']<-'well96'
colnames(out2)[colnames(out2)=='Plate']<-'plate96'
out2$in131<-paste(out2$well48,out2$plate48) %in% paste(out$well48,out$plate48)
write.csv(out2[,c('well48','plate48','Note','well96','plate96','in131')],'out/2020-02-04_positiveSpotsConverted.csv')
out2$date<-'2020/02/03'

spots3<-read.csv('ice/2020-01-24_positiveSpots.csv',stringsAsFactors=FALSE)
convert3<-convert96To48(spots3$Well)
out3<-cbind(spots3,convert3)
out3$plate48<-mapply('[',strsplit(out3$Plate,'/'),out3$plate)
out3$well48<-paste(out3$row,out3$col,sep='')
colnames(out3)[colnames(out3)=='Well']<-'well96'
colnames(out3)[colnames(out3)=='Plate']<-'plate96'
out3$date<-'2020/01/24'

spots4<-read.csv('ice/2020-02-10_positiveSpots.csv',stringsAsFactors=FALSE)
spots4$bak<-spots4$Well
spots4$Well[spots4$Plate=='9-10/Repeat9-10']<-sub('[78]','10',sub('(9|10)','11',sub('(11|12)','12',spots4$Well[spots4$Plate=='9-10/Repeat9-10'])))
convert4<-convert96To48(spots4$Well)
out4<-cbind(spots4,convert4)
out4$plate48<-mapply('[',strsplit(out4$Plate,'/'),out4$plate)
out4$well48<-paste(out4$row,out4$col,sep='')
colnames(out4)[colnames(out4)=='Well']<-'well96'
colnames(out4)[colnames(out4)=='Plate']<-'plate96'
out4[out4$plate48=='Repeat9-10','plate48']<-'9-10'
out4$date<-'2020/02/10'
out4$previous<-paste(out4$plate48,out4$well48) %in% paste(c(out$plate48,out2$plate48,out3$plate48),c(out$well48,out2$well48,out3$well48))
out4<-out4[order(out4$plate48,out4$row,out4$col),]
out4$num<-NA
out4$num[!out4$previous]<-80+1:sum(!out4$previous)
write.csv(out4[,c('num','plate48','well48','previous')],'out/positiveSpots_2020-02-10.csv')

spots5<-read.csv('ice/2020-02-17_positiveSpots.csv',stringsAsFactors=FALSE)
spots5$bak<-spots5$Well
convert5<-convert96To48(spots5$Well)
out5<-cbind(spots5,convert5)
out5$plate48<-mapply('[',strsplit(out5$Plate,'/'),out5$plate)
out5$well48<-paste(out5$row,out5$col,sep='')
colnames(out5)[colnames(out5)=='Well']<-'well96'
colnames(out5)[colnames(out5)=='Plate']<-'plate96'
out5$date<-'2020/02/17'
out5$previous<-paste(out5$plate48,out5$well48) %in% paste(c(out$plate48,out2$plate48,out3$plate48,out4$plate48),c(out$well48,out2$well48,out3$well48,out4$plate48))
out5$num<-NA
out5$num[!out5$previous]<-100+1:sum(!out5$previous)
write.csv(out5[,c('num','plate48','well48','previous')],'out/positiveSpots_2020-02-17.csv')


#order is important here. out wasn't filtered for previous out3 but rest were
comboSpots<-rbind(out[,colnames(out3)],out3,out2[,colnames(out3)],out4[,colnames(out3)],out5[,colnames(out3)])
comboSpots$animal<-apply(comboSpots[,c('plate48','row')],1,function(xx){out<-plateLookup[plateLookup$plate==xx[1]&grepl(xx[2],plateLookup$row),'animal'];ifelse(length(out)==1,out,NA)})
comboSpots$day<-apply(comboSpots[,c('plate48','row')],1,function(xx){out<-plateLookup[plateLookup$plate==xx[1]&grepl(xx[2],plateLookup$row),'day'];ifelse(length(out)==1,out,NA)})
comboSpots$type<-apply(comboSpots[,c('plate48','row')],1,function(xx){out<-plateLookup[plateLookup$plate==xx[1]&grepl(xx[2],plateLookup$row),'type'];ifelse(length(out)==1,out,NA)})
comboSpots$virus<-apply(comboSpots[,c('plate48','row')],1,function(xx){out<-plateLookup[plateLookup$plate==xx[1]&grepl(xx[2],plateLookup$row),'virus'];ifelse(length(out)==1,out,NA)})
comboSpots$vl<-apply(comboSpots[,c('plate48','row')],1,function(xx){out<-plateLookup[plateLookup$plate==xx[1]&grepl(xx[2],plateLookup$row),'vl'];ifelse(length(out)==1,out,NA)})
bak<-comboSpots
comboSpots<-comboSpots[!duplicated(paste(comboSpots$plate48,comboSpots$well48)),]
comboSpots$flaskId
comboSpots[comboSpots$date=='2020/01/31','flaskId']<-1:48
if(comboSpots[comboSpots$well48=='E5'&comboSpots$plate48=='5-6','flaskId']!=22)stop('Flask ID spot check failed')
if(comboSpots[comboSpots$well48=='A7'&comboSpots$plate48=='C','flaskId']!=39)stop('Flask ID spot check failed')
if(comboSpots[comboSpots$well48=='F7'&comboSpots$plate48=='D','flaskId']!=48)stop('Flask ID spot check failed')
comboSpots[comboSpots$date=='2020/02/03','flaskId']<-c(49:76,78:79)
if(comboSpots[comboSpots$well48=='D3'&comboSpots$plate48=='9-10','flaskId']!=57)stop('Flask ID spot check failed')
if(comboSpots[comboSpots$well48=='A1'&comboSpots$plate48=='C','flaskId']!=65)stop('Flask ID spot check failed')
if(comboSpots[comboSpots$well48=='E4'&comboSpots$plate48=='B','flaskId']!=79)stop('Flask ID spot check failed')
comboSpots[comboSpots$date=='2020/02/10','flaskId']<-81:93
comboSpots[comboSpots$date=='2020/02/10'&comboSpots$well48=='D4'&comboSpots$plate48=='9-10','flaskId']<-82
comboSpots[comboSpots$date=='2020/02/10'&comboSpots$well48=='D2'&comboSpots$plate48=='9-10','flaskId']<-83
if(comboSpots[comboSpots$well48=='D4'&comboSpots$plate48=='9-10','flaskId']!=82)stop('Flask ID spot check failed')
if(comboSpots[comboSpots$well48=='F4'&comboSpots$plate48=='9-10','flaskId']!=91)stop('Flask ID spot check failed')
if(comboSpots[comboSpots$well48=='C5'&comboSpots$plate48=='A','flaskId']!=93)stop('Flask ID spot check failed')
if(comboSpots[comboSpots$well48=='D5'&comboSpots$plate48=='3-4','flaskId']!=81)stop('Flask ID spot check failed')
comboSpots$name<-sprintf('%s.d%s.%s.%s',comboSpots$animal,comboSpots$day,ifelse(comboSpots$type=='plasma','PLAS','PBMC'),comboSpots$well48)
comboSpots[comboSpots$name=='RM10N011.d942.PLAS.E5'&comboSpots$plate48=='D','name']<-'RM10N011.d942.PLAS.1E5'
comboSpots[comboSpots$date=='2020/02/17','flaskId']<-100:(100+sum(comboSpots$date=='2020/02/17')-1)
if(any(table(comboSpots$name)>1))stop('Duplicated name')
write.csv(comboSpots[c('flaskId','name','date','animal','day','type','virus','Note','plate48','well48','plate96','well96')],'out/combinedMacaqueIds.csv',row.names=FALSE)
rt<-read.csv('ice/RT list for macaque isolates_20200218.csv',skip=1)[,-1]
rt$corrected<-sub('6563.d1085','6563.d1065',sub('VOA','PBMC',sub('^.*_','',rt$Isolates.ID)))
rownames(rt)<-rt$corrected
if(any(!rt$corrected %in% comboSpots$name))stop('Name mismatch')
comboSpots$rt<-rt[comboSpots$name,'RT.ng.ul']

table(comboSpots$animal,comboSpots$day)
attempted<-unique(plateLookup[,c('animal','day','vl','type')])
attempted$positives<-sapply(paste(attempted$animal,attempted$day,attempted$type),function(xx)sum(xx == paste(comboSpots$animal,comboSpots$day,comboSpots$type)))
attempted$positiveMoreThan1Spot<-sapply(paste(attempted$animal,attempted$day,attempted$type),function(xx)withAs(comboSpots=comboSpots[!grepl('spot',comboSpots$Note),],sum(xx == paste(comboSpots$animal,comboSpots$day,comboSpots$type))))
attempted$confirmedPositives<-sapply(paste(attempted$animal,attempted$day,attempted$type),function(xx)withAs(comboSpots=comboSpots[!is.na(comboSpots$rt)&comboSpots$rt>0.004,],sum(xx == paste(comboSpots$animal,comboSpots$day,comboSpots$type))))
attempted$confirmedNegatives<-sapply(paste(attempted$animal,attempted$day,attempted$type),function(xx)withAs(comboSpots=comboSpots[!is.na(comboSpots$rt)&comboSpots$rt<0.004,],sum(xx == paste(comboSpots$animal,comboSpots$day,comboSpots$type))))
attempted$awaitingHarvest<-sapply(paste(attempted$animal,attempted$day,attempted$type),function(xx)withAs(comboSpots=comboSpots[is.na(comboSpots$rt)&!grepl('spot',comboSpots$Note),],sum(xx == paste(comboSpots$animal,comboSpots$day,comboSpots$type))))
attempted$awaitingHarvest1Spot<-sapply(paste(attempted$animal,attempted$day,attempted$type),function(xx)withAs(comboSpots=comboSpots[is.na(comboSpots$rt)&grepl('spot',comboSpots$Note),],sum(xx == paste(comboSpots$animal,comboSpots$day,comboSpots$type))))
attempted$positivesRetro<-sapply(paste(attempted$animal,attempted$day,attempted$type),function(xx)sum(xx == paste(comboSpots$animal,comboSpots$day,comboSpots$type)&grepl('[AB]',comboSpots$plate48)))
attempted$positivesRetroEFC<-sapply(paste(attempted$animal,attempted$day,attempted$type),function(xx)sum(xx == paste(comboSpots$animal,comboSpots$day,comboSpots$type)&grepl('[CD]',comboSpots$plate48)))
attempted<-attempted[order(attempted$animal,suppressWarnings(as.numeric(attempted$day))),]
write.csv(attempted,'out/positive_vs_attempted_macaque.csv',row.names=FALSE)
tmp<-bak[bak$date=='2020/01/24',c('animal','well48','day','type')]
tmp<-tmp[order(tmp$animal,tmp$well48),]
rownames(tmp)<-1:nrow(tmp)
write.csv(tmp,'out/positives_2020-01-24.csv',row.names=FALSE)


