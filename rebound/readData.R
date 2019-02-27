library(vipor)
library(dnar)
dat<-read.csv('A09 and A06 data tables for scott.csv',stringsAsFactors=FALSE)
colnames(dat)[colnames(dat)=='Virus']<-'virus'
dat$pat<-ifelse(grepl('621803B',dat$virus),'A06',
  ifelse(grepl('621818B|A09-|Rebound ',dat$virus),'A09',
    ifelse(grepl('^B[0-9]+_',dat$virus),sub('_.*','',dat$virus),'')))
dat$type<-ifelse(grepl('pre',dat$virus),'Pre-ATI',
  ifelse(grepl('post',dat$virus),'Post-ATI', 
    ifelse(grepl('Rebound|A09-',dat$virus),'Rebound',
      ifelse(grepl('^B[0-9]+_',dat$virus),'Lorenzi et al.',
        ifelse(grepl('UK61|MM34',dat$virus),'Control',NA)))))
qvoaGroups<-c('Pre-ATI','Post-ATI','Lorenzi et al.')
dat$class<-ifelse(dat$type %in% qvoaGroups,'QVOA',dat$type)
dat$label<-sprintf('%s%s%s',dat$type,ifelse(dat$pat!='',' ',''),dat$pat)
dat$ic50<-apply(dat[,grepl('IFNA2',colnames(dat))],1,mean)
dat$source<-'stephanie'

#Add Shilpa acute vs chronic, add Marvin nadir/first, add Marvin QVOA
mm<-read.csv('firstNadir.csv',stringsAsFactors=FALSE)
mm<-mm[mm$isFirst|mm$qvoa|mm$isSix|mm$isLast|mm$isNadir,]
mm$class<-ifelse(mm$qvoa,'QVOA',ifelse(mm$isNadir,'Nadir',ifelse(mm$isSix,'6 Month',ifelse(mm$isLast,'Last','Acute'))))
mm$type<-'MM cohort'
#mm$label<-sprintf('MM %s',ifelse(mm$class=='QVOA',sprintf('Outgrowth %s',mm$pat),mm$class))
mm$label<-sprintf('%s',ifelse(mm$class=='QVOA',sprintf('Outgrowth %s',mm$pat),mm$class))
mm<-mm[!is.na(mm$ic50),]
mm[,colnames(dat)[!colnames(dat) %in% colnames(mm)]]<-NA
mm$source<-'marvin'

#write.csv(hiv[hiv$select=='UT'&hiv$fluid=='PL',c('sample','donor','IFNa2.Pooled.Donor.cells.IC50..pg..ml','IFNbeta.Pooled.Donor.cells.IC50..pg.ml')],'out/donorRecipient.csv')
pair<-read.csv('donorRecipient.csv',stringsAsFactors=FALSE)
pair$class<-ifelse(pair$donor,'Donor','Recipient')
pair$type<-'CHAVI cohort'
#pair$label<-sprintf('CHAVI %s',pair$class)
pair$label<-sprintf('%s %s',ifelse(pair$class=='Donor','Chronic','Acute'),pair$class)
pair$ic50<-pair$IFNa2.Pooled.Donor.cells.IC50..pg..ml
pair[pair$class=='Recipient','class']<-'Acute'
pair$pat<-pair$sample
pair[,colnames(dat)[!colnames(dat) %in% colnames(pair)]]<-NA
pair$source<-'shilpa'

mont<-read.csv('../out/montaner_ic50.csv',row.names=1)
mont<-mont[grepl('Alpha|alpha',rownames(mont))&!grepl('UK61',rownames(mont)),]
mont$class<-'Rebound'
mont$virus<-sub(' \\([^)]+\\) - ?','-',sub(' \\(Alpha.*','',rownames(mont)))
mont$type<-sprintf('Montaner %s',sapply(strsplit(rownames(mont),'[()]'),'[',2))
mont$pat<-sub(' .*','',rownames(mont))
mont$label<-sprintf('Rebound %s',mont$pat)
mont$source<-'Montaner'
mont[,colnames(dat)[!colnames(dat) %in% colnames(mont)]]<-NA

newRebound<-rbind(
  read.csv('../out/reboundIc50.csv',row.names=1),
  read.csv('../out/rebound_2019-01-30_Ic50.csv',row.names=1),
  read.csv('../out/voa_2019-01-10.csv',row.names=1)
)
newRebound<-newRebound[grepl('Alpha|alpha',rownames(newRebound))&!grepl('UK61|SG3|CH40|WEAU',rownames(newRebound)),]
newRebound$virus<-sub(' \\([^)]+\\) - ?','-',sub(' \\(plate 1 - Alpha.*','',rownames(newRebound)))
newRebound$class<-ifelse(grepl('Rebound|S-22|BEAT',newRebound$virus),'Rebound','QVOA')
newRebound$pat<-sub('[ _.].*','',sub('QVOA ','',sub('Rebound ','',rownames(newRebound))))
newRebound$type<-ifelse(grepl('Rebound|S-22|BEAT',newRebound$virus),'Rebound',ifelse(grepl('MM',newRebound$virus),'Outgrowth','Pre-ATI'))
newRebound$label<-sprintf('%s %s',newRebound$type,newRebound$pat)
newRebound$source<-'newRebound'
newRebound[,colnames(dat)[!colnames(dat) %in% colnames(newRebound)]]<-NA
newRebound<-newRebound[!grepl('MM[0-9]+\\.',rownames(newRebound)),]

combo<-rbind(dat,mm[,colnames(dat)],pair[,colnames(dat)],mont[,colnames(dat)],newRebound[,colnames(dat)])
combo<-combo[combo$class!='Control',]
combo$study<-ifelse(grepl('^A[0-9]+$',combo$pat),'VRC01',
  ifelse(grepl('^B[0-9]+$',combo$pat),'Reservoir',
    ifelse(grepl('^BEAT-',combo$pat),'BEAT',
      ifelse(grepl('^92[0-9][0-9]|60[0-9]',combo$pat),'3BNC117/10-1074',
        ifelse(grepl('^MM[0-9]+|WEAU',combo$pat),'MM',
          ifelse(grepl('^S-[0-9]+',combo$pat),'ATI',
            ifelse(grepl('^Donor|Recipient',combo$pat),'Transmission',
                'UNKNOWN'
)))))))
if(any(combo$study=='UNKNOWN'))stop('Unknown study')

write.csv(combo[!combo$source %in% c('marvin','shilpa'),c('virus','pat','type','class','study','source','ic50')],'comboAlpha.csv',row.names=FALSE)

ordering<-c('Rebound 601','Pre-ATI 9201','Rebound 9201','Rebound 9202','Pre-ATI 9203','Rebound 9203',
  'Rebound 9207','Rebound BEAT-004','Rebound BEAT-030','Rebound BEAT-044','Rebound S-22','Rebound S-30','Pre-ATI A06','Post-ATI A06','Pre-ATI A09','Rebound A09','Post-ATI A09','Rebound A08','Lorenzi et al. B106','Lorenzi et al. B199','Outgrowth MM14','Outgrowth MM15','Outgrowth MM23','Outgrowth MM34','Outgrowth MM40','Acute','6 Month','Nadir','Last','Acute Recipient','Chronic Donor')
#typeOrder<-c('Pre-ATI','Post-ATI','Rebound','Nussenzweig','Control','MM cohort','CHAVI cohort')
#classOrder<-c('QVOA','Rebound','Acute','Donor','Nadir')
pos<-structure(1:length(unique(combo$label)),.Names=unique(combo$label[orderIn(combo$label,ordering)]))
deemphasize<-c('Acute','6 Month','Donor','Nadir','Last')
uniqClass<-unique(combo$class)
uniqClass<-c(uniqClass[uniqClass!='Rebound'],'Rebound')
uniqClass[uniqClass %in% deemphasize]<-deemphasize[deemphasize %in% uniqClass]
classCols<-structure(rainbow.lab(length(uniqClass),alpha=ifelse(uniqClass %in% deemphasize,.7,.9),lightMultiple=ifelse(uniqClass %in% deemphasize,.65,.75)),.Names=uniqClass)
#classCols[deemphasize]<-sub('..$','88',classCols[deemphasize])

#classCols2<-structure(rainbow.lab(length(unique(combo$class)),alpha=.4),.Names=c(unique(combo$class[combo$class!='Rebound']),'Rebound'))



ranges<-do.call(rbind,tapply(combo$ic50,combo$label,range))
rangeClass<-sapply(rownames(ranges),function(xx)combo[combo$label==xx,'class'][1])


subs<-c('A06'='Patient A06','A09'='Patient A09','B106'='Patient B106','B199'='Patient B199','Recipient'='Recipients','Donor'='Donors','Lorenzi et al.'='Outgrowth','S-30'='Patient S-30','BEAT-044'='Patient BEAT-044','BEAT-030'='Patient BEAT-030','MM23'='Patient MM23','MM34'='Patient MM34','ATI'='ATI Outgrowth')
newNames<-sapply(names(pos),function(xx){for(ii in names(subs))xx<-sub(ii,subs[ii],xx);xx})
pos<-pos+cumsum(names(pos)=='Acute')*.5+cumsum(names(pos)=='Acute Recipient')*.5 #+cumsum(names(pos)=='Lorenzi et al. B106')*.5+cumsum(names(pos)=='Outgrowth MM23')*.5
subsets<-rev(list('mm'=names(pos)[!grepl('ATI|Rebound|Lorenzi',names(pos))],'lorenzi'=names(pos)[!grepl('ATI|Rebound',names(pos))],'prePost'=names(pos)[!grepl('Rebound',names(pos))],'noMont'=names(pos)[!grepl('BEAT|S-[0-9]',names(pos))],'all'=names(pos)))
pdf('qvoa_compare.pdf',width=9,height=6)
ylim=range(combo$ic50)
spread<-offsetX(log10(combo$ic50),combo$label,width=.25)
for(ii in names(subsets)){
  message(ii)
  selector<-combo$label %in% subsets[[ii]]
  if(ii=='all')marSpace<-0
  else marSpace<-marSpaces[ii]
  if(ii=='prePost')marSpace<-marSpaces['noMont']
  par(mar=c(9.8,4+marSpace,.1,3))
  plot(pos[combo$label[selector]]+spread[selector],combo$ic50[selector],log='y',yaxt='n',ylab='IFNa2 IC50 (pg/ml)',xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim)
  if(ii=='all')marSpaces<-sapply(subsets,function(xx)diff(convertUserToLine(1:0,2))*sum(!names(pos) %in% xx))
  slantAxis(1,pos[names(pos) %in% subsets[[ii]]],newNames[names(pos) %in% subsets[[ii]]],srt=-45)
  logAxis(las=1)
  #segments(pos[rownames(ranges)],ranges[,1],pos[rownames(ranges)],ranges[,2],col=classCols2[rangeClass],lwd=2)
  #,cex=ifelse(combo$class %in% deemphasize,1.5)
  points(pos[combo$label[selector]]+spread[selector],combo$ic50[selector],pch=21,bg=classCols[combo$class[selector]],col=ifelse(combo$class[selector] %in% deemphasize,'#00000066','#000000CC'),lwd=ifelse(combo$class[selector]=='Rebound',1.75,1.5),cex=2)
  abline(v=pos['Acute Recipient']-.75,lty=2,col='#00000099')
  abline(v=pos['Acute']-.75,lty=2,col='#00000099')
  #abline(v=pos['Outgrowth MM23']-.75,lty=2,col='#00000099')
  #abline(v=pos['Lorenzi et al. B106']-.75,lty=2,col='#00000099')
}
dev.off()



controls<-list('A09'=c('UK61.1-P2A3','UK61.13-P2C3'),'A06'=c('MM34.12.21D1','MM34.15.11D3'))
basePos<-structure(1:4,.Names=c('Pre-ATI','Rebound','Post-ATI','Control'))
cols<-c('Pre-ATI'='#008000','Post-ATI'='#1CFF1C','Rebound'='#FF00FF','Control'='#FF8000')
pdf('A06A09.pdf',width=5,height=4)
par(mar=c(2,4,1,.1))
for(ii in names(controls)){
  for(addLines in c(FALSE,TRUE)){
  thisDat<-dat[(dat$pat==ii&!is.na(dat$pat))|dat$virus %in% controls[[ii]],]
  pos<-basePos[names(basePos) %in% thisDat$type]
  pos[1:length(pos)]<-1:length(pos)
  spread<-offsetX(log10(thisDat$ic50),thisDat$label,width=.35,varwidth=TRUE)
  spread[thisDat$type!='Rebound']<-withAs(thisDat=thisDat[thisDat$type!='Rebound',],ave(thisDat$Identical,thisDat$label,FUN=function(xx){out<-ifelse(is.na(xx),0,.2*(1+xx-min(xx,na.rm=TRUE)));out<-out-max(out)/2;out}))
  plot(pos[thisDat$type]+spread,thisDat$ic50,log='y',yaxt='n',ylab='IFNa2 IC50 (pg/ml)',xlab='',xaxt='n',type='n',main=ii,ylim=range(dat$ic50,dat$Marvin.s.Value,na.rm=TRUE),xlim=range(pos)+c(-.5,.5))
  if(any(thisDat$type=='Control')){
    segments(pos['Control']-.125,thisDat$Marvin.s.Value,pos['Control']+.125,thisDat$Marvin.s.Value,pch=21,cex=2,col=cols['Control'],lwd=2)
    segments(pos['Control'],thisDat$Marvin.s.Value,pos['Control'],thisDat$ic50,pch=21,cex=2,col=cols['Control'])
    #points(rep(pos['Control'],nrow(thisDat)),thisDat$Marvin.s.Value,cex=2,col=cols[thisDat$type])
  }
  #if(any(thisDat$type=='Control'))segments(pos['Control'],thisDat$Marvin.s.Value,pos['Control'],thisDat$ic50,pch=21,cex=2,bg=cols[thisDat$type],col=cols['Control'])
  for(jj in unique(thisDat$Identical[!is.na(thisDat$Identical)])){
    selector<-thisDat$Identical==jj&!is.na(thisDat$Identical)
    thisX<-pos[thisDat[selector,'type']]+spread[selector]
    meanX<-tapply(thisX,paste(thisDat[selector,'Identical'],thisDat[selector,'label']),mean)
    meanY<-tapply(thisDat[selector,'ic50'],paste(thisDat[selector,'Identical'],thisDat[selector,'label']),mean)
    #lines(thisX[order(thisX)],thisDat[selector,'ic50'][order(thisX)],col='#00000011',lwd=4)
    if(addLines)lines(meanX,meanY,col='#00000022',lwd=4)
  }
  points(pos[thisDat$type]+spread,thisDat$ic50,cex=2.5,bg=sprintf('%sCC',cols[thisDat$type]),pch=20+ifelse(is.na(thisDat$Identical),1,thisDat$Identical+1))
  logAxis(las=1)
  #slantAxis(1,pos,names(pos))
  axis(1,pos,names(pos),las=1)
  }
}
dev.off()


spread<-offsetX(log10(combo$ic50),combo$label,width=.4,varwidth=TRUE)
pdf('studies.pdf',height=2.5,width=6)
for(ii in sort(unique(combo$study))){
  message(ii)
  selector<-combo$study==ii
  par(mar=c(3,4,.1,.1))
  pos<-structure(1:length(unique(combo$label[selector])),.Names=unique(combo$label[selector][orderIn(combo$label[selector],ordering)]))
  plot(pos[combo$label[selector]]+spread[selector],combo$ic50[selector],log='y',yaxt='n',ylab='IFNa2 IC50 (pg/ml)',xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim,xlim=c(.5,length(pos)+.5))
  par(lheight=.75)
  prettyNames<-sub(' MM','\nMM',sub('ATI ','ATI\n',sub('Rebound ','Rebound\n',sub('outgrowth ','outgrowth\n',sub('Outgrowth ','Outgrowth\n',newNames[names(pos)])))))
  twoLines<-grepl('\n',prettyNames)
  for(jj in 1:length(pos))axis(1,pos[jj],prettyNames[jj],mgp=c(3,0,0),padj=1)
  logAxis(las=1)
  points(pos[combo$label[selector]]+spread[selector],combo$ic50[selector],pch=21,bg=classCols[combo$class[selector]],col=ifelse(combo$class[selector] %in% deemphasize,'#00000066','#000000CC'),lwd=ifelse(combo$class[selector]=='Rebound',1.75,1.5),cex=2)
}
dev.off()


