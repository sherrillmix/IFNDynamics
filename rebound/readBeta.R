
beta<-read.csv('stephanieBeta.csv',stringsAsFactors=FALSE)
beta$pat<-ifelse(grepl('621803B',beta$virus),'A06',
  ifelse(grepl('621818B|A09-|Rebound ',beta$virus),'A09',
    ifelse(grepl('^B[0-9]+_',beta$virus),sub('_.*','',beta$virus),'')))
beta$type<-ifelse(grepl('pre',beta$virus),'Pre-ATI',
  ifelse(grepl('post',beta$virus),'Post-ATI', 
    ifelse(grepl('Rebound|A09-',beta$virus),'Rebound',
      ifelse(grepl('^B[0-9]+_',beta$virus),'Lorenzi et al.',
        ifelse(grepl('UK61|MM34',beta$virus),'Control',NA)))))
qvoaGroups<-c('Pre-ATI','Post-ATI','Lorenzi et al.')
beta$class<-ifelse(beta$type %in% qvoaGroups,'QVOA',beta$type)
beta$label<-sprintf('%s%s%s',beta$type,ifelse(beta$pat!='',' ',''),beta$pat)
beta$beta<-apply(beta[,grepl('ic50',colnames(beta)),drop=FALSE],1,mean)
beta$source<-'stephanie'
beta<-beta[colnames(beta)!='ic50',]

mm<-read.csv('firstNadir.csv',stringsAsFactors=FALSE)
mm$class<-ifelse(mm$qvoa,'QVOA',ifelse(mm$isNadir,'Nadir',ifelse(mm$isSix,'6 Month',ifelse(mm$isLast,'Last','Acute'))))
mm$type<-'MM cohort'
#mm$label<-sprintf('MM %s',ifelse(mm$class=='QVOA',sprintf('Outgrowth %s',mm$pat),mm$class))
mm$label<-sprintf('%s',ifelse(mm$class=='QVOA',sprintf('Outgrowth %s',mm$pat),mm$class))
mm<-mm[!is.na(mm$beta),]
mm[,colnames(beta)[!colnames(beta) %in% colnames(mm)]]<-NA
mm$source<-'marvin'

pair<-read.csv('donorRecipient.csv',stringsAsFactors=FALSE)
pair$class<-ifelse(pair$donor,'Donor','Recipient')
pair$type<-'CHAVI cohort'
#pair$label<-sprintf('CHAVI %s',pair$class)
pair$label<-sprintf('%s %s',ifelse(pair$class=='Donor','Chronic','Acute'),pair$class)
pair$beta<-pair$IFNbeta.Pooled.Donor.cells.IC50..pg.ml
#adjust for super IFN
pair$beta<-pair$beta*6386
pair[pair$class=='Recipient','class']<-'Acute'
pair[,colnames(beta)[!colnames(beta) %in% colnames(pair)]]<-NA
pair$source<-'shilpa'

mont<-read.csv('../out/montaner_ic50.csv',row.names=1)
mont<-mont[grepl('Beta|beta',rownames(mont))&!grepl('UK61',rownames(mont)),]
mont$class<-'Rebound'
mont$Virus<-sub(' \\([^)]+\\) - ?','-',sub(' \\(Alpha.*','',rownames(mont)))
mont$type<-sprintf('Montaner %s',sapply(strsplit(rownames(mont),'[()]'),'[',2))
mont$pat<-sub(' .*','',rownames(mont))
mont$label<-sprintf('Rebound %s',mont$pat)
mont$source<-'Montaner'
mont$beta<-mont$ic50
mont[,colnames(beta)[!colnames(beta) %in% colnames(mont)]]<-NA

combo<-rbind(beta,mm[,colnames(beta)],pair[,colnames(beta)],mont[,colnames(beta)])
combo<-combo[combo$class!='Control',]
write.csv(combo,'comboBeta.csv')

ordering<-c('Rebound BEAT-030','Rebound BEAT-044','Rebound S-30','Pre-ATI A06','Post-ATI A06','Pre-ATI A09','Rebound A09','Post-ATI A09','Lorenzi et al. B106','Lorenzi et al. B199','Outgrowth MM23','Outgrowth MM34','Acute','6 Month','Nadir','Last','Acute Recipient','Chronic Donor')
pos<-structure(1:length(unique(combo$label)),.Names=unique(combo$label[orderIn(combo$label,ordering)]))
deemphasize<-c('Acute','6 Month','Donor','Nadir','Last')
uniqClass<-unique(combo$class)
uniqClass<-c(uniqClass[uniqClass!='Rebound'],'Rebound')
uniqClass[uniqClass %in% deemphasize]<-deemphasize[deemphasize %in% uniqClass]
classCols<-structure(rainbow.lab(length(uniqClass),alpha=ifelse(uniqClass %in% deemphasize,.7,.9),lightMultiple=ifelse(uniqClass %in% deemphasize,.65,.75)),.Names=uniqClass)

ranges<-do.call(rbind,tapply(combo$beta,combo$label,range))
rangeClass<-sapply(rownames(ranges),function(xx)combo[combo$label==xx,'class'][1])

subs<-c('A06'='Patient 1','A09'='Patient 2','B106'='Patient 1','B199'='Patient 2','Recipient'='Recipients','Donor'='Donors','Lorenzi et al.'='Nussenzweig')
newNames<-sapply(names(pos),function(xx){for(ii in names(subs))xx<-sub(ii,subs[ii],xx);xx})
pos<-pos+cumsum(names(pos)=='Acute')*.5+cumsum(names(pos)=='Acute Recipient')*.5 #+cumsum(names(pos)=='Lorenzi et al. B106')*.5+cumsum(names(pos)=='Outgrowth MM23')*.5
subsets<-rev(list('mm'=names(pos)[!grepl('ATI|Rebound|Lorenzi',names(pos))],'lorenzi'=names(pos)[!grepl('ATI|Rebound',names(pos))],'prePost'=names(pos)[!grepl('Rebound',names(pos))],'noMont'=names(pos)[!grepl('BEAT|S-[0-9]',names(pos))],'all'=names(pos)))
pdf('beta_qvoa_compare.pdf',width=8,height=5)
ylim=range(combo$beta)
spread<-offsetX(log10(combo$beta),combo$label,width=.25)
for(ii in names(subsets)){
  message(ii)
  selector<-combo$label %in% subsets[[ii]]
  if(ii=='all')marSpace<-0
  else marSpace<-marSpaces[ii]
  if(ii=='prePost')marSpace<-marSpaces['noMont']
  par(mar=c(7.3,4+marSpace,.1,3))
  plot(pos[combo$label[selector]]+spread[selector],combo$beta[selector],log='y',yaxt='n',ylab='IFNb IC50 (pg/ml)',xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim)
  if(ii=='all')marSpaces<-sapply(subsets,function(xx)diff(convertUserToLine(1:0,2))*sum(!names(pos) %in% xx))
  slantAxis(1,pos[names(pos) %in% subsets[[ii]]],newNames[names(pos) %in% subsets[[ii]]],srt=-45)
  logAxis(las=1)
  #segments(pos[rownames(ranges)],ranges[,1],pos[rownames(ranges)],ranges[,2],col=classCols2[rangeClass],lwd=2)
  #,cex=ifelse(combo$class %in% deemphasize,1.5)
  points(pos[combo$label[selector]]+spread[selector],combo$beta[selector],pch=21,bg=classCols[combo$class[selector]],col=ifelse(combo$class[selector] %in% deemphasize,'#00000066','#000000CC'),lwd=ifelse(combo$class[selector]=='Rebound',1.75,1.5),cex=2)
  abline(v=pos['Acute Recipient']-.75,lty=2,col='#00000099')
  abline(v=pos['Acute']-.75,lty=2,col='#00000099')
  #abline(v=pos['Outgrowth MM23']-.75,lty=2,col='#00000099')
  #abline(v=pos['Lorenzi et al. B106']-.75,lty=2,col='#00000099')
}
dev.off()


