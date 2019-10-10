library(vipor)
library(dnar)
source('functions.R')

beta<-read.csv('stephanieBeta.csv',stringsAsFactors=FALSE)
beta2<-read.csv('stephanieBeta2.csv',stringsAsFactors=FALSE)
beta2<-beta2[,c('Sample','IFN.b')]
colnames(beta2)<-c('virus','ic50')
extraBeta<-read.csv('QVOA for A06 + Nussenz. - Beta.csv')
beta<-rbind(beta,extraBeta[,c('virus','ic50')],beta2)
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
beta<-beta[,colnames(beta)!='ic50']
repCap<-read.csv('stephanieRepCap.csv',stringsAsFactors=FALSE,row.names=1)
beta$repCap<-repCap[beta$virus,'Rep.Cap.']

mm<-read.csv('firstNadir.csv',stringsAsFactors=FALSE)
mm<-mm[mm$isFirst|mm$qvoa|mm$isSix|mm$isLast|mm$isBetaNadir,]
mm$isNadir<-mm$isBetaNadir
mm$class<-ifelse(mm$qvoa,'QVOA',ifelse(mm$isNadir,'Nadir',ifelse(mm$isSix,'6 Month',ifelse(mm$isLast,'Last','Acute'))))
mm$type<-'MM cohort'
#mm$label<-sprintf('MM %s',ifelse(mm$class=='QVOA',sprintf('Outgrowth %s',mm$pat),mm$class))
mm$label<-sprintf('%s',ifelse(mm$class=='QVOA',sprintf('Outgrowth %s',mm$pat),mm$class))
mm<-mm[!is.na(mm$beta),]
mm$repCap<-mm$replication
mm[,colnames(beta)[!colnames(beta) %in% colnames(mm)]]<-NA
mm$source<-'marvin'
mm$virus<-mm[,1]


pair<-read.csv('donorRecipient.csv',stringsAsFactors=FALSE)
pair$class<-ifelse(pair$donor,'Donor','Recipient')
pair$type<-'CHAVI cohort'
#pair$label<-sprintf('CHAVI %s',pair$class)
pair$label<-sprintf('%s %s',ifelse(pair$class=='Donor','Chronic','Acute'),pair$class)
pair$beta<-pair$IFNbeta.Pooled.Donor.cells.IC50..pg.ml
#adjust for super IFN
pair$beta<-pair$beta*6386
pair[pair$class=='Recipient','class']<-'Acute'
pair$pat<-pair$sample
pair$repCap<-pair$Replicative.capacity.Pooled.Donor.cells.p24.d7
pair[,colnames(beta)[!colnames(beta) %in% colnames(pair)]]<-NA
pair$source<-'shilpa'

mont<-read.csv('../out/montaner_ic50.csv',row.names=1)
mont<-mont[grepl('Beta|beta',rownames(mont))&!grepl('UK61',rownames(mont)),]
mont$class<-'Rebound'
mont$virus<-sub(' \\([^)]+\\) - ?','-',sub(' \\(Alpha.*','',rownames(mont)))
mont$type<-sprintf('Montaner %s',sapply(strsplit(rownames(mont),'[()]'),'[',2))
mont$pat<-sub(' .*','',rownames(mont))
mont$label<-sprintf('Rebound %s',mont$pat)
mont$source<-'Montaner'
mont$beta<-mont$ic50
mont[,colnames(beta)[!colnames(beta) %in% colnames(mont)]]<-NA

newRebound<-rbind(
  read.csv('../out/reboundIc50.csv',row.names=1),
  read.csv('../out/rebound_2019-01-30_Ic50.csv',row.names=1),
  read.csv('../out/voa_2019-01-10.csv',row.names=1)
)
tmp<-read.csv('../out/pendingCombined_20190815_ic50.csv',row.names=1)
tmp<-tmp[!grepl('UK61|A08.21_P2F4',row.names(tmp)),!colnames(tmp) %in% c('source','id')]
rownames(tmp)<-sub('alpha and beta','',rownames(tmp))
rownames(tmp)<-sub('(QVOA|Rebound) (.*)','\\2 \\1',rownames(tmp))
newRebound<-rbind(newRebound,tmp)
newRebound<-newRebound[grepl('Beta|beta',rownames(newRebound))&!grepl('UK61|SG3|CH40|WEAU',rownames(newRebound)),]
newRebound$virus<-sub(' \\([^)]+\\)$','',sub(' (Alpha|Beta)?$','',rownames(newRebound)))
newRebound$class<-ifelse(grepl('Rebound|S-22|BEAT',newRebound$virus),'Rebound','QVOA')
newRebound$type<-ifelse(grepl('Rebound|S-22|BEAT',newRebound$virus),'Rebound',ifelse(grepl('MM',newRebound$virus),'Outgrowth',ifelse(grepl('Post',newRebound$virus),'Post-ATI','Pre-ATI')))
newRebound$pat<-sub('[ _.].*','',sub('QVOA ','',sub('Rebound ','',rownames(newRebound))))
newRebound$label<-sprintf('%s %s',newRebound$type,newRebound$pat)
newRebound$source<-'newRebound'
newRebound$beta<-newRebound$ic50
newRebound<-newRebound[!is.na(newRebound$ic50),]
newRebound[,colnames(beta)[!colnames(beta) %in% colnames(newRebound)]]<-NA
newRebound<-newRebound[!grepl('MM[0-9]+\\.',rownames(newRebound)),]

combo<-rbind(beta,mm[,colnames(beta)],pair[,colnames(beta)],mont[,colnames(beta)],newRebound[,colnames(beta)])
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

standardRegex<-'MM14|MM23|MM33|MM34|MM39|MM40'
fastRegex<-'MM15|WEAU'
slowRegex<-'MM55|MM62'
combo$speed<-ifelse(grepl(fastRegex,combo$virus),'Fast',ifelse(grepl(slowRegex,combo$virus),'Slow',ifelse(grepl(standardRegex,combo$virus),'Standard','Other')))
#combo$speed[combo$class=='QVOA']<-'Standard'



write.csv(combo[!combo$source %in% c('marvin','shilpa'),c('virus','pat','type','class','study','source','beta')],'out/comboBeta.csv',row.names=FALSE)

ordering<-c('Rebound S-22','Rebound S-23','Rebound S-30','Pre-ATI A06','Post-ATI A06','Pre-ATI A08','Rebound A08','Post-ATI A08','Pre-ATI A09','Rebound A09','Post-ATI A09','Rebound A08','Rebound 601','Pre-ATI 9201','Rebound 9201','Pre-ATI 9202','Rebound 9202','Pre-ATI 9203','Rebound 9203','Pre-ATI 9207','Rebound 9207','Rebound BEAT-004','Rebound BEAT-030','Rebound BEAT-044','Lorenzi et al. B106','Lorenzi et al. B199','Outgrowth MM14','Outgrowth MM15','Outgrowth MM23','Outgrowth MM34','Outgrowth MM40','Outgrowth MM34','Acute','6 Month','Nadir','Last','Acute Recipient','Chronic Donor')
pos<-structure(1:length(unique(combo$label)),.Names=unique(combo$label[orderIn(combo$label,ordering)]))
posStudy<-sapply(names(pos),function(xx)combo[combo$label==xx,'study'][1])
studySpace<-1.5
pos<-pos+cumsum(c(0,posStudy[-length(posStudy)]!=posStudy[-1]))*studySpace
deemphasize<-c('Acute','6 Month','Donor','Nadir','Last')
uniqClass<-unique(combo$class)
uniqClass<-c(uniqClass[uniqClass!='Rebound'],'Rebound')
uniqClass[uniqClass %in% deemphasize]<-deemphasize[deemphasize %in% uniqClass]
classCols<-structure(rainbow.lab(length(uniqClass),alpha=ifelse(uniqClass %in% deemphasize,.7,.9),lightMultiple=ifelse(uniqClass %in% deemphasize,.65,.75)),.Names=uniqClass)

ranges<-do.call(rbind,tapply(combo$beta,combo$label,range))
rangeClass<-sapply(rownames(ranges),function(xx)combo[combo$label==xx,'class'][1])

subs<-c('A06'='Patient A06','A09'='Patient A09','B106'='Patient B106','B199'='Patient B199','Recipient'='Recipients','Donor'='Donors','Lorenzi et al.'='Outgrowth','S-30'='Patient S-30','BEAT-044'='Patient BEAT-044','BEAT-030'='Patient BEAT-030','MM23'='Patient MM23','MM34'='Patient MM34','ATI'='ATI Outgrowth')
newNames<-sapply(names(pos),function(xx){for(ii in names(subs))xx<-sub(ii,subs[ii],xx);xx})
acuteSpace<-.5
pos<-pos+cumsum(names(pos)=='Acute')*acuteSpace+cumsum(names(pos)=='Acute Recipient')*acuteSpace #+cumsum(names(pos)=='Lorenzi et al. B106')*.5+cumsum(names(pos)=='Outgrowth MM23')*.5
extraSpace<-1
pos<-pos+cumsum(grepl('Acute|6 Month|Nadir|Last|Chronic',names(pos)))*extraSpace
subsets<-rev(list('mm'=names(pos)[!grepl('ATI|Rebound|Lorenzi',names(pos))],'lorenzi'=names(pos)[!grepl('ATI|Rebound',names(pos))],'prePost'=names(pos)[!grepl('Rebound',names(pos))],'noMont'=names(pos)[!grepl('BEAT|S-[0-9]',names(pos))],'all'=names(pos)))

pdf('out/beta_qvoa_compare.pdf',width=9,height=6)
ylim=range(combo$beta)
spread<-offsetX(log10(combo$beta),combo$label,width=.25)
for(ii in names(subsets)){
  message(ii)
  selector<-combo$label %in% subsets[[ii]]
  if(ii=='all')marSpace<-0
  else marSpace<-marSpaces[ii]
  if(ii=='prePost')marSpace<-marSpaces['noMont']
  par(mar=c(9.8,4+marSpace,.1,3))
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

pdf('out/beta_qvoa_compare2.pdf',width=9,height=3.5)
plotQvoa2(combo$beta,combo$label,pos,combo$class,combo$study,combo$speed,ylab='IFNb IC50 (pg/ml)')
dev.off()

pdf('out/studiesBeta.pdf',height=2.5,width=7)
plotStudies(combo$study,combo$label,combo$beta,combo$pat,combo$speed,combo$class,ylab='IFNb IC50 (pg/ml)')
dev.off()

pdf('out/A09beta.pdf',height=3,width=7)
  selector<-grepl('A09',combo$label)
  plotStudies(combo$study[selector],combo$label[selector],combo$beta[selector],combo$pat[selector],combo$speed[selector],combo$class[selector],ylab='IFNb IC50 (pg/ml)')
dev.off()


pdf('out/beta_repCap_compare.pdf',width=8,height=5.5)
  withAs(combo=combo[!is.na(combo$repCap),],plotQvoa2(combo$repCap,combo$label,pos,combo$class,combo$study,combo$speed,ylab='Replication capacity'))
dev.off()


acuteRebound<-combo[combo$class %in% c("Acute","Rebound") & combo$source != "shilpa",]
pdf('out/ifnb_acuteRebound.pdf',width=4,height=4)
par(mar=c(2,4,.1,.1))
suppressWarnings(vipor::vpPlot(acuteRebound$class,acuteRebound$beta,ylab='IFNb IC50',pch=21,bg=classCols[acuteRebound$class],las=1,log='y',yaxt='n'))
dnar::logAxis(las=1)
dev.off()
10^(t.test(log10(acuteRebound[acuteRebound$class=='Rebound','beta']),log10(acuteRebound[acuteRebound$class=='Acute','beta']))$conf.int)
10^(diff(t.test(log10(acuteRebound[acuteRebound$class=='Acute','beta']),log10(acuteRebound[acuteRebound$class=='Rebound','beta']))$estimate))
means<-tapply(log10(acuteRebound$beta),list(acuteRebound$class,acuteRebound$pat),mean)
t.test(means['Acute',],means['Rebound',])
