source('functions.R')
concBetaMac<-c(0,4400*10^(-4:0))

p27<-read.csv('ice/02.19.2020 - IC50 Beta Mac.Isol. - Macaque CD4 cells.csv',stringsAsFactors=FALSE)
p27$row<-c(LETTERS[c(1:3,1:3)],NA,LETTERS[c(1:6,1,1:6,1)])
serial<-p27[grepl('serial',p27[,2]),]
p27<-p27[!grepl('serial',p27[,2]),]
tmp<-p27[,11:16]; colnames(tmp)<-colnames(p27)[5:10]
p27Split<-rbind(p27[,5:10],tmp)
p27Split[,]<-as.numeric(sub('[<>]','',unlist(p27Split)))
p27Split$virus<-c(p27[,3],p27[,4])
p27Split$dil<-as.numeric(sub('2/3','1.5',sub(' ?x','',c(p27[,2],p27[,2]))))
p27Split$plate<-c(p27[,1],p27[,1])
table(p27Split$virus)
p27Split<-p27Split[order(p27Split$virus,p27Split$plate,p27Split$dil),]
p27Split[,1:6][p27Split[,1:6]==0]<-10
pdf('out/macIc50_20200226_ic50.pdf',width=20,height=6)
  for(ii in unique(p27Split$virus)){
    par(mfrow=c(2,5))
    for(jj in 1:sum(p27Split$virus==ii)){
      thisDat<-p27Split[p27Split$virus==ii,,drop=FALSE][jj,]
      thisDat$sample<-sprintf('%s %s %sx',thisDat$virus,thisDat$plate,thisDat$dil)
      plotIfns(thisDat,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(p27Split[,1:6]),ylab='p27 concentration (pg/ml)')
    }
  }
dev.off()

p27Split[,c('ic50','vres','max')]<-t(apply(p27Split[,1:6],1,function(xx)calculateBasicIc50(concBetaMac,xx,means=matrix(xx,nrow=1))))[,c(1,3,4)]
p27Split[p27Split$vres>50&is.na(p27Split$ic50),'ic50']<-max(concBetaMac)
p27Split[p27Split$vres<50&is.na(p27Split$ic50),'ic50']<-min(concBetaMac[concBetaMac>0])
p27Split[p27Split[,1]>4000,'ic50']<-NA
p27Split[p27Split[,1]<200,'ic50']<-NA
p27Split[p27Split[,'vres']>100,'ic50']<-NA
p27Split$virus<-sub(' \\(.*','',p27Split$virus)
p27Split[p27Split[,1]<p27Split[,2]/1.5,'ic50']<-NA
#p27Split[(p27Split$virus=='RM10N011.d56.PLAS.C7'&p27Split$plate=='Fresh'&p27Split$X1==449.038),'ic50']<-NA #weird blip between 1st 2nd dose 
#p27Split[(p27Split$virus=='RM10N011.d56.PLAS.C7'&p27Split$plate=='Frozen'&p27Split$X1==872.392),'ic50']<-NA #weird blip between 1st 2nd dose 
p27Split[!is.na(p27Split$ic50),c('virus','dil','plate','ic50','vres')]

pdf('out/macIc50_20200226_ic50_pretty.pdf',width=20,height=6)
  for(ii in unique(p27Split$virus)){
    par(mfrow=c(2,5))
    for(jj in which(!is.na(p27Split[p27Split$virus==ii,'ic50']))){
      thisDat<-p27Split[p27Split$virus==ii,,drop=FALSE][jj,]
      thisDat$sample<-sprintf('%s %s %sx',thisDat$virus,thisDat$plate,thisDat$dil)
      plotIfns(thisDat,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(p27Split[,1:6]),ylab='p27 concentration (pg/ml)')
    }
  }
dev.off()



rowLookup<-data.frame('plate'=rep(c('A','B'),c(8,7)),'row'=LETTERS[c(1:8,1:7)],'origPlate'=rep(c('Frozen','Human','Fresh'),c(6,6,3)),'origRow'=LETTERS[c(1:6,1:6,1:3)],stringsAsFactors=FALSE)
rownames(rowLookup)<-paste(rowLookup$plate,rowLookup$row)
virusLookup<-unique(p27[,c('Plate','row','X','X.1')])
rownames(virusLookup)<-paste(virusLookup$Plate,virusLookup$row)
tzm<-read.csv('out/tzmblMacIc50_20200226.csv',stringsAsFactors=FALSE)
tzm<-rbind(tzm,read.csv('out/tzmblMacIc50_20200228.csv',stringsAsFactors=FALSE))
tzm$ab<-toupper(substring(tzm$plate,2))
tzm$day<-substring(tzm$plate,1,1)
tzm[,c('origPlate','origRow')]<-rowLookup[paste(tzm$ab,tzm$row),c('origPlate','origRow')]
tzm[,c('virus','virus2')]<-virusLookup[paste(tzm$origPlate,tzm$origRow),c('X','X.1')]
tzm<-tzm[!is.na(tzm$virus),]
tmp<-list(tzm[,c(colnames(tzm)[1:6],'virus','origPlate','day')],tzm[,c(colnames(tzm)[7:12],'virus2','origPlate','day')])
colnames(tmp[[2]])<-colnames(tmp[[1]])
tzmStack<-do.call(rbind,tmp)
tzmStack<-tzmStack[order(tzmStack$virus,tzmStack$origPlate,tzmStack$day),]
tzmStack[,c('ic50','vres','max')]<-t(apply(tzmStack[,1:6],1,function(xx)calculateBasicIc50(concBetaMac,xx,means=matrix(xx,nrow=1))))[,c(1,3,4)]
tzmStack[tzmStack$vres>50&is.na(tzmStack$ic50),'ic50']<-max(concBetaMac)
tzmStack[tzmStack$vres<50&is.na(tzmStack$ic50),'ic50']<-min(concBetaMac[concBetaMac>0])
tzmStack$virus<-sub(' \\(.*','',tzmStack$virus)
tzmStack[,c('virus','origPlate','day','ic50','vres')]

pdf('out/macIc50_20200226Tzmbl.pdf',width=20,height=3)
  for(ii in unique(tzmStack$virus)){
    for(day in unique(tzmStack$day)){
      par(mfrow=c(1,6))
      for(jj in 1:sum(tzmStack$virus==ii&tzmStack$day==day)){
        thisDat<-tzmStack[tzmStack$virus==ii&tzmStack$day==day,,drop=FALSE][jj,]
        thisDat$sample<-sprintf('%s %s d%s',thisDat$virus,thisDat$origPlate,thisDat$day)
        plotIfns(thisDat,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(tzmStack[,1:6]),ylab='IU/ul')
      }
    }
  }
dev.off()

tmp<-p27
tmp[,5:16]<-suppressWarnings(apply(tmp[,5:16],2,as.numeric))
tmp[,5:16][tmp[,5:16]<200]<-NA
tmp[,5:16][tmp[,5:16]>6000]<-NA
tmp$dil<-as.numeric(sub('2/3','1.5',sub(' ?x','',tmp[,2])))
tmp[,5:16]<-apply(tmp[,5:16],2,function(xx)xx*tmp$dil)
tmp2<-do.call(rbind,by(tmp,paste(tmp$Plate,tmp$row),function(xx){xx[1,5:16]<-apply(xx[,5:16],2,mean,na.rm=TRUE);xx[1,]}))
tmp3<-withAs(tzm=tzm[tzm$day==7,],tzm[order(tzm$origPlate,tzm$origRow),c(1:12,17:18)])
tmp2[,5:16]/tmp3[,1:12]



plates<-unique(tzmStack$origPlate)
days<-unique(tzmStack$day)
cols<-structure(dnar::rainbow.lab(length(plates),alpha=.7),.Names=plates)
offset<-structure(seq(-.1,.1,length.out=length(plates)),.Names=plates)
dayOffset<-structure(seq(-.3,.1,length.out=length(days)),.Names=days)
cols<-structure(dnar::rainbow.lab(length(plates),alpha=.7),.Names=unique(tzmStack$origPlate))
pdf('out/tzmMacIc50Summary.pdf',width=12)
  par(mar=c(9,4,.1,1))
  pos<-structure(1:length(unique(tzmStack$virus)),.Names=unique(tzmStack$virus))
  plot(1,1,xlim=range(pos)+c(-.5,.5),ylim=range(concBetaMac[concBetaMac>0]),ylab='IFNb IC50',xlab='',xaxt='n',log='y',type='n',yaxt='n')
  abline(v=0:max(pos)+.5,col='#00000033')
  abline(v=0:max(pos)+.3,col='#00000033',lty=3)
  dnar::logAxis(las=1)
  slantAxis(1,pos,names(pos))
  points(pos[tzmStack$virus]+offset[tzmStack$origPlate]+dayOffset[tzmStack$day],tzmStack$ic50,pch=21,bg=cols[tzmStack$origPlate])
  points(pos[p27Split$virus]+offset[p27Split$plate]+.3,p27Split$ic50,pch=22,bg=cols[p27Split$plate])
  legend('bottomleft',names(cols),pch=21,pt.bg=cols,inset=c(-.05,-.25),xpd=NA)
  for(ii in names(cols)){
    plot(1,1,xlim=range(pos)+c(-.5,.5),ylim=range(concBetaMac[concBetaMac>0]),ylab='IFNb IC50',xlab='',xaxt='n',log='y',type='n',yaxt='n')
    abline(v=0:max(pos)+.5,col='#00000033')
    abline(v=1:max(pos)+.3,col='#00000033',lty=3)
    dnar::logAxis(las=1)
    slantAxis(1,pos,names(pos))
    withAs(tzmStack=tzmStack[tzmStack$origPlate==ii,],points(pos[tzmStack$virus]+dayOffset[tzmStack$day],tzmStack$ic50,pch=21,bg=cols[tzmStack$origPlate]))
    withAs(p27Split=p27Split[p27Split$plate==ii,],points(pos[p27Split$virus]+offset[p27Split$plate]+.3,p27Split$ic50,pch=22,bg=cols[p27Split$plate]))
    legend('bottomleft',names(cols),pch=21,pt.bg=cols,inset=c(-.05,-.25),xpd=NA)
  }
dev.off()

pdf('out/tzmMacIc50Summary_d5.pdf')
par(mar=c(9,4,.1,1))
pos<-structure(1:length(unique(tzmStack$virus)),.Names=unique(tzmStack$virus))
for(ii in names(cols)){
plot(1,1,xlim=range(pos)+c(-.5,.5),ylim=range(concBetaMac[concBetaMac>0]),ylab='IFNb IC50',xlab='',xaxt='n',log='y',type='n',yaxt='n')
abline(v=0:max(pos)+.5,col='#00000033')
dnar::logAxis(las=1)
slantAxis(1,pos,names(pos))
withAs(tzmStack=tzmStack[tzmStack$day=='5'&tzmStack$origPlate==ii,],points(pos[tzmStack$virus]+offset[tzmStack$origPlate],tzmStack$ic50,pch=21,bg=cols[tzmStack$origPlate],main=ii))
}
legend('bottomleft',names(cols),pch=21,pt.bg=cols,inset=c(-.12,-.25),xpd=NA)
dev.off()

pdf('out/tzmMacIc50Summary_d7.pdf')
par(mar=c(9,4,.1,1))
pos<-structure(1:length(unique(tzmStack$virus)),.Names=unique(tzmStack$virus))
for(ii in names(cols)){
plot(1,1,xlim=range(pos)+c(-.5,.5),ylim=range(concBetaMac[concBetaMac>0]),ylab='IFNb IC50',xlab='',xaxt='n',log='y',type='n',yaxt='n')
abline(v=0:max(pos)+.5,col='#00000033')
dnar::logAxis(las=1)
slantAxis(1,pos,names(pos))
withAs(tzmStack=tzmStack[tzmStack$day=='7'&tzmStack$origPlate==ii,],points(pos[tzmStack$virus]+offset[tzmStack$origPlate],tzmStack$ic50,pch=21,bg=cols[tzmStack$origPlate],main=ii))
}
legend('bottomleft',names(cols),pch=21,pt.bg=cols,inset=c(-.12,-.25),xpd=NA)
dev.off()



pdf('test.pdf');plot(2*2^(0:11),as.numeric(sub('>','',serial[,5:16])),log='yx',xlab='Dilution',ylab='p27',xaxt='n',yaxt='n');logAxis(las=1);logAxis(1);dev.off()

tzm<-read.csv('out/tzmblMacIc50_20200304.csv',stringsAsFactors=FALSE)
#tzm[,c(1,3,7,9)]<-tzm[,c(3,1,9,7)]
tzm<-tzm[tzm$row!='H',]
tzm<-tzm[order(tzm$plate,tzm$row),]
vir<-read.csv('ice/2020-03-04-macIc50Virus.csv',stringsAsFactors=FALSE)
rownames(vir)<-paste(vir$plate,vir$row)
tzm$virus<-vir[paste(tzm$plate,tzm$row),'virus']
tzm$virPlate<-paste(tzm$plate,tzm$row,tzm$virus)
tzm[,c('ic50','vres','max')]<-t(apply(tzm[,1:12],1,function(xx){yy<-convertIfn6(matrix(xx,ncol=12));calculateBasicIc50(concBetaMac,yy,means=yy)})[c(1,3,4),])
tzm[tzm$vres>40&is.na(tzm$ic50),'ic50']<-max(concBetaMac)
tzm$d<-sub('d','',sapply(strsplit(tzm$virus,'\\.'),'[',2))

tzmConvert<-convertIfn6(tzm[,1:12])
tzmConvert$sample<-rep(sprintf('%s',tzm$virus),each=2)
pdf('out/macIc50_20200304Tzmbl.pdf',width=6,height=6)
        #thisDat<-tzmStack[tzmStack$virus==ii&tzmStack$day==day,,drop=FALSE][jj,]
        #thisDat$sample<-sprintf('%s %s d%s',thisDat$virus,thisDat$origPlate,thisDat$day)
        #plotIfns(thisDat,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(tzmStack[,1:6]),ylab='IU/ul')
for(ii in unique(tzm$virPlate)){
  thisDat<-convertIfn6(tzm[tzm$virPlate==ii,1:12,drop=FALSE])
  thisDat$sample<-ii
  #thisDat[,1:6]<-t(apply(thisDat[,1:6],1,function(xx)xx/xx[1]))
  plotIfns(thisDat,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=calculateBasicIc50,nRep=1,ylab='IU/ul')
}
dev.off()

p27<-read.csv('macaque/p27 ResMac beta IC50 03.06.2020 - Wenge cells.csv',stringsAsFactors=FALSE)
tmpP27<-p27
tmpP27[,c(3,5,9,11)]<-tmpP27[,c(5,3,11,9)]
tmpP272<-p27
p27[,c('ic50','vres','max')]<-t(apply(p27[,3:14],1,function(xx){yy<-convertIfn6(matrix(xx,ncol=12));calculateBasicIc50(concBetaMac,yy,means=yy)})[c(1,3,4),])
tmpP27[,c('ic50','vres','max')]<-t(apply(tmpP27[,3:14],1,function(xx){yy<-convertIfn6(matrix(xx,ncol=12));calculateBasicIc50(concBetaMac,yy,means=yy)})[c(1,3,4),])
t(apply(p27[,c(3:4,6:10,12:14)],1,function(xx){yy<-convertIfn6(matrix(xx,ncol=10));calculateBasicIc50(concBetaMac[-3],yy,means=yy)})[c(1,3,4),])

rt<-read.csv('macaque/RT ResMac beta IC50 03.06.2020 - Wenge cells.csv',stringsAsFactors=FALSE)
colnames(rt)[2]<-'virus'
pdf('out/compare_RT_p27_TZMbl.pdf')
plot(unlist(p27[,3:14]),unlist(tzm[tzm$plate=='p2',1:12]),log='xy',xlab='p27 concentration',ylab='IU/ul')
plot(unlist(rt[,3:14]),unlist(tzm[tzm$plate=='p2',1:12]),log='xy',xlab='RT concentration',ylab='IU/ul')
plot(unlist(p27[,3:14]),unlist(rt[,3:14]),log='xy',xlab='p27 concentration',ylab='RT concentration')
summary(lm(unlist(p27[,3:14]*50/1000)~unlist(rt[,3:14])+0))
summary(lm(unlist(p27[,3:14]*50/1000)~unlist(tzm[tzm$plate=='p2',1:12])+0))
summary(lm(unlist(rt[,3:14])~unlist(tzm[tzm$plate=='p2',1:12])+0))
dev.off()

pdf('out/macIc50_20200306_ic50_all.pdf',width=8,height=6)
  for(ii in unique(p27$virus)){
      thisDat<-p27[p27$virus==ii,c(3:14,1:2,15:ncol(p27)),drop=FALSE]
      thisDat$sample<-sprintf('%s',thisDat$virus)
      tmp<-convertIfn6(thisDat[,1:12])
      tmp$sample<-rep(thisDat$sample,each=2)
      plotIfns(tmp,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(p27[,3:14]),ylab='p27 concentration (pg/ml)')
      thisDat<-tmpP27[tmpP27$virus==ii,c(3:14,1:2,15:ncol(tmpP27)),drop=FALSE]
      thisDat$sample<-sprintf('%s (switch 1 and 3)',thisDat$virus)
      tmp<-convertIfn6(thisDat[,1:12])
      tmp$sample<-rep(thisDat$sample,each=2)
      plotIfns(tmp,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(p27[,3:14]),ylab='p27 concentration (pg/ml)')
      thisDat<-rt[rt$virus==ii,c(3:14,1:2,15:ncol(rt)),drop=FALSE]
      thisDat$sample<-sprintf('%s',thisDat$virus)
      thisDat[,1:12]<-thisDat[,1:12]+.0001
      tmp<-convertIfn6(thisDat[,1:12])
      tmp$sample<-rep(thisDat$sample,each=2)
      plotIfns(tmp,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(rt[,3:14]+.0001),ylab='RT concentration (pg/ml)')
    }
dev.off()

p27_13<-read.csv('macaque/p27 ResMac beta IC50 03.09.2020 - Plate 1 and 3 - Wenge cells.csv',stringsAsFactors=FALSE)
p27_13[,3:14]<-apply(p27_13[,3:14],2,function(xx)as.numeric(sub('>','',xx,xx)))
p27_13[,c('ic50','vres','max')]<-t(apply(p27_13[,3:14],1,function(xx){yy<-convertIfn6(matrix(xx,ncol=12));calculateBasicIc50(concBetaMac,yy,means=yy)})[c(1,3,4),])
p27<-rbind(p27,p27_13)
p27$virus<-ifelse(p27$virus=='SHIV.D.375M (human CD4 viral stock)',paste(p27$virus,'plate',p27$plate),p27$virus)
tzm$virus<-ifelse(tzm$virus=='SHIV.D.375M (human CD4 viral stock)',paste(tzm$virus,'plate',sub('p','',tzm$plate)),tzm$virus)
tmpP27$virus<-ifelse(tmpP27$virus=='SHIV.D.375M (human CD4 viral stock)',paste(tmpP27$virus,'plate',sub('p','',tmpP27$plate)),tmpP27$virus)
rt$virus<-ifelse(rt$virus=='SHIV.D.375M (human CD4 viral stock)',paste(rt$virus,'plate',2),rt$virus)

rm10<-p27[!grepl('SHIV.D|10/17/2016',p27$virus),c(3:ncol(p27),1:2)]
rm10$day<-as.numeric(sub('d','',sapply(strsplit(rm10$virus,'\\.'),'[',2)))
rm10$monkey<-sub('d','',sapply(strsplit(rm10$virus,'\\.'),'[',1))
rm10$source<-sub('d','',sapply(strsplit(rm10$virus,'\\.'),'[',3))
sourceCol<-c('PBMC'=NA,'PLAS'='#FF000099')
sourceCex<-c('PBMC'=1,'PLAS'=2)
vir<-unique(read.csv('ice/2020-02-01_plateIds.csv',stringsAsFactors=FALSE)[,c('virus','animal')])
vir$virus[vir$virus=='MT145.Q171K']<-'MT145'
rownames(vir)<-vir$animal
pdf('out/macIfnOptions.pdf',width=15,height=6)
par(mfrow=c(2,3)) 
par(mar=c(5,4,1,1))
for(ii in 2:6){
  thisVres<-apply(rm10[,c(ii,ii+6)]/rm10[,c(1,7)],1,mean)
  pos<-structure(1:length(unique(rm10$day)),.Names=sort(unique(rm10$day)))
  plot(pos[as.character(rm10$day)],thisVres,xaxt='n',xlab='',ylab='Proportion of untreated p27',main=sprintf('%s ng/ml IFNb',format(concBetaMac[ii])),las=1,xlim=range(pos)+c(-.5,.5),ylim=c(0,2),bg=sourceCol[rm10$source],pch=21,cex=sourceCex[rm10$source])
  posInfo<-sapply(names(pos),function(xx){
    monkeys<-unique(rm10$monkey[rm10$day==xx])
    viruses<-vir[monkeys,'virus']
    paste(monkeys,viruses,sep='\n',collapse='\n')
  })
  axis(1,pos,paste(names(pos),'\n',posInfo),padj=1,mgp=c(3,.1,0))
  mtext('Day:',1,line=.7,at=dnar::convertLineToUser(-1.5,2),cex=.8,adj=1)
  mtext('Animal:',1,line=1.7,at=dnar::convertLineToUser(-1.5,2),cex=.8,adj=1)
  mtext('Virus:',1,line=2.7,at=dnar::convertLineToUser(-1.5,2),cex=.8,adj=1)
  abline(h=1,lty=2)
}
dev.off()



pdf('out/macIc50_20200306_ic50.pdf',width=15,height=6)
  for(ii in unique(p27$virus)){
      par(mfrow=c(1,3))
      thisDat<-p27[p27$virus==ii,c(3:14,1:2,15:ncol(p27)),drop=FALSE]
      thisDat$sample<-sprintf('%s',thisDat$virus)
      #thisDat[,1:12]<-as.numeric(sub('1 to ','',thisDat$Dilution))*thisDat[,1:12]/1000
      #tmp<-thisDat[,7:12,drop=FALSE]
      #tmp$sample<-thisDat$sample
      tmp<-convertIfn6(thisDat[,1:12])
      tmp<-as.numeric(sub('1 to ','',thisDat$Dilution))*tmp/1000
      tmp$sample<-sprintf('%s',thisDat$sample)
      plotIfns(tmp,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(p27[,3:14]*as.numeric(sub('1 to ','',p27$Dilution))/1000),ylab='p27 concentration (ng/ml)',log='xy',showFit=FALSE,showPercent=FALSE)
      rect(10^par('usr')[1],50*7000/1000,10^par('usr')[2],10^par('usr')[4],col='#FF000012',border=NA)
      rect(10^par('usr')[1],50*50/1000,10^par('usr')[2],10^par('usr')[1],col='#FF000012',border=NA)
      thisDat[,1:6]<-t(apply(thisDat[,1:6],1,function(xx)xx/xx[1]))*100
      thisDat[,7:12]<-t(apply(thisDat[,7:12],1,function(xx)xx/xx[1]))*100
      #thisDat[,1:12]<-as.numeric(sub('1 to ','',thisDat$Dilution))*thisDat[,1:12]/1000
      tmp<-convertIfn6(thisDat[,1:12])
      tmp$sample<-rep(thisDat$sample,each=2)
      plotIfns(tmp,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=c(0,450),ylab='Percent of untreated p27 concentration',log='x',showFit=FALSE)
      thisDat<-tzm[tzm$virus==ii,,drop=FALSE]
      thisDat$sample<-sprintf('%s',thisDat$virus)
      #thisDat[,1:12]<-as.numeric(sub('1 to ','',thisDat$Dilution))*thisDat[,1:12]/1000
      thisDat[,1:6]<-t(apply(thisDat[,1:6],1,function(xx)xx/xx[1]))*100
      thisDat[,7:12]<-t(apply(thisDat[,7:12],1,function(xx)xx/xx[1]))*100
      #thisDat[,1:12]<-as.numeric(sub('1 to ','',thisDat$Dilution))*thisDat[,1:12]/1000
      #tmp<-thisDat[,7:12,drop=FALSE]
      #tmp$sample<-thisDat$sample
      tmp<-convertIfn6(thisDat[,1:12])
      tmp$sample<-rep(thisDat$sample,each=2)
      #plotIfns(tmp,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(p27[,3:14]*as.numeric(sub('1 to ','',p27$Dilution))/1000),ylab='p27 concentration (ng/ml)',log='x')
      plotIfns(tmp,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=c(0,450),ylab='Percent of untreated TZMbl IU',log='x',showFit=FALSE)
    }
dev.off()

tzm<-read.csv('out/tzmblIc50_20200330.csv',stringsAsFactors=FALSE)
#tzm[,c(1,3,7,9)]<-tzm[,c(3,1,9,7)]
tzm<-tzm[!grepl('^mac',tzm$plate),]
tzm<-tzm[order(tzm$plate,tzm$row),]
tzm[,1:12]<-tzm[,12:1]
vir<-read.csv('ice/2020-03-30-ic50Virus.csv')
rownames(vir)<-paste(vir$plate,vir$row)
tzm$virus<-vir[paste(tzm$plate,tzm$row),'virus']
tzm$virPlate<-paste(tzm$plate,tzm$row,tzm$virus)
tzm[,c('ic50','vres','max')]<-t(apply(tzm[,1:12],1,function(xx){yy<-convertIfn6(matrix(xx,ncol=12));calculateBasicIc50(concBeta6,yy,means=yy)})[c(1,3,4),])
tzm[tzm$vres>40&is.na(tzm$ic50),'ic50']<-max(concBeta6)
tzmConvert<-convertIfn6(tzm[,1:12])
tzmConvert$sample<-rep(sprintf('%s',tzm$virus),each=2)
pdf('out/ic50_20200330Tzmbl.pdf',width=6,height=6)
        #thisDat<-tzmStack[tzmStack$virus==ii&tzmStack$day==day,,drop=FALSE][jj,]
        #thisDat$sample<-sprintf('%s %s d%s',thisDat$virus,thisDat$origPlate,thisDat$day)
        #plotIfns(thisDat,concBetaMac,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=function(concs,p24s)calculateBasicIc50(concs,p24s,means=p24s),nRep=1,ylims=range(tzmStack[,1:6]),ylab='IU/ul')
for(ii in unique(tzm$virPlate)){
  thisDat<-convertIfn6(tzm[tzm$virPlate==ii,1:12,drop=FALSE])
  thisDat$sample<-ii
  #thisDat[,1:6]<-t(apply(thisDat[,1:6],1,function(xx)xx/xx[1]))
  plotIfns(thisDat,concBeta6,'IFNb concentration (pg/ml)',condenseTechs=FALSE,findVresIc50=findVresIc50,nRep=1,ylab='IU/ul')
}
dev.off()
#looks like garbage


p27<-read.csv('macaque/03.19.2020 - IC50 Beta Mac.Isol. - Macaque CD4 cells. csv .csv',stringsAsFactors=FALSE,check.names=FALSE)
p27<-p27[,!grepl('Standard',colnames(p27))]
p27[,grep('UT|untreat',colnames(p27))] #"not so good" cells sucked. ignoring
p27Cols<-colnames(p27)[grep('IFN|untreat',colnames(p27))]
p27$dil<-as.numeric(sub('1 to ','',p27$dilution))
p27Stack<-data.frame(
  'virus'=rep(p27$Virus,length(p27Cols)),
  'plate'=rep(p27$Plate,length(p27Cols)),
  'ifn'=rep(p27Cols,each=nrow(p27)),
  'col'=rep(1:length(p27Cols),each=nrow(p27)),
  'dil'=rep(p27$dil,length(p27Cols)),
  'orig'=unlist(p27[,p27Cols]),
  stringsAsFactors=FALSE
)
p27Stack<-p27Stack[order(p27Stack$virus,p27Stack$plate,p27Stack$col,p27Stack$dil),]
p27Stack$p27<-as.numeric(sub('[><]','',p27Stack$orig))
p27Stack$p27Convert<-p27Stack$p27*(p27Stack$dil)#*(1+p27Stack$dil)
tmp<-paste(p27Stack$virus,p27Stack$ifn)
if(any(tmp[seq(1,length(tmp),2)]!=tmp[seq(2,length(tmp),2)]))stop('Mismatch in p27')
if(any(round(p27Stack[seq(1,length(tmp),2),'dil'],1)!=round(p27Stack[seq(2,length(tmp),2),'dil']/10,1)))stop('Mismatch in p27')
p27Merge<-merge(p27Stack[p27Stack$dil==2.5,],p27Stack[p27Stack$dil==25,],all.x=TRUE,all.y=TRUE,by=c('virus','plate','col','ifn'),suffixes=c('_2.5','_25'))
if(nrow(p27Merge)!=nrow(p27Stack)/2)stop('Merge problem')
p27Merge[grepl('SL92b',p27Merge$virus),'virus']<-apply(p27Merge[grepl('SL92b',p27Merge$virus),c('virus','plate')],1,paste,collapse=' ')
p27Merge$p27Best<-apply(cbind(p27Merge[,c('p27_2.5','p27_25','p27Convert_2.5','p27Convert_25')]),1,function(xx){
  isHigh<-xx[1:2]>5000
  isLow<-xx[1:2]<200
  if(all(isLow))return(xx[3])
  if(all(isHigh)){warning('Both high');return(xx[4])}
  if(all(isHigh|isLow)){warning('Both high/low');return(mean(xx[3:4]))}
  return(mean(xx[3:4][!isHigh&!isLow]))
})

virusIfn<-tapply(p27Merge$p27Best,p27Merge[,c('virus','ifn')],c)
virusIfn[virusIfn<50]<-50
simpleNames<-sub('\\.[0-9]+$','',colnames(virusIfn))
conc<-as.numeric(ifelse(grepl('untreated',colnames(virusIfn)),0,sub(',','',sub('pg.*','',colnames(virusIfn)))))
isAlpha<-grepl('IFNa2',colnames(virusIfn))
colPos<-structure(1:length(unique(simpleNames)),.Names=unique(simpleNames[order(isAlpha,conc)]))
pdf('out/macIc50_20200330.pdf')
  for(ii in 1:nrow(virusIfn)){
    plot(1,1,type='n',xlab='',xaxt='n',log='y',yaxt='n',xlim=range(colPos)+c(-.5,.5),ylim=range(virusIfn),main=rownames(virusIfn)[ii],ylab='p27 concentration',xaxs='i')
    logAxis(las=1)
    points(colPos[simpleNames],virusIfn[ii,],pch=21,cex=1.2,bg='#00000033')
    axis(1,colPos,sub(' ','\n',names(colPos)),padj=1,mgp=c(3,.2,0))
    abline(v=colPos[c(1,length(colPos))]+c(.5,-.5),lty=1,col='#00000066')
    axis(4,exp(mean(log(virusIfn[ii,conc==0])))/c(1,2),c('UT','50%'),las=1,mgp=c(3,.2,0),tcl=0)
    abline(h=exp(mean(log(virusIfn[ii,conc==0])))/c(1,2),lty=2)
  }
dev.off()

pdf('out/macDilCompare.pdf')
#plot(p27Stack[p27Stack$dil==2.5,'p27Convert'],p27Stack[p27Stack$dil==25,'p27Convert'],xlab='2.5x dilution p27',ylab='25x dilution p27')
#abline(0,1)
plot(p27Merge[,'p27Convert_2.5'],p27Merge[,'p27Convert_25'],xlab='2.5x dilution p27',ylab='25x dilution p27',main='1 to 2.5 = 3.5x')
abline(0,1)
#plot(p27Merge[,'p27Convert_2.5']/3.5*2.5,p27Merge[,'p27Convert_25']/26*25,xlab='2.5x dilution p27',ylab='25x dilution p27',main='1 to 2.5 = 2.5x')
#abline(0,1)
dev.off()



virusIfn<-tapply(p27Merge$p27Best,list(p27Merge$virus,sub('\\.[0-9]+$','',p27Merge$ifn)),function(xx)exp(mean(log(xx))))
virusIfn[virusIfn<50]<-50
virusIfn<-virusIfn[,order(colnames(virusIfn)!='untreated',!grepl('IFNb',colnames(virusIfn)),as.numeric(gsub('[a-zA-Z, ]+','',colnames(virusIfn))))]
virusIfn<-virusIfn[virusIfn[,1]>50,]
virusIfn<-virusIfn[!grepl('4x less',rownames(virusIfn)),]
rm<-data.frame('virus'=rownames(virusIfn),stringsAsFactors=FALSE)
rownames(rm)<-rm$virus
rm$control<-!grepl('^RM',rm$virus)
rm$day<-as.numeric(sub('d','',ifelse(rm$control,NA,sapply(strsplit(rm$virus,'\\.'),'[',2))))
rm$monkey<-ifelse(rm$control,NA,sapply(strsplit(rm$virus,'\\.'),'[',1))
rm$source<-ifelse(rm$control,'control',sapply(strsplit(rm$virus,'\\.'),'[',3))
rm$short<-sub('\\.375[A-Z]','',sub('SHIV\\.A?\\.?','',sub(' \\(human CD4 viral stock\\)','',rm$virus)))
vir<-unique(read.csv('ice/2020-02-01_plateIds.csv',stringsAsFactors=FALSE)[,c('virus','animal')])
vir$virus[vir$virus=='MT145.Q171K']<-'MT145'
rownames(vir)<-vir$animal
rm$input<-vir[rm$monkey,'virus']
rm$input[grepl('BG505',rm$virus)]<-'BG505'
rm$input[grepl('1176',rm$virus)]<-'1176'
rm$input[is.na(rm$input)]<-sub(' \\(.*','',rm$virus[is.na(rm$input)])
rm$note<-trimws(sub('^[^ ]+','',rm$short))
virusIfn<-virusIfn[order(rm$input,!rm$control,rm$monkey,rm$day),]
rm<-rm[rownames(virusIfn),]
rm$label<-paste(ifelse(is.na(rm$day),'',rm$day),ifelse(is.na(rm$monkey),'',rm$monkey),rm$input,sep='\n')
rm$pos<-as.numeric(factor(paste(rm$monkey,rm$day,rm$input),levels=unique(paste(rm$monkey,rm$day,rm$input))))



#sourceCol<-c('PBMC'=NA,'PLAS'='#FF000099','control'=NA)
sourceCol<-c('PBMC'='#00000033','PLAS'='#FF000033','control'='#00000033')
#sourceCex<-c('PBMC'=1,'PLAS'=2,'control'=1)
sourceCex<-c('PBMC'=1.5,'PLAS'=1.5,'control'=1.5)
pdf('out/macIfns_20200406.pdf',width=18,height=6,useDingbats=FALSE)
par(mfrow=c(2,3)) 
par(mar=c(5,4,1,1))
for(ii in 2:ncol(virusIfn)){
  thisVres<-virusIfn[,ii]/virusIfn[,1]
  plot(rm$pos,thisVres,xaxt='n',xlab='',ylab='Proportion of untreated p27',main=colnames(virusIfn)[ii],las=1,xlim=range(rm$pos)+c(-.5,.5),ylim=c(0,2.5),bg=sourceCol[rm$source],pch=21,cex=sourceCex[rm$source])
  sapply(which(!duplicated(rm$pos)),function(xx)axis(1,rm$pos[xx],rm$label[xx],padj=1,mgp=c(3,.1,0)))
  mtext('Day:',1,line=.7,at=dnar::convertLineToUser(-1.25,2),cex=.8,adj=1)
  mtext('Animal:',1,line=1.7,at=dnar::convertLineToUser(-1.25,2),cex=.8,adj=1)
  mtext('Virus:',1,line=2.7,at=dnar::convertLineToUser(-1.25,2),cex=.8,adj=1)
  abline(h=1,lty=2)
  abline(v=c(max(rm[rm$input=='BG505','pos']),max(rm[rm$input=='1176','pos']))+.5,lty=1)
}
dev.off()


virus<-c( "SL92b", "SHIV.A.BG505.375Y", "RM10N011.d56.PLAS.C7", "RM10N011.d942.PLAS.E5", "RM10N011.d942.PLAS.F4", "RM6563.d1065.PBMC.D7", "RM6563.d1456.PBMC.C7")
p27<-read.csv('data/p27 ResMac beta IC50 06.19.2020.csv',stringsAsFactors=FALSE,check.names=FALSE,nrow=7)
ifns<-c( "untreated", "0.044pg IFNb", "4.4pg IFNb", "4,400pg IFNb", "440,000pg IFNb", "5.5pg IFNa2")
p27Stack<-data.frame(
  'virus'=rep(virus,12),
  'ifn'=rep(rep(ifns,2),each=nrow(p27)),
  'orig'=unlist(p27[,2:13]),
  stringsAsFactors=FALSE
)
p27Stack$p27<-as.numeric(sub('[><]','',p27Stack$orig))
p27Stack$censor<-p27Stack$p27
p27Stack$censor[p27Stack$censor<50]<-50
colPos<-structure(1:length(unique(ifns)),.Names=ifns)
pdf('out/macIc50_20200624.pdf',height=8,width=18)
par(mfrow=c(2,4))
  for(ii in virus){
    thisDat<-p27Stack[p27Stack$virus==ii,]
    thisDat<-thisDat[order(thisDat$ifn),]
    repCap<-exp(mean(log(thisDat[thisDat$ifn=='untreated'&thisDat$censor>50,'p27'])))*200/1000
    plot(1,1,type='n',xlab='',xaxt='n',log='y',yaxt='n',xlim=range(colPos)+c(-.5,.5),ylim=range(p27Stack$censor),main=sprintf('%s\nRep cap: %0.1f',ii,repCap),ylab='p27 concentration',xaxs='i')
    dnar::logAxis(las=1)
    points(colPos[thisDat$ifn]+c(-.05,.05),thisDat$censor,pch=21,cex=1.2,bg='#00000033')
    axis(1,colPos,sub(' ','\n',names(colPos)),padj=1,mgp=c(3,.2,0))
    abline(v=colPos[c(1,length(colPos))]+c(.5,-.5),lty=1,col='#00000066')
    axis(4,exp(mean(log(thisDat[thisDat$ifn=='untreated'&thisDat$censor>50,'p27'])))/c(1,2),c('UT','50%'),las=1,mgp=c(3,.2,0),tcl=0)
    abline(h=exp(mean(log(thisDat[thisDat$ifn=='untreated'&thisDat$censor>50,'p27'])))/c(1,2),lty=2)
    abline(h=50,lty=3)
  }
dev.off()


