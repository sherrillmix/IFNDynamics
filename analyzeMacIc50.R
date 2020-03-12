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

