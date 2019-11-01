speedPch<-c('Fast'=8,'Slow'=21,'Standard'=21,'Other'=21)
speedCex<-c('Fast'=1,'Slow'=2,'Standard'=2,'Other'=2)

plotQvoa2<-function(ic50,label,pos,class,study,speed,ylab='IFNa2 IC50 (pg/ml)',mar=c(6.9,4,.1,3.9),cex.axis=1.25,startDown=FALSE){
  spread<-offsetX(log10(ic50),label,width=.25)
  ylim=range(ic50)
  spread<-offsetX(log10(ic50),label,width=.25)
  ii<-'all'
  marSpace<-0
  selector<-label %in% subsets[[ii]]
  par(mar=mar)
  plot(pos[label[selector]]+spread[selector],ic50[selector],log='y',yaxt='n',ylab=ylab,xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim,mgp=c(2.5,1,0),xlim=range(pos)+c(-1,1),xaxs='i')
  if(ii=='all')marSpaces<-sapply(subsets,function(xx)diff(convertUserToLine(1:0,2))*sum(!names(pos) %in% xx))
  slantSelect<-!grepl('Outgrowth',names(pos)) &names(pos) %in% subsets[['mm']]
  textOffsets<-c('Acute'=-.3,'6 Month'=-.15,'Nadir'=.15,'Last'=.3,'Acute Recipients'=-.3,'Chronic Donors'=.3)[newNames[slantSelect]]
  textOffsets[is.na(textOffsets)]<-0
  slantAxis(1,pos[slantSelect],newNames[slantSelect],srt=-45,cex=cex.axis,location=.8,xpd=NA,textOffsets=textOffsets)
  logAxis(las=1,mgp=c(3,.7,0))
  abline(v=pos,col='#00000055',lty=3)
  points(pos[label[selector]]+spread[selector],ic50[selector],pch=speedPch[speed[selector]],bg=classCols[class[selector]],col=ifelse(speedPch[speed[selector]]>20,ifelse(class[selector] %in% deemphasize,'#00000066','#000000CC'),classCols[class[selector]]),lwd=ifelse(class[selector]=='Rebound',1.75,1.5),cex=ifelse(class[selector]=='QVOA',max(speedCex),speedCex[speed[selector]]))
  abline(v=pos['Acute Recipient']-(pos['Acute Recipient']-pos[which(names(pos)=='Acute Recipient')-1])/2,lty=2,col='#00000099')
  abline(v=pos['Outgrowth MM14']-(pos['Outgrowth MM14']-pos[which(names(pos)=='Outgrowth MM14')-1])/2,lty=1,col='#00000099')
  cols<-c('#00000033','#00000000')
  studyOrder<-unique(sapply(names(pos),function(xx)study[label==xx][1]))
  studyOrder<-studyOrder[studyOrder!='Transmission'&!is.na(studyOrder)]
  counter<-1
    for(ii in studyOrder){
      if(ii=='MM'){
        minPos<-pos[min(which(names(pos) %in% label[study==ii&class=='QVOA']))]
        maxPos<-pos[max(which(names(pos) %in% label[study==ii&class=='QVOA']))]
      }else{
        minPos<-pos[min(which(names(pos) %in% label[study==ii]))]
        maxPos<-pos[max(which(names(pos) %in% label[study==ii]))]
        if(ii!='Reservoir')abline(v=maxPos+.5+studySpace/2,lty=2)
      }
      #rect(minPos-.5,10^par('usr')[3],maxPos+.5,10^par('usr')[4],col=cols[counter%%2+1],border=NA)
      counter<-counter+1
      axis(1,mean(c(minPos,maxPos)),sub('MM','Outgrowth',sub('ATI','Interrupt',sub('BEAT','IFNa2',sub('/','/\n',ii)))),padj=1,mgp=c(3,.2+2*(startDown+counter)%%2,0),tcl=-.7+-2*(startDown+counter)%%2,cex.axis=cex.axis)
    }
  return(list(ylim=ylim))
}

plotStudies<-function(study,label,ic50,pat,speed,class,ylab='IFNa2 IC50 (pg/ml)'){
spread<-offsetX(log10(ic50),label,width=.4,varwidth=TRUE)
for(ii in sort(unique(study))){
  message(ii)
  selector<-study==ii
  par(mar=c(3.1,4,.1,.1))
  thisPos<-structure(1:length(unique(label[selector])),.Names=unique(label[selector][orderIn(label[selector],ordering)]))
  plot(thisPos[label[selector]]+spread[selector],ic50[selector],log='y',yaxt='n',ylab=ylab,xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim,xlim=c(.5,max(thisPos)+.5),xaxs='i')
  par(lheight=.75)
  prettyNames<-sub('\n\n','\n',sub(' ?Patient ?','',sub(' MM','\nMM',sub('ATI ','ATI\n',sub('Rebound ','Rebound\n',sub('outgrowth ','outgrowth\n',sub('Outgrowth (Patient )?','Outgrow\n',newNames[names(thisPos)])))))))
  twoLines<-grepl('\n',prettyNames)
  counter<-1
  nLines<-sapply(gregexpr('\n',prettyNames),length)
  maxLines<-max(nLines)
  for(jj in 1:length(thisPos)){
    if(length(thisPos)>7)suppressWarnings(axis(1,thisPos[jj],sub('growth','grow',prettyNames[jj]),mgp=c(3,0,0),padj=1,mgp=c(3,-.4+counter%%2,0),tcl=-.25+-1*counter%%2,cex.axis=1.2))
    else suppressWarnings(axis(1,thisPos[jj],sprintf('%s%s',ifelse(nLines[jj]<maxLines,'\n',''),prettyNames[jj]),mgp=c(3,0,0),padj=1,mgp=c(3,-.4,0),tcl=-.25,cex.axis=1.2))
    counter<-counter+1
  }
  logAxis(las=1)
  if(any(class[selector]=='Rebound')){
    for(reb in unique(pat[selector&class=='Rebound'])){
      if(any(pat==reb&class=='QVOA')){
        rect(thisPos[min(grep(reb,prettyNames))]-.2,10^par('usr')[3],thisPos[max(grep(reb,prettyNames))]+.2,10^par('usr')[4],col='#00000022',border=NA)
      }
    }
  }
  points(thisPos[label[selector]]+spread[selector],ic50[selector],pch=speedPch[speed[selector]],bg=classCols[class[selector]],col=ifelse(speedPch[speed[selector]]>20,ifelse(class[selector] %in% deemphasize,'#00000066','#000000CC'),classCols[class[selector]]),lwd=ifelse(class[selector]=='Rebound',1.75,1.5),cex=ifelse(class[selector]=='QVOA',max(speedCex),speedCex[speed[selector]]))
}
}
