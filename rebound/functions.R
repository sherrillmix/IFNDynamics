plotQvoa2<-function(ic50,label,pos,class,study,speed,ylab='IFNa2 IC50 (pg/ml)',mar=c(6.9,4,.1,3.9),cex.axis=1.25,startDown=FALSE,pats=NULL,classCols,labelXAxis=TRUE,log='y',ylim=range(ic50)){
  speedPch<-c('Fast'=22,'Slow'=21,'Standard'=21,'Other'=21)
  speedCex<-c('Fast'=1,'Slow'=1,'Standard'=1,'Other'=1)
  deemphasize<-c("Acute","6 Month","Donor","Nadir","Last","Acute Recipients","Chronic Donors","Acute Recipient","Chronic Donor",'1 Year','Acute','Chronic')
  spread<-vipor::offsetX(log10(ic50),label,width=.45)
  ii<-'all'
  marSpace<-0
  if(ii!='all')selector<-label %in% subsets[[ii]]
  else selector<-rep(TRUE,length(label))
  par(mar=mar)
  plot(pos[label[selector]]+spread[selector],ic50[selector],log=log,yaxt='n',ylab=ylab,xlab='',xaxt='n',type='n',cex.lab=1.2,ylim=ylim,mgp=c(2.4,1,0),xlim=range(pos)+c(-1,1),xaxs='i')
  spread2<-ave(ic50,label,FUN=function(xx)beeswarm::swarmx(rep(0,length(xx)),xx,cex=.5,priority=ifelse(length(xx)>2,'density','ascending'))$x)
  spread<-ifelse(ave(ic50,label,FUN=length)>5,spread,spread2)
  if(ii=='all'&exists('subsets'))marSpaces<-sapply(subsets,function(xx)diff(convertUserToLine(1:0,2))*sum(!names(pos) %in% xx))
  if(!is.null(pats)){
    ranges<-tapply(pos[label],pats,range)
    #print(ranges[sapply(ranges,diff)>0.1])
    axFunc<-if(log=='y')function(xx)10^xx else function(xx)xx
    lapply(ranges[!is.na(sapply(ranges,diff))&sapply(ranges,diff)>0.1&!grepl('MM|^$',names(ranges))],function(xx)rect(min(xx)-.4,axFunc(par('usr')[3]),max(xx)+.4,axFunc(par('usr')[4]),col='#00000014',border=NA))
  }
  slantSelect<-!grepl('Outgrowth|VOA',names(pos)) &grepl('VOA MM|Outgrowth MM|Acute|6 Month|1 Year|Nadir|Last|Chronic',names(pos))
  #textOffsets<-c('Acute'=-.2,'6 Month'=-.2,'Nadir'=0,'Last'=.2,'Acute Recipients'=-.3,'Chronic Donors'=.3)[names(pos)[slantSelect]]
  if(any(slantSelect)){
    textOffsets<-c('Acute'=0.3,'6 Month'=0,'Nadir'=0,'Last'=0,'Chronic'=.3,'Acute Recipients'=.1,'Chronic Donors'=.1)[names(pos)[slantSelect]]
    textOffsets[is.na(textOffsets)]<-0
    slantAxis(1,pos[slantSelect],names(pos)[slantSelect],srt=-45,cex=cex.axis,location=.8,xpd=NA,textOffsets=textOffsets)
  }
  if(log=='y'){
    if(diff(par('usr')[3:4])>1)logAxis(las=1,mgp=c(1,.7,0))
    else axis(2,las=1,mgp=c(3,.5,0),tcl=-.3)
  }else{
    axis(2,las=1,mgp=c(3,.5,0),tcl=-.3)
  }
  abline(v=pos,col='#00000055',lty=3)
  points(pos[label[selector]]+spread[selector],ic50[selector],pch=speedPch[speed[selector]],bg=classCols[class[selector]],col=ifelse(speedPch[speed[selector]]>20,ifelse(class[selector] %in% deemphasize,'#00000066','#000000CC'),classCols[class[selector]]),lwd=ifelse(class[selector]=='Rebound',1,1),cex=ifelse(class[selector]=='QVOA',max(speedCex),speedCex[speed[selector]]))
  abline(v=pos['Acute Recipients']-(pos['Acute Recipients']-pos[which(names(pos)=='Acute Recipients')-1])/2,lty=1,col='#000000')
  #abline(v=pos['Acute']-(pos['Acute']-pos[which(names(pos)=='Acute')-1])/2,lty=1,col='#000000')
  #abline(v=pos['Outgrowth MM14']-(pos['Outgrowth MM14']-pos[which(names(pos)=='Outgrowth MM14')-1])/2,lty=1,col='#00000099')
  #abline(v=pos['Outgrowth MM14']-(pos['Outgrowth MM14']-pos[which(names(pos)=='Outgrowth MM14')-1])/2,lty=1,col='#00000099')
  #abline(v=pos['VOA MM14']-(pos['VOA MM14']-pos[which(names(pos)=='VOA MM14')-1])/2,lty=1,col='#000000')
  cols<-c('#00000033','#00000000')
  studyOrder<-unique(sapply(names(pos),function(xx)study[label==xx][1]))
  studyOrder<-studyOrder[studyOrder!='Transmission'&!is.na(studyOrder)]
  counter<-1
  for(ii in studyOrder){
    if(ii=='MM'){
      minPos<-pos[min(which(names(pos) %in% label[study==ii&class %in% c('VOA','QVOA','Outgrowth')]))]
      maxPos<-pos[max(which(names(pos) %in% label[study==ii&class %in% c('VOA','QVOA','Outgrowth')]))]
    }else{
      minPos<-pos[min(which(names(pos) %in% label[study==ii]))]
      maxPos<-pos[max(which(names(pos) %in% label[study==ii]))]
      nextPos<-pos[max(which(names(pos) %in% label[study==ii]))+1]
      if(ii!='Reservoir')abline(v=(maxPos+nextPos)/2,lty=2)
    }
    #rect(minPos-.5,10^par('usr')[3],maxPos+.5,10^par('usr')[4],col=cols[counter%%2+1],border=NA)
    counter<-counter+1
    if(is.na(minPos))browser()
    if(labelXAxis)axis(1,mean(c(minPos,maxPos)),sub('MM','Outgrowth',sub('ATI','Interrupt',sub('BEAT','IFNa2',sub('/','/\n',ii)))),padj=1,mgp=c(3,.2+2*(startDown+counter)%%2,0),tcl=-.7+-2*(startDown+counter)%%2,cex.axis=cex.axis)
  }
  return(list(ylim=ylim,pos=pos))
}

plotStudies<-function(study,label,ic50,pat,speed,class,ylab='IFNa2 IC50 (pg/ml)'){
spread<-vipor::offsetX(log10(ic50),label,width=.4,varwidth=TRUE)
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
