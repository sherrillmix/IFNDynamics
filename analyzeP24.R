library(xlsx)
library(abind)
readP24<-function(file,exclude=c(),firstCol=2){
  wb <- loadWorkbook(file)
  counter<-0
  p24s<-lapply(getSheets(wb),function(sheet){
    counter<<-counter+1
    if(counter %in% exclude)return(NULL)
    message(' Sheet ',counter)
    #can get weird error otherwise
    #rows<-lapply(1:40,function(xx)tryCatch(getCells(getRows(sheet)[xx])[[1]],error=function(e)return(NULL)))
    if(length(getRows(sheet))==0)return(NULL)
    rows<-getCells(getRows(sheet,1:40),1:50,simplify=FALSE)
    vals<-lapply(rows,function(row){
      tmp<-sapply(row,function(xx)ifelse(is.null(xx),NA,getCellValue(xx)))
      out<-rep(NA,50)
      names(out)<-1:50
      out[names(tmp)[names(tmp)!='']]<-tmp[names(tmp)!='']
      return(out)
    })
    if(any(sapply(vals[1:10],function(xx)is.null(xx[[3]])))){message('  Not found');return(NULL)}
    #assuming col 3, first row contains 1
    firstRow<-min(c(Inf,which(!is.na(sapply(vals[1:10],'[[',firstCol))&sapply(vals[1:10],'[[',firstCol)=='1')))+1
    if(firstRow==Inf){message('  Not found');return(NULL)}
    lastRow<-firstRow+min(c(Inf,which(sapply(vals[firstRow:40],function(xx)is.null(xx)||all(is.na(xx[firstCol+0:5])|is.null(xx[firstCol+0:5]))))))-1-1
    if(lastRow==Inf)return(NULL)
    #if((lastRow-firstRow+1) %%2 !=0)warning('Number of rows not a multiple of 2 on sheet ',counter,' of ',file)
    dat<-as.data.frame(apply(do.call(rbind,lapply(vals[firstRow:lastRow],function(zz)zz[firstCol+0:23])),2,function(xx){as.numeric(ifelse(xx==''|is.na(xx)|grepl('\\?\\?',xx),NA,sub('[><]','',xx)))}),stringsAsFactors=FALSE)
    colnames(dat)<-1:ncol(dat)
    rownames(dat)<-LETTERS[1:nrow(dat)]
    return(dat)
  })
  names(p24s)<-names(getSheets(wb))
  p24s
}

processP24<-function(p24){
  reps<-abind(t(p24[,seq(1,ncol(p24),2)]),t(p24[,seq(2,ncol(p24),2)]),along=3)
  nCol<-ncol(reps)/2
  nRow<-nrow(reps)/2
  quads<-expand.grid(row=c(1,nRow+1),col=c(1,nCol+1))
  splits<-lapply(1:nrow(quads),function(xx){
    colRow<-unlist(quads[xx,])
    out<-reps[colRow['row']:(colRow['row']+nRow-1),colRow['col']:(colRow['col']+nCol-1),][,nCol:1,]
    colnames(out)<-1:ncol(out)
    rownames(out)<-LETTERS[1:nrow(out)]
    return(out)
  })
  return(splits)
}

condenseP24<- function(xx,name=NULL){
  out<-data.frame('row'=rep(rownames(xx),ncol(xx)),'col'=rep(colnames(xx),each=nrow(xx)),'p24_1'=as.vector(xx[,,1]),'p24_2'=as.vector(xx[,,2]),stringsAsFactors=FALSE)
  if(!is.null(name)){
    out$plate<-name
    out<-out[,c('plate',colnames(out)[colnames(out)!='plate'])]
  }
  return(out)
}

p24s<-readP24('data/p24 9.6.2018 isolation rebounds.xlsx')
splits<-lapply(p24s,processP24)
nameSplits<-strsplit(sub('12wells 24wells','12wells_24wells',names(p24s)),' ')
allPlates<-unlist(splits,recursive=FALSE)
names(allPlates)<-unlist(nameSplits)
plateDf<-do.call(rbind,mapply(condenseP24,allPlates,names(allPlates),SIMPLIFY=FALSE))
plateDf$meanP24<-(plateDf$p24_1+plateDf$p24_2)/2
plateDf$topBottom<-ifelse(plateDf$row %in% LETTERS[1:3],'top',ifelse(plateDf$row %in% LETTERS[4:6],'bottom',NA))
plateDf$study<-NA
plateDf$sample<-NA


dat829<-data.frame('Plate'=rep(1:9,each=2),'TopBottom'=rep(c('Top','Bottom'),9),'Study'=rep(c('Katie rebound','Neussensweig'),c(6,12)),'Sample'=c('A01.30','A02.17','A07.15','A08.21','A13.26','A14.2','601','602','609','610','9201','9202','9203','9207','9208','9209','601','9203'),stringsAsFactors=FALSE)
dat830<-read.csv('data/2018-08-30_isolation.csv',skip=1,stringsAsFactors=FALSE)
for(ii in 1:nrow(dat830)){
  selector<-plateDf$topBottom==tolower(dat830[ii,'TopBottom']) & plateDf$plate==sprintf('30-%d',dat830[ii,'Plate'])
  plateDf[selector,c('study','sample')]<-dat830[ii,c('Study','Sample')]
}
for(ii in 1:nrow(dat829)){
  selector<-plateDf$topBottom==tolower(dat829[ii,'TopBottom']) & plateDf$plate==sprintf('29-%d',dat829[ii,'Plate'])
  plateDf[selector,c('study','sample')]<-dat829[ii,c('Study','Sample')]
}
plateDf<-plateDf[order(structure(1:length(allPlates),.Names=names(allPlates))[plateDf$plate],plateDf$topBottom=='bottom',plateDf$row,plateDf$col),]
write.csv(plateDf,'out/p24_2018-09-06.csv',row.names=FALSE)


p24s2<-readP24('data/p24 9.14.2018 isolation rebounds.xlsx')
splits<-lapply(p24s2,processP24)
nameSplits<-strsplit(sub('12wells 24wells','12wells_24wells',names(p24s2)),' ')
allPlates<-unlist(splits,recursive=FALSE)
names(allPlates)<-unlist(nameSplits)
plateDf2<-do.call(rbind,mapply(condenseP24,allPlates,names(allPlates),SIMPLIFY=FALSE))
plateDf2$meanP24<-(plateDf2$p24_1+plateDf2$p24_2)/2
plateDf2$topBottom<-ifelse(plateDf2$row %in% LETTERS[1:3],'top',ifelse(plateDf2$row %in% LETTERS[4:6],'bottom',NA))
plateDf2$study<-NA
plateDf2$sample<-NA


for(ii in 1:nrow(dat830)){
  selector<-plateDf2$topBottom==tolower(dat830[ii,'TopBottom']) & plateDf2$plate==sprintf('30-%d',dat830[ii,'Plate'])
  plateDf2[selector,c('study','sample')]<-dat830[ii,c('Study','Sample')]
}
for(ii in 1:nrow(dat829)){
  selector<-plateDf2$topBottom==tolower(dat829[ii,'TopBottom']) & plateDf2$plate==sprintf('29-%d',dat829[ii,'Plate'])
  plateDf2[selector,c('study','sample')]<-dat829[ii,c('Study','Sample')]
}
plateDf2<-plateDf2[order(structure(1:length(allPlates),.Names=names(allPlates))[plateDf2$plate],plateDf2$topBottom=='bottom',plateDf2$row,plateDf2$col),]
write.csv(plateDf2,'out/p24_2018-09-14.csv',row.names=FALSE)

p24s3<-readP24('data/p24 9.21.2018 isolation rebounds.xlsx')
splits<-lapply(p24s3,processP24)
nameSplits<-strsplit(sub('12wells 24wells','12wells_24wells',names(p24s3)),' ')
allPlates<-unlist(splits,recursive=FALSE)
names(allPlates)<-unlist(nameSplits)
plateDf3<-do.call(rbind,mapply(condenseP24,allPlates,names(allPlates),SIMPLIFY=FALSE))
plateDf3$meanP24<-(plateDf3$p24_1+plateDf3$p24_2)/2
plateDf3$topBottom<-ifelse(plateDf3$row %in% LETTERS[1:3],'top',ifelse(plateDf3$row %in% LETTERS[4:6],'bottom',NA))
plateDf3$study<-NA
plateDf3$sample<-NA
for(ii in 1:nrow(dat830)){
  selector<-plateDf3$topBottom==tolower(dat830[ii,'TopBottom']) & plateDf3$plate==sprintf('30-%d',dat830[ii,'Plate'])
  plateDf3[selector,c('study','sample')]<-dat830[ii,c('Study','Sample')]
}
for(ii in 1:nrow(dat829)){
  selector<-plateDf3$topBottom==tolower(dat829[ii,'TopBottom']) & plateDf3$plate==sprintf('29-%d',dat829[ii,'Plate'])
  plateDf3[selector,c('study','sample')]<-dat829[ii,c('Study','Sample')]
}
plateDf3<-plateDf3[order(structure(1:length(allPlates),.Names=names(allPlates))[plateDf3$plate],plateDf3$topBottom=='bottom',plateDf3$row,plateDf3$col),]
plateDf3<-merge(plateDf3,expand[,c('plate','row','col','expansion')],all.x=TRUE,all.y=TRUE)
write.csv(plateDf3,'out/p24_2018-09-21.csv',row.names=FALSE)
write.csv(plateDf3[plateDf3$meanP24>2000&is.na(plateDf3$expansion),],'expand.csv')




expand<-read.csv('data/Expansions.csv',stringsAsFactors=FALSE,skip=1,header=TRUE)
expand<-expand[expand$row %in% LETTERS,]
colnames(expand)[1]<-'expansion'

merged<-merge(plateDf,plateDf2,by=c('plate','row','col','study','sample','topBottom'),suffixes=c('_9-6','_9-14'),all.x=TRUE,all.y=TRUE)
tmp<-plateDf3
colnames(tmp)[colnames(tmp)==c('p24_1','p24_2','meanP24')]<-c('p24_1_9-21','p24_2_9-21','meanP24_9-21')
merged<-merge(merged,tmp,by=c('plate','row','col','study','sample','topBottom'),all.x=TRUE,all.y=TRUE)
write.csv(merged,'out/p24_merged.csv')

p24s4<-readP24('data/p24 9.26.2018 isolation rebounds.xlsx')
#sheets not named this time
names(p24s4)<-names(p24s3)
splits<-lapply(p24s4,processP24)
nameSplits<-strsplit(sub('12wells 24wells','12wells_24wells',names(p24s4)),' ')
allPlates<-unlist(splits,recursive=FALSE)
names(allPlates)<-unlist(nameSplits)
plateDf4<-do.call(rbind,mapply(condenseP24,allPlates,names(allPlates),SIMPLIFY=FALSE))
plateDf4$meanP24<-(plateDf4$p24_1+plateDf4$p24_2)/2
plateDf4$topBottom<-ifelse(plateDf4$row %in% LETTERS[1:3],'top',ifelse(plateDf4$row %in% LETTERS[4:6],'bottom',NA))
plateDf4$study<-NA
plateDf4$sample<-NA
for(ii in 1:nrow(dat830)){
  selector<-plateDf4$topBottom==tolower(dat830[ii,'TopBottom']) & plateDf4$plate==sprintf('30-%d',dat830[ii,'Plate'])
  plateDf4[selector,c('study','sample')]<-dat830[ii,c('Study','Sample')]
}
for(ii in 1:nrow(dat829)){
  selector<-plateDf4$topBottom==tolower(dat829[ii,'TopBottom']) & plateDf4$plate==sprintf('29-%d',dat829[ii,'Plate'])
  plateDf4[selector,c('study','sample')]<-dat829[ii,c('Study','Sample')]
}
plateDf4<-plateDf4[order(structure(1:length(allPlates),.Names=names(allPlates))[plateDf4$plate],plateDf4$topBottom=='bottom',plateDf4$row,plateDf4$col),]
write.csv(plateDf4,'out/p24_2018-09-26.csv',row.names=FALSE)
expands<-read.csv('data/Expansions2.csv',stringsAsFactors=FALSE)
targets<-plateDf4[plateDf4$meanP24>2000,]
targets<-merge(targets,expands[,c('plate','row','col','id')],all.x=TRUE)
write.csv(targets,'out/expansion_9-27.csv',row.names=FALSE)

