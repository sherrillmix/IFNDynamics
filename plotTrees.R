library(ggtree)
library(ggplot2)
library(phangorn)
library(dnar)
library(ape)
library(grid)

if(!exists('meta'))source('readNewData.R')


pdf('out/trees.pdf',width=3,height=6)
for(ii in list.files('trees','ph$',full.name=TRUE)){
  message(ii)
  tree<-read.tree(ii)
  if(grepl('MM23',ii))tree<-midpoint(tree)
  times<-sapply(strsplit(tree$tip.label,'\\.'),'[[',2)
  datIds<-sapply(strsplit(tree$tip.label,'\\.'),function(xx){
    if(!grepl('isol',xx[3])&&length(xx)!=3)return(NA) 
    id<-which(dat$pat==xx[1]&as.numeric(dat$visit)==as.numeric(xx[2])&dat$virusId==xx[length(xx)])
    if(length(id)!=1){
      warning('Problem finding ',paste(xx,collapse='.'))
      return(nrow(dat)+9999)
    }
    return(id)
  })
  datIds[is.na(datIds)]<-nrow(dat)+9999
  ifns<-dat[datIds,c('ic50','beta')]
  rownames(ifns)<-tree$tip.label
  ids<-sapply(strsplit(tree$tip.label,'\\.'),function(xx)paste(xx[1:2],collapse='.'))
  days<-as.numeric(meta[ids,'DFOSx'])
  subject<-sapply(strsplit(tree$tip.label,'\\.'),'[[',1)
  timeCols<-rainbow.lab(length(unique(times)))
  names(timeCols)<-unique(times)[order(as.numeric(gsub('[^0-9]','',unique(times))))]
  #,layout='unrooted'
  out<-ggtree(tree)+
    geom_tippoint(color=timeCols[times],show.legend=TRUE,size=1.5)+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')))+
    geom_treescale(offset=-5,fontsize=3,x=.01,y=-8)+
    scale_color_manual('',breaks=names(timeCols),values=timeCols)
  if(!grepl('MM23',ii))out<-out+scale_y_reverse()
  print(out)
  out<-ggtree(tree)
  out$data$ifn<-ifns[out$data$label,'ic50']
  out<-out+geom_tippoint(aes(color=ifn),show.legend=TRUE,size=1.5)+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')))+
    geom_treescale(offset=-5,fontsize=3,x=.01,y=-8)+
    scale_color_gradient(low = "blue", high = "red",na.value='#FFFFFF33')
  if(!grepl('MM23',ii))out<-out+scale_y_reverse()
  print(out)
}
dev.off()


pdf('out/voa_trees.pdf',width=3,height=6)
for(ii in list.files('trees/qvoa','ph$',full.name=TRUE)){
  message(ii)
  tree<-ape::read.tree(ii)
  days<-tree$tip.label
  regex<-'^[^.]+\\.d0*([1-9][0-9]+)\\..*$'
  isNum<-grepl(regex,days)
  days[isNum]<-sub(regex,'\\1',days[isNum])
  days[grepl('pbmc',days)]<-'PBMC'
  days[grepl('VOA|MM23.18',days)]<-'VOA'
  nDays<-max(as.numeric(days[isNum]))
  #colorFunc<-colorRampPalette(c('black','red'),space='Lab')
  colorFunc<-rainbow.lab
  timeCols<-c(colorFunc(nDays),'grey','purple')
  names(timeCols)<-c(1:nDays,'PBMC','VOA')
  subject<-sub('\\..*','',basename(ii))
  out<-ggtree(tree)+
    #magic numbers to get dash aligned with tip
    geom_tiplab(color=timeCols[days],size=7.5,label='-',vjust=.35)+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')))+
    geom_treescale(width=.01,offset=-5,fontsize=3,x=0,y=length(days)*.9)
    #scale_color_gradientn(colours = rainbow(5))
    #scale_fill_continuous(guide = guide_colorbar(direction = "horizontal"))
    #scale_color_manual('',breaks=names(timeCols),values=timeCols)
  vp<-viewport()
  plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  print(out,vp=vp)
  insetScale(c(1:nDays-.5,nDays+.5),timeCols[1:nDays],main='Days after onset of symptoms',insetPos = c(0.025, 0.3, 0.04, 0.95),at=c(1,500,1000,1500,2000))
  #plot.phylo(tree,show.tip.label=FALSE,tip.color=timeCols[days])
}
dev.off()

pdf('out/update_trees.pdf',width=2.,height=6)
for(ii in list.files('trees/updated2','ph$',full.name=TRUE)){
  message(ii)
  tree<-read.tree(ii)
  isHxb2<-grepl('HXB2',tree$tip.label)
  if(any(isHxb2)){
    message(ii,' has HXB2')
    tree<-root(tree,outgroup=tree$tip.label[isHxb2],resolve.root=TRUE)
    tree<-drop.tip(tree,tree$tip.label[isHxb2])
  }
  times<-sapply(strsplit(tree$tip.label,'\\.'),'[[',2)
  #datIds<-sapply(strsplit(tree$tip.label,'\\.'),function(xx){
    #if(!grepl('isol',xx[3])&&length(xx)!=3)return(NA) 
    #id<-which(dat$pat==xx[1]&as.numeric(dat$visit)==as.numeric(xx[2])&dat$virusId==xx[length(xx)])
    ##if(length(id)!=1){
      #warning('Problem finding ',paste(xx,collapse='.'))
      #return(nrow(dat)+9999)
    #}
    #return(id)
  #})
  #datIds[is.na(datIds)]<-nrow(dat)+9999
  ids<-sub('^(MM[0-9]+)[.-]([0-9]+)-.*','\\1.\\2',sub('^(Marvin_)?(VOA_)?','',sapply(strsplit(tree$tip.label,'\\.'),function(xx)paste(xx[1:2],collapse='.'))))
  #days<-as.numeric(meta[ids,'DFOSx'])
  days<-sapply(ids,function(xx)ifelse(any(compiledMeta$sample==xx),
      compiledMeta[compiledMeta$sample==xx,'time'],
      ifelse(grepl('MM[0-9]+[.-]d([0-9]+)$',xx),as.numeric(sub('MM[0-9]+[.-]d([0-9]+)$','\\1',xx)),NA)
  ))
  if(any(is.na(days))){
    warning('Removing tips ',paste(tree$tip.label[is.na(days)],collapse=', '))
    tree<-drop.tip(tree,tree$tip.label[is.na(days)])
    times<-times[!is.na(days)]
    days<-days[!is.na(days)]
  }
  isVoa<-grepl('VOA',tree$tip.label)
  isPbmc<-grepl('pbmc|PBMC',tree$tip.label)
  fakeDays<-min(c(0,days)):max(days[!isVoa&!isPbmc])
  cols<-rainbow.lab(length(fakeDays))
  timeCols<-cols[fakeDays %in% days[!isVoa&!isPbmc]]
  names(timeCols)<-sort(unique(days[!isVoa&!isPbmc]))
  voaDays<-unique(days[isVoa])
  pbmcDays<-unique(days[isPbmc])
  timeCols<-c(timeCols,c('voa'='purple','pbmc'='grey'))
  times<-as.character(ifelse(isPbmc,'pbmc',ifelse(isVoa,'voa',days)))
  subject<-dnar::mostAbundant(sapply(strsplit(tree$tip.label,'[-.]'),'[[',1))
  #timeCols<-rainbow.lab(length(unique(times)))
  #names(timeCols)<-unique(times)[order(as.numeric(gsub('[^0-9]','',unique(times))))]
  #,layout='unrooted'
  firsts<-which(days==min(days))
  firsts<-firsts[order(!grepl('consensus',tree$tip.label[firsts]),!grepl('[Vv]1',tree$tip.label[firsts]))]
  print(tree$tip.label[firsts])
  tree<-root(tree,tree$tip.label[firsts[1]])
  out<-ggtree(tree)+
    #geom_tippoint(color=timeCols[times],show.legend=TRUE,size=5/log(length(days)))+
    geom_tiplab(color=timeCols[times],size=7.5,label='-',vjust=.35)+
    geom_tiplab(color=timeCols[times],size=7.5,label=ifelse(times=='voa','  *',NA),vjust=.8)+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')),plot.margin=margin(2,0,20,0))+
    geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=0)+
    scale_x_continuous(limits = c(0,max(node.depth.edgelength(tree))*1.2))
    #scale_color_manual('',breaks=names(timeCols),values=timeCols)
  vp<-viewport()
  par(mar=c(0,0,0,0))
  plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  print(out,vp=vp)
  ticks<-pretty(fakeDays,n=3)
  par(lheight=.75)
  insetScale(c(fakeDays-.5,tail(fakeDays,1)+.5),cols,main='Days after onset\nof symptoms',insetPos = c(0.025, 0.1, 0.04, 0.9),at=ticks[ticks<max(fakeDays)],cex=.7)
  out<-ggtree(tree)+
    #geom_tippoint(color=timeCols[times],show.legend=TRUE,size=5/log(length(days)))+
    geom_tiplab(color=timeCols[times],size=1,vjust=.35)+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12),plot.margin=margin(2,0,20,0))+
    geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=-1)+
    scale_x_continuous(limits = c(0,max(node.depth.edgelength(tree))*2))
    #scale_color_manual('',breaks=names(timeCols),values=timeCols)
  print(out)
}
dev.off()

recomb<-read.csv('lindsey/recombinant.csv',stringsAsFactors=FALSE)
reroot<-c('MM23'="MM23.PLAS.ISO.00009.002.01",'MM14'='MM14.PLAS.SGS.00040.007.01')
pdf('out/lindsey_trees.pdf',width=2.,height=6)
for(ii in list.files('lindsey/rooted/','MM[0-9]+.bs100.rooted.nex$',full.name=TRUE)){
  message(ii)
  tree<-read.nexus(ii)
  tree$tip.label<-gsub("'",'',tree$tip.label)
  days<-as.numeric(sapply(strsplit(tree$tip.label,'\\.'),'[[',4))
  ids<-sapply(strsplit(tree$tip.label,'\\.'),function(xx)paste(xx[2:3],collapse='.'))
  print(unique(ids))
  isVoa<-grepl('VOA',ids)
  isSgs<-grepl('SGS',ids)
  isPbmc<-grepl('PBMC',ids)&!isVoa
  fakeDays<-min(c(0,days)):max(days[!isVoa&!isPbmc])
  cols<-rainbow.lab(length(fakeDays))
  timeCols<-cols[fakeDays %in% days[!isVoa&!isPbmc]]
  names(timeCols)<-sort(unique(days[!isVoa&!isPbmc]))
  timeCols<-c(timeCols,c('voa'='purple','pbmc'='grey'))
  times<-as.character(ifelse(isPbmc,'pbmc',ifelse(isVoa,'voa',days)))
  subject<-sub("'",'',dnar::mostAbundant(sapply(strsplit(tree$tip.label,'[-.]'),'[[',1)))
  firsts<-which(days==min(days))
  firsts<-firsts[order(!grepl('consensus',tree$tip.label[firsts]),!grepl('[Vv]1',tree$tip.label[firsts]))]
  #lindsey order
  if(subject %in% names(reroot)){
    tree<-root(tree,reroot[subject])
  }
  out<-ggtree(tree)+
    #geom_tippoint(color=timeCols[times],show.legend=TRUE,size=5/log(length(days)))+
    geom_tiplab(color=timeCols[times],size=7.5,label='-',vjust=1.5)+
    geom_tiplab(color=timeCols[times],size=7.5,label=ifelse(times=='voa','  *',NA),vjust=.8)+
    geom_tiplab(color=timeCols[times],size=3,label=ifelse(tree$tip.label %in% recomb$Sequence,'         r',NA))+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')),plot.margin=margin(2,0,20,0))+
    geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=0)+
    scale_x_continuous(limits = c(0,max(node.depth.edgelength(tree))*1.2))
    #scale_color_manual('',breaks=names(timeCols),values=timeCols)
  vp<-viewport()
  par(mar=c(0,0,0,0))
  plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  print(out,vp=vp)
  ticks<-pretty(fakeDays,n=3)
  par(lheight=.75)
  insetScale(c(fakeDays-.5,tail(fakeDays,1)+.5),cols,main='Days after onset\nof symptoms',insetPos = c(0.025, 0.1, 0.04, 0.9),at=ticks[ticks<max(fakeDays)],cex=.7)
  out<-ggtree(tree)+
    #geom_tippoint(color=timeCols[times],show.legend=TRUE,size=5/log(length(days)))+
    geom_tiplab(color=timeCols[times],size=1,vjust=.35)+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12),plot.margin=margin(2,0,20,0))+
    geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=-1)+
    scale_x_continuous(limits = c(0,max(node.depth.edgelength(tree))*2))
    #scale_color_manual('',breaks=names(timeCols),values=timeCols)
  print(out)
}
dev.off()

censorLength<-.15
noCensor<-c('MM23')
#recomb<-data.frame('Sequence'=NULL)
recomb<-read.csv('lindsey/recombinant.csv',stringsAsFactors=FALSE)
reroot<-c('MM14'='MM14.PLAS.ISO.00040.001.01','MM15'='MM15.PLAS.SGS.00027.001.01')
#rotate<-list('MM39'=c(99,76))
maxDepth<-apply(do.call(rbind,lapply(list.files('trees/10trees','\\.txt$',full.name=TRUE),function(xx){
    tree<-read.tree(xx)
    out<-rep(NA,2)
    out[1]<-max(node.depth.edgelength(tree))
    tree$edge.length[tree$edge.length>censorLength]<-censorLength
    out[2]<-max(node.depth.edgelength(tree))
    return(out)
}
)),2,max)
print(maxDepth)
makeLegend<-function(out,tree,insetPos=c(0.025, 0.1, 0.04, 0.7)){
  vp<-viewport()
  par(mar=c(0,0,0,0))
  plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  print(out,vp=vp)
  ticks<-pretty(fakeDays,n=3)
  par(lheight=.75)
  insetScale(c(fakeDays-.5,tail(fakeDays,1)+.5),cols,main='Days after onset\nof symptoms',insetPos = insetPos,at=ticks[ticks<max(fakeDays)],cex=.7)
  if(any(grepl('PBMC',tree$tip.label))){
    selector<-c(any(grepl('PBMC\\.VOA',tree$tip.label)),any(grepl('PBMC\\.SGS',tree$tip.label)))
    legend(grconvertX(insetPos[4]+.02,'npc','user'),grconvertY(mean(c(.025,.04)),'npc','user'),c('VOA isolate','PBMC sequence')[selector],fill=c('purple','grey')[selector],cex=.65,yjust=.4,bty='n',title='Post-ART')
  }
}
pdf('out/lindsey_trees2.pdf',width=4,height=6)
for(ii in list.files('trees/10trees','\\.txt$',full.name=TRUE)){
  message(ii)
  tree<-read.tree(ii)
  #if(grepl('MM34',ii))tree<-drop.tip(tree,which(node.depth.edgelength(tree)>.5))
  #which(node.depth.edgelength(tree)>.5)
  #which(tree$edge.length>.5)
  tree$tip.label<-gsub("'",'',tree$tip.label)
  days<-as.numeric(sapply(strsplit(tree$tip.label,'\\.'),'[[',4))
  ids<-sapply(strsplit(tree$tip.label,'\\.'),function(xx)paste(xx[2:3],collapse='.'))
  print(unique(ids))
  isVoa<-grepl('VOA',ids)
  isSgs<-grepl('SGS',ids)
  isPbmc<-grepl('PBMC',ids)&!isVoa
  fakeDays<-min(c(0,days)):max(days[!isVoa&!isPbmc])
  cols<-rainbow.lab(length(fakeDays))
  timeCols<-cols[fakeDays %in% days[!isVoa&!isPbmc]]
  names(timeCols)<-sort(unique(days[!isVoa&!isPbmc]))
  timeCols<-c(timeCols,c('voa'='purple','pbmc'='grey'))
  times<-as.character(ifelse(isPbmc,'pbmc',ifelse(isVoa,'voa',days)))
  subject<-sub("'",'',dnar::mostAbundant(sapply(strsplit(tree$tip.label,'[-.]'),'[[',1)))
  firsts<-which(days==min(days))
  firsts<-firsts[order(!grepl('consensus',tree$tip.label[firsts]),!grepl('[Vv]1',tree$tip.label[firsts]))]
  #lindsey order
  if(subject %in% names(reroot)){
    tree<-root(tree,reroot[subject])
  }
  #if(subject %in% names(rotate)){
  #  tree<-reorder(ladderize(tree,FALSE))
  #  tree<-rotate(tree,rotate[[subject]])
  #}
  edge<-data.frame(tree$edge,edge_num=1:nrow(tree$edge))
  colnames(edge)<-c('parent','node','edge_num')
  lengthBak<-tree$edge.length
  if(!unique(subject) %in% noCensor)tree$edge.length[tree$edge.length>censorLength]<-censorLength
  if(any(tree$tip.label %in% recomb$Sequence))message('RECOMB ',paste(tree$tip.label[tree$tip.label %in% recomb$Sequence],collapse=', '))
  out<-ggtree(tree,size=.2)+#,ladderize=FALSE)+ 
    #geom_tippoint(color=timeCols[times],show.legend=TRUE,size=5/log(length(days)))+
    geom_tiplab(color=timeCols[times],size=min(7.5,1500/length(tree$tip.label)),label='-',vjust=.325,hjust=.1)+
    #geom_tiplab(color=timeCols[times],size=7.5,label=ifelse(times=='voa','  *',NA),vjust=.8)+
    #geom_tiplab(color=timeCols[times],size=3,label=ifelse(tree$tip.label %in% recomb$Sequence,'         r',NA))+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')),plot.margin=margin(2,0,20,0))+
    geom_treescale(offset=-length(days)/55,fontsize=2.5,x=.01,y=-.01*length(tree$tip.label),width=.05)+
    scale_x_continuous(limits = c(0,maxDepth[2]*1.2))
    #geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=-.01*length(tree$tip.label),width=.05)+
  censoredEdges<-out$data[out$data$branch.length==censorLength,]
  if(nrow(censoredEdges)>0){
    message('CENSORING')
    slash<-.005
    left<-censoredEdges$x-censorLength/20-censorLength/2
    right<-censoredEdges$x+censorLength/20-censorLength/2
    out<-out+
      annotate('rect',xmin=left,xmax=right,ymin=censoredEdges$y-.4,ymax=censoredEdges$y+.4,fill='white')+
      annotate('segment',x=c(left,right)-slash,xend=c(left,right)+slash,y=rep(censoredEdges$y,2)-.7,yend=rep(censoredEdges$y,2)+.7,size=.2)
  }
  #out<-out %<+% edge + geom_point(aes(x=branch,y=edge_num))
    #scale_x_continuous(limits = c(0,max(node.depth.edgelength(tree))*1.2))
    #scale_color_manual('',breaks=names(timeCols),values=timeCols)
  #out$data[out$data$branch.length>censorLength,'x']<-out$data[out$data$branch.length>censorLength,'x']-out$data[out$data$branch.length>censorLength,'branch.length']+censorLength
  makeLegend(out,tree,insetPos=c(0.035, 0.05, 0.05, 0.45))
  #tree$edge.length<-lengthBak
  out<-ggtree(tree,size=.2)+
    #geom_tippoint(color=timeCols[times],show.legend=TRUE,size=5/log(length(days)))+
    #1.4/(max(1,length(tree$tip.label)/17))^.4
    geom_tiplab(color=timeCols[times],size=min(1.5,100/length(tree$tip.label)),label=sprintf('%s%s',tree$tip.label,ifelse(tree$tip.label %in% recomb$Sequence,'','')),vjust=.5)+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12),plot.margin=margin(2,0,20,0))+
    #geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=-1)+
    geom_tiplab(color=timeCols[times],size=3,label=ifelse(tree$tip.label %in% recomb$Sequence,'*',NA),vjust=.8,offset=ifelse(unique(subject)=='MM14',0.15,ifelse(unique(subject)=='MM23',.073,0))*strwidth(tree$tip.label)/strwidth(tree$tip.label[1]))+ #magic numbers
    geom_treescale(offset=-length(days)/55,fontsize=2.5,x=.01,y=-.01*length(tree$tip.label),width=.05)+
    scale_x_continuous(limits = c(0,maxDepth[2]*1.2))
    if(nrow(censoredEdges)>0){
      message('CENSORING')
      out<-out+
        annotate('rect',xmin=left,xmax=right,ymin=censoredEdges$y-.6,ymax=censoredEdges$y+.6,fill='white')+
        annotate('segment',x=c(left,right),xend=c(left,right),y=rep(censoredEdges$y,2)-.6,yend=rep(censoredEdges$y,2)+.6,size=.2)
    }
    #scale_x_continuous(limits = c(0,max(node.depth.edgelength(tree))*2))+
    #geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=-1,width=.05)
    #geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=5)
    #scale_color_manual('',breaks=names(timeCols),values=timeCols)
  makeLegend(out,tree)
  #display tree with equal length branches
  if(FALSE){
    tmp<-tree
    tmp$edge.length[]<-1
    tmp$node.label<-1:length(tmp$node.label)
    out<-ggtree(tmp)+
      geom_tippoint(color=timeCols[times],show.legend=TRUE,size=5/log(length(days)))+
      geom_tiplab(color=timeCols[times],size=1,vjust=.35)+
      ggtitle(unique(subject))+
      theme(plot.title = element_text(hjust = .5,size=12),plot.margin=margin(2,0,20,0))+
      geom_nodelab(size=2,col='red')
    print(out)
  }
}
dev.off()
system('pdftk out/lindsey_trees2.pdf cat 1 4 7 13 19 output tmp.pdf')
system('pdfjam tmp.pdf --nup 5x1 --landscape --outfile out/voaStackedTrees.pdf')
system('pdftk out/lindsey_trees2.pdf cat 2 4 6 8 10 12 14 16 18 20 output out/separateTrees.pdf')
system('pdfjam out/separateTrees.pdf --nup 5x2 --landscape --outfile out/allStackedTrees.pdf')
#file.copy('out/allStackedTrees.pdf','out/Fig. S2.pdf')
file.copy('out/separateTrees.pdf','out/Fig. S2.pdf',overwrite=TRUE)
#system('pdftk out/lindsey_trees2.pdf cat 9 output out/Fig._4.pdf') #tweaked in inkscape

zz<-lapply(list.files('trees/10trees','\\.txt$',full.name=TRUE),function(ii)read.tree(ii)$tip.label)
lapply(zz,function(zz)zz[grep('PBMC',zz)])
#note that 430 is VOA not SGS
#1512,641,1152,997,647,736,430

pdf('out/treeMM34.pdf',width=3.5,height=5)
  ii<-'trees/10trees/MM34.bs100.txt'
  message(ii)
  tree<-read.tree(ii)
  tree$tip.label<-gsub("'",'',tree$tip.label)
  days<-as.numeric(sapply(strsplit(tree$tip.label,'\\.'),'[[',4))
  ids<-sapply(strsplit(tree$tip.label,'\\.'),function(xx)paste(xx[2:3],collapse='.'))
  print(unique(ids))
  isVoa<-grepl('VOA',ids)
  isSgs<-grepl('SGS',ids)
  isPbmc<-grepl('PBMC',ids)&!isVoa
  fakeDays<-min(c(0,days)):max(days[!isVoa&!isPbmc])
  cols<-rainbow.lab(length(fakeDays))
  timeCols<-cols[fakeDays %in% days[!isVoa&!isPbmc]]
  names(timeCols)<-sort(unique(days[!isVoa&!isPbmc]))
  timeCols<-c(timeCols,c('voa'='purple','pbmc'='grey'))
  times<-as.character(ifelse(isPbmc,'pbmc',ifelse(isVoa,'voa',days)))
  subject<-sub("'",'',dnar::mostAbundant(sapply(strsplit(tree$tip.label,'[-.]'),'[[',1)))
  firsts<-which(days==min(days))
  firsts<-firsts[order(!grepl('consensus',tree$tip.label[firsts]),!grepl('[Vv]1',tree$tip.label[firsts]))]
  #lindsey order
  if(subject %in% names(reroot)){
    tree<-root(tree,reroot[subject])
  }
  edge<-data.frame(tree$edge,edge_num=1:nrow(tree$edge))
  colnames(edge)<-c('parent','node','edge_num')
  lengthBak<-tree$edge.length
  if(!unique(subject) %in% noCensor)tree$edge.length[tree$edge.length>censorLength]<-censorLength
  out<-ggtree(tree,size=.2)+#,ladderize=FALSE)+ 
    #geom_tippoint(color=timeCols[times],show.legend=TRUE,size=5/log(length(days)))+
    geom_tiplab(color=timeCols[times],size=min(7.5,1500/length(tree$tip.label)),label='-',vjust=.35)+
    #geom_tiplab(color=timeCols[times],size=7.5,label=ifelse(times=='voa','  *',NA),vjust=.8)+
    geom_tiplab(color=timeCols[times],size=3,label=ifelse(tree$tip.label %in% recomb$Sequence,'         r',NA))+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')),plot.margin=margin(2,0,20,0))+
    geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=-.01*length(tree$tip.label),width=.05)#+
    #scale_x_continuous(limits = c(0,maxDepth[2]*1.2))
    #geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=-.01*length(tree$tip.label),width=.05)+
  censoredEdges<-out$data[out$data$branch.length==censorLength,]
  if(nrow(censoredEdges)>0){
    message('CENSORING')
    left<-censoredEdges$x-censorLength/20-censorLength/2
    right<-censoredEdges$x+censorLength/20-censorLength/2
    slant<-.003
    out<-out+
      annotate('rect',xmin=left,xmax=right,ymin=censoredEdges$y-.6,ymax=censoredEdges$y+.6,fill='white')+
      #annotate('segment',x=c(left,right),xend=c(left,right),y=rep(censoredEdges$y,2)-.6,yend=rep(censoredEdges$y,2)+.6,size=.2)
      annotate('segment',x=c(left-slant,right-slant),xend=c(left+slant,right+slant),y=rep(censoredEdges$y,2)-.8,yend=rep(censoredEdges$y,2)+.8,size=.2)
  }
  makeLegend(out,tree,insetPos=c(0.025, 0.1, 0.04, 0.6))
dev.off()
#file.copy('out/treeMM34.pdf','out/Fig._4.pdf',TRUE)

pdf('out/lindsey_treesCircular.pdf',width=4,height=6)
for(ii in list.files('trees/10trees','\\.txt$',full.name=TRUE)){
  message(ii)
  tree<-read.tree(ii)
  #if(grepl('MM34',ii))tree<-drop.tip(tree,which(node.depth.edgelength(tree)>.5))
  #which(node.depth.edgelength(tree)>.5)
  #which(tree$edge.length>.5)
  tree$tip.label<-gsub("'",'',tree$tip.label)
  days<-as.numeric(sapply(strsplit(tree$tip.label,'\\.'),'[[',4))
  ids<-sapply(strsplit(tree$tip.label,'\\.'),function(xx)paste(xx[2:3],collapse='.'))
  print(unique(ids))
  isVoa<-grepl('VOA',ids)
  isSgs<-grepl('SGS',ids)
  isPbmc<-grepl('PBMC',ids)&!isVoa
  fakeDays<-min(c(0,days)):max(days[!isVoa&!isPbmc])
  cols<-rainbow.lab(length(fakeDays))
  timeCols<-cols[fakeDays %in% days[!isVoa&!isPbmc]]
  names(timeCols)<-sort(unique(days[!isVoa&!isPbmc]))
  timeCols<-c(timeCols,c('voa'='purple','pbmc'='grey'))
  times<-as.character(ifelse(isPbmc,'pbmc',ifelse(isVoa,'voa',days)))
  subject<-sub("'",'',dnar::mostAbundant(sapply(strsplit(tree$tip.label,'[-.]'),'[[',1)))
  firsts<-which(days==min(days))
  firsts<-firsts[order(!grepl('consensus',tree$tip.label[firsts]),!grepl('[Vv]1',tree$tip.label[firsts]))]
  #lindsey order
  if(subject %in% names(reroot)){
    tree<-root(tree,reroot[subject])
  }
  #if(subject %in% names(rotate)){
  #  tree<-reorder(ladderize(tree,FALSE))
  #  tree<-rotate(tree,rotate[[subject]])
  #}
  edge<-data.frame(tree$edge,edge_num=1:nrow(tree$edge))
  colnames(edge)<-c('parent','node','edge_num')
  lengthBak<-tree$edge.length
  if(!unique(subject) %in% noCensor)tree$edge.length[tree$edge.length>censorLength]<-censorLength
  if(any(tree$tip.label %in% recomb$Sequence))message('RECOMB ',paste(tree$tip.label[tree$tip.label %in% recomb$Sequence],collapse=', '))
  out<-ggtree(tree,layout='circular',size=.2)+#,ladderize=FALSE)+ 
    geom_tippoint(color=timeCols[times],show.legend=TRUE,size=.1/log(length(days)))+
    #geom_tiplab(color=timeCols[times],size=min(7.5,1500/length(tree$tip.label)),label='-',vjust=.325,hjust=.1)+
    #geom_tiplab(color=timeCols[times],size=7.5,label=ifelse(times=='voa','  *',NA),vjust=.8)+
    #geom_tiplab(color=timeCols[times],size=3,label=ifelse(tree$tip.label %in% recomb$Sequence,'         r',NA))+
    ggtitle(unique(subject))+
    theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')),plot.margin=margin(2,0,2,0))
    #geom_treescale(offset=-length(days)/55,fontsize=2.5,x=.01,y=-.01*length(tree$tip.label),width=.05)+
    #scale_x_continuous(limits = c(0,maxDepth[2]*1.2))
    #geom_treescale(offset=-length(days)/90,fontsize=2.5,x=.01,y=-.01*length(tree$tip.label),width=.05)+
  censoredEdges<-out$data[out$data$branch.length==censorLength,]
  if(nrow(censoredEdges)>0){
    message('CENSORING')
    slash<-.005
    left<-censoredEdges$x-censorLength/20-censorLength/2
    right<-censoredEdges$x+censorLength/20-censorLength/2
    out<-out+
      annotate('rect',xmin=left,xmax=right,ymin=censoredEdges$y-.4,ymax=censoredEdges$y+.4,fill='white')+
      annotate('segment',x=c(left,right)-slash,xend=c(left,right)+slash,y=rep(censoredEdges$y,2)-.7,yend=rep(censoredEdges$y,2)+.7,size=.2)
  }
  #out<-out %<+% edge + geom_point(aes(x=branch,y=edge_num))
    #scale_x_continuous(limits = c(0,max(node.depth.edgelength(tree))*1.2))
    #scale_color_manual('',breaks=names(timeCols),values=timeCols)
  #out$data[out$data$branch.length>censorLength,'x']<-out$data[out$data$branch.length>censorLength,'x']-out$data[out$data$branch.length>censorLength,'branch.length']+censorLength
  print(out)
  #makeLegend(out,tree,insetPos=c(0.035, 0.05, 0.05, 0.45))
}
dev.off()

