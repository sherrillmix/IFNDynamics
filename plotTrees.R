library(ggtree)
library(phangorn)
library(dnar)

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


library(ape)
library(grid)
#pdf('out/voa_trees.pdf',width=3,height=6)
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
#dev.off()
