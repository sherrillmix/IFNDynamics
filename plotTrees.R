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
