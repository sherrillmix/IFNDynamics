library(ggtree)
library(phangorn)
library(dnar)

if(!exists('meta'))source('readNewData.R')


pdf('out/trees.pdf',width=3,height=6)
for(ii in list.files('trees','ph$',full.name=TRUE)){
tree<-read.tree(ii)
if(grepl('MM23',ii))tree<-midpoint(tree)
times<-sapply(strsplit(tree$tip.label,'\\.'),'[[',2)
ids<-sapply(strsplit(tree$tip.label,'\\.'),function(xx)paste(xx[1:2],collapse='.'))
days<-as.numeric(meta[ids,'DFOSx'])
subject<-sapply(strsplit(tree$tip.label,'\\.'),'[[',1)
timeCols<-rainbow.lab(length(unique(times)))
names(timeCols)<-unique(times)[order(as.numeric(gsub('[^0-9]','',unique(times))))]
#,layout='unrooted'
print(ggtree(tree)+geom_tippoint(color=timeCols[times],show.legend=TRUE,size=1.5)+ggtitle(unique(subject))+theme(plot.title = element_text(hjust = .5,size=12,margin=margin(b=-10,unit='pt')))+geom_treescale(offset=-5,fontsize=3,x=.01,y=-8)+scale_color_manual('',breaks=names(timeCols),values=timeCols))
}
dev.off()
