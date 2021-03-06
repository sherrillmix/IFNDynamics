source('functions.R')
x<-read.csv('data/Test old vs new beta and 7 dose concentrations.csv',skip=44,nrows=16,stringsAsFactors=FALSE,check.names=FALSE)
y<-read.csv('data/Test old vs new beta and 7 dose concentrations.csv',skip=67,nrows=16,stringsAsFactors=FALSE,check.names=FALSE)
meltP24s<-function(xx){
  concs<-x[,2]
  xx<-xx[,-1:-2]
  out<-data.frame('concText'=rep(concs,ncol(xx)),'virus'=sub('\\.[0-9]$','',rep(colnames(xx),each=nrow(xx))),'beta'=rep(c('old','new'),each=nrow(xx)*ncol(xx)/2),'cells'=rep(c('old','new','old','new'),each=nrow(xx)*ncol(xx)/4),'techRep'=1:2,'p24'=unlist(xx),stringsAsFactors=FALSE)
  out$conc<-out$concText
  out[out$conc=='UT','conc']<-0
  out$conc<-as.numeric(out$conc)
  out
}
p24s<-rbind(cbind(meltP24s(x),'bioRep'=1),cbind(meltP24s(y),'bioRep'=2))

pdf('out/oldNewIc50.pdf')
ic50s<-do.call(rbind,lapply(unique(p24s$virus),function(virus){
  out<-do.call(rbind,lapply(unique(p24s$beta),function(beta){
    out<-do.call(rbind,lapply(unique(p24s$cells),function(cells){
      xx<-p24s[p24s$virus==virus&p24s$beta==beta&p24s$cells==cells,]
      fit<-plotIfnStacked(xx$conc,xx$p24,xx$bioRep,xx$techRep,xlab='IFNb concentration (pg/ml)',main=sprintf('%s %s cells %s beta',virus,cells,beta))
      means<-tapply(xx$p24,list(xx$bioRep,xx$conc),mean)
      out<-as.data.frame(t(calculateBasicIc50(as.numeric(colnames(means)),means=means)))
      out$cells<-cells
      out$lineInflect<-fit[1]
      out$line50<-fit[2]
      out
    }))
    out$beta<-beta
    out
  }))
  out$virus<-virus
  out
}))
dev.off()
print(ic50s[,c('virus','cells','beta','ic50','percVres','line50','lineInflect')],row.names=FALSE,digits=3)

pdf('out/oldNewIc50_lessConc.pdf',width=10)
par(mfrow=c(1,2))
ic50sLess<-do.call(rbind,lapply(unique(p24s$virus),function(virus){
  out<-do.call(rbind,lapply(unique(p24s$beta),function(beta){
    out<-do.call(rbind,lapply(unique(p24s$cells),function(cells){
      xx<-p24s[p24s$virus==virus&p24s$beta==beta&p24s$cells==cells& !p24s$concText %in% c('0.028','11'),]
      xx2<-p24s[p24s$virus==virus&p24s$beta==beta&p24s$cells==cells,]
      fit<-plotIfnStacked(xx$conc,xx$p24,xx$bioRep,xx$techRep,xlab='IFNb concentration (pg/ml)',main=sprintf('%s %s cells %s beta',virus,cells,beta))
      plotIfnStacked(xx2$conc,xx2$p24,xx2$bioRep,xx2$techRep,xlab='IFNb concentration (pg/ml)',main=sprintf('%s %s cells %s beta',virus,cells,beta))
      means<-tapply(xx$p24,list(xx$bioRep,xx$conc),mean)
      out<-as.data.frame(t(calculateBasicIc50(as.numeric(colnames(means)),means=means)))
      out$cells<-cells
      out$lineInflect<-fit[1]
      out$line50<-fit[2]
      out
    }))
    out$beta<-beta
    out
  }))
  out$virus<-virus
  out
}))
dev.off()
print(ic50sLess[,c('virus','cells','beta','ic50','percVres','line50','lineInflect')],row.names=FALSE,digits=3)
