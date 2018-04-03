facsFiles<-list.files('2018-04-02-iceTropism/jltr_2018-4-2/','CSV$',full.names=TRUE)
counts<-do.call(rbind,lapply(facsFiles,function(xx){out<-read.csv(xx,stringsAsFactors=FALSE);if(nrow(out)==0)return(NULL);out$file<-xx;out}))
counts$row<-sub('.*\\.([A-H])[01][0-9].CSV$','\\1',counts$file)
counts$col<-sub('.*\\.[A-H]([01][0-9]).CSV$','\\1',counts$file)
hist(counts[grep('H0[1-5].CSV',counts$file),'GRN.HLog'],breaks=100)
round(tapply(counts$GRN.HLog>3,list(counts$row,counts$col),function(xx)sum(xx)/length(xx))*1000,3)

