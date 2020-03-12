if(!exists('compiledMeta'))source('readNewData.R')
out<-compiledMeta[,c('mm','rDate','vl','cd4','time')]
out$time<-as.numeric(out$time)
dat$hasIc50<-!is.na(dat$ic50)|!is.na(dat$beta)
noBulkN<-withAs(dat=dat[!dat$bulk&!dat$qvoa&dat$hasIc50,],tapply(dat$id,paste(dat$pat,dat$time),length))
bulkN<-withAs(dat=dat[dat$bulk&!dat$qvoa&dat$hasIc50,],tapply(dat$id,paste(dat$pat,dat$time),length))
voaN<-withAs(dat=dat[dat$qvoa&dat$hasIc50,],tapply(dat$id,paste(dat$pat,dat$time),length))

out$singleIsolate<-noBulkN[paste(out$mm,out$time)]
out$bulkIsolate<-bulkN[paste(out$mm,out$time)]
out$voaIsolate<-voaN[paste(out$mm,out$time)]
if(sum(out$singleIsolate,na.rm=TRUE)!=sum(noBulkN))stop('Problem assigning single')
if(sum(out$bulkIsolate,na.rm=TRUE)!=sum(bulkN))stop('Problem assigning bulk')
if(sum(out$voaIsolate,na.rm=TRUE)!=sum(voaN))stop('Problem assigning VOA')


seqNames<-unlist(lapply(list.files('trees/10trees','\\.txt$',full.name=TRUE),function(ii)ape::read.tree(ii)$tip.label))
seqInfo<-do.call(rbind,lapply(strsplit(seqNames,'\\.'),'[',1:4))

plasSgsN<-withAs(seqInfo=seqInfo[seqInfo[,3]=='SGS'&seqInfo[,2]=='PLAS',],tapply(seqInfo[,1],paste(seqInfo[,1],as.numeric(seqInfo[,4])),length))
pbmcSgsN<-withAs(seqInfo=seqInfo[seqInfo[,3]=='SGS'&seqInfo[,2]=='PBMC',],tapply(seqInfo[,1],paste(seqInfo[,1],as.numeric(seqInfo[,4])),length))
isoSeqN<-withAs(seqInfo=seqInfo[seqInfo[,3]=='ISO',],tapply(seqInfo[,1],paste(seqInfo[,1],as.numeric(seqInfo[,4])),length))


missingWeau<-names(plasSgsN)[!names(plasSgsN) %in% paste(out$mm,out$time)&grepl('WEAU',names(plasSgsN))]
tmp<-out[1,]
tmp[,]<-NA
tmp<-tmp[rep(1,length(missingWeau)),]
tmp$mm<-"WEAU"
tmp$time<-as.numeric(sub('WEAU ','',missingWeau))
tmp$rDate<-weauSymptomDate+tmp$time
out<-rbind(out,tmp)

out$plasmaSga<-plasSgsN[paste(out$mm,out$time)]
out$pbmcSga<-pbmcSgsN[paste(out$mm,out$time)]
if(sum(out$plasmaSga,na.rm=TRUE)!=sum(plasSgsN))stop('Problem assigning plasma SGS')
if(sum(out$pbmcSga,na.rm=TRUE)!=sum(pbmcSgsN))stop('Problem assigning PBMC SGS')

out$isolateSequence<-isoSeqN[paste(out$mm,out$time)]
if(sum(out$isolateSequence,na.rm=TRUE)!=sum(isoSeqN))stop('Problem assigning isolate seqs')

if(any((!is.na(out$pbmcSga)|!is.na(out$voaIsolate))&(!is.na(out$plasmaSga)|!is.na(out$singleIsolate)|!is.na(out$bulkIsolate))))stop('Plasma and pbmc mixed')
out$sampleType<-ifelse(!is.na(out$pbmcSga)|!is.na(out$voaIsolate),'PBMC',ifelse((!is.na(out$plasmaSga)|!is.na(out$singleIsolate)|!is.na(out$bulkIsolate)),'Plasma',NA))

out<-out[apply(!is.na(out[,c('vl','cd4','singleIsolate','bulkIsolate','voaIsolate','plasmaSga','pbmcSga')]),1,any),]

#Add iso sequence column

write.csv(out,'out/compiledMeta.csv',row.names=FALSE)


seqs<-dnar::read.fa('out/mmLongitudinalSequences.fa')
#seqs<-read.fa('out/mmLongitudinalSequences.fa.gz')
aligns<-dnar::readFaDir('trees/10trees/aln/alns/nomask/')
if(any(!aligns$name %in% seqs$name))stop('Name mismatch')
rownames(seqs)<-seqs$name
idents<-sapply(1:nrow(aligns),function(ii)grepl(dnar::degap(aligns$seq[ii]),seqs[aligns$name[ii],'seq'],perl=TRUE))
if(any(!idents))stop('Sequence mismatch')
idConvert<-read.csv('out/newNames.csv',stringsAsFactors=FALSE)
rownames(idConvert)<-idConvert$newName
if(any(!aligns$name %in% idConvert$newName))stop('Conversion mismatch')
idents<-sapply(1:nrow(aligns),function(ii)grepl(dnar::degap(aligns$seq[ii]),idConvert[aligns$name[ii],'seq'],perl=TRUE))
if(any(!idents))stop('Sequence mismatch')

#dists<-levenR::leven(aligns$seq,seqs[aligns$name,'seq'],oneToOne=TRUE,substring2=TRUE)

isMin<-as.numeric(seqInfo[,4])==ave(as.numeric(seqInfo[,4]),seqInfo[,1],FUN=min)
firsts<-seqNames[isMin]
table(seqInfo[isMin,1])
firstSeqs<-seqs[seqs$name %in% firsts,]
firstSeqs<-firstSeqs[order(firstSeqs$name),]
firstAligns<-aligns[aligns$name %in% firsts,]
firstAligns<-firstAligns[order(firstAligns$name),]
if(nrow(firstSeqs)!=length(firsts)|nrow(firstAligns)!=length(firsts))stop('Missing sequences')
write.fa(firstSeqs$name,firstSeqs$seq,'out/firstSeqs.fa')
write.fa(firstAligns$name,firstAligns$seq,'out/firstAligns.fa')

isMinSeq<-withAs(mms=sapply(strsplit(seqs$name,'\\.'),'[[',1),seqDay=as.numeric(sapply(strsplit(seqs$name,'\\.'),'[[',4)),seqDay==ave(seqDay,mms,FUN=min))
firstAllSeqs<-seqs[isMinSeq,]
table(sapply(strsplit(firstAllSeqs$name,'\\.'),'[',1))
write.fa(firstAllSeqs$name,firstAllSeqs$seq,'out/firstAllSeqs.fa')



xx<-read.fa('data/Jesse_MM14.fasta')

xx<-read.csv('data/List of bulks for name changing.csv',header=FALSE,stringsAsFactors=FALSE)
colnames(xx)<-'old'
xx$newName<-sapply(xx$old,function(xxx){
  out<-idConvert[idConvert$isolate==xxx&!is.na(idConvert$isolate),'newName']
  if(length(out)==0)return(NA)
  if(length(out)>1)stop('Multiple isolate match')
  return(out)
})
xx$seqName<-sapply(xx$old,function(xxx){
  out<-idConvert[idConvert$isolate==xxx&!is.na(idConvert$isolate),'name']
  if(length(out)==0)return(NA)
  if(length(out)>1)stop('Multiple isolate match')
  return(out)
})
write.csv(xx,'out/listOfBulksForNameChanging.csv')
isolatesMissingSeqs<-dat[(!is.na(dat$ic50)|!is.na(dat$beta))&!dat$id %in% idConvert$isolate,c('id','Patient.Original.ID','ic50','beta','Sequence..Full.Length.or.3.half.')]
write.csv(isolatesMissingSeqs,'out/isolatesMissingSeqs.csv')


if(any(!idConvert[!is.na(idConvert$isolate),'isolate'] %in% dat$id))stop('Unidentified isolate sequence')
zz<-dat[dat$id %in% idConvert[!is.na(idConvert$isolate),'isolate'],c('ic50','beta','id')]
zz$newName<-sapply(zz$id,function(xx)idConvert[idConvert$isolate==xx&!is.na(idConvert$isolate),'newName'])
zz$oldName<-sapply(zz$id,function(xx)idConvert[idConvert$isolate==xx&!is.na(idConvert$isolate),'name'])
write.csv(zz[apply(is.na(zz[,c('ic50','beta')]),1,sum)==2,],'out/isolatesWithSequenceButNoIC50.csv',row.names=FALSE)

align<-read.fa('out/mmLongitudinalSequencesAlign/AlignSeq.nt.fasta.gz')
#alignIso$firstPos<-regexpr('[ACGT]',alignIso$seq)
#alignIso$lastPos<-regexpr('[ACGT][^ACGT]*$',alignIso$seq)
align$firstPos<-gap2NoGap(align[1,'seq'],sapply(align$seq,noGap2Gap,50))
align$lastPos<-gap2NoGap(align[1,'seq'],mapply(noGap2Gap,align$seq,nchar(degap(align$seq))-50))
align$oldName<-sapply(align$name,function(xx)ifelse(xx %in% idConvert$newName,idConvert[idConvert$newName==xx,'isolate'],NA))
align$seqName<-sapply(align$name,function(xx)ifelse(xx %in% idConvert$newName,idConvert[idConvert$newName==xx,'name'],NA))
align$type<-ifelse(align$firstPos<3000,
                 ifelse(align$lastPos>8500,'Full',"5'"),
                 ifelse(align$lastPos>9000,"3'","Env")
               )
alignIso<-align[align$name %in% c(align$name[1],idConvert[!is.na(idConvert$isolate),'newName']),]
write.csv(alignIso[,c('oldName','name','seqName','type')],'out/isolateSequencingType.csv',row.names=FALSE)

png('out/isolateAlign.png',width=5000,height=3000,res=300);withAs(xx=alignIso[order(alignIso$firstPos,alignIso$lastPos),],dnaplotr::plotDNA(xx$seq,refSeq=xx$seq[1],xlab='HXB2 coordinates',groups=xx$type));dev.off()
png('test.png',width=5000,height=3000,res=300);dnaplotr::plotDNA(alignIso$seq[c(1,which(alignIso$lastPos>10846&alignIso$lastPos<10870))],refSeq=alignIso$seq[1]);dev.off()

#https://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html for RT coords
pol<-read.fa('out/mmLongitudinalSequencesAlign/Pol.nt.fasta.gz')
rtStart<-regexpr('PISPIETVPV',dna2aa(degap(pol[1,'seq'])))
rtEnd<-regexpr('KEPIVGAETF',dna2aa(degap(pol[1,'seq'])))+9
substring(dna2aa(degap(pol[1,'seq'])),rtStart,rtEnd)
pol$rt<-substring(pol$seq,noGap2Gap(pol[1,'seq'],(rtStart-1)*3+1),noGap2Gap(pol[1,'seq'],rtEnd*3))
pol$rtProt<-dna2aa(degap(pol$rt))
pol[grep('WEAU',pol$name),'rtProt']

dat$nIfna2Ic50<-apply(!is.na(dat[,ifna2_ic50]),1,sum)
dat$nIfnbIc50<-apply(!is.na(dat[,ifnb_ic50]),1,sum)
rownames(dat)<-dat$id
s2<-read.csv('data/Table S2.nearfinal.csv',skip=2,stringsAsFactors=FALSE,check.names=FALSE)
if(any(!s2[,'Isolate ID'] %in% dat$id & !is.na(s2[,'Isolate ID']) &s2[,'Isolate ID']!=''&!grepl('IMC',s2[,'Isolate ID'])))stop('Problem finding S2 virus')
s2$nIfna2Ic50<-dat[s2[,'Isolate ID'],'nIfna2Ic50']
s2$nIfnbIc50<-dat[s2[,'Isolate ID'],'nIfnbIc50']
seqTypes<-read.csv('out/isolateSequencingType.csv',stringsAsFactors=FALSE)
seqTypes<-seqTypes[seqTypes$name!='B.FR.83.HXB2_LAI_IIIB_BRU_K03455',]
#if(any(!s2[,'Isolate ID'] %in% seqTypes$oldName))stop('Missing seq')
rownames(seqTypes)<-seqTypes$oldName
s2$seqType<-seqTypes[s2[,'Isolate ID'],'type']
convert<-c('MM15.2._P2A4'='MM15.02.2A4.bulk','MM15.02.1C2'='MM15.02.1C2.bulk','MM15.3_P2B1'='MM15.03.2B1.bulk','MM15.14.2C4'='MM15.14.2C4.bulk')
trop<-read.csv('out/mmTropism.csv',stringsAsFactors=FALSE,row.names=1)
trop[trop$id %in% names(convert),'id']<-convert[trop[trop$id %in% names(convert),'id']]
trop[trop$id %in% dat$ID.for.Publications&!trop$id %in% s2[,'Isolate ID'],'id']<-sapply(trop[trop$id %in% dat$ID.for.Publications&!trop$id %in% s2[,'Isolate ID'],'id'],function(xx)dat[dat$ID.for.Publications==xx,'id'])
rownames(trop)<-trop$id
trop<-trop[trop$tropism!='Problem',]
s2$tropism<-trop[s2[,'Isolate ID'],'tropism']
write.csv(s2,'out/updatedS2.csv')


imcFiles<-c(
  'out/fred_20200224_ic50.csv',
  'out/imc_20191105_ic50.csv',
  'out/imc_20191115_ic50.csv',
  'out/zapRedo_20190712.csv'
)
imcIc50s<-lapply(imcFiles,read.csv,stringsAsFactors=FALSE,row.names=1)
imcIc50s<-lapply(imcIc50s,function(xx)xx[grep('MM33',rownames(xx)),])
desiredCols<-c('ic50.IFNa2','ic50.IFNb','repCap.IFNa2','repCap.IFNb','vres.IFNa2','vres.IFNb')
standardized<-do.call(rbind,mapply(function(xx,yy){
  if(yy=='out/imc_20191105_ic50.csv'){
    colnames(xx)<-sprintf('%s.IFNb',colnames(xx))
    xx<-xx[grepl('50$',rownames(xx)),]
  }
  if(yy=='out/zapRedo_20190712.csv')colnames(xx)[colnames(xx)=='IFNa2_repCap.1']<-'IFNb_repCap'
  colnames(xx)<-sub('Vres','vres',sub('IC50','ic50',sub('(IFN[^_.]+)[_.](.*)','\\2.\\1',colnames(xx))))
  xx<-xx[!grepl('MM33.17.2A4|MM33.1.13C1',rownames(xx)),]
  xx$file<-yy
  xx$virus<-sub('^.*MM33[ ._](TF|13|17).*$','MM33.\\1',rownames(xx))
  xx[,desiredCols[!desiredCols %in% colnames(xx)]]<-NA
  xx[,c(desiredCols,'virus','file')]
},imcIc50s,imcFiles,SIMPLIFY=FALSE))
out<-do.call(rbind,by(standardized,standardized$virus,function(xx){
  out<-xx[1,]
  out[,desiredCols]<-NA
  out[,desiredCols[!grepl('ic50',desiredCols)]]<-apply(xx[,desiredCols[!grepl('ic50',desiredCols)]],2,mean,na.rm=TRUE)
  out[,desiredCols[grepl('ic50',desiredCols)]]<-apply(xx[,desiredCols[grepl('ic50',desiredCols)]],2,function(xx)exp(mean(log(xx),na.rm=TRUE)))
  out[,desiredCols]
}))
standardized<-standardized[order(standardized$virus,standardized$file),]
write.csv(standardized,'out/mm33_imc_allIc50.csv')
write.csv(out,'out/mm33_imc_meanIc50.csv')



