rebQ<-read.csv('../data/Rebounds and QVOAs .csv',stringsAsFactors=FALSE)
rebQ$Study<-dnar::fillDown(rebQ$Study)

source('readData.R')
comboRQ<-combo[combo$class %in% c('QVOA','Rebound'),]
comboRQ$simple<-trimws(sub(' pre-ATI','',sub(' - ','-',sub('\\.VOA','',sub(' \\([^)]+\\)','',sub('  +',' ',sub('Rebound ','',sub(' Alpha.*','',comboRQ$virus))))))))
rebQ$simple<-trimws(sub(' - ','-',sub(' \\([^)]+\\)','',sub('  +',' ',sub('Rebound |QVOA | pre-ATI|Rb ','',rebQ$REBOUNDS.ID)))))
if(any(!comboRQ$simple %in% rebQ$simple)||any(!rebQ$simple %in% comboRQ$simple))stop('Mismatch in viruses')
rownames(comboRQ)<-comboRQ$simple
rebQ$ic50_IFNa2<-comboRQ[rebQ$simple,'ic50']
rebQ$repCap<-comboRQ[rebQ$simple,'repCap']
rebQ$vres_IFNa2<-comboRQ[rebQ$simple,'vres']


source('readBeta.R')
comboRQ<-combo[combo$class %in% c('QVOA','Rebound'),]
comboRQ$simple<-trimws(sub('A09 ','',sub(' [12]0$','',sub(' - ','-',sub('\\.VOA','',sub(' \\([^)]+\\)','',sub('  +',' ',sub('Rebound | Pre ATI| pre-ATI','',sub(' Beta.*','',comboRQ$virus)))))))))
rebQ$simple<-trimws(sub(' Pre ATI','',sub(' - ','-',sub(' \\([^)]+\\)','',sub('  +',' ',sub('Rebound |QVOA | pre-ATI|Rb ','',rebQ$REBOUNDS.ID))))))
if(any(!comboRQ$simple %in% rebQ$simple)||any(!rebQ$simple %in% comboRQ$simple))stop('Mismatch in viruses')
rownames(comboRQ)<-comboRQ$simple
rebQ$ic50_IFNb<-comboRQ[rebQ$simple,'beta']
rebQ$vres_IFNb<-comboRQ[rebQ$simple,'vres']

dat<-read.csv('../out/allLongitudinal.csv',stringsAsFactors=FALSE)
mmVoa<-dat[dat$qvoa,]
rownames(mmVoa)<-sub('\\.VOA','',mmVoa$id)
if(any(!rownames(mmVoa) %in% rebQ$simple))stop('Extra MM')
if(any(!rebQ[grep('^MM',rebQ$simple),'simple'] %in% rownames(mmVoa)))stop('Missing MM')
rebQ[grep('^MM',rebQ$simple),'vres_IFNa2']<-mmVoa[rebQ[grep('^MM',rebQ$simple),'simple'],'vres']
rebQ[grep('^MM',rebQ$simple),'vres_IFNb']<-mmVoa[rebQ[grep('^MM',rebQ$simple),'simple'],'betaVres']


colnames(rebQ)[colnames(rebQ)=='REBOUNDS.ID']<-'ID'
write.csv(rebQ[,!colnames(rebQ) %in% c('Has.IC50','X0.5ng','X500ul')],'out/rebQvoaMaster.csv',row.names=FALSE)



