##########################################################################
##########  convert fixedStep average.wig file to variableStep  ##########
########## average.wig from averaging input subtracted wig files #########
##########################################################################
#make sure that the fixedStep file has a 1bp step (look at the header)

args <- commandArgs(TRUE)

wig<-read.delim(file=args[1],header=F)
chrI<-as.data.frame(wig[2:15072425,])
chrII<-as.data.frame(wig[15072426:30351771,])
chrIII<-as.data.frame(wig[30351772:44135472,])
chrIV<-as.data.frame(wig[44135473:61629266,])
chrMtDNA<-as.data.frame(wig[61629267:61643061,])
chrV<-as.data.frame(wig[61643062:82567211,])
chrX<-as.data.frame(wig[82567212:100286078,])
rm(wig)
chrX$base<-c(0:17718866)
colnames(chrX)<-c('signal','base')
chrXvariable<-chrX[which(!chrX$signal == " "),]
chrV$base<-c(0:20924149)
colnames(chrV)<-c('signal','base')
chrVvariable<-chrV[which(!chrV$signal == " "),]
chrIV$base<-c(0:17493793)
colnames(chrIV)<-c('signal','base')
chrIVvariable<-chrIV[which(!chrIV$signal == " "),]
chrIVvariable<-chrIVvariable[,c(2,1)]
chrVvariable<-chrVvariable[,c(2,1)]
chrXvariable<-chrXvariable[,c(2,1)]
chrIII$base<-c(0:13783700)
colnames(chrIII)<-c('signal','base')
chrIIIvariable<-chrIII[which(!chrIII$signal == " "),]
chrIIIvariable<-chrIIIvariable[,c(2,1)]
chrII$base<-c(0:15279345)
colnames(chrII)<-c('signal','base')
chrIIvariable<-chrII[which(!chrII$signal == " "),]
chrIIvariable<-chrIIvariable[,c(2,1)]
chrI$base<-c(0:15072423)
colnames(chrI)<-c('signal','base')
chrIvariable<-chrI[which(!chrI$signal == " "),]
chrIvariable<-chrIvariable[,c(2,1)]
chrMtDNA$base<-c(0:13794)
colnames(chrMtDNA)<-c('signal','base')
chrMtDNAvariable<-chrMtDNA[which(!chrMtDNA$signal == " "),]
chrMtDNAvariable<-chrMtDNAvariable[,c(2,1)]
write.table(chrIvariable,file=args[2],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(chrIIvariable,file=args[3],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(chrIIIvariable,file=args[4],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(chrIVvariable,file=args[5],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(chrVvariable,file=args[6],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(chrMtDNAvariable,file=args[7],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(chrXvariable,file=args[8],col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)

