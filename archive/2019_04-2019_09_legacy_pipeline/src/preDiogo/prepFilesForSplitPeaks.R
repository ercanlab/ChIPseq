########################################################################################
##this script is used to prepare files to be input into getFinalSplitPeaksAndSummits.R##
########################################################################################


args <- commandArgs(TRUE)

###get subpeaks.overlap.bed file###
final.peaks<-read.delim(file=args[1],header=F)
split<-read.delim(file='temp.sorted.bed',header=F)
nm<-c("V2","V3")
split[nm] <- lapply(nm, function(x) final.peaks[[x]][match(split$V4, final.peaks$V4)])
write.table(split,file='temp.sorted_1.bed',row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

###get summits that match peaks_final.bed###
summits<-read.delim(file=args[2],header=F)
summits.final<-summits[which(summits[,4] %in% final.peaks[,4]),] #get summits of original peaks
write.table(summits.final,file=args[3],row.names=F,col.names=F,quote=F,sep='\t')

