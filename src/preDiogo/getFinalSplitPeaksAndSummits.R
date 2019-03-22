####################################################################
##this script gets the final split peak and split summit bed files##
####################################################################


args <- commandArgs(TRUE)

# read original peaks: final_peaks.bed file (merged bam run through macs filtered with e-5 peak files from individual reps)
fr <- read.delim(file=args[1], header=F)
all.peaks <- as.matrix(fr)

fr <- read.delim(file=args[2], header=F)
all.summits <- as.matrix(fr)

# read peaks that were splitted: peaks_final_splitted.bed
fr <- read.delim(file=args[3], header=F)
splitted.peaks <- as.matrix(fr)

# read splitted peak information: subpeaks.overlap.bed file
fr <- read.delim(file=args[4], header=F)
splits <- as.matrix(fr)

###
#get split_peaks_final.bed
final.peaks <- all.peaks[which(!all.peaks[,4] %in% splitted.peaks[,4]),] #get all original peaks that were not split
res <- splits[which(splits[,9] %in% splitted.peaks[,4]),] #get the splits for the peaks that were split
peaks.to.split <- unique(res[,9]) #get the MACS_peak_#s of the original peaks that were split

new.peaks <- NULL

for (j in 1:length(peaks.to.split)) { #for each original peak that was split...
  cur <- res[which(res[,9]==peaks.to.split[j]),] #get the split peaks
  nlen <- nrow(cur) #get the number of new peaks under the original peak
  temp <- cur
  temp[1,2] <- cur[1,7] #these next four lines split up the original peak based on the locations of the split peaks (the new peaks will cover the same total coordinates as the original peak)
  temp[1,3] <- round(as.numeric(cur[1,3]) + ((as.numeric(cur[2,2]) - as.numeric(cur[1,3])) / 2))
  temp[nlen,2] <- round(as.numeric(cur[nlen,2]) - ((as.numeric(cur[nlen,2]) - as.numeric(cur[(nlen-1),3])) / 2))+1
  temp[nlen,3] <- cur[nlen,8]
  if (nlen > 2) { #if there are more than 2 split peaks in the orignial peak, these lines will do the same splitting as above but for a larger number of split peaks
    for (i in 2:(nlen-1)) {
      temp[i,2] <- round(as.numeric(cur[i,2]) - ((as.numeric(cur[i,2]) - as.numeric(cur[(i-1),3])) / 2))+1
      temp[i,3] <- round(as.numeric(cur[i,3]) + ((as.numeric(cur[(i+1),2]) - as.numeric(cur[i,3])) / 2))
    }
  }
  new.peaks <- rbind(new.peaks, temp) #add the new split peak coordinates to the object new.peaks
}

final.peaks <- rbind(final.peaks[,1:4], new.peaks[,c(1:3,9)]) #get the chr, start, stop, and MACS_peak_# for old peaks that were not split and new split peaks
final.peaks[,2:3] <- as.numeric(final.peaks[,2:3])
dim(final.peaks)
sorted <- order(as.numeric(sub("MACS_peak_","",final.peaks[,4])))
peaks <- unique(final.peaks[sorted,4])
for (peak in peaks) {
  ind <- which(final.peaks[sorted,4] %in% peak)
  if (length(ind) > 1) {
    for (i in 1:length(ind)) {
      final.peaks[sorted[ind[i]],4] <- paste(final.peaks[sorted[ind[i]],4], i, sep=".") #rename the split peaks MACS_peak_#.1, MACS_peak_#.2, etc.
    }
  }
}

write.table(final.peaks[sorted,], file=args[5], col.names=F, row.names=F, quote=F, sep="\t") #write the new bed file with the final peaks after splitting

###
#get split_summits_final.bed
final.summits <- all.summits[which(!all.summits[,4] %in% splitted.peaks[,4]),] #get the summits for all original peaks that were not split
res <- splits[which(splits[,9] %in% splitted.peaks[,4]),c(1,5,9)] #get the summits of the split peaks
res <- cbind(res[,1:2], as.numeric(res[,2])+1, res[,3]) #make a summit stop coordinate 1bp after the start coordinate
final.summits <- rbind(final.summits[,1:4], res)
final.summits[,2:3] <- as.numeric(final.summits[,2:3])
dim(final.summits)
sorted <- order(as.numeric(sub("MACS_peak_","",final.summits[,4])))
summits <- unique(final.summits[sorted,4])
for (summit in summits) {
  ind <- which(final.summits[sorted,4] %in% summit)
  if (length(ind) > 1) {
    for (i in 1:length(ind)) {
      final.summits[sorted[ind[i]],4] <- paste(final.summits[sorted[ind[i]],4], i, sep=".") #rename the split peaks MACS_peak_#.1, MACS_peak_#.2, etc.
    }
  }
}

write.table(final.summits[sorted,], file=args[6], col.names=F, row.names=F, quote=F, sep="\t") #write the new bed file with the final summits after splitting

