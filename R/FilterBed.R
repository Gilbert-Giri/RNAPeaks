FilterBed<-function(bed,chr,start,end,strand,omit=c(),collapse){
  bed<-bed[which(as.character(bed$seqnames)==chr & bed$start>=start & bed$end<=end & bed$strand==strand),]

  #Removing Targets if given
  bed<-bed[which(!(bed$target %in% omit)),]

  #Converting into Granges list to merge
  bed_gr<-GenomicRanges::makeGRangesListFromDataFrame(bed,keep.extra.columns = T,split.field="target")
  sp<-list("Human"="Homo_sapiens","Mouse"="Mus_musculus")
  bed_gr<-GenomicRanges::reduce(bed_gr,min.gapwidth=collapse)

  #converting back to df
  red_bed<-data.frame(bed_gr)
  if (nrow(red_bed)==0){
    stop("No Peaks Found!")
  }
  return(red_bed)
}
