
PlotRegion<-function(Chr=NULL,
                     Start=NULL,
                     End=NULL,
                     Strand=NULL,
                     geneID=NA,
                     gtf=NULL,
                     bed=NULL,
                     species="Human",
                     TxID=NA,
                     Target_col=NULL,
                     omit=c(),
                     order_by="Count",
                     order_in=NULL,
                     merge=0,
                     peaks_width=0.3,
                     utr_col="dark gray",
                     peak_col="Blue",
                     exon_width=0.5,
                     utr_width=0.3,
                     exon_col="black"){

  if (is.null(Chr) | is.null(Start) | is.null(End) | is.null(Strand)){
    stop("Need to provide Chr, Start, Stop and Strand parameters to visualize region.")
  }

  Range<-GenomicRanges::GRanges(seqnames = Chr,strand = Strand,ranges = IRanges(start = Start,end=End))
  print(Range)
  GenomeInfoDb::seqlevelsStyle(Range)<-"NCBI"

  #Loading a genome each time takes long so giving an option to load it once and input it
  Region<-GetGene(geneID=geneID,species=species,TxID=TxID,gtf=gtf)

  Region_gr<-GenomicRanges::makeGRangesFromDataFrame(Region,keep.extra.columns = T)
  Region_gr<-data.frame(Region_gr[queryHits(IRanges::findOverlaps(Region_gr,Range))])
  print(Region_gr)

  #check the bed file
  colnames(bed)[which(colnames(bed)==Target_col)]<-"target"
  bed<-checkBed(bed)
  bed<-FilterBed(bed=bed,chr=as.character(seqnames(Range)),start=Start,end=End,strand=Strand,omit=omit,collapse=merge)
  print(bed)

  #Option to give an order to rank proteins in
  if (!is.null(order_in)){
    Target_rank<-order_in
  } else {
    Target_rank<-OrderPeak(bed=bed,order_by=order_by)
  }

  plot<-Get_Plot_by_Region(Chr=Chr,
                           Start=Start,
                           End=End,
                           Strand=Strand,
                           bed=bed,
                           gene=Region_gr,
                           rank_=Target_rank,
                           utr_col=utr_col,
                           peak_col=peak_col,
                           peaks_width=peaks_width,
                           exon_width=exon_width,
                           utr_width=utr_width,
                           exon_col=exon_col)

  return(plot)
}

