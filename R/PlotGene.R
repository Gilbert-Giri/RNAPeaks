PlotGene<-function(bed=NULL,
                   geneID=NA,
                   gtf=NULL,
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

  #check the bed file
  colnames(bed)[which(colnames(bed)==Target_col)]<-"target"
  bed<-checkBed(bed)

  #Loading a genome each time takes long so giving an option to load it once and input it
  Region<-GetGene(geneID=geneID,species=species,TxID=TxID,gtf=gtf)

  bed<-FilterBed(bed=bed,chr=Region$seqnames[1],start=min(Region$start),end=max(Region$end),strand=Region$strand[1],omit=omit,collapse=merge)

  #Option to give an order to rank proteins in
  if (!is.null(order_in)){
    Target_rank<-order_in
  } else {
    Target_rank<-OrderPeak(bed=bed,order_by=order_by)
  }

  #Plotting the peak plot
  Plot<-Get_Plot_by_Gene(bed=bed,gene =Region ,rank_=Target_rank,utr_col=utr_col,peak_col=peak_col,peaks_width=peaks_width,
                         exon_width=exon_width,utr_width=utr_width,exon_col=exon_col)
  return(Plot)
}
