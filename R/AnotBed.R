AnotBed<-function(peak=NULL,enDb=NULL,species="Human"){

  if (is.null(peak)){stop("Provide peak bed file.")}
  if(is.null(enDb)){
    ah<-AnnotationHub::AnnotationHub()
    if (species=="Human"){
      #query(ah, pattern = c("Homo sapiens", "EnsDb"))
      enDb<-ah[["AH109606"]]
    } else if (species=="Mouse"){
      #query(ah, pattern = c("Mus musculus", "EnsDb"))
      enDb<-ah[["AH89211"]]
    } else (stop ("The species is not supported."))
  }

  if (is.data.frame(peak)){
    peak<-makeGRangesFromDataFrame(peak,
                                   seqnames.field = colnames(peak)[1],
                                   start.field = colnames(peak)[2],
                                   end.field = colnames(peak)[3],
                                   strand.field=colnames(peak)[6],
                                   starts.in.df.are.0based=T)
  }
  seqlevelsStyle(peak)<-"NCBI"

  #Getting regions
  miRNA<-unlist(exonsBy(enDb,by="gene")[genes(enDb)$gene_id[which(genes(enDb)$gene_biotype=="miRNA")]])
  Coding_Exons<-unlist(cdsBy(enDb,by="gene"))
  Introns<-unlist(intronsByTranscript(enDb))
  Utr5<-unlist(fiveUTRsByTranscript(enDb))
  Utr3<-unlist(threeUTRsByTranscript(enDb))
  NCexons<-unlist(exonsBy(enDb,by="gene")[genes(enDb)$gene_id[which(genes(enDb)$gene_biotype=="lncRNA")]])

  #Getting Splice Site annotation
  all_cds<-cdsBy(enDb,by="gene")
  all_ncexons<-exonsBy(enDb,by="gene")[genes(enDb)$gene_id[which(genes(enDb)$gene_biotype=="lncRNA")]]

  all_cds<-all_cds[-which(lapply(all_cds, length)==1)]
  all_ncexons<-all_ncexons[-which(lapply(all_ncexons, length)==1)]

  #Getting 5' Splice site
  all_cds_fivess<-unlist(lapply(all_cds,Get_Junctions_five))
  all_ncexons_fivess<-unlist(lapply(all_ncexons,Get_Junctions_five))

  #Getting 3' Splice site
  all_cds_threes<-lapply(all_cds,Get_Junctions_three)
  all_ncexons_threes<-lapply(all_ncexons,Get_Junctions_three)

  #converting them to GRanges
  all_cds_fivess<-unlist(GRangesList(all_cds_fivess))
  all_ncexons_fivess<-unlist(GRangesList(all_ncexons_fivess))
  all_cds_threes<-unlist(GRangesList(all_cds_threes))
  all_ncexons_threes<-unlist(GRangesList(all_ncexons_threes))

  splicesite5<-c(all_cds_fivess,all_ncexons_fivess)
  splicesite3<-c(all_cds_threes,all_ncexons_threes)


  #Getting the peaks

  #miRNA
  #length(peak)
  miRNA_peaks<-peak[Get_overlaps(peak = peak,anot = miRNA)]
  miRNA_peaks$Region<-"miRNA"
  peak<-peak[-Get_overlaps(peak = peak,anot = miRNA)]
  #length(peak)
  Coding_Exons_peaks<-peak[Get_overlaps(peak = peak,anot = Coding_Exons)]
  Coding_Exons_peaks$Region<-"CDS"
  peak<-peak[-Get_overlaps(peak = peak,anot = Coding_Exons)]
  #length(peak)
  Utr5_peaks<-peak[Get_overlaps(peak = peak,anot = Utr5)]
  Utr5_peaks$Region<-"UTR5"
  peak<-peak[-Get_overlaps(peak = peak,anot = Utr5_peaks)]
  #length(peak)
  Utr3_peaks<-peak[Get_overlaps(peak = peak,anot = Utr3)]
  Utr3_peaks$Region<-"UTR3"
  peak<-peak[-Get_overlaps(peak = peak,anot = Utr3)]
  #length(peak)
  NCExons_peaks<-peak[Get_overlaps(peak = peak,anot = NCexons)]
  NCExons_peaks$Region<-"NCExon"
  peak<-peak[-Get_overlaps(peak = peak,anot = NCexons)]
  #length(peak)
  ss3_peaks<-peak[Get_overlaps(peak = peak,anot = splicesite3)]
  ss3_peaks$Region<-"SpliceSite3"
  peak<-peak[-Get_overlaps(peak = peak,anot = splicesite3)]
  #length(peak)
  ss5_peaks<-peak[Get_overlaps(peak = peak,anot = splicesite5)]
  ss5_peaks$Region<-"SpliceSite5"
  peak<-peak[-Get_overlaps(peak = peak,anot = splicesite5)]
  #length(peak)
  Introns_peaks<-peak[Get_overlaps(peak = peak,anot = Introns)]
  Introns_peaks$Region<-"Intron"
  peak<-peak[-Get_overlaps(peak = peak,anot = Introns)]
  #length(peak)
  peak$Region<-"Others"

  comb<-c(miRNA_peaks,Coding_Exons_peaks,Utr5_peaks,Utr3_peaks,NCExons_peaks,ss3_peaks,ss5_peaks,Introns_peaks,peak)
  return(data.frame(comb))

}

Get_overlaps<-function(peak,anot){
  ov <- findOverlaps(peak, anot)
  overlaps <- pintersect(peak[queryHits(ov)], anot[subjectHits(ov)])
  percentOverlap <- width(overlaps) / width(peak[queryHits(ov)])
  m<-queryHits(ov)[which(percentOverlap>0.5)]
  return(unique(m))
}


Get_Junctions_five<-function(x){
  if (length(x)==1){return(NULL)}
  strand<-strand(x)@values[1]
  if (strand=="+"){
    x_resize<-GenomicRanges::flank(x[1:(length(x)-1)],width=50,start=F,both = T)
  } else if (strand=="-") {
    x_resize<-GenomicRanges::flank(x[-1],width=50,start=F,both = T)
  }
  return(x_resize)
}


Get_Junctions_three<-function(x){
  if (length(x)==1){return(NULL)}
  strand<-strand(x)@values[1]
  if (strand=="+"){
    x_resize<-GenomicRanges::flank(x[-1],width=50,start=T,both = T)
  } else if (strand=="-") {
    x_resize<-GenomicRanges::flank(x[1:(length(x)-1)],width=50,start=T,both = T)
  }
  return(x_resize)
}


