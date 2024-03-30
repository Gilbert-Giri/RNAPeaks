GetGene<-function(geneID,species,TxID,gtf){
  if (is.null(geneID) & is.null(TxID)){
    stop("Enter Gene or Transcript ID!")
  }

  if (is.null(gtf) & species=="Human"){
    #AnnotationHub::query(ah, pattern = c("Homo Sapiens","gtf"))
    ah<-AnnotationHub::AnnotationHub()
    gtf<-data.frame(ah[["AH110867"]])
  } else if (is.null(gtf) & species=="Mouse"){
    #query(ah, pattern = c("Mus musculus", "gtf"))[[1]]
    ah<-AnnotationHub::AnnotationHub()
    gtf<-data.frame(ah[["AH47076"]])
  }

  if (!(geneID %in% gtf$gene_id | geneID %in% gtf$gene_name | (!is.na(TxID) & (TxID %in% gtf$transcript_id)))){
    stop("Gene cannot be found!")
  }

  if (!is.na(TxID)){
    gtf<-gtf[which(gtf$transcript_id==TxID),]
  } else if (is.na(TxID)){
    if (grepl("ENSG",geneID)){
      gtf<-gtf[which(gtf$gene_id==geneID),]
    } else {
      gtf<-gtf[which(gtf$gene_name==geneID),]
    }
    transcripts_all<-unique(gtf[which(gtf$type=="transcript"),c("transcript_id","width")])
    TxID<-transcripts_all$transcript_id[which.max(transcripts_all$width)]
    gtf<-gtf[which(gtf$transcript_id==TxID),]
  }
  return(gtf)
}
