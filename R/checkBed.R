checkBed<-function(df){
  colnames(df)<-tolower(colnames(df))

  if (!(all(c("start","end","strand","chr") %in% colnames(df)))) {
    stop("Must Contain Chr,Start,End and Strand columns!")
  }

  if(is.character(df$chr) &
             all(as.numeric(df$end)>=as.numeric(df$start)) &
             is.character(df$strand) &
             all(df$strand %in% c("+","-"))){

    #converting the version of the seqnames
    df_gr<-GenomicRanges::makeGRangesFromDataFrame(df,keep.extra.columns = T)
    GenomeInfoDb::seqlevelsStyle(df_gr)<-"NCBI"
    df<-data.frame(df_gr)
    return(df)
  } else {
    stop("Check Bed File!")
  }
}
