SplicingMapSE<-function(splice_file=NULL,
                        peak_file=NULL,
                        psi=0.05,
                        pval=0.05,
                        ctrl_size=2,
                        FDR=0.05,
                        in_intron=300,
                        in_exon=50,
                        min_count=30,
                        ctrl_grp=1,
                        skipdf=NULL,
                        retdf=NULL,
                        ctrldf=NULL,
                        Process=T,
                        exon_width=0.3,
                        intron_col="gray",
                        SE_exon_col="orange",
                        Surround_exon_col="black",
                        Retained_col="red",
                        Skipped_col="blue",
                        Control_col="black",
                        Gap_offset=10){

  #checking parameters
  if (isTRUE(Process)){
    if (is.null(splice_file)){
      stop("Please provide bedtools path and rMATS SE Junction Counts file.")
    }
  }
  if (!isTRUE(Process)){
    if (is.null(skipdf) | is.null(retdf) | is.null(ctrldf)){
      stop("Please provide skipped, retained and control exons SE Junction Counts file in rMATS format.")
    }
  }

  if (is.null(peak_file)){
    stop("Missing peak bed file.")
  }

  #preparing splicing events
  if (isTRUE(Process)){
    splice_file_filtered<-Remove_Overlaps(df=splice_file,min_count=min_count)
    grouped_events<-Group_events(df=splice_file_filtered,pval=pval,FDR=FDR,psi=psi,ctrl_grp=ctrl_grp,ctrl_size=ctrl_size)
  } else {
    grouped_events<-list("Retained"=retdf,"Skipped"=skipdf,"Control"=ctrldf)
  }

  #converting peak file to Granges
  peak_file_gr<-makeGRangesFromDataFrame(peak_file,
                                      keep.extra.columns = T,
                                      seqnames.field = colnames(peak_file)[1],
                                      start.field = colnames(peak_file)[2],
                                      end.field = colnames(peak_file)[3],
                                      strand.field =colnames(peak_file)[6])

  #Creating ranges for exons
  grouped_events_Upstream<-lapply(grouped_events, Upstream_Ranges,int=in_intron,ex=in_exon)
  grouped_events_Downstream<-lapply(grouped_events, Downstream_Ranges,int=in_intron,ex=in_exon)
  grouped_events_SE5<-lapply(grouped_events, SE5_Ranges,int=in_intron,ex=in_exon)
  grouped_events_SE3<-lapply(grouped_events, SE53_Ranges,int=in_intron,ex=in_exon)

  #making a new list to apply functions to all of these
  comb<-list("grouped_events_Upstream"=grouped_events_Upstream,
             "grouped_events_Downstream"=grouped_events_Downstream,
             "grouped_events_SE5"=grouped_events_SE5,
             "grouped_events_SE3"=grouped_events_SE3)

  comb<-lapply(comb, function(x){lapply(x,tile,width=1)})
  comb<-lapply(comb, function(x){lapply(x,function(x){lapply(x,function(x){mcols(x)<-paste0("C",1:(in_intron+in_exon));return(x)})})})
  comb<-lapply(comb, function(x){lapply(x,function(x){lapply(x,Count_ov,peak_file_gr)})})

  comb_df<-lapply(comb,function(x){lapply(x,function(x){data.frame(unlist(GRangesList(x)))})})
  comb_df<-lapply(comb_df,function(x){lapply(x,function(x){reshape2::dcast(x[c("Identifier","X","Ov")],Identifier~X)})})
  comb_df<-lapply(comb_df,function(x){lapply(x,function(x){rownames(x)<-x$Identifier;return(x[,-1])})})

  cols<-paste0("C",1:(in_intron+in_exon))
  comb_df<-lapply(comb_df,function(x){lapply(x,function(x){return(x[cols])})})

  #Combining the dataframes
  Retained<-cbind(comb_df$grouped_events_Upstream$Retained[rownames(comb_df$grouped_events_Upstream$Retained),],
                  comb_df$grouped_events_SE5$Retained[rownames(comb_df$grouped_events_Upstream$Retained),],
                  comb_df$grouped_events_SE3$Retained[rownames(comb_df$grouped_events_Upstream$Retained),],
                  comb_df$grouped_events_Downstream$Retained[rownames(comb_df$grouped_events_Upstream$Retained),])
  Skipped<-cbind(comb_df$grouped_events_Upstream$Skipped[rownames(comb_df$grouped_events_Upstream$Skipped),],
                  comb_df$grouped_events_SE5$Skipped[rownames(comb_df$grouped_events_Upstream$Skipped),],
                  comb_df$grouped_events_SE3$Skipped[rownames(comb_df$grouped_events_Upstream$Skipped),],
                  comb_df$grouped_events_Downstream$Skipped[rownames(comb_df$grouped_events_Upstream$Skipped),])
  Control<-cbind(comb_df$grouped_events_Upstream$Control[rownames(comb_df$grouped_events_Upstream$Control),],
                  comb_df$grouped_events_SE5$Control[rownames(comb_df$grouped_events_Upstream$Control),],
                  comb_df$grouped_events_SE3$Control[rownames(comb_df$grouped_events_Upstream$Control),],
                  comb_df$grouped_events_Downstream$Control[rownames(comb_df$grouped_events_Upstream$Control),])

  #Accounting for strand
  ret_negs<-which(sapply(rownames(Retained),function(x){strsplit(x,"_")[[1]][5]})=="-")
  ski_negs<-which(sapply(rownames(Skipped),function(x){strsplit(x,"_")[[1]][5]})=="-")
  ctrl_negs<-which(sapply(rownames(Control),function(x){strsplit(x,"_")[[1]][5]})=="-")

  Retained_neg<-rev(Retained[ret_negs,])
  Skipped_neg<-rev(Skipped[ski_negs,])
  Control_neg<-rev(Control[ctrl_negs,])

  Retained_pos<-Retained[-ret_negs,]
  Skipped_pos<-Skipped[-ski_negs,]
  Control_pos<-Control[-ctrl_negs,]

  colnames(Retained_neg)<-colnames(Retained_pos)
  colnames(Skipped_neg)<-colnames(Skipped_pos)
  colnames(Control_neg)<-colnames(Control_pos)

  Retained<-rbind(Retained_pos,Retained_neg)
  Skipped<-rbind(Skipped_pos,Skipped_neg)
  Control<-rbind(Control_pos,Control_neg)

  Retained_means<-colMeans(Retained)
  Skipped_means<-colMeans(Skipped)
  Control_means<-colMeans(Control)

  Region<-c(rep("Upsteam",(in_intron+in_exon)),
            rep("SE5",(in_intron+in_exon)),
            rep("SE3",(in_intron+in_exon)),
            rep("Downstream",(in_intron+in_exon)))

  Combined<-data.frame(rbind(cbind("Mean"=Retained_means,"Type"="Retained","Region"=Region,"Pos"=1:length(Retained_means)),
                             cbind("Mean"=Skipped_means,"Type"="Skipped","Region"=Region,"Pos"=1:length(Skipped_means)),
                             cbind("Mean"=Control_means,"Type"="Control","Region"=Region,"Pos"=1:length(Control_means))))

  Combined$Mean<-as.numeric(Combined$Mean)
  Combined$Pos<-as.numeric(Combined$Pos)
  Combined$Region<-factor(Combined$Region,levels=c("Upsteam","SE5","SE3","Downstream"))



  #Making the tables ready for plotting
  Combined2<-Combined
  Combined2$Pos[which(Combined2$Region=="SE5")]<-Combined2$Pos[which(Combined2$Region=="SE5")]+Gap_offset
  Combined2$Pos[which(Combined2$Region=="SE3")]<-Combined2$Pos[which(Combined2$Region=="SE3")]+Gap_offset*2
  Combined2$Pos[which(Combined2$Region=="Downstream")]<-Combined2$Pos[which(Combined2$Region=="Downstream")]+Gap_offset*3
  Combined2<-Combined2 %>% mutate(Color =
                                        case_when(Type == "Retained" ~ Retained_col,
                                                  Type == "Skipped" ~ Skipped_col,
                                                  Type == "Control" ~ Control_col))


  exon_s<-Get_Exon_st(df=Combined2,
                      exon_width=exon_width,
                      in_intron=in_intron,
                      in_exon=in_exon,
                      intron_col=intron_col,
                      SE_exon_col=SE_exon_col,
                      Surround_exon_col=Surround_exon_col,
                      Retained_col=Retained_col,
                      Skipped_col=Skipped_col,
                      Control_col=Control_col,
                      Gap=Gap_offset)

  labs<-data.frame("X"=c(1,in_exon,(in_exon+in_intron),(in_exon+in_intron+gap),(in_exon+in_intron*2+gap),(in_exon*2+in_intron*2+gap)),
                      "Y"=rep(-exon_width*1.02,6),
                      "Label"=c(-in_exon,1,in_intron,-in_intron,1,in_exon),
                      "HJust"=c(1,0,1,0,1,0),
                      "VJust"=c(1,1,1,1,1,1))


  g<-ggplot(Combined2, aes(Pos,Mean,color=Color)) +
    geom_line()+
    scale_color_identity()+
    geom_rect(inherit.aes = F,data=exon_s,mapping=aes(xmin=Start,xmax=End,ymin=Ystart,ymax=Yend,fill=Color),size = 1,show.legend = F)+
    scale_fill_identity()+
    theme_classic()+
    geom_segment(aes(y=0.0001, yend=max(Combined2$Mean), x=in_exon, xend=in_exon), linetype="dotted",col="Black")+
    geom_segment(aes(y=0.0001, yend=max(Combined2$Mean), x=(in_exon+in_intron*2+gap), xend=(in_exon+in_intron*2+gap)), linetype="dotted",col="Black")+
    geom_segment(aes(y=0.0001, yend=max(Combined2$Mean), x=(in_exon*3+in_intron*2+gap*2), xend=(in_exon*3+in_intron*2+gap*2)), linetype="dotted",col="Black")+
    geom_segment(aes(y=0.0001, yend=max(Combined2$Mean), x=(in_exon*3+in_intron*4+gap*3), xend=(in_exon*3+in_intron*4+gap*3)), linetype="dotted",col="Black")+
    geom_text(inherit.aes = F,data=labs, aes(x=X, y=Y,label = Label, hjust=HJust, vjust=VJust),size=5)+
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "bottom",
          legend.text=element_text(size=15),
          legend.title=element_blank())
  return(g)
}

Get_Exon_st<-function(df,exon_width,in_intron,in_exon,intron_col,SE_exon_col,Surround_exon_col,Retained_col,Skipped_col,Control_col,Gap){
  tot<-in_intron+in_exon
  Start<-c(1,in_exon,tot+gap,tot+gap+in_intron,tot*2+gap*2,tot*2+gap*2+in_exon,tot*3+gap*3,tot*3+gap*3+in_intron,tot,tot*2+gap,tot*3+gap*2)
  End<-c(in_exon,tot,tot+gap+in_intron,tot*2+gap,tot*2+gap*2+in_exon,tot*3+gap*2,tot*3+gap*3+in_intron,tot*4+gap*3,tot+gap,tot*2+gap*2,tot*3+gap*3)
  St<-c("Exon","Intron","Intron","Exon","Exon","Intron","Intron","Exon","Gap","Gap","Gap")
  comb<-data.frame(Start,End,St)
  comb$Ystart<-c(-exon_width,-exon_width*0.75,-exon_width*0.75,-exon_width,-exon_width,-exon_width*0.75,-exon_width*0.75,-exon_width,-exon_width,-exon_width,-exon_width)
  comb$Yend<-c(0,-exon_width*0.25,-exon_width*0.25,0,0,-exon_width*0.25,-exon_width*0.25,0,max(df$Mean),max(df$Mean),max(df$Mean))
  comb$Color<-c(Surround_exon_col,intron_col,intron_col,SE_exon_col,SE_exon_col,intron_col,intron_col,Surround_exon_col,"white","white","white")
  return(comb)
}




#Count Overlaps function without removing ranges
Count_ov<-function(gr,peak){
  ov<-findOverlaps(gr,peak)
  mcols(gr)$Ov<-0
  mcols(gr)$Ov[unique(queryHits(ov))]<-1
  mcols(gr)$Identifier<-names(gr)
  #if (sum(mcols(gr)$Ov)>0){print(names(gr)[1])}
  return(gr)
}


#Functions to return Granges
Upstream_Ranges<-function(df,int,ex){
  df$NS<-df$upstreamEE-(ex-1)
  df$NE<-df$upstreamEE+int
  df_gr<-makeGRangesFromDataFrame(df[c("seqnames","NS","NE","strand","Identifier","GeneID","geneSymbol")],seqnames.field = "seqnames",start.field = "NS",end.field = "NE",strand.field = "strand",keep.extra.columns=T)
  names(df_gr)<-df_gr$Identifier
  return(df_gr)
}

Downstream_Ranges<-function(df,int,ex){
  df$NS<-df$downstreamES-int
  df$NE<-df$downstreamES+(ex-1)
  df_gr<-makeGRangesFromDataFrame(df[c("seqnames","NS","NE","strand","Identifier","GeneID","geneSymbol")],seqnames.field = "seqnames",start.field = "NS",end.field = "NE",strand.field = "strand",keep.extra.columns=T)
  names(df_gr)<-df_gr$Identifier
  return(df_gr)
}

SE5_Ranges<-function(df,int,ex){
  df$NS<-df$start-int
  df$NE<-df$start+(ex-1)
  df_gr<-makeGRangesFromDataFrame(df[c("seqnames","NS","NE","strand","Identifier","GeneID","geneSymbol")],seqnames.field = "seqnames",start.field = "NS",end.field = "NE",strand.field = "strand",keep.extra.columns=T)
  names(df_gr)<-df_gr$Identifier
  return(df_gr)
}

SE53_Ranges<-function(df,int,ex){
  df$NS<-df$end-(ex-1)
  df$NE<-df$end+int
  df_gr<-makeGRangesFromDataFrame(df[c("seqnames","NS","NE","strand","Identifier","GeneID","geneSymbol")],seqnames.field = "seqnames",start.field = "NS",end.field = "NE",strand.field = "strand",keep.extra.columns=T)
  names(df_gr)<-df_gr$Identifier
  return(df_gr)
}

#Grouping Events
Group_events<-function(df,pval,FDR,psi,ctrl_grp,ctrl_size){
  if (ctrl_grp==1){

    skipdf=df[which(df$FDR<pval & df$PValue<FDR & df$IncLevelDifference>psi),]
    retdf=df[which(df$FDR<pval & df$PValue<FDR & df$IncLevelDifference<psi*(-1)),]
    ctrldf=df[which(df$Mean_Inc_Level_1>0.1 & df$Mean_Inc_Level_1<0.9 & df$FDR>FDR & df$PValue>pval),]

    if (nrow(ctrldf)>ctrl_size*(nrow(skipdf)+nrow(retdf))){
      ctrldf<-ctrldf[sample(1:nrow(ctrldf),ctrl_size*(nrow(skipdf)+nrow(retdf)),replace = FALSE),]
    }

  } else if (ctrl_grp==2){

    skipdf=df[which(df$FDR<pval &df$PValue<FDR & df$IncLevelDifference<(psi*(-1))),]
    retdf=df[which(df$FDR<pval & df$PValue<FDR & df$IncLevelDifference>psi),]
    ctrldf=df[which(df$Mean_Inc_Level_1>0.1 & df$Mean_Inc_Level_1<0.9 & df$FDR>FDR & df$PValue>pval),]

    if (nrow(ctrldf)>ctrl_size*(nrow(skipdf)+nrow(retdf))){
      ctrldf<-ctrldf[sample(1:nrow(ctrldf),ctrl_size*(nrow(skipdf)+nrow(retdf)),replace = FALSE),]
    }
  }
  event_list<-list("Retained"=retdf,"Skipped"=skipdf,"Control"=ctrldf)
  return(event_list)
}


Remove_Overlaps<-function(df,min_count,bedtools_path=bedtools_path){
  df$Identifier<-paste(df$chr,"_",df$exonStart_0base,"-",df$exonEnd,
                       "_",df$upstreamES,"-",df$upstreamEE,"_",
                       df$downstreamES,"-",df$downstreamEE,"_",df$strand,sep="")

  if (length(unique(df$Identifier)) != nrow(df)){
    stop ("Please remove duplicated events!")
  }

  df<-Sum_Counts(df)
  df<-df[which(df$Total_Count>min_count),]

  #back to genomic ranges
  df_gr<-makeGRangesFromDataFrame(df,
                                  start.field = "exonStart_0base",
                                  end.field = "exonEnd",
                                  keep.extra.columns = T)
  df_gr$cluster <- subjectHits(findOverlaps(df_gr, reduce(df_gr, min.gapwidth = 0L)))
  df_gr<-data.frame(df_gr)

  df_gr_trim<-df_gr %>%
    group_by(cluster) %>%
    slice_max(Total_IJC,n = 1,with_ties = F) %>%
    data.frame()

  return(df_gr_trim)
}

#function to calculate total counts from IJC and SJC columns
Sum_JC<-function(x){sum(as.numeric(strsplit(x,split = ",")[[1]]))}
Sum_Counts<-function(df){
  df$IJC_Sum_1<-sapply(df$IJC_SAMPLE_1,Sum_JC)
  df$IJC_Sum_2<-sapply(df$IJC_SAMPLE_2,Sum_JC)
  df$SJC_Sum_1<-sapply(df$SJC_SAMPLE_1,Sum_JC)
  df$SJC_Sum_2<-sapply(df$SJC_SAMPLE_2,Sum_JC)

  df$Total_IJC<-df$IJC_Sum_1+df$IJC_Sum_2
  df$Total_SJC<-df$SJC_Sum_1+df$SJC_Sum_2

  df$Mean_Inc_Level_1<-sapply(df$IncLevel1,function(x){mean(as.numeric(strsplit(x,",")[[1]]),na.rm=T)})
  df$Mean_Inc_Level_2<-sapply(df$IncLevel2,function(x){mean(as.numeric(strsplit(x,",")[[1]]),na.rm=T)})

  df$Total_Count_1<-df$IJC_Sum_1+df$SJC_Sum_1
  df$Total_Count_2<-df$IJC_Sum_2+df$SJC_Sum_2
  df$Total_Count<-df$Total_Count_1+df$Total_Count_2
  return(df)
}


