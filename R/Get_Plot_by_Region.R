
Get_Plot_by_Region<-function(Chr,
                             Start,
                             End,
                             Strand,
                             bed,
                             gene,
                             rank_,
                             density=F,
                             utr_col="dark gray",
                             peak_col="Blue",
                             peaks_width=0.3,
                             exon_width=0.5,
                             utr_width=0.3,
                             exon_col="black"){

  bed$rank<-match(bed$group_name,rev(rank_))

  #adding coordinates for the peaks
  bed$y_start<-peaks_width*bed$rank-peaks_width
  bed$y_end<-peaks_width*bed$rank
  bed$col<-peak_col

  #Plotting gene structure
  if (gene$gene_biotype[1]=="protein_coding"){
    Gene_s<-Gene_strcuture_prot(gtf=gene,levels=max(bed$rank),peaks_width=peaks_width,
                                exon_width=exon_width,utr_width=utr_width,exon_col=exon_col,utr_col=utr_col)
  } else {
    Gene_s<-Gene_strcuture_nc(gtf=gene,levels=max(bed$rank),peaks_width=peaks_width,
                              exon_width=exon_width,utr_width=utr_width,exon_col=exon_col,utr_col=utr_col)
  }

  #Setting the Start and End for the plot area
  Start_min<-min(c(Start,min(Gene_s$start)))
  End_max<-max(c(End,max(Gene_s$end)))

  #Background for making sure the peaks are clearly separated
  Bg_tab<-Background_Table(df=bed,Start=Start_min,End=End_max)

  Intron_s<-Gene_s[which(Gene_s$type=="intron"),]
  Gene_s<-Gene_s[which(Gene_s$type!="intron"),]
  Intron_s$start<-Start_min
  Intron_s$end<-End_max

  bed$xpos<-Start_min-50

  #Get a df for labeling 5' and 3'
  if (Gene_s$strand[1]=="+"){
    strand_anot_left<-data.frame("Label"="5'","X"=Start_min-50,"Y"=c(min(Gene_s$y_start)))
    strand_anot_right<-data.frame("Label"="3'","X"=End_max+50,"Y"=c(min(Gene_s$y_start)))
  } else if (Gene_s$strand[1]=="-"){
    strand_anot_left<-data.frame("Label"="3'","X"=Start_min-50,"Y"=c(min(Gene_s$y_start)))
    strand_anot_right<-data.frame("Label"="5'","X"=End_max+50,"Y"=c(min(Gene_s$y_start)))
  }

  area_rect<-data.frame("xstart"=Start,"xend"=End,"ystart"=min(Bg_tab$y_start),"yend"=max(Gene_s$y_end)+0.1)

  #Plotting Gene
  g<-ggplot(Gene_s)+
    scale_x_continuous(name="x", limits = c(Start_min-500,End_max+500))+
    theme_classic()+
    geom_rect(inherit.aes = F,data=Intron_s,aes(xmin=start,xmax=end,ymin=y_start,ymax=y_end),fill=Intron_s$col)+
    geom_rect(data=Gene_s,aes(xmin=start-1,xmax=end+1,ymin=y_start,ymax=y_end),fill=Gene_s$col)+
    #Gene structure done
    theme(strip.background = element_blank()) +
    theme(legend.position="bottom", legend.direction = "horizontal", legend.text = element_text(size=5))+
    ggtitle(paste(gene$gene_name[1],"(",gene$transcript_id[1],")",sep=""))+
    geom_rect(data=Bg_tab,mapping=aes(xmin=x_start,xmax=x_end,ymin=y_start,ymax=y_end),fill=Bg_tab$col,color="black")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          legend.text = element_text(size=10, angle = 90),
          panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
          plot.background = element_rect(fill = "transparent",colour = NA),
          plot.title = element_text(hjust=0.5,size=30))+
    geom_rect(inherit.aes = FALSE,data=bed,mapping=aes(xmin=start,xmax=end,ymin=y_start,ymax=y_end),fill=bed$col,alpha=0.5)+
    geom_text(data = bed[!duplicated(bed$group_name),],aes(label=group_name,x=xpos,y=y_start+peaks_width/2,hjust=1),size=5)+
    geom_text(data = strand_anot_left,aes(label=Label,x=X,y=Y,hjust=1),size=5)+
    geom_text(data = strand_anot_right,aes(label=Label,x=X,y=Y,hjust=0),size=5)+
    ggplot2::geom_vline(xintercept = c(Start,End),linetype="dotted")
    #geom_line(inherit.aes = FALSE,data=area_rect,mapping=aes(xmin=xstart,xmax=xend,ymin=ystart,ymax=yend),)
  return(g)
}

#Gene Strcuture-protein coding genes
Gene_strcuture_prot<-function(gtf,levels,peaks_width,exon_width,utr_width,exon_col,utr_col){
  gtf<-gtf[which(gtf$type %in% c("CDS","five_prime_utr","three_prime_utr")),]
  exon_center<-levels*peaks_width+0.2+exon_width/2

  gtf$y_start<-ifelse(gtf$type=="CDS",exon_center-exon_width/2,exon_center-utr_width/2)
  gtf$y_end<-ifelse(gtf$type=="CDS",exon_center+exon_width/2,exon_center+utr_width/2)
  gtf$col<-ifelse(gtf$type=="CDS",exon_col,utr_col)
  gtf<-gtf[c("seqnames","type","start","end","strand","y_start","y_end","col")]
  Intron_s<-data.frame("seqnames"=gtf$seqnames[1],"type"="intron","start"=min(gtf$start),
                       "end"=max(gtf$end),"strand"=gtf$strand[1],
                       "y_start"=exon_center-utr_width/4,"y_end"=exon_center+utr_width/4,"col"="light gray")
  gtf<-rbind(gtf,Intron_s)
  return(gtf)
}


#Gene Strcuture-non coding genes
Gene_strcuture_nc<-function(gtf,levels,peaks_width,exon_width,utr_width,exon_col,utr_col){
  gtf<-gtf[which(gtf$type %in% c("exon")),]
  exon_center<-levels*peaks_width+0.2+exon_width/2
  gtf$y_start<-exon_center-exon_width/2
  gtf$y_end<-exon_center+exon_width/2
  gtf$col<-exon_col
  gtf<-gtf[c("seqnames","type","start","end","strand","y_start","y_end","col")]
  Intron_s<-data.frame("seqnames"=gtf$seqnames[1],"type"="intron","start"=min(gtf$start),
                       "end"=max(gtf$end),"strand"=gtf$strand[1],
                       "y_start"=exon_center-utr_width/4,"y_end"=exon_center+utr_width/4,"col"="light gray")
  gtf<-rbind(gtf,Intron_s)
  return(gtf)
}



#Caclualte Density of Peaks
Density<-function(Chr,Start,End,df){
  peaks_gr<-makeGRangesFromDataFrame(df,keep.extra.columns = T,
                                     start.field = "Start",
                                     end.field = "End",
                                     seqnames.field = "seqnames")
  bins<-GRanges(seqnames=Chr,strand="*",
                ranges = IRanges(start = seq(Start,End,by=10),width = 10))
  Overlaps<-countOverlaps(bins,peaks_gr)
  values(bins)<-DataFrame(Density=Overlaps)
  bins<-as.data.frame(bins)
  return(bins)
}


#Preapre Background
Background_Table<-function(df,Start,End){
  Background<-df[c("y_start","y_end")]
  Background<-Background[!duplicated(Background),]
  Background<-Background[order(Background$y_start),]
  Background$x_start<-Start
  Background$x_end<-End
  Background$col<-rep(c("#33333333","#11111111"),length.out = nrow(Background))
  return(Background)
}


