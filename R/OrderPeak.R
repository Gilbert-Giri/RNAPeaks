
OrderPeak<-function(bed,order_by){
  if (order_by=="Region"){
    rank<-By_region(bed)
  } else if (order_by=="Target"){
    rank<-By_name(bed)
  } else if (order_by=="Count") {
    rank<-By_count(bed)
  } else {
    stop("Please provide acceptable feature to rank by: Target,Count")
  }
  return(rank)
}


#Order by name
By_name<-function(bed){
  return(sort(unique(bed$group_name)))
}


#Order by number of peaks
By_count<-function(bed){
  df<-names(sort(table(bed$group_name),decreasing = T))
  return(df)
}



#Cluster by most similar
By_region<-function(bed){
  bins<-IRanges(start = min(bed$start):max(bed$end),width=1,names = min(bed$start):max(bed$end))
  regions_by_target<-list()
  for (i in 1:nrow(bed)){
    ran<-IRanges(start = bed$start[i]:bed$end[i],width=1)
    if (!(bed$group_name[i] %in% names(regions_by_target))){
      regions_by_target[[bed$group_name[i]]]<-ran
    } else if (bed$group_name[i] %in% names(regions_by_target)){
      regions_by_target[[bed$group_name[i]]]<-c(regions_by_target[[bed$group_name[i]]],ran)
    }
  bins_df<-data.frame("Region"=names(bins))
  for (i in 1:length(names(regions_by_target))){
    m<-rep(0,length(bins))
    m[queryHits(findOverlaps(bins,regions_by_target[[names(regions_by_target)[i]]]))]<-1
    bins_df<-cbind(bins_df,m)
  }
  rownames(bins_df)<-bins_df$Region
  bins_df<-bins_df[,-1]
  colnames(bins_df)<-names(regions_by_target)
  }
}



