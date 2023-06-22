# The following function keeps rows where at least X samples have at least
# Y reads and outputs the processed data.frame adding the taxid column

filter_reads<- function(data, min.sample, min.reads) {
  require(dplyr)
  data$taxid<-rownames(data)
  data$filter<-0
  for (i in 1:dim(data)[1]){
    if(sum(data[i,1:length(data)]>=min.reads, na.rm = TRUE)<min.sample){
      data$filter[i] <- 1}   
  }
  data<-filter(data, filter==0)
  data$filter<-NULL
  return(data)
}

### Cut counts table at given taxonomic rank (sum reads below given level)
tax_level<-function(counts, rank, samples_column, show_parent=F) {
  if (show_parent==F)
    counts=counts
  else counts[,rank]<-paste(counts[,which(colnames(counts)==rank )], " (", counts[,which(colnames(counts)==rank )-1], ")", sep="")
  filtered<-data.frame(lapply(counts[,samples_column],
                              function(x) tapply(x, counts[,rank], sum, na.rm=TRUE)))
  filtered<-filtered[rowSums(is.na(filtered)) != ncol(filtered),]
  return(filtered)
}

### Traspose, order for reads number and select N top-taxa and add the "Other" category. taxa as row & samples as column table as input
reduce_taxa<-function(df, N, method="sum", keep_rownames=F) {
  require(reshape2)
  df_t<-t(df)
  {
  if (method=="sum")
  df_t_ord<-df_t[,order(colSums(df_t), decreasing=TRUE)]
  if (method=="mean")
  df_t_ord<-df_t[,order(colMeans(df_t), decreasing=TRUE)]
  }
  taxa_list<-colnames(df_t_ord)[1:N]
  df_reduced<-data.frame(df_t_ord[,colnames(df_t_ord) %in% taxa_list], Others=rowSums(df_t_ord[,!colnames(df_t_ord) %in% taxa_list], na.rm = T))
  if (keep_rownames==F) 
   return(df_reduced)
  else colnames(df_reduced)[1:N]<-colnames(df_t_ord)[1:N]
  return(df_reduced)
}
