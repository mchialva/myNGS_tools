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
LCA_level<-function(counts, rank, samples_column) {
  filtered<-data.frame(lapply(counts[,samples_column],
                              function(x) tapply(x, counts[,rank], sum, na.rm=TRUE)))
  filtered<-filtered[rowSums(is.na(filtered)) != ncol(filtered),]
  return(filtered)
}

### Traspose, order for reads number and select N top-taxa and add the "Other" category. taxa as row & samples as column table as input
reduce_taxa<-function(df, N) {
  require(reshape2)
  df_t<-t(df)
  df_t_ord<-df_t[,order(colSums(df_t), decreasing=TRUE)]
  taxa_list<-colnames(df_t_ord)[1:N]
  df_reduced<-data.frame(df_t_ord[,colnames(df_t_ord) %in% taxa_list], Others=rowSums(df_t_ord[,!colnames(df_t_ord) %in% taxa_list], na.rm = T))
  return(df_reduced)
}
