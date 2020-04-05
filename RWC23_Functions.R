# matrix_select() subsets a matrix based on vactor of gene names.
# Pleasure ensure gene naming is consistent between files.
# It is best to use WBGeneIDs for this process

matrix_select <- function(count_matrix, gene_subset_vector){
  nameselect <- rownames(count_matrix) %in% gene_subset_vector
  nameMatrix <- count_matrix[nameselect,]
  return(nameMatrix)
}

# row_normalize_matrix() performs a row normalization by subtracting the row mean from each value in the row.

row_normalize_matrix <- function(count_matrix){
  namevarRowNormalized <- count_matrix - rowMeans(count_matrix)
  return(namevarRowNormalized)
}

# row_normalize_matrix_cutoff() subsets the count matrix based on value in variance_cutoff
# It then performs a row normalization by subtracting the row mean from each value in the row.

row_normalize_matrix_cutoff <- function(count_matrix, variance_cutoff){
  varselect <- rowVars(count_matrix) >= variance_cutoff
  count_matrix <- count_matrix[varselect, ]
  namevarRowNormalized <- count_matrix - rowMeans(count_matrix)
  return(namevarRowNormalized)
}

# row_zscore_matrix() assigns a row Z score to each row value. 

row_zscore_matrix <- function(count_matrix){
  namevarZscore <- (count_matrix - rowMeans(count_matrix))/rowSds(count_matrix)
  return(namevarZscore)
  
}

# row_zscore_matrix_cutoff() assigns a row Z score to each row value. 
# Then uses variance to subset data.

row_zscore_matrix_cutoff <- function(count_matrix, variance_cutoff){
  varselect <- rowVars(count_matrix) >= variance_cutoff
  count_matrix <- count_matrix[varselect, ]
  namevarZscore <- (count_matrix - rowMeans(count_matrix))/rowSds(count_matrix)
  return(namevarZscore)
  
}

# id2name() replaces WBGeneIDs in The input count matrix to 
# gene names for easier interpretation of genes on a heatmap.
# This function requires access to the internet for proper funcitoning.

id2name <- function(count_matrix){
  id2name.df <- suppressMessages(biomaRt::getBM(mart = paramart, 
                                                filter=c("wbps_gene_id"), 
                                                value=rownames(count_matrix), 
                                                attributes = c('wbps_gene_id','external_gene_id')
  ))
  rownames(count_matrix) <- id2name.df$external_gene_id[ 
    match(rownames(count_matrix), id2name.df$wbps_gene_id)
    ]
  return(count_matrix)
}

# myPheatmap() is a custom pheatmap plotting function

myPheatmap <- function(count_matrix){
  pheatmap(count_matrix, 
           cluster_cols = FALSE, 
           cluster_rows = TRUE, 
           show_rownames = FALSE, 
           border_color = NA)
}
