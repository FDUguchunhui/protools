#' @export
#' @title Plot heatmap for SpC_list object with customized setting
#' @description This is an wrapper function that use pheatmap to plot heatmap
#' for the SpC_list object and provide additional data processing and plotting setting that
#' is useful for proteomics spectral count data
#' @param x
#' @param upper_limit
#' @param log_scale
#' @param row_keep
#' @param show_rownames 
#' @param show_colnames, show_colnames, cluster_cols, cluster_rows, ... pheatmap 
#' parameters that will be used when calling pheatmap function


plot_heat_map <- function(x, 
                          upper_limit=NULL, 
                          log_scale=TRUE, 
                          row_keep=NULL, 
                          show_rownames = F,
                          show_colnames = F,
                          cluster_cols =  F,
                          cluster_rows = T,
                          ...) {
  
  
  data_matrix <- x$matrix
  
  if (!is.null(row_keep)) {
    data_matrix <- data_matrix[row_keep, ]
  }
  
  if (log_scale == TRUE) {
    data_matrix <- log2(data_matrix + 1)
  }
 
  
  # if there is a upper limit for heatmap intensity 
  # intensity higher than the upper_limit will be the rightmost color in the color palette
  if (!is.null(upper_limit)) {
    breaksList = seq(0, upper_limit, by = 0.05)
  } else {
    breaksList = seq(0, max(data_matrix), by = 0.05)
  }
  
  pheatmap::pheatmap(data_matrix,
                     color = colorRampPalette(colors = c("black", "yellow", "darkorange"))(length(breaksList)),
                     breaks = breaksList,
                     #fontsize_col = 5, 
                     #fontsize_row = 4,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     cluster_cols =  cluster_cols,
                     cluster_rows = cluster_rows,
                     annotation_col=y$annotation,
                     ...)
  
} 




basic_describe <- function(x, boxplot_main='Distribution', piechart_main='Group', plot_boxplot=TRUE, plot_piechart=TRUE, plot_heatmap=FALSE, ...) {
  
  # boxplot for spectral count distribution in each replicate
  if (plot_boxplot) {
    boxplot(y$matrix, names=FALSE, main=boxplot_main)
  }
 
  if (plot_piechart) {
    # piechart for proportion of each cancer type in the data
    pie_data <- x$annotation %>% group_by(disease) %>% summarise(count=n()) %>% arrange(desc(count))
    pie(pie_data$count, labels = paste0(pie_data$disease, '(', pie_data$count, ')'), main=piechart_main)
  }
  
  
  # heatmap for the entire dataset
  if (plot_heatmap) {
    plot_heat_map(x, ...)
  }
 
  
}












