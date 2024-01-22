

remove_NA_rows <- function(data, fraction_NAs = 0.5){
  # remove all rows in a matrix with fraction of NAs greater than fraction_NAs (default is 0.5)
  
  na_fraction <- apply(data, 1, function(x) {sum(is.na(x))/length(x)})
  row.names <- rownames(data)
  new_matrix <- data[which(na_fraction <= fraction_NAs), ,drop=FALSE]
  rownames(new_matrix) <- unlist(row.names[which(na_fraction <= fraction_NAs)])

  return(new_matrix)
  
}




replace_all_NAs_with_row_means <- function(matrix){
  # replaces all NAs in a matrix with the row means
  output_matrix <- t(apply(matrix, 1, function(row){row[is.na(row)] <- mean(row, na.rm = TRUE) ; return(row)}))
  
  return(as.data.frame(output_matrix))
}





plot_PCA <- function(matrix, dataset_name, labels, label_name, pcs= c(1,2)){
  # plots a PCA of the data coloured by labels
  pca_data <- list(BiocSingular::runSVD(t(matrix),k=10,center=T))[[1]]
  
  percentage  <- pca_data$d ^ 2 / sum(pca_data$d ^ 2) * 100
  percentage  <- sapply(seq_along(percentage), function(i)
    {round(percentage [i], 1)})
  
  data <- data.frame(pca_data$u)
  colnames(data) <- paste0('PC',c(1:10))
  data$label <- labels
  
  first_pc <- sym(paste0('PC',pcs[1]))
  second_pc <- sym(paste0('PC',pcs[2]))
  
  xlab <- paste0(first_pc,' (',percentage[pcs[1]],'% variance)')
  ylab <- paste0(second_pc,' (',percentage[pcs[2]],'% variance)')
  
  
  p <- ggplot(data, mapping= aes(x = !!first_pc, y = !!second_pc, fill= factor(label))) +
    geom_point(aes(fill=label), color = 'grey55', alpha=0.72 , shape=21) +
    scale_fill_discrete(name=label_name) +
    labs(title = paste0("PCA of ",dataset_name," grouped by ",label_name), x = xlab, y = ylab) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme_minimal() +
    theme(legend.background = element_blank(), 
          aspect.ratio = 1, 
          legend.key = element_blank(),
          legend.position = "right",
          panel.background = element_blank(), 
          axis.line = element_line(colour = "grey35",linewidth = 1.1, lineend='round'),
          panel.grid.major = element_line(color = "grey96"))
  
  p <- p + theme(legend.position = "right")
  
  xdens <- axis_canvas(p, axis = "x")+
    geom_density(data, mapping = aes(x = !!first_pc, fill = labels), color= 'grey45', alpha = 0.50, size = 0.2) +
    theme(legend.position = "none")
  
  ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
    geom_density(data, mapping = aes(x = !!second_pc, fill = labels), color= 'grey45', alpha = 0.50, size = 0.2) +
    theme(legend.position = "none")+
    coord_flip()
  
  p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  
  pList <- ggdraw(p2)

  return(pList)
  
}






remove_all_low_number_rows <- function(matrix, number){
  # removes all rows in a matrix with max value less than the number
  max_values <- apply(matrix, 1, max)  # Calculate the max of each row
  
  # Subset the matrix to keep only rows with a max value of less than 3
  filtered_matrix <- matrix[max_values > number, ]
  
  return(as.data.frame(filtered_matrix))
}





lineup_samples <- function(matrix_list, direction = 'columns'){
  # taking a list of matrices and direction ('rows' or 'columns'), find common samples in all matrices, discard all 
  # others, and line up all the common samples in the same order in all matrices.
  # outputs a list of the processed matrices.
  
  if( direction == 'columns'){
    common_samples <- Reduce(intersect, lapply(matrix_list,colnames))
    matrix_list_aligned <- lapply(matrix_list, function(matrix) {matrix[, common_samples, drop = FALSE]})
  }
  
  if( direction == 'rows'){
    common_samples <- Reduce(intersect, lapply(matrix_list,rownames))
    matrix_list_aligned <- lapply(matrix_list, function(matrix) {matrix[common_samples, , drop = FALSE]})
  }
    
  else{stop("The direction you entered is invalid. You can only enter 'rows' or 'columns'.")}
  
  return(matrix_list_aligned)
}






remove_all_zero_variance_rows <- function(matrix){
  # taking a matrix, removes all rows with 0 variance
  variances <- apply(matrix, 1, var)
  
  # Identify rows with zero variance
  rows_to_keep <- variances != 0
  
  # Subset the matrix to keep only rows with non-zero variance
  return(matrix[rows_to_keep, , drop = FALSE])
  
}






compute_log_transformed_RLE <- function(matrix){
  # taking log transformed matrix, output a matrix of RLE values.
  # try to use rowMedians() instead
  medians <- apply(matrix, 1, median)
  
  RLE_scores <- sweep(matrix, 1, medians, FUN='-')
  
  return(RLE_scores)
  
}






plot_RLE <- function(RLE_matrix, batch_info, ylimit=c(-10,10)){
  # taking a matrix of RLE scores, and a vector of batch information (in the same order as your samples), and y axis limits as a vector of 2 numbers, output an RLE graph ordered by batch
  
  RLE_long <- reshape2::melt(RLE_matrix)
  
  RLE_long$Batch <- factor(rep(batch_info, each= nrow(RLE_matrix)))
  
  RLE_long <- RLE_long %>%
    arrange(Batch, variable) %>%
    mutate(variable = factor(variable, levels = unique(variable)))
  
  ggplot(RLE_long, aes(x = variable, y = value, fill = Batch)) +
    geom_boxplot(outlier.shape = NA,alpha=0.9, linewidth=0.5) +  # Set coef to Inf to extend whiskers to max/min + 
    theme(axis.text.x = element_blank())+
    labs(x = "Sample", y = "RLE Score", fill = "Batch", color = 'Batch')+
    coord_cartesian(ylim=ylimit)+
    theme_minimal() +
    theme(panel.border = element_rect(colour = "grey85", fill=NA, size=1.1),axis.text.x = element_blank(),
          panel.grid.major = element_line(color = "grey96"),)
  
}






plot_linear_correlation <- function(matrix_list, dataset_names, variable_vector, variable_name, num_pcs=10){
  # taking a list of matrices, a vector of the names of the omics, a vector of numerical values (in order) of the values you want to compare, and number of PCs, plot the linear correlation R^2 values

  if (length(matrix_list)>30) {stop('Your omics matrices must be a list. Did you accidentally use c() instead of list()?')}
  
  if (is.list(dataset_names)) {stop('Your dataset names cannot be a list. Has to be a vector. use c() instead of list()!')}
  
  pcas <- lapply(matrix_list, function(m) {prcomp(t(m), scale = F, center=T)})
  
  r_squared_results <- sapply(pcas, function(matrix, variable, num) {
    sapply(1:num_pcs, function(y) {
      lm_model <- summary(lm(variable ~ matrix$x[,1:y]))$r.squared})
  }, variable = variable_vector, num=num_pcs)
  
  
  PCs <- rep(paste0('PC1:', 1:num_pcs), length(matrix_list))
  Datasets <- rep(dataset_names, each = num_pcs)
  R_squared <- as.vector(r_squared_results)
  
  pc_correlations <- data.frame(PCs, R_squared, Datasets)
  pc_correlations$PCs <- factor(pc_correlations$PCs, levels = paste0('PC1:', 1:num_pcs))
  
  ggplot(pc_correlations, aes(x = PCs, y = R_squared, color = Datasets, group = Datasets)) +
    geom_line(size=0.5,alpha=0.8) +
    geom_point(alpha=0.8) +
    labs(x = 'Principal Components', y = 'R-squared', color = 'Dataset') +
    ylim(0,1) +
    ggtitle(paste0('Correlation between PC1 and the first ',num_pcs,' PCs with ',variable_name))+
    theme_minimal() +
    theme(axis.line = element_line(colour = "grey83", linewidth=1.1),
          panel.border = element_rect(colour = "grey90", fill=NA, size=0.7),
          panel.grid.major = element_line(color = "grey96"),
          aspect.ratio = 1/1.2)
}

  
  
  


plot_vector_correlation <- function(pca_list, dataset_names, variable_vector, variable_name, num_pcs=10){
  # taking a list of pca objects, a vector of categorical values (in order) of the values you want to compare, and number of PCs, plot the correlation R^2 values
  
  if (length(pca_list)>30) {stop('Your pcas must be a list. Did you accidentally use c() instead of list()?')}
  
  if (is.list(dataset_names)) {stop('Your dataset names cannot be a list. Has to be a vector. use c() instead of list()!')}
  
  dummies <- fastDummies::dummy_cols(variable_vector)[,-1]
  
  cancor_scores <- sapply(pca_list, function(m) {lapply(1:num_pcs, function(y) {cca <- stats::cancor(x= m$x[,1:y, drop=F], 
                      y= dummies)  
  1 - prod(1 - cca$cor^2)})
  })
  
  PCs <- rep(paste0('PC1:', 1:num_pcs), length(pca_list))
  Datasets <- rep(dataset_names, each = num_pcs)
  cor <- unlist(as.vector(cancor_scores))
  
  pc_correlations <- data.frame(PCs, cor, Datasets)
  pc_correlations$PCs <- factor(pc_correlations$PCs, levels = paste0('PC1:', 1:10))
  
  
  ggplot(pc_correlations, aes(x = PCs, y = cor, color = Datasets, group = Datasets)) +
    geom_line(size=0.5,alpha=0.8) +
    geom_point(alpha=0.8) +
    labs(x = 'Principal Components', y = 'Correlation', color = 'Dataset') +
    ylim(0,1) +
    theme_minimal() +
    theme(axis.line = element_line(colour = "grey83", linewidth = 1.1),
          panel.border = element_rect(colour = "grey90", fill=NA, size=0.7),
          panel.grid.major = element_line(color = "grey96"),
          aspect.ratio = 1/1.2) +
    ggtitle(paste0('Vector correlation between PC1 and the first ',num_pcs,' PCs with ',variable_name))
  
}





plot_boxplot_categorical <- function(data_vector, category_vector, names){
  # taking a vector of numerical scores and a vector of categorical variables and their names as vector, plot boxplot grouped by the categorical variables.
  data <- as.data.frame(cbind(data_vector, category_vector))
  colnames(data) <- names
  
  ggplot(data, aes(x=data[,2], y=as.numeric(data[,1]), fill=data[,2]))+
    geom_boxplot()+
    theme(axis.text.x = element_blank())+
    labs(x = names[2], y = names[1], fill=names[2] )+
    ggtitle(paste0('Boxplot of ',names[1],' grouped by ',names[2]))
  
  
}





get_high_correlation_features <- function(matrix, variable_vector, threshold, method = 'pearson'){
  # taking a matrix and a vector (numerical), find features (rows) that have absolute correlation score (both positive and negative) greater than your set threshold value with your variable vector. 
  correlations <- lapply(1:nrow(matrix), function(x) {cor(as.vector(as.numeric(matrix[x,])),as.vector(variable_vector), method=method)})
  
  correlations <- as.matrix(unlist(correlations))
  
  rownames(correlations) <- rownames(as.matrix(matrix))
  
  high_corr_features <- correlations[which(abs(correlations) > threshold), ]
  
  return(high_corr_features)
  
}






ftest <- function(matrix, variable, is.log=T, n.cores=8){
  # taking a matrix of data and a vector of categorical variables, calculate F and P values for each feature.
  # Note that your data is first log transformed. If you set is.log =F, we'll transform it first and then do the ANOVA.
  
  data <- matrix
  if (!is.matrix(data)){stop(paste0("We detect your data is a '",type(matrix),"' not a matrix. Please double check"))}
  if (is.log) data <- data
  else {data <- log(data + 1)}
  average.exp <- rowMeans(data)
  f.test <- parallel::mclapply( 1:nrow(data), function(x) 
    {MASS::dropterm(lm(data[x , ] ~ variable), test = 'F')[c(5:6)]}, mc.cores = n.cores)
  
  f.test <- data.frame(
    FValue = round(unlist(lapply(f.test, function(x) {x$`F Value`[2]})), digits = 4) ,
    PValue = unlist(lapply(f.test, function(x) {x$`Pr(F)`[2]})),
    Adj.PValue = p.adjust(unlist(lapply(f.test, function(x) {x$`Pr(F)`[2]})), method = 'BH'),
    Mean = round(average.exp, digits = 2))
  return(f.test)
}





ANOVA_2way <- function(matrix, variable_1, variable_2, is.log=T, n.cores=8, pval_adj_method = 'BH', sort_by = NULL) {
  # Taking a matrix and two categorical variable vectors, conducts 2 way ANOVA for every row
  # returns F score, P-value, adjusted P-value and mean of each row. 
  
  #Check data
  if (!is.matrix(matrix)){stop(paste0("We detect your data is a '",class(matrix),"' not a matrix. Please double check"))}
  if (!is.log) {matrix <- log(matrix + 1)}
  
  data <- as.data.frame(matrix)
  
  # Applying the two_way_anova function across all genes
  results <- parallel::mclapply(1:nrow(data), function(x) {fit <- aov(as.numeric(data[x,]) ~ variable_1*variable_2)
  return(summary(fit)[[1]])}, mc.cores = n.cores)
  
  # Format the results
  results <- data.frame(
    Interaction_F_score = round(unlist(lapply(results, function(x) {x[['F value']][3]})), digits = 4),
    Interaction_PValue = unlist(lapply(results, function(x) {x[['Pr(>F)']][3]})),
    Adj.PValue = p.adjust(unlist(lapply(results, function(x) {x[['Pr(>F)']][3]})), method = pval_adj_method),
    Mean = rowMeans(matrix))
  
  if (!is.null(sort_by)) {
    if (!sort_by %in% colnames(results)){stop("The column you want to sort by doesn't exist. You can choose from:'Interaction_F_score','Interaction_PValue','Adj.PValue', 'Mean'.")}
    else {results <- results[order(results[[sort_by]]),]}
  }
  
  return(results)
}






wilcoxon_test <- function(matrix, variable, is.log=T, n.cores = 8){
  # 
  # 
  
  if (!is.matrix(matrix)){stop('Inevitable in retrospect. Please uninstall R and quit bioinformatics.')}#paste0("We detect your data is a '",type(matrix),"' not a matrix. Please double check"))}
  if(is.log){expr.data <- matrix}
  else{expr.data <- log(matrix + 1)}
  pval <- parallel::mclapply(row.names(expr.data), function(x) stats::wilcox.test(expr.data[x ,] ~ variable)[[3]], mc.cores = n.cores)
  results <- data.frame(genes = row.names(expr.data),
                        Pvalue = unlist(pval),
                        Adj.pvalue = p.adjust(p = unlist(pval), method = 'BH'))
  
  
  return(results)
}
  
  





kruskal_wallis_test <- function(matrix, variable, is.log=T, n.cores=8, pval_adj_method = 'BH', sort_by = NULL){
  # Taking a matrix of data and a vector of categorical variables, calculate Kruskal Wallis chisq value and P values for each feature.
  # Note that if you put is.log=False, we'll transform it first and then do the ANOVA.
  
  data <- matrix
  if (!is.matrix(data)){
    if (!is.data.frame(data)){ stop(paste0("We detect your data is a '",type(matrix),"' not a matrix or dataframe. Please double check"))}}
  if (is.log) data <- data
  else {data <- log(data + 1)}
  means <- rowMeans(data)
  
  k.test <- parallel::mclapply(1:nrow(data), function(x) {kruskal.test(as.numeric(data[x,]), variable)[c(1,3)]})
  
  k.test <- data.frame(
    kruskal_wallis_chisq = round(unlist(lapply(k.test, function(x) {x[1]})), digits = 4) ,
    PValue = unlist(lapply(k.test, function(x) {x[2]})))
  
  k.test$Adj.PValue <- p.adjust(k.test$PValue, method = pval_adj_method)
  k.test$Mean <- round(means, digits = 2)
  
  rownames(k.test) <- rownames(matrix)
  
  if (!is.null(sort_by)) {
    if (!sort_by %in% colnames(k.test)){stop("The column you want to sort by doesn't exist. You can choose from:
                                              'kruskal_wallis_chisq','PValue','Adj.PValue','Mean'.")}
    else {k.test <- k.test[order(k.test[[sort_by]]),]}
    }
  
  return(k.test)
}



ScheirerRayHare_test <- function(matrix, variable_1, variable_2, is.log=T, n.cores=8, pval_adj_method = 'BH', sort_by = NULL){
  # Non-parametric version of the two-way ANOVA. 
  
  
  
  
  
}





plot_UMAP <- function(matrix, metadata_vector, matrix_name, metadata_name){
  # taking a matrix and a vector of metadata, plot UMAP of the matrix colored by the groups in the metadata vector
  
  umap_result = umap(t(matrix))
  df = data.frame(umap_result$layout)
  colnames(df) = c("UMAP1", "UMAP2")
  df$metadata = metadata_vector  # Replace with your batch or metadata vector
  
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = metadata)) +
    geom_point() +
    ggtitle(paste0('UMAP of ',matrix_name,' grouped by ',metadata_name))
  
}



plot_tSNE <- function(matrix, metadata_vector, matrix_name, metadata_name){
  # taking a matrix and a vector of metadata, plot tSNE of the matrix colored by the groups in the metadata vector
  
  tsne_result = Rtsne(t(matrix))
  df = data.frame(tsne_result$Y)
  colnames(df) = c("tSNE1", "tSNE2")
  df$metadata = metadata_vector  # Replace with your batch or metadata vector
  
  ggplot(df, aes(x = tSNE1, y = tSNE2, color = metadata)) +
    geom_point() +
    ggtitle(paste0('tSNE of ',matrix_name,' grouped by ',metadata_name))
  
}











basic_count_normalization <- function(matrix, type='', log.transform = F){
  # taking a matrix of raw counts, convert to FPKM, FPKM_UQ or TPM. 
  # if log.transform is set to True the output matrix will also be log transformed for you (log base 2 + 1 pseudo count)
  
  
  
  
}



