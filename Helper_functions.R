

remove_NA_rows <- function(data, fraction_NAs = 0.5){
  # remove all rows in a matrix with fraction of NAs greater than fraction_NAs (default is 0.5)
  
  na_fraction <- apply(data, 1, function(x) {sum(is.na(x))/length(x)})
  new_matrix <- data[which(na_fraction <= fraction_NAs),]

  return(new_matrix)
  
}




replace_all_NAs_with_row_means <- function(matrix){
  # replaces all NAs in a matrix with the row means
  #rowmeans <- rowMeans(matrix, na.rm = T)
  output_matrix <-apply(matrix, 1, function(row){row[is.na(row)] <- mean(row, na.rm = TRUE); return(row)})
  #output_list <- parallel::mclapply(1:nrow(matrix), function(i) {
  #  row <- matrix[i, ]
  #  row[is.na(row)] <- rowmeans[i]
  #  ;return(row)}, mc.cores = 6)
  
  #output_matrix <- do.call(rbind, output_list)
  #print(dim(output_matrix))
  #rownames(output_matrix) <- rownames(matrix)
  
  return(t(output_matrix))
}





run_PCA <- function(matrix, pcs= 10){
  # compute PCA for a set number of PCs using BiocSingular SVD. 
  pca_obj <- list(BiocSingular::runSVD(t(matrix),k=pcs,center=T, BSPARAM =BiocSingular::bsparam()))[[1]]
  colnames(pca_obj$u) <- paste0('PC',1:pcs)
  
  return(pca_obj)
}




plot_PCA <- function(matrix, dataset_name, labels, label_name, pcs= c(1,2), is_pca_obj= F, is_continuous = F, pcs_to_compute = 10){
  # plots a PCA of the data coloured by labels
  
  if(is_pca_obj){pca_data <- matrix}
  else{pca_data <- run_PCA(matrix, pcs = pcs_to_compute)}
  
  percentage  <- pca_data$d ^ 2 / sum(pca_data$d ^ 2) * 100
  percentage  <- sapply(seq_along(percentage), function(i)
    {round(percentage [i], 1)})
  
  data <- data.frame(pca_data$u)
  colnames(data) <- paste0('PC',c(1:pcs_to_compute))
  data$label <- labels
  
  first_pc <- sym(paste0('PC',pcs[1]))
  second_pc <- sym(paste0('PC',pcs[2]))
  
  xlab <- paste0(first_pc,' (',percentage[pcs[1]],'% variance)')
  ylab <- paste0(second_pc,' (',percentage[pcs[2]],'% variance)')
  
  if(is_continuous){
    p <- ggplot(data, mapping= aes(x = !!first_pc, y = !!second_pc, fill= label)) +
      geom_point(aes(color =label), alpha=0.8 , shape=21) +
      scale_fill_viridis(name = label_name) +
      scale_color_viridis(name = label_name) +
      labs(title = paste0("PCA of ",dataset_name," grouped by ",label_name), x = xlab, y = ylab) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      theme_minimal() +
      theme(legend.background = element_blank(), 
            aspect.ratio = 1, 
            legend.key = element_blank(),
            legend.position = "right",
            panel.background = element_blank(), 
            axis.line = element_line(colour = "grey43",linewidth = 1.1, lineend='round'),
            panel.grid.major = element_line(color = "grey96"))
  }
  
  if(!is_continuous){
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
          axis.line = element_line(colour = "grey43",linewidth = 1.1, lineend='round'),
          panel.grid.major = element_line(color = "grey96"))
  }
  p <- p + theme(legend.position = "right")
  
  if(!is_continuous){
  xdens <- axis_canvas(p, axis = "x")+
    geom_density(data, mapping = aes(x = !!first_pc, fill = labels), color= 'grey45', alpha = 0.50, size = 0.2) +
    theme(legend.position = "none")
  
  ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
    geom_density(data, mapping = aes(x = !!second_pc, fill = labels), color= 'grey45', alpha = 0.50, size = 0.2) +
    theme(legend.position = "none")+
    coord_flip()
  }
  
  if(is_continuous){
    min <- min(labels)
    max <- max(labels)
    weights <- (labels - min)/(max- min)
    xdens <- axis_canvas(p, axis = "x")+
      geom_density(data, mapping = aes(x = !!first_pc), fill= '#1f9a8eff', 
                   color= 'grey45', alpha = 0.45, size = 0.2) +
      geom_density(data, mapping = aes(x = !!first_pc, weight = rank(labels)^2), fill= '#8fd744ff', 
                   color= 'grey45', alpha = 0.45, size = 0.2, bw = "nrd0") +
      scale_fill_manual(name = 'Densities') +
      theme(legend.position = 'right')
    
    ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
      geom_density(data, mapping = aes(x = !!second_pc), fill= '#1f9a8eff', 
                   color= 'grey45', alpha = 0.45, size = 0.2,) +
      geom_density(data, mapping = aes(x = !!second_pc, weight = rank(labels)^2), fill= '#8fd744ff', 
                   color= 'grey45', alpha = 0.45, size = 0.2, bw = "nrd0") +
      theme(legend.position = "none") +
      coord_flip()
  }
  
  p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  
  pList <- ggdraw(p2)

  return(pList)
  
}







plot_3D_PCA <- function(matrix, label_vector = NULL, title = "3D PCA Plot" , is_PCA_obj = F){
  # Plot an interactive 3D PCA using plotly of the first 3 PCs. 
  
  if(is_PCA_obj){data = as.data.frame(matrix$u)}
  else{data = as.data.frame(run_PCA(matrix, pcs=10)$u)}
  
  plot_ly(data = data, 
          x = ~PC1, y = ~PC2, z = ~PC3, color = label_vector, 
          type = 'scatter3d', mode = 'markers',
          marker = list(size = 4, opacity = 0.8)) %>%
    layout(title = title,
           scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))
  
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
  
  if(direction == 'rows'){
    common_samples <- Reduce(intersect, lapply(matrix_list,rownames))
    matrix_list_aligned <- lapply(matrix_list, function(matrix) {matrix[common_samples, , drop = FALSE]})
  }
  
  if(direction == 'columns'){
    common_samples <- Reduce(intersect, lapply(matrix_list,colnames))
    matrix_list_aligned <- lapply(matrix_list, function(matrix) {matrix[, common_samples, drop = FALSE]})
  }
    
  if(!(direction %in% c('columns','rows'))){stop("The direction you entered is invalid. You can only enter 'rows' or 'columns'.")}
  
  return(matrix_list_aligned)
}






remove_all_zero_variance_rows <- function(matrix){
  # taking a matrix, removes all rows with 0 variance
  variances <- apply(matrix, 1, var)
  rows_to_keep <- variances != 0
  return(matrix[rows_to_keep, , drop = FALSE])
  
}






compute_log_transformed_RLE <- function(matrix){
  # taking log transformed matrix, output a matrix of RLE values.
  # try to use rowMedians() instead
  medians <- apply(matrix, 1, median)
  
  RLE_scores <- sweep(matrix, 1, medians, FUN='-')
  
  return(RLE_scores)
  
}






plot_RLE <- function(matrix, batch_info, ylimit=c(-4,4), is_continuous = F, variable_name = NULL){
  # taking a matrix of RLE scores, and a vector of batch information (in the same order as your samples), and y axis limits as a vector of 2 numbers, output an RLE graph ordered by batch
  
  RLE_matrix <- as.data.frame(compute_log_transformed_RLE(matrix))
  
  RLE_long <- reshape2::melt(RLE_matrix)
  
  if(!is_continuous){RLE_long$Batch <- factor(rep(batch_info, each= nrow(RLE_matrix)))}
  if(is_continuous){RLE_long$Batch <- rep(batch_info, each= nrow(RLE_matrix))}
  
  RLE_long <- RLE_long %>%
    arrange(Batch, variable) %>%
    mutate(variable = factor(variable, levels = unique(variable)))
  
  if(!is_continuous){return(
  ggplot(RLE_long, aes(x = variable, y = value, fill = Batch)) +
    geom_boxplot(outlier.shape = NA,alpha=0.9, linewidth=0.5) +  # Set coef to Inf to extend whiskers to max/min + 
    labs(x = "Sample", y = "RLE Score", fill = variable_name, color = 'Batch')+
    coord_cartesian(ylim=ylimit)+
    theme_minimal() +
    theme(panel.border = element_rect(colour = "grey85", fill=NA, size=1.1),
          axis.text.x = element_blank(),
          panel.grid.major = element_line(color = "grey96"),)
  )
  }
  if(is_continuous){return(
    ggplot(RLE_long, aes(x = variable, y = value, fill = Batch)) +
      geom_boxplot(outlier.shape = NA,alpha=0.9, linewidth=0.5) +  # This will use the viridis color scale for the continuous variable
      scale_fill_viridis() +  
      theme_minimal() +
      labs(x = "Sample", y = "RLE Score", fill = variable_name, color = 'Batch')+
      coord_cartesian(ylim=ylimit)+
      labs(x = "Sample", y = "RLE Score", color = "Continuous Variable") +
      theme(panel.border = element_rect(colour = "grey85", fill=NA, size=1.1),
            axis.text.x = element_blank(),
            panel.grid.major = element_line(color = "grey96"),)
  )
  }
}






plot_linear_correlation <- function(matrix_list, dataset_names, variable_vector, variable_name, num_pcs=10, is.pca.obj = T){
  # taking a list of matrices, a vector of the names of the omics, a vector of numerical values (in order) of the values you want to compare, and number of PCs, plot the linear correlation R^2 values

  if(length(matrix_list)>30) {stop('Your omics matrices must be a list. Did you accidentally use c() instead of list()?')}
  
  if(is.list(dataset_names)) {stop('Your dataset names cannot be a list. Has to be a vector. use c() instead of list()!')}
  
  if(!is.pca.obj){pcas <- lapply(matrix_list, function(x) {run_PCA(x, pcs = num_pcs)})}
  if(is.pca.obj){pcas <- matrix_list}
  
  r_squared_results <- sapply(pcas, function(matrix, variable, num) {
    sapply(1:num_pcs, function(y) {
      lm_model <- summary(lm(variable ~ matrix$u[,1:y]))$r.squared})
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

  
  
  


plot_vector_correlation <- function(data_list, dataset_names, variable_vector, variable_name, num_pcs=10, is.pca.obj = T){
  # taking a list of pca objects, a vector of categorical values (in order) of the values you want to compare, and number of PCs, plot the correlation R^2 values
  
  if (length(data_list)>30) {stop('Your data must be a list. Did you accidentally use c() instead of list()?')}
  if (!class(variable_vector) %in% c('factor', 'character')){stop('Your variable is not a factor or character.')}
  if (is.list(dataset_names)) {stop('Your dataset names cannot be a list. Has to be a vector. use c() instead of list()!')}
  
  if(!is.pca.obj){pca_list <- lapply(data_list, function(x) {run_PCA(x,pcs = num_pcs)})}
  if(is.pca.obj){pca_list <- data_list}
  
  dummies <- fastDummies::dummy_cols(variable_vector)[,-1]
  
  cancor_scores <- sapply(pca_list, function(m) {lapply(1:num_pcs, function(y) {cca <- stats::cancor(x= m$u[,1:y, drop=F], 
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





plot_boxplot_categorical <- function(data_vector, category_vector, names, aspect_ratio=1.3){
  # taking a vector of numerical scores and a vector of categorical variables and their names as vector, plot boxplot grouped by the categorical variables.
  data <- as.data.frame(cbind(data_vector, category_vector))
  colnames(data) <- names
  
  ggplot(data, aes(x=data[,2], y=as.numeric(data[,1]), fill=data[,2]))+
    geom_boxplot()+
    theme(axis.text.x = element_blank())+
    labs(x = names[2], y = names[1], fill=names[2] )+
    ggtitle(paste0('Boxplot of ',names[1],' grouped by ',names[2]))+
    theme_minimal() +
    theme(panel.border=element_rect(colour = "grey80", fill=NA, size=0.8),
          aspect.ratio = 1/aspect_ratio,
          axis.line = element_line(colour = "grey75", linewidth = 1.1),
          panel.grid.major = element_line(color = "grey96"),)
  
  
}





plot_violin <- function(data_vector, category_vector, names, aspect_ratio=1.3, overlay_type = 'boxplot'){
  # taking a vector of numerical scores and a vector of categorical variables and their names as vector, plot violin plot & swarm plot grouped by the categorical variables.
  # requires the ggbeeswarm package.
  if(!(overlay_type %in% c('boxplot','swarm'))){stop("The overlay type you entered doesn't exist. You can choose from 'boxplot' and 'swarm'.")}
  
  data <- as.data.frame(cbind(data_vector, category_vector))
  colnames(data) <- names
  
  p <- ggplot(data, aes(x=data[,2], y=as.numeric(data[,1]), fill=data[,2]))+
    geom_violin()+
    labs(x = names[2], y = names[1], fill=names[2] )+
    ggtitle(paste0('Violin plot of ',names[1],' grouped by ',names[2]))+
    theme_minimal() +
    theme(panel.border=element_rect(colour = "grey80", fill=NA, size=0.8),
          aspect.ratio = 1/aspect_ratio,
          axis.line = element_line(colour = "grey75", linewidth = 1.1),
          panel.grid.major = element_line(color = "grey96"),
          axis.text.x = element_text(size = 10,angle = 45,hjust = 1))
  
  if (overlay_type == 'swarm') {
    final_plot <- p + geom_beeswarm()
  } else if (overlay_type == 'boxplot') {
    final_plot <- p + geom_boxplot(width=0.1, alpha=0.5, outlier.shape =NA)
  }
  
  return(final_plot)
  
}




plot_ridgeplot <- function(matrix, label_vector, title = 'Ridgeplot', aspect_ratio = 1/1){
  # taking a matrix and a vector of categorical labels, plot ridgeplot for each category.
  require('ggridges')
  data <- as.data.frame(t(matrix))
  data <- tidyr::pivot_longer(data, cols = everything())
  data$label <- rep(label_vector, each = nrow(matrix))
  
  p <- ggplot(data, aes(x=value, y=label, fill = label))+
    geom_density_ridges2() +
    labs()+
    ggtitle(title)+
    theme_minimal() +
    theme(panel.border=element_rect(colour = "grey80", fill=NA, size=0.8),
          aspect.ratio = aspect_ratio,
          axis.line = element_line(colour = "grey75", linewidth = 1.1),
          panel.grid.major = element_line(color = "grey96"),
          axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
          legend.position = 'none')
  
  return(p)
}







get_high_correlation_features <- function(matrix, variable_vector, threshold, method = 'pearson'){
  # taking a matrix and a vector (numerical), find features (rows) that have absolute correlation score (both positive and negative) greater than your set threshold value with your variable vector. 
  correlations <- lapply(1:nrow(matrix), function(x) {cor(as.vector(as.numeric(matrix[x,])),as.vector(variable_vector), method=method)})
  
  correlations <- as.matrix(unlist(correlations))
  
  rownames(correlations) <- rownames(as.matrix(matrix))
  
  high_corr_features <- correlations[which(abs(correlations) > threshold), ]
  
  return(high_corr_features)
  
}






ftest <- function(matrix, variable, is.log=T, n.cores=8, sort_by = NULL){
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
  
  if (!is.null(sort_by)) {
    if (!sort_by %in% colnames(f.test)){
      stop("The column you want to sort by doesn't exist. You can choose from:'FValue', 'PValue', 'Adj.PValue', 'Mean'.")}
    else {f.test <- f.test[order(f.test[[sort_by]]),]}
  }
  
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
    Var1_F_score = round(unlist(lapply(results, function(x) {x[['F value']][1]})), digits = 4),
    Var1_Adj_PValue = p.adjust(unlist(lapply(results, function(x) {x[['Pr(>F)']][1]})), method = pval_adj_method),
    Var2_F_score = round(unlist(lapply(results, function(x) {x[['F value']][2]})), digits = 4),
    Var2_Adj_PValue = p.adjust(unlist(lapply(results, function(x) {x[['Pr(>F)']][2]})), method = pval_adj_method),
    Interaction_F_score = round(unlist(lapply(results, function(x) {x[['F value']][3]})), digits = 4),
    Interaction_Adj_PValue = p.adjust(unlist(lapply(results, function(x) {x[['Pr(>F)']][3]})), method = pval_adj_method),
    Mean = rowMeans(matrix))
  
  if (!is.null(sort_by)) {
    if (!sort_by %in% colnames(results)){
      stop("The column you want to sort by doesn't exist. You can choose from:'Var1_F_score', 'Var1_Adj_PValue', 
      'Var2_F_score', 'Var2_Adj_PValue', 'Interaction_F_score','Interaction_Adj_PValue', 'Mean'.")}
    else {results <- results[order(results[[sort_by]]),]}
  }
  
  return(results)
}




wilcoxon_test <- function(matrix, variable, is.log=T, n.cores = 8){
  # Taking a matrix of data and a vector of categorical variables, do wilcox test and calculate P values & adjusted Pvalues for each feature.
  # Note that if you put is.log=False, we'll transform it first and then do the ANOVA.
  variable <- as.factor(variable)
  
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
  
  if (!is.null(sort_by)) {
    if (!sort_by %in% colnames(k.test)){
      stop("The column you want to sort by doesn't exist. You can choose from:
       'kruskal_wallis_chisq','PValue','Adj.PValue','Mean'.")
    }}
  
  data <- matrix
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
  
  if (!is.null(sort_by)){k.test <- k.test[order(k.test[[sort_by]]),]}
  
  return(k.test)
}




plot_UMAP <- function(matrix, metadata_vector, title = 'UMAP', aspect_ratio = 1/1, run_umap = T, label_is_continuous = F, 
                      continuous_var_upper_lim = NULL, alpha = 1){
  # taking a matrix and a vector of metadata, plot UMAP of the matrix colored by the groups in the metadata vector
  
  if(run_umap){
    umap_result = umap(t(matrix))
    df = data.frame(umap_result$layout)}
  if(!run_umap){df <- as.data.frame(matrix)}
  colnames(df) = c("UMAP1", "UMAP2")
  
  if(!is.null(continuous_var_upper_lim)){
    if(class(continuous_var_upper_lim) != 'numeric'){stop("You didn't enter a numerical value for the continous variable uppe limit. Please only enter numbers.")}
    else(metadata_vector <- lapply(metadata_vector,
                                   function(x) ifelse(x > continuous_var_upper_lim,
                                                      continuous_var_upper_lim,
                                                      x)) %>% unlist())}
  
  df$metadata = metadata_vector 
  
  if(!label_is_continuous){
    centroids <- aggregate(cbind(UMAP1, UMAP2) ~ metadata, df, mean)
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = metadata)) +
      geom_point(size = 0.07, alpha = alpha) +
      ggtitle(title) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "grey50", linewidth = 0.9),
            panel.border = element_blank(),  #element_rect(colour = "grey90", fill=NA, size=0.7),
            panel.grid.major = element_blank(),  #element_line(color = "grey96"),
            panel.grid.minor = element_blank(),
            aspect.ratio = aspect_ratio,
            legend.position = "none")+
      geom_text(data = centroids, aes(label = metadata), size = 3, color = "black", hjust = 0.5, vjust = 0.5)
  }
  
  if(label_is_continuous){
    ggplot(df, aes(x = UMAP1, y = UMAP2, color = metadata)) +
      geom_point(size = 0.07, alpha = alpha) +
      scale_fill_viridis() +
      scale_color_viridis() +
      ggtitle(title) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "grey50", linewidth = 0.9),
            panel.border = element_rect(colour = "grey90", fill=NA, size=0.7),
            panel.grid.major = element_line(color = "grey96"),
            aspect.ratio = 1/1.05)  
  }
  
  if(!label_is_continuous){
    xdens <- axis_canvas(p, axis = "x")+
      geom_density(df, mapping = aes(x = UMAP1, fill = metadata_vector), color= 'grey55', alpha = 0.50, size = 0.2) +
      theme(legend.position = "none")
    
    ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
      geom_density(df, mapping = aes(x = UMAP2, fill = metadata_vector), color= 'grey55', alpha = 0.50, size = 0.2) +
      theme(legend.position = "none")+
      coord_flip()
  }
  
  p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  pList <- ggdraw(p2)
  
  return(pList)
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





plot_count_map <- function(variable1, variable2, color_threshold = NULL){
  #
  #
  
  map <- data.frame(V1 = variable1, V2 = variable2) %>%
    dplyr::count(V1, V2)
  
  if(!is.null(color_threshold)){
  map$use <- ifelse(map$n > color_threshold, 'yes','no')
  return(
  ggplot(data = map, aes(x = V1, y = V2)) +
    geom_count(aes(color = use)) +
    geom_text(aes(label = n, hjust = 0.5, vjust = 0.5)) +
    theme_bw() +
    theme(axis.line = element_line(colour = 'grey35', size = 1.1),
          axis.title.x = element_text(size = 0),
          axis.title.y = element_text(size = 0),
          axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
          axis.text.y = element_text(size = 12, angle = 45, hjust = 1),
          legend.position = 'none',
          panel.grid.major = element_line(color = "grey95")))}
  if(is.null(color_threshold)){
    return(
    ggplot(data = map, aes(x = V1, y = V2)) +
      geom_count(aes(color = '')) +
      geom_text(aes(label = n, hjust = 0.5, vjust = 0.5)) +
      theme_bw() +
      theme(axis.line = element_line(colour = 'grey35', size = 1.1),
            axis.title.x = element_text(size = 0),
            axis.title.y = element_text(size = 0),
            axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
            axis.text.y = element_text(size = 12, angle = 45, hjust = 1),
            legend.position = 'none',
            panel.grid.major = element_line(color = "grey95")))}
}





plot_GO_enrichment <- function(gene_name_vector, gene_id_type = 'ENSEMBL', ontology = 'BP', pval_cutoff = 0.2){
  # Taking a vector of gene names, plot GO enrichment dotplot. Default ontology is biological process. 
  
  if(!(gene_id_type %in% c('ENSEMBL','ENTREZID'))){stop("The gene id type you entered doesn't exist. You can choose from 'ENTREZID' and 'ENSEMBL'.")}
  if(!(ontology %in% c('BP','CC','MF'))){stop("The ontology you entered doesn't exist. You can choose from 'BP', 'CC' and 'MF'.")}
  
  enrich_obj <- enrichGO(gene = gene_name_vector,
                         OrgDb = 'org.Hs.eg.db', 
                         keyType = gene_id_type, 
                         ont = ontology, 
                         pAdjustMethod = "BH",
                         pvalueCutoff = pval_cutoff)
  
  enrichplot::dotplot(enrich_obj, 
                      label_format = 100)
  
}



compute_silhouette <- function(matrix, label_vector, run_dim_reduction = NULL, num_pcs = 3, result_format = 'average'){
  # taking a matrix of either raw data (rows are featues and columns are samples) 
  # or dimensionally reduced data (PCA, umap or t-sne) (rows are samples and columns are components), calculate silhouette coefficient for each cluster.
  # if you input raw data, dimensionality reduction will be done for you.
  
  if(!(result_format %in% c('average','per_cluster', 'plot'))){stop("Invalid results format. You can choose 
                                                                    between 'average', 'per_cluster' or 'plot'")}
  if(is.null(run_dim_reduction)){
  if(nrow(matrix) != length(label_vector)){stop("Your matrix has unequal number of rows than the length of your label vector.
                                                If you're inputting already dimensionally reduced data, please double 
                                                check that your rows are samples and columns are scores.")}}
  if(!is.null(run_dim_reduction)){
    if(!(run_dim_reduction %in% c('pca','umap','tsne'))){stop("You inputted an invalid dimensionality reduction method. 
                                                              You can choose from 'pca', 'umap' or 'tsne'.")}
    if(run_dim_reduction == 'pca'){matrix <- run_PCA(matrix, pcs = num_pcs)$u}
    if(run_dim_reduction == 'umap'){matrix <- umap(t(matrix))$layout}
    if(run_dim_reduction == 'tsne'){matrix <- Rtsne(t(matrix))$Y}
  }
  
  silhouette_result <- summary(cluster::silhouette(as.numeric(as.factor(label_vector)), 
                                                   stats::dist(matrix)))
  
  if(result_format == 'average'){return(round(silhouette_result$avg.width, digits = 4))}
  if(result_format == 'per_cluster'){results <- data.frame(labels = levels(factor(label_vector)),
                                                           n = as.vector(table(label_vector)),
                                                           silhouette_score = round(silhouette_result$clus.avg.widths, digits=4))
  return(results)
  }
  if(result_format == 'plot'){results <- data.frame(labels = levels(factor(label_vector)),
                                                    n = as.vector(table(label_vector)),
                                                    silhouette_score = round(silhouette_result$clus.avg.widths, digits=4))
  
    return(ggplot(results, aes(x=labels, y= silhouette_score, fill= labels))+
             geom_bar(stat= 'Identity')+
             theme_minimal()+
             theme(axis.line = element_line(colour = "grey88", linewidth=1.1),
                   panel.border = element_rect(colour = "grey90", fill=NA, size=0.7),
                   panel.grid.major = element_line(color = "grey96"),
                   axis.text.x = element_text(angle = 45, hjust = 1)))}
}



compute_ARI <- function(matrix, label_vector, run_dim_reduction = NULL, num_pcs = 3){
  # taking a matrix of either raw data or dimensionally reduced data (PCA, umap or t-sne), calculate adjusted Rand index for each cluster.
  # if you input raw data, dimensionality reduction will be done for you.
  
  if(is.null(run_dim_reduction)){
    if(nrow(matrix) != length(label_vector)){stop("Your matrix has unequal number of rows than the length of your label vector.
                                                If you're inputting already dimensionally reduced data, please double 
                                                check that your rows are samples and columns are scores.")}}
  if(!is.null(run_dim_reduction)){
    if(!(run_dim_reduction %in% c('pca','umap','tsne'))){stop("You inputted an invalid dimensionality reduction method. 
                                                              You can choose from 'pca', 'umap' or 'tsne'.")}
    if(run_dim_reduction == 'pca'){matrix <- run_PCA(matrix, pcs = num_pcs)$u}
    if(run_dim_reduction == 'umap'){matrix <- umap(t(matrix))$layout}
    if(run_dim_reduction == 'tsne'){matrix <- Rtsne(t(matrix))$Y}
  }
  
  computed_clusters <- mclust::Mclust(data = (matrix), 
                                      x= mclust::mclustBIC(matrix))$classification
  
  ARI_score <- mclust::adjustedRandIndex(computed_clusters, label_vector)
  
  return(ARI_score)
  
}



ScheirerRayHare_test <- function(data, variable_1, variable_2, is.log = TRUE, pval_adj_method = 'BH', sort_by = NULL) {
  # Non-parametric version of the two-way ANOVA. Requires the rcompanion package.
  
  if (!requireNamespace("rcompanion", quietly = TRUE)) {
    stop("The rcompanion package is needed for this function to work. Please install it.", call. = FALSE)
  }
  
  if (!is.log) {
    data <- log2(data + 1)
  }
  
  means <- rowMeans(data)
  results <- data.frame(
    Var1_chisq = numeric(),
    Var1_PValue = numeric(),
    Var1_Adj_PValue = numeric(),
    Var2_chisq = numeric(),
    Var2_PValue = numeric(),
    Var2_Adj_PValue = numeric(),
    Interaction_chisq = numeric(),
    Interaction_PValue = numeric(),
    Interaction_Adj_PValue = numeric(),
    Mean = numeric()
  )
  
  if (!is.null(sort_by)) {
    if (!(sort_by %in% colnames(results))) {
      stop("The column you want to sort by doesn't exist. You can choose from: 'Var1_chisq', 'Var1_PValue', 'Var1_Adj_PValue', 'Var2_chisq', 'Var2_PValue', 'Var2_Adj_PValue', 'Interaction_chisq', 'Interaction_PValue', 'Interaction_Adj_PValue', 'Mean'.")
    }}
  
  for (x in 1:nrow(data)) {
    ranks <- rank(as.numeric(data[x, ]))
    data_frame <- data.frame(ranks = ranks, variable_1 = variable_1, variable_2 = variable_2)
    model <- rcompanion::scheirerRayHare(ranks ~ variable_1 * variable_2, data = data_frame)
    
    results[x, "Var1_chisq"] <- round(model$H[1], digits = 4)
    results[x, "Var1_PValue"] <- model$p.value[1]
    results[x, "Var2_chisq"] <- round(model$H[2], digits = 4)
    results[x, "Var2_PValue"] <- model$p.value[2]
    results[x, "Interaction_chisq"] <- round(model$H[3], digits = 4)
    results[x, "Interaction_PValue"] <- model$p.value[3]
    results[x, "Mean"] <- round(means[x], digits = 2)
  }
  
  results$Var1_Adj_PValue <- p.adjust(results$Var1_PValue, method = pval_adj_method)
  results$Var2_Adj_PValue <- p.adjust(results$Var2_PValue, method = pval_adj_method)
  results$Interaction_Adj_PValue <- p.adjust(results$Interaction_PValue, method = pval_adj_method)
  rownames(results) <- rownames(data)
  
  results <- results[order(results[[sort_by]]), ]
  return(results)
}









