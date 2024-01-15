BiocManager::install("TCGAbiolinks")

BiocManager::install("MOFA2")
BiocManager::install('sesameData')
BiocManager::install('sesame')
BiocManager::install('GenomicRanges',force = TRUE)
library(GenomicRanges)
library(sesameData)
library(sesame)
library(TCGAbiolinks)
library(MOFA2)



# Download the Data (transcriptome, methylation, proteome, micro-RNA)

query_transcriptome <- GDCquery(project = "TCGA-PRAD", 
                                data.category = "Transcriptome Profiling",
                                data.type = "Gene Expression Quantification")


GDCdownload(query_transcriptome)

data_transcriptome <- GDCprepare(query_transcriptome)

query_methylation <- GDCquery(project = "TCGA-PRAD", 
                              data.category = "DNA Methylation",
                              data.type = "Methylation Beta Value")

GDCdownload(query_methylation)

data_methylation <- GDCprepare(query_methylation)

query_proteome <- GDCquery(project = "TCGA-PRAD", 
                              data.category = "Proteome Profiling",
                              data.type = "Protein Expression Quantification")

GDCdownload(query_proteome)

data_proteome <- GDCprepare(query_proteome)



query_mir <- GDCquery(project = "TCGA-PRAD", 
                                 data.category = "Transcriptome Profiling",
                                 data.type = "miRNA Expression Quantification")

GDCdownload(query_mir)

data_mir <- GDCprepare(query_mir)


# Data prep for mofa

# Make the data into matrices with rows being features and columns being patients

proteome_matrix <- data_proteome[,5:ncol(data_proteome)]

selected_columns <- seq(3, ncol(data_mir), by=3)

mir <- data_mir[,1]

mir_matrix <- cbind(mir,data_mir[,selected_columns])

transcriptome_matrix <- assay(data_transcriptome, 'fpkm_uq_unstrand')

methylation_matrix <- assay(data_methylation)

methylation_matrix <- methylation_matrix[!apply(is.na(methylation_matrix), 1, all), ]


# find patients that are in multiple omic data types 


# Extract patient IDs from column names
get_patient_id <- function(barcode) {
  substr(barcode, 1, 12)
}

get_mir_patient_id <- function(barcode) {
  substr(barcode, 32, 43)
}

proteome_patients <- sapply(colnames(proteome_matrix), get_patient_id)
mir_patients <- sapply(colnames(mir_matrix), get_mir_patient_id)
transcriptome_patients <- sapply(colnames(transcriptome_matrix),get_patient_id)
methylation_patients <- sapply(colnames(methylation_matrix),get_patient_id)



# Identify common patients
common_patients <- Reduce(intersect, list(proteome_patients,mir_patients,transcriptome_patients,methylation_patients))

# Filter matrices to include only common patients
filtered_proteome_matrix <- proteome_matrix[, sapply(colnames(proteome_matrix), get_patient_id) %in% common_patients]
filtered_mir_matrix <- mir_matrix[, sapply(colnames(mir_matrix), get_mir_patient_id) %in% common_patients]
filtered_transcriptome_matrix <- transcriptome_matrix[, sapply(colnames(transcriptome_matrix), get_patient_id) %in% common_patients]
filtered_methylation_matrix <- methylation_matrix[, sapply(colnames(methylation_matrix), get_patient_id) %in% common_patients]

dim(filtered_transcriptome_matrix)
dim(filtered_methylation_matrix)
dim(filtered_mir_matrix)
dim(filtered_proteome_matrix)



patient_ids_transcriptome <- sapply(colnames(transcriptome_matrix), get_patient_id)

# Tabulate frequency of each patient ID
patient_freq <- table(patient_ids_transcriptome)

# Identify patient IDs with frequency greater than 1
duplicated_patients <- sort(names(patient_freq[patient_freq > 1]))

print(duplicated_patients)


duplicated_barcodes <- colnames(filtered_transcriptome_matrix)[sapply(colnames(filtered_transcriptome_matrix), get_patient_id) %in% duplicated_patients]

# Sort and print the duplicated barcodes
sorted_duplicated_barcodes <- sort(duplicated_barcodes)
print(sorted_duplicated_barcodes)


get_sample_type <- function(barcode) {
  substr(barcode, 14, 15)
}
get_mir_sample_type <- function(barcode) {
  substr(barcode, 45, 46)
}
# Filter out solid tissue normal samples
filtered_transcriptome_matrix <- filtered_transcriptome_matrix[, get_sample_type(colnames(filtered_transcriptome_matrix)) != "11"]
filtered_methylation_matrix <- filtered_methylation_matrix[, get_sample_type(colnames(filtered_methylation_matrix)) != "11"]
filtered_mir_matrix <- filtered_mir_matrix[, get_mir_sample_type(colnames(filtered_mir_matrix)) != "11"]
filtered_proteome_matrix <- filtered_proteome_matrix[, get_sample_type(colnames(filtered_proteome_matrix)) != "11"]


patient_ids_transcriptome <- sapply(colnames(filtered_transcriptome_matrix), get_patient_id)
patient_freq <- table(patient_ids_transcriptome)

# Identify patient IDs with frequency greater than 1
duplicated_patients <- sort(names(patient_freq[patient_freq > 1]))

print(duplicated_patients)


duplicated_barcodes <- colnames(filtered_transcriptome_matrix)[sapply(colnames(filtered_transcriptome_matrix), get_patient_id) %in% duplicated_patients]


print(duplicated_barcodes)


samples_to_remove <- c("TCGA-V1-A9O5-06A", "TCGA-HC-7740-01B")


extract_first_16 <- function(barcode) {
  substr(barcode, 1, 16)
}
extract_mir <- function(barcode) {
  substr(barcode, 32, 47)
}

# Filter out columns from the transcriptome matrix that match the specified prefixes
filtered_transcriptome_matrix <- filtered_transcriptome_matrix[, !extract_first_16(colnames(filtered_transcriptome_matrix)) %in% samples_to_remove]

# Repeat for the methylation matrix
filtered_methylation_matrix <- filtered_methylation_matrix[, !extract_first_16(colnames(filtered_methylation_matrix)) %in% samples_to_remove]

# Repeat for the miRNA matrix
filtered_mir_matrix <- filtered_mir_matrix[, !extract_mir(colnames(filtered_mir_matrix)) %in% samples_to_remove]

# Repeat for the proteome matrix
filtered_proteome_matrix <- filtered_proteome_matrix[, !extract_first_16(colnames(filtered_proteome_matrix)) %in% samples_to_remove]


# check if all our omic matrices have the same samples
dim(filtered_transcriptome_matrix)
dim(filtered_methylation_matrix)
dim(filtered_mir_matrix)
dim(filtered_proteome_matrix)



# Sort matrices by patient ID

ordered_patient_ids <- sapply(colnames(filtered_transcriptome_matrix), get_patient_id)

# Sort the transcriptome matrix
filtered_transcriptome_matrix <- filtered_transcriptome_matrix[, match(ordered_patient_ids, sapply(colnames(filtered_transcriptome_matrix), get_patient_id))]

# Sort the methylation matrix
filtered_methylation_matrix <- filtered_methylation_matrix[, match(ordered_patient_ids, sapply(colnames(filtered_methylation_matrix), get_patient_id))]

# Sort the miRNA matrix
filtered_mir_matrix <- filtered_mir_matrix[, match(ordered_patient_ids, sapply(colnames(filtered_mir_matrix), get_mir_patient_id))]

# Sort the proteome matrix
filtered_proteome_matrix <- filtered_proteome_matrix[, match(ordered_patient_ids, sapply(colnames(filtered_proteome_matrix), get_patient_id))]


# Now for the actual MOFA

is.matrix(filtered_transcriptome_matrix)
is.matrix(filtered_mir_matrix)
is.matrix(filtered_proteome_matrix)
is.matrix(filtered_methylation_matrix)

filtered_mir_matrix <- as.matrix(filtered_mir_matrix)
filtered_proteome_matrix <- as.matrix(filtered_proteome_matrix)



colnames(filtered_transcriptome_matrix) <- substr(colnames(filtered_transcriptome_matrix), 1, 12)
colnames(filtered_methylation_matrix) <- substr(colnames(filtered_methylation_matrix), 1, 12)
colnames(filtered_mir_matrix) <- substr(colnames(filtered_mir_matrix), 32, 43)
colnames(filtered_proteome_matrix) <- substr(colnames(filtered_proteome_matrix), 1, 12)

original_mir_names <- mir_matrix[,1]
original_protein_names <- proteome_matrix[,1]

rownames(filtered_mir_matrix) <- original_mir_names
rownames(filtered_proteome_matrix) <- original_protein_names


data_list <- list(
  filtered_transcriptome_matrix,
  filtered_methylation_matrix,
  filtered_mir_matrix,
  filtered_proteome_matrix
)


mofa_object <- create_mofa(data_list)
mofa_object <- prepare_mofa(mofa_object)
trained_mofa <- run_mofa(mofa_object, use_basilisk = TRUE)


plot_data_overview(mofa_object)

plot_variance_explained(trained_mofa, x="view", y="factor")

plot_factor(trained_mofa, factor= 1:15)

plot_data_heatmap(trained_mofa, view = 'view_2', factor = 1)

plot_data_heatmap(trained_mofa, view = 'view_2', factor = 2)

plot_data_heatmap(trained_mofa, view = 'view_1', factor = 1)

plot_data_heatmap(trained_mofa, view = 'view_2', factor = 3)

plot_data_scatter(trained_mofa,
                  view = "view_2",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,  )

plot_data_scatter(trained_mofa,
                  view = "view_4",         # view of interest
                  factor = 4,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,  )







# examine the data types more to see batch effects
# get batch information as metadata and do PCA to see if theres batch effect
# try complex heat map??
# examine methylation and RNA (get biology subtypes)
# try to do feature engineering (remove useless features) and identify outliers
# PURITY as metadata to see what's happening

# outliers: add up all the beta values for each patient (global methylation of the patient)
# and see if any factors picks up the global methylation
# library size: add up all the transcripts
# use as metadata to see what is happening







# Examining batch effects by PCA

library(ggplot2)

batch_info <- data_transcriptome@colData$paper_Batch
transcriptome_matrix1 <- assay(data_transcriptome, "fpkm_uq_unstrand")

zero_rows <- apply(transcriptome_matrix1, 1, function(row) all(row == 0))

filtered_transcriptome_matrix1 <- transcriptome_matrix1[!zero_rows, ]

percentage_zeros <- apply(filtered_transcriptome_matrix1, 1, function(row) {
  sum(row == 0) / length(row) * 100
})

# Identify rows with more than 90% zeros
rows_to_remove <- percentage_zeros > 80

print(rows_to_remove)

filtered_transcriptome_matrix2 <- filtered_transcriptome_matrix1[!rows_to_remove, ]

selected_batches <- c(91, 357, 184, 370, 389)
selected_samples <- batch_info %in% selected_batches

# Subset the SummarizedExperiment object
filtered_transcriptome_matrix3 <- filtered_transcriptome_matrix2[, selected_samples]

pca_result <- prcomp(t(filtered_transcriptome_matrix3), center = TRUE, scale. = TRUE)
pca_result <- prcomp(t(filtered_transcriptome_matrix3), center = TRUE, scale. = FALSE)

pca_scores <- as.data.frame(pca_result$x)

selected_batch_info_transcriptome <- batch_info[selected_samples]

pca_scores$Batch <- as.factor(selected_batch_info_transcriptome)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  labs(title = "PCA Plot of transcriptome", 
       x = paste("PC1: ", round(pca_result$sdev[1]^2/sum(pca_result$sdev^2) * 100, 2), "% variance"),
       y = paste("PC2: ", round(pca_result$sdev[2]^2/sum(pca_result$sdev^2) * 100, 2), "% variance")) +
xlim(c(-60,60))



# for methylation
batch_info_methylation <- data_methylation@colData$paper_Batch

selected_batches <- c(91, 370, 389)
selected_samples <- batch_info_methylation %in% selected_batches

# Subset the SummarizedExperiment object
filtered_methylation_matrix1 <- methylation_matrix[, selected_samples]

# Compute the mean for each probe excluding NAs
na_fraction <- apply(filtered_methylation_matrix1, 1, function(x) sum(is.na(x))/length(x))

# Identify rows where the fraction of NAs is greater than 90%
rows_to_remove <- which(na_fraction > 0.5)

print(rows_to_remove)

# Remove these rows
filtered_methylation_matrix2 <- filtered_methylation_matrix1[-rows_to_remove, ]


# Compute the mean for each probe excluding NAs
probe_means <- apply(filtered_methylation_matrix2, 1, function(x) mean(x, na.rm = TRUE))

# Replace NAs with the computed means for each probe
for(i in 1:nrow(filtered_methylation_matrix2)) {
  filtered_methylation_matrix2[i, is.na(filtered_methylation_matrix2[i,])] <- probe_means[i]
}


pca_methylation <- prcomp(t(filtered_methylation_matrix2), center = TRUE, scale. = TRUE)

pca_methylation <- prcomp(t(filtered_methylation_matrix2), center = TRUE, scale. = FALSE)

pca_methylation_scores <- as.data.frame(pca_methylation$x)

selected_batch_info_methylation <- batch_info_methylation[selected_samples]

pca_methylation_scores$Batch <-as.factor(selected_batch_info_methylation)

ggplot(pca_methylation_scores, aes(x = PC1, y = PC2, color= Batch)) +
  geom_point() +
  labs(title = "PCA Plot of methylation", 
       x = paste("PC1: ", round(pca_result$sdev[1]^2/sum(pca_result$sdev^2) * 100, 2), "% variance"),
       y = paste("PC2: ", round(pca_result$sdev[2]^2/sum(pca_result$sdev^2) * 100, 2), "% variance")
)



# for miRNA:
na_fraction <- apply(mir_matrix, 1, function(x) sum(is.na(x))/length(x))

# Identify rows where the fraction of NAs is greater than 90%
rows_to_remove <- which(na_fraction > 0.5)

print(rows_to_remove)

zeros_fraction <- apply(mir_matrix, 1, function(x) mean(x==0))

print(zeros_fraction)

zero_rows2 <- which(zeros_fraction < 0.1)

print(zero_rows2)

mir_matrix1 <- mir_matrix[-zero_rows2,]

print(table(data_proteome$set_id))








