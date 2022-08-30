

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("DEqMS")

library(DEqMS)
library(tidyverse)
library(pheatmap)

#==============================================================================
#====================================Reading data==============================
#==============================================================================

df_missing_normalized <- read_tsv("data/CellLineCompartment_TCE_SumTickNormalization.tsv", skip=2)
df_missing_normalized_all <- read.csv("data/CellLineCompartment_TCE_SumTickNormalization.tsv",sep = "\t", header = F)

## annotation column
col_annotation <- data.frame("Cancer type"= t(df_missing_normalized_all[1,2:dim(df_missing_normalized_all)[2]]))

col_annotation$X2 <- as.character(col_annotation$X1)
col_annotation$X2[grep("Breast",col_annotation$X1)] = ("Breast")
col_annotation$X2[grep("Leukemia",col_annotation$X1)] = ("Leukemia")
col_annotation <- data.frame(col_annotation[,-c(1)])
rownames(col_annotation) <- colnames(df_missing_normalized[, -1])
colnames(col_annotation) <- c('simplified_cancer_type')



#==============================================================================
#====================================Exploration==============================
#==============================================================================

df_missing_normalized_matrix <- as.matrix(df_missing_normalized[-1])
rownames(df_missing_normalized_matrix) <-  pull(df_missing_normalized[1])

# check median normalization
# the data is not median normalized
# normalization should be done after log2-transformation
# normalization should be done after filtering missing protein




#===========================check distribution ==================================
# df_missing_normalized_matrix[df_missing_normalized_matrix == 0] <- NA
hist(apply(df_missing_normalized_matrix[, -1], 2, median, na.rm=T))
boxplot(df_missing_normalized_matrix[, 1:20])


#============================missing data handling============================
df_missing_normalized_matrix[is.na(df_missing_normalized_matrix)] <-  0
df_missing_normalized_matrix <- log10(df_missing_normalized_matrix + 1)




#==========================heatmap============================================
# annotation_col is a dataframe with rownames corresponding to the colnames in the feeded matrix
pheatmap(df_missing_normalized_matrix,
         #fontsize_col = 5, 
         #fontsize_row = 4,
         show_rownames = F,
         show_colnames = F,
         main = "Surface", 
         cluster_cols = F,
         cluster_rows = F,
         annotation_col = col_annotation)


#============================Data filtering====================================
# ??
# only leave those spread apart? 
# df_missing_normalized_matrix_frequestnidentified <- 
#   df_missing_normalized_matrix[apply(df_missing_normalized_matrix, 1, function(c) sd(c)>0.1), ]








# test case

# 
# # col_annotation <- data.frame(col_annotation[2], row.names = col_annotation$full_name)
# breast_leukemia_dat <- contrast_dependent_missing_handling(df_missing_normalized_matrix, col_annotation, 'Breast', 'Leukemia', 2)
# 
# 
# #====================================simple t-tests============================
# apply(breast_leukemia_dat, 1, function(x) t.test(x)$p.value)
# 
# fish_exact_test(df_missing_normalized_matrix, 'Breast', col_annotation$simplified_cancer_type)
# 


