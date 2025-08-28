#Hierarchical clustering for meta-transcriptomics

# Libraries
library(pheatmap)
library(dendextend)
library(readr)
library(dplyr)
library(tidyverse)

# Load data
data <- read_csv("M7 alternate/Metatranscriptomics_reduced.csv")

data <- data[,-c(1, 2, 3)]

data <- scale(data)

data <- t(data)

data_names <- rownames(data)
data <- as.matrix(data)

# Calculate distances with euclidean
euclidean_dist_rows <- dist(data, method = "euclidean")

# Perform hierarchical clustering for rows
complete_clusters_euclidean_rows <- hclust(euclidean_dist_rows, method = "complete")

# Calculate distances for columns
euclidean_dist_cols <- dist(t(data), method = "euclidean")

# Perform hierarchical clustering for columns
complete_clusters_euclidean_cols <- hclust(euclidean_dist_cols, method = "complete")


#Euclidean distance heatmap
pheatmap(as.matrix(data),
         cluster_rows = complete_clusters_euclidean_rows,
         cluster_cols = complete_clusters_euclidean_cols,
         main = "Meta-transcriptomics Euclidean Distance Heatmap",
         fontsize = 12)