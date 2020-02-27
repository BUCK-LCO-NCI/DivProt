#heatmap
library("pheatmap")

args = commandArgs(trailingOnly=TRUE)
align_file = read.csv(args[1], header = TRUE, row.names = 1, check.names = FALSE)

pdf(file = "heatmap_Rplots1.pdf")

#algnscore_matrix <- read.csv("./algnscore_matrix.csv") #remeber this matrix needs to be fully done w/ R script before this vis
algnscore_matrix <-  align_file #var passed from the python wrapper script
matrix_confirm <- as.matrix(algnscore_matrix)
mode(matrix_confirm) <- "numeric"
pheatmap(matrix_confirm, cex = 0.5, main = "Heatmap")

dev.off()
