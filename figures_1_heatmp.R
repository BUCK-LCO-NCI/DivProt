#heatmap
library("pheatmap")

pdf(file = "heatmap_Rplots1.pdf")

algnscore_matrix <- read.csv("./algnscore_matrix.csv") #remeber this matrix needs to be fully done w/ R script before this vis

matrix_confirm <- as.matrix(algnscore_matrix)
mode(matrix_confirm) <- "numeric"
pheatmap(shaped_algn, cex = 0.7, main = "Heatmap")

dev.off()
