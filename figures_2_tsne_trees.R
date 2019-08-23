library(Matrix)
library(Rtsne)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
align_file = read.csv(args[1], header = TRUE, row.names = 1)
var_score_to_pass = args[2]

train <- tril(as.matrix(align_file))
train[is.infinite(train)] <- 0

#this tsne code is from https://www.r-bloggers.com/playing-with-dimensions-from-clustering-pca-t-sne-to-carl-sagan/   !!!!!!
tsne_model_1 = Rtsne(as.matrix(train), check_duplicates=FALSE, pca=TRUE, perplexity=20, theta=0.5, dims=2) #perplexity may need to change with more complicated data

d_tsne_1 = as.data.frame(tsne_model_1$Y)

ggplot(d_tsne_1, aes(x=V1, y=V2)) +
  geom_point(size=0.25) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")


## keeping original data
d_tsne_1_original=d_tsne_1

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1), var_score_to_pass) #this var will be imported from py wrapper
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)


# plot the t-SNE map.
pdf(file = "tsne_clust_trees_Rplots.pdf")

plot_cluster=function(data, var_cluster, palette)
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=1.0) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("T-SNE clusters on k-means") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.box = "horizontal") +
    scale_colour_brewer(palette = palette)
}

plot_cluster(d_tsne_1_original, "cl_kmeans", "Accent")


#now let's extract the clusters and build those trees
trainlist = as.matrix(train[,0])

dtsne1clkmeans = as.matrix(d_tsne_1_original$cl_kmeans)
input_and_clust <- data.frame(trainlist, d_tsne_1_original$cl_kmeans) #has ids and inde and vector length of [1]
#this is, in fact, necessary
input_and_clust2 <- data.frame(names = row.names(input_and_clust), input_and_clust)
rownames(input_and_clust2) <- NULL


#algnscore_matrix <- read.csv("./algnscore_matrix.csv") nope this needs to be dyanmic b/c someone could use main 3-mat csv, weighted csv, or even one of the 3 method mat csvs if they want
algnscore_matrix <- align_file #variable from python wrapper script, the input csv

library("stats")
library("ape")
#dynamic subsetting + trees from clust subsets subsetting algnscore
#TODO make this less ugly, also include which cluster in title name of trees
for(i in unique(input_and_clust2$d_tsne_1_original.cl_kmeans)) {
  i_clusts <- input_and_clust2[input_and_clust2$d_tsne_1_original.cl_kmeans==i, 'names']
  sub_mat <- algnscore_matrix[i_clusts,i_clusts]
  for_hclust <- dist(sub_mat, method = "canberra", diag = FALSE, upper = FALSE, p = 2) #method can be changed. From our tests, the one chosen here produced the most resonable trees
  hc <- hclust(for_hclust, method = "average", members = NULL)
  plot(hc, hang = -0.5, cex = 0.3, edge.width=0.4, main = "Cluster dendorgram", cex.main = 1.0, xlab = "ID", ylab = "Height") ##reg tree
  myphylo <- as.phylo.hclust(hc)
  plot(myphylo, type = "unrooted", cex = 0.4, label.offset = 0.5, edge.width=0.4, lab4ut="axial", no.margin = FALSE, main = "Unrooted cluster dendorgram", cex.main = 1.0) ##unrooted
  #plot(myphylo, type = "fan", cex = 0.3, edge.width=0.4) ##Fan/circular
  }

dev.off()

#hey so don't worry about the 'There were __ warnings...' message you're going to get out to termainal, it's all good
