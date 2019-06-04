library(Matrix)
library(Rtsne)


train <- tril(as.matrix(algnscore_matrix))
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
fit_cluster_kmeans=kmeans(scale(d_tsne_1), variable_to_pass) #this var will be imported from py wrapper
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)


# plot the t-SNE map.
pdf(file = "tane_clust_trees_Rplots.pdf")

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
trainlist = list(train[,0])
input_and_clust <- data.frame(trainlist, d_tsne_1_original$cl_kmeans)


algnscore_matrix <- read.csv("./algnscore_matrix.csv")


library("stats")
library("ape")
#dynamic subsetting + trees from clust subsets subsetting algnscore
#TODO make this less ugly
for(i in unique(input_and_clust$d_tsne_1_original.cl_kmeans)) {
  i_clusts <- input_and_clust[input_and_clust$d_tsne_1_original.cl_kmeans==i, 'names']
  sub_mat <- algnscore_matrix[i_clusts,i_clusts]
  for_hclust <- dist(sub_mat, method = "canberra", diag = FALSE, upper = FALSE, p = 2) #method can be changed. From our tests, the one chosen here produced the most resonable trees
  hc <- hclust(for_hclust, method = "average", members = NULL)
  plot(hc, hang = -0.5, cex = 0.3, edge.width=0.4, main = "Cluster dendorgram", cex.main = 1.0, xlab = "ID", ylab = "Height") ##reg tree
  myphylo <- as.phylo.hclust(hc)
  plot(myphylo, type = "unrooted", cex = 0.4, label.offset = 0.5, edge.width=0.4, lab4ut="axial", no.margin = FALSE, main = "Unrooted cluster dendorgram", cex.main = 1.0) ##unrooted
  #plot(myphylo, type = "fan", cex = 0.3, edge.width=0.4) ##Fan/circular
  }

dev.off()