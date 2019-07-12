#THIS IS A DRAFT, THRESHOLDING METHOD IS ON IT'S WAY

args = commandArgs(trailingOnly=TRUE)
matrix confirm = args[1]
var_score_to_pass = args[2]

 #csv var passed from python wrapper script

#actual igraph, two for different layout options
library("igraph")
ig <- graph.adjacency(matrix_confirm, mode="undirected", weighted=TRUE, diag = FALSE, add.colnames = NULL)

community_clustering <- fastgreedy.community(ig)
cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)

l <- layout <- layout.reingold.tilford(ig, circular=T)
ll <- layout.fruchterman.reingold(ig, niter=10000)
#layout_nicel  with the same params = same graph. This seems to be the best network layout for this data. If you want to play around with layouts, we highly recommend a high number of iterations like we've used, as it makes a significant difference in clustering. Higher niter = MUCH more accurate
#circlar plot provided b/c sometimes the characteristics of data can make a typical network hard to read. It's usually not going to be espacially valuable.

pdf(file = "networks_Rplots3.pdf")

plot(ig, layout=l,
     edge.arrow.size=0.5,
     vertex.label.cex=0.5,
     vertex.label.family="Helvetica",
     vertex.label.font=1.5,
     vertex.shape="circle",
     vertex.size=3,
     vertex.color=cluster_colours[membership(community_clustering)],
     vertex.label.color="black",
     edge.width=0.3)

plot(ig, layout=ll,
     edge.arrow.size=0.5,
     vertex.label.cex=0.5,
     vertex.label.family="Helvetica",
     vertex.label.font=1.5,
     vertex.shape="circle",
     vertex.size=3,
     vertex.color=cluster_colours[membership(community_clustering)],
     vertex.label.color="black",
     edge.width=0.3)

dev.off()
