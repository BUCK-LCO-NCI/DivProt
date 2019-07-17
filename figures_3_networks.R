#THIS IS A DRAFT, THRESHOLDING METHOD IS ON IT'S WAY

args = commandArgs(trailingOnly=TRUE)
matrix_confirm_temp = args[1]
var_score_to_pass = args[2]

#csv var passed from python wrapper script

#formatting
matrix_confirm <- as.matrix(matrix_confirm_temp)
mode(matrix_confirm) <- "numeric"   

#actual igraph, two for different layout options
library("igraph")
ig <- graph.adjacency(matrix_confirm, mode="undirected", weighted=TRUE, diag = FALSE, add.colnames = NULL)

community_clustering <- fastgreedy.community(ig)
cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)

l <- layout <- layout.reingold.tilford(ig, circular=T)
ll <- layout_with_graphopt(ig, niter=5000, charge = 0.01)
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
     main="[title]",
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
     main="[title]",
     edge.width=0.3)


#For edgeweight introduction
te <- delete_edges(ig, E(ig)[weight<1976])
community_clustering <- fastgreedy.community(te)
cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)
te2 <- layout.fruchterman.reingold(te, niter=5000)

plot(te,
     layout=te2,
     edge.arrow.size=0.5,
     vertex.label.cex=0.5,
     vertex.label.family="Helvetica",
     vertex.label.font=1.5,
     vertex.shape="circle",
     vertex.size=3,
     vertex.color=cluster_colours[membership(community_clustering)],
     vertex.label.color="black",
     edge.width=0.3,
     main="[title]")

#For delete zero edge nodes of edgeweight-adjusted
#gives more spaced + readable clusters but totherwise same as above
dv <- delete.vertices(te, V(te)[degree(te)==0])

community_clustering <- multilevel.community(dv)
cluster_colors <- rainbow(max(membership(community_clustering)), alpha = 0.5)

lld <- layout_with_graphopt(dv, niter=5000,charge = 0.01) #this is BY FAR the best layout

plot(dv,
     layout=lld, 
     edge.arrow.size=0.5, 
     vertex.label.cex=0.5, 
     vertex.label.family="Helvetica",
     vertex.label.font=1.5,
     vertex.shape="circle", 
     vertex.size=3,
     vertex.color=cluster_colors[membership(community_clustering)], 
     vertex.label.color="black", 
     edge.width=0.4, 
     main="[title]")



dev.off()
