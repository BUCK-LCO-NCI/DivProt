#!/usr/bin/env r

#produce matrices

library("tidyr")
library("dplyr")
library("reshape")

##subject file
struct_align_score_table <- read.table("struct_score_table.csv", sep = ",", header = FALSE)

colnames(struct_align_score_table) <- c("qseqid", "sseqid", "align_score")

struct_align_score_table$align_score = as.numeric(gsub("Score=", "", struct_align_score_table$align_score)) #gets rid of Score=

####align score
shaped_algn_score <- reshape(struct_align_score_table,idvar = "qseqid", timevar = "sseqid", direction = "wide")
names(shaped_algn_score) <- substring(names(shaped_algn_score),13,100) 


#this creates rownames (seq ids are in a  column/row now)
tmpb <- shaped_algn_score[,-1]
rownames(tmpb) <- shaped_algn_score[,1]
shaped_algn <- tmpb

algnscore_matrix <- shaped_algn

#Below replaces NA with 0, was necessary for amino acid script as perfect alignments gave us NA out. 
#Here perfect alignment across the diagonal is a real , and variable,value. Not a problem for the tees and networks as we can ignore diagonal 
#it is visible on the heatmap though. I don't think it's a problem, as the diagonal always = best score, which one can see. I can edit though to produce a prefct diagonal if we decide we want that though
#algnscore_matrix[is.na(algnscore_matrix)] <- 0 #matrix with 0s to replace Nas

write.csv(algnscore_matrix, "algnscore_matrix.csv") 

#dendrogram / tree   
library("stats")
pdf(file = "Rplots.pdf")
for_hclust <- dist(algnscore_matrix, method = "canberra", diag = FALSE, upper = FALSE, p = 2)
hc <- hclust(for_hclust, method = "average", members = NULL)
plot(hc, hang = -0.5, cex = 0.6) 

#more trees
library("ape")
##unrooted
myphylo <- as.phylo.hclust(hc)
plot(as.phylo.hclust(hc), type = "unrooted", cex = 0.6,
     no.margin = TRUE)

##Fan/circular
plot(as.phylo.hclust(hc), type = "fan")

#heatmap
library("pheatmap")
matrix_confirm <- as.matrix(algnscore_matrix)
mode(matrix_confirm) <- "numeric"
plot(matrix_confirm, cex = 0.6) 


#actual igrpah, two for different layout options 
ig <- graph.adjacency(matrix_confirm, mode="undirected", weighted=TRUE, diag = TRUE)

community_clustering <- fastgreedy.community(ig)
cluster_colors <- rainbow(max(membership(community_clustering)))

l <- layout <- layout.reingold.tilford(ig, circular=T)
ll <- layout.fruchterman.reingold(ig, niter=10000)

plot(ig, layout=l, 
     edge.arrow.size=0.5, 
     vertex.label.cex=0.5, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.shape="circle", 
     vertex.size=3,
     vertex.color=colors[membership(fc)], 
     vertex.label.color="black", 
     edge.width=0.4)

plot(ig, layout=ll, 
     edge.arrow.size=0.5, 
     vertex.label.cex=0.5, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.shape="circle", 
     vertex.size=3,
     vertex.color=colors[membership(fc)], 
     vertex.label.color="black", 
     edge.width=0.4)

dev.off()

