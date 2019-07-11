#!/usr/bin/env r

#produce matrices

library("tidyr")
library("dplyr")
library("reshape")

##subject file
PSI_table <- read.table("PSIBLAST_results.txt", sep = "\t", header = FALSE)

colnames(PSI_table) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

PSI_evalue <- PSI_table[,c(1:2,11)]
PSI_bitscore <- PSI_table[,c(1:2,12)]
PSI_pident <- PSI_table[,c(1:3)]

#####evalue 
shaped_eval <- reshape(PSI_evalue,idvar = "qseqid", timevar = "sseqid", direction = "wide")
names(shaped_eval) <- substring(names(shaped_eval),8,100) #the 8 and 100 are the characters I want to keep, need the names of columns to be the same as rows

tmp <- shaped_eval[,-1]
rownames(tmp) <- shaped_eval[,1]
shaped_eval <- tmp

evalue_matrix <- shaped_eval[,rownames(shaped_eval)]

write.csv(evalue_matrix, "evalue_matrix.csv")

#####%identity

shaped_pident <- reshape(PSI_pident,idvar = "qseqid", timevar = "sseqid", direction = "wide")
names(shaped_pident) <- substring(names(shaped_pident),8,100)

tmpi <- shaped_pident[,-1]
rownames(tmpi) <- shaped_pident[,1]
shaped_pident <- tmpi


pident_matrix <- shaped_pident[,rownames(shaped_pident)]

write.csv(pident_matrix, "pident_matrix.csv")

#####bitscore
shaped_bit <- reshape(PSI_bitscore,idvar = "qseqid", timevar = "sseqid", direction = "wide")
names(shaped_bit) <- substring(names(shaped_bit),10,100)

tmpb <- shaped_bit[,-1]
rownames(tmpb) <- shaped_bit[,1]
shaped_bit <- tmpb


bitscore_matrix <- shaped_bit[,rownames(shaped_bit)]

bitscore_matrix[is.na(bitscore_matrix)] <- 0 #matrix with 0s to replace Nas
#bitscore_matrix[bitscore_matrix == 0] <- 1 #for a matrix that undergo the transformation, we need to replace all 0s with 1s because of weighting. We can't weight with 0s. 1s are just fine, and don't influence the actual scores (they're always much much higher)

A1 = bitscore_matrix/apply(bitscore_matrix,1,max) ###normalise all values in column to highest value (different bitscore maximums because of prot length, but this is not good - max 100% identity needs to have the same value)
A2 = t((bitscore_matrix)/apply(bitscore_matrix,2,max))
A1 <- as.matrix(A1)
A2 <- as.matrix(A2)
result = ifelse(A1>A2,A1,A2)
bitscore_matrix_2 <- as.data.frame(result) 
##doing the transformation gives you different values on different sides of the matrix from dividing by different maximums
#we have three options for dealing with this: choosing the smallest, choosing the largest, or avging them
#after discussing we're taking the largest, since that value is the one from the smaller max bitscore division, which 
#equals the greater the "confidence" in the bitscore from a longer protein (even though all are 100% similar,
#longer prots have more weight = smaller = sortabetter bitscore)

bitscore_matrix_2 <- bitscore_matrix_2 * 10000 #more info on this value in the paper, but basically we need to scale the matrix up post-normalisation to actually have any effect in weighting, what we're assuming you will use it for

write.csv(bitscore_matrix, "bitscore_matrix_original.csv")
write.csv(bitscore_matrix_2, "bitscore_matrix_diag_normd_scaled.csv") 

#TO-DO: CHANGE THIS(?) (choose one) WHEN WE DECIDE WHICH IS BEST!!

#dendrogram / tree
library("stats")
pdf(file = "Rplots.pdf")
for_hclust <- dist(bitscore_matrix, method = "canberra", diag = FALSE, upper = FALSE, p = 2)
hc <- hclust(for_hclust, method = "average", members = NULL)
plot(hc, hang = -0.5, cex = 0.5) 

#more trees
library("ape")
##unrooted
myphylo <- as.phylo.hclust(hc)
plot(as.phylo.hclust(hc), type = "unrooted", cex = 0.5,
     no.margin = TRUE)

#heatmap
library("pheatmap")
ph <- pheatmap(bitscore_matrix, width = 10, height = 10)
plot(ph, cex = 0.5)

#igraph
matrix_confirm <- as.matrix(bitscore_matrix)
mode(matrix_confirm) <- "numeric"


ig <- graph.adjacency(matrix_confirm, mode="undirected", weighted=TRUE, diag = FALSE)

community_clustering <- fastgreedy.community(ig)
cluster_colours <- rainbow(max(membership(community_clustering)))

l <- layout <- layout.reingold.tilford(ig, circular=T)
ll <- layout.fruchterman.reingold(ig, niter=10000)


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

