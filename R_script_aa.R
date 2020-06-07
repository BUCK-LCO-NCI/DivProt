#!/usr/bin/env r

#produce matrices

library("tidyr")
library("dplyr")
library("reshape")

#create sub-directories for edge-adjusted nets
dir.create("./networks")
dir.create("./networks/adjusted_algn_score_csvs")
dir.create("./networks/community_csvs")

##subject file
PSI_table <- read.csv("PSIBLAST_results.txt", sep = "\t", header = FALSE)

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

write.csv(evalue_matrix, "aa_align_evalue_matrix_original.csv")

#####%identity

shaped_pident <- reshape(PSI_pident,idvar = "qseqid", timevar = "sseqid", direction = "wide")
names(shaped_pident) <- substring(names(shaped_pident),8,100)

tmpi <- shaped_pident[,-1]
rownames(tmpi) <- shaped_pident[,1]
shaped_pident <- tmpi


pident_matrix <- shaped_pident[,rownames(shaped_pident)]

write.csv(pident_matrix, "aa_align_pident_matrix_original.csv")

#####bitscore
shaped_bit <- reshape(PSI_bitscore,idvar = "qseqid", timevar = "sseqid", direction = "wide")
names(shaped_bit) <- substring(names(shaped_bit),10,100)

tmpb <- shaped_bit[,-1]
rownames(tmpb) <- shaped_bit[,1]
shaped_bit <- tmpb


bitscore_matrix <- shaped_bit[,rownames(shaped_bit)]

bitscore_matrix[is.na(bitscore_matrix)] <- 0 #matrix with 0s to replace Nas
bitscore_matrix[bitscore_matrix == 0] <- 1 #for a matrix that undergo the transformation, we need to replace all 0s with 1s because of weighting. We can't weight with 0s. 1s are just fine, and don't influence the actual scores (they're always much much higher)

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

write.csv(bitscore_matrix, "aa_align_bitscore_original.csv")
write.csv(bitscore_matrix_2, "aa_align_bitscore_diag_normd_scaled.csv") 

#Figures
pdf(file = "Rplots_heatmap_trees.pdf")

#hmap
library("pheatmap")

matrix_confirm <- as.matrix(bitscore_matrix_2)
mode(matrix_confirm) <- "numeric"
pheatmap(matrix_confirm, cex = 0.3, border_color = NA, main = "Heatmap: Just AA, transformed, scaled")

#dendrogram / tree
library("stats")
for_hclust <- dist(bitscore_matrix_2, method = "canberra", diag = FALSE, upper = FALSE, p = 2)
hc <- hclust(for_hclust, method = "average", members = NULL)
plot(hc, hang = -0.5, cex = 0.23) 

#unrooted tree
library("ape")
myphylo <- as.phylo.hclust(hc)
plot(as.phylo.hclust(hc), type = "unrooted", cex = 0.23,
     no.margin = TRUE)


dev.off()

#maths for upper iqr for net edge-adjustments (***NOTE: it's effectivness at trimming networks is based ss scoring, here it is looking to be less effective, but is better than nothing. Future optimisation here)
library(Matrix)
library(stats)

#for upper half of each matrix, excluding the diagonal
group1 <- matrix_confirm
group1[lower.tri(group1,diag=TRUE)] <- 0

#2a. n = number of observations (half of mat, excluding diag)
n <- ((((length(group1) /2 ) * (length(group1) /2)) / 2) - (length(group1) /2))

#2b. x = mean
x <- mean(group1[group1!=0])

#2c. s = standard deviation
s <- sd(group1[group1!=0])

#2d. z = z-score constant = 1.960 (95 int) / 1.15 (75 int)
z <- 1.960  

#2e. m = median (just for fun)
m <- median(group1[group1!=0])

#2f. r = range (just for fun)
rmin <- min(group1[group1!=0])
rmax <- max(group1[group1!=0])

#2g. 95% confidence interval = x +- (z*(s/sqrt(n)))
# + 
conf_95_pos <- x + (z*(s/sqrt(n)))
conf_95_neg <- x - (z*(s/sqrt(n)))

#3.IQR and upper IQR
g1 <-as.numeric(unlist(group1[group1!=0]))

iqr <- IQR(g1)
iqr_up <- quantile(g1, 3/4)

############
#write input fastasome fo the above stats to a .txt file for usrs. not so useful here as with secondary structure, but may be of interest so some

txt_var <- file("./aa_align_stats.txt")
writeLines(c(paste("Number of observations (1/2 of matrix, excluding diagonal self-self scoring): ", n),
             paste("Align score mean: ", x), paste("Align score standard deviation: ", s), paste("Align score median: ", m),
             paste("Align score min: ", rmin), paste("Align score max: ", rmax), paste("Align score 95% confidence interval (upper): ",conf_95_pos),
             paste("Align score 95% confidence interval (lower): ",conf_95_neg), paste("Align score IQR: ", iqr), 
             paste("Align score upper IQR value (Q3): ", iqr_up)), txt_var)
close(txt_var)

############
#``NOTE DEVELOPMENT HERE``#
###################################
#this is from querying PDB structures pinned at different AA lengths and running through DP 
#(eq = slope of the line of plotting upper IQR (y) against avg. AA len (x) -- see paper for more info)
#T_base <- ((3.657* avg_AA_len) + 260.1)

#this is what users will change if desired -- see figures_3_networks_custom.R and use --networks_custom to use this!
#Tem_custom <- ((3.657* avg_AA_len) + 260.1) + (((3.657* avg_AA_len) + 260.1) * Ts) 
#10-50% scale
t_scale_array <- as.array(c(0.1,0.2,0.3,0.4,0.5))
new_array <- ((3.657*386) + 260.1) + (((3.657* 386) + 260.1) * t_scale_array)

####################################
#4. Now back to the regular network construction
#formatting
matrix_confirm <- as.matrix(matrix_confirm_temp)
mode(matrix_confirm) <- "numeric"   

#actual igraph, two for different layout options
library("igraph")
ig <- graph.adjacency(matrix_confirm, mode="undirected", weighted=TRUE, diag = FALSE, add.colnames = NULL)

community_clustering <- fastgreedy.community(ig)
cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)


#l <- layout <- layout.reingold.tilford(ig, circular=T) #circ is not useful for connections, but is for name reading / ref for IDs / colour grouping
ll <- layout.fruchterman.reingold(ig, niter=1000) #don't change this to the better layout_with_graphopt below - it's not actually better when the data is contains so many edges like this original will before transformations below
#layout_nicely = same as fruchterman reingold in my tests   


####################################
#. Apply edgeweight cutoff
#define pdf outside fo the fucntion
pdf(file = "./networks/networks_Rplots3_all_nodes.pdf")

#FUNction 1 - all nodes, edge-adjusted
Run_net_full <- lapply(new_array, function(x) {

te <- delete_edges(ig, E(ig)[weight< x]) #x = Tem just fyi
community_clustering <- fastgreedy.community(te)
cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)
te2 <- layout.fruchterman.reingold(te, niter=1000)

#plot with edge weight applied
print(plot(te,
     layout=te2,
     edge.arrow.size=0.4,
     vertex.label.cex=0.38,
     vertex.label.family="Helvetica",
     vertex.label.font=1.5,
     vertex.shape="circle",
     vertex.size=2,
     vertex.color=cluster_colours[membership(community_clustering)],
     vertex.label.color="black",
     edge.width=0.3,
     main= sprintf("All nodes: Network at edge-weight cutoff of %.2f",x))) #maybe change this to include 10-50%

#to make those cluster community lists at each cutoff. File title gives the numerical value at each increasing 10%

clu <- components(te) 
clu_who <- groups(clu)
lapply(clu_who, function(clust) write.table(as.data.frame(clust), file=sprintf("./networks/community_csvs/all_nodes_community_clust_at_edge_cutoff_%.2f.csv",x), append= T, sep=',' ))

#to make variations of the 'Final_3_...' matrix for potential downstream applications       
matrix_confirm[matrix_confirm < x] <- 0
write.csv(matrix_confirm, file = sprintf("./networks/adjusted_algn_score_csvs/adjust_net_Final_3_align_matrix_at_%.2f.csv", x), row.names = TRUE)

#maybe change the above two file sets to include mentioning 10-50% scaling in the title...

})
dev.off()

####################################

#6. delete zero edge nodes of edgeweight-adjusted
#gives more spaced + readable clusters but totherwise same as above

pdf(file = "./networks/networks_Rplots3_zero_nodes_deleted.pdf")

Run_net_edge_del <- lapply(new_array, function(x) {
  
te <- delete_edges(ig, E(ig)[weight< x])
dv <- delete.vertices(te, V(te)[degree(te)==0])
community_clustering <- multilevel.community(dv)
cluster_colors <- rainbow(max(membership(community_clustering)), alpha = 0.5)
lld <- layout_with_graphopt(dv, niter=1000,charge = 0.01) #this is BY FAR the best layout here

#plot with edge weight applied and zero edge weight nodes deleted

print(plot(dv,
     layout=lld, 
     edge.arrow.size=0.4, 
     vertex.label.cex=0.38, 
     vertex.label.family="Helvetica",
     vertex.label.font=1.5,
     vertex.shape="circle", 
     vertex.size=2,
     vertex.color=cluster_colors[membership(community_clustering)], 
     vertex.label.color="black", 
     edge.width=0.3, 
     main= sprintf("Zero edge nodes deleted: Network at edge-weight cutoff of %.2f",x)))

#to make those cluster community lists at each cutoff. File title gives the numerical value at each increasing 10%
clu <- components(dv) 
clu_who <- groups(clu)
lapply(clu_who, function(cluster) write.table(as.data.frame(cluster), file=sprintf("./networks/community_csvs/zero_nodes_del_community_clust_at_edge_cutoff_%.2f.csv",x), append= T, sep=',' ))

       
#we don't need to make the 'variations of the 'Final_3_...' matrix for potential downstream applications' again (5 total, not 10)        

})

dev.off()

