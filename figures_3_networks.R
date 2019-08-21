#STILL A DRAFT

args = commandArgs(trailingOnly=TRUE)
matrix_confirm_temp = args[1]
orig_fasta = arg[2]


dir.create("./networks")
dir.create("./networks/adjusted_algn_score_csvs")
dir.create("./networks/community_csvs")

####################################
#1. calculate orig_fasta average aa len info
library(Biostrings)
for_avg_len <- fasta.seqlengths(orig_fasta)
avg_AA_len <- mean(for_avg_len)
min_AA_len <- min(for_avg_len)
max_AA_len <- max(for_avg_len)

####################################
#2. this is the math to compute the threshold expect modified equation value. See the paper for further information.
library(Matrix)
library(stats)

#for upper half of each matrix, excluding the diagonal
group <- matrix_confirm
group[lower.tri(group1,diag=TRUE)] <- 0

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
library(stats)
g1 <-as.numeric(unlist(group1[group1!=0]))

iqr <- IQR(g1)
iqr_up <- quantile(g1, 3/4)

####################################
#write input fasta stats to a .txt file for people to view

txt_var <- file("./networks/input_fa_and_align_stats.txt")
writeLines(c(paste("Average AA length: ", avg_AA_len), paste("Range AA length lower: ", min_AA_len), paste("Range AA length upper: ", max_AA_len),
             paste("Number of observations (1/2 of matrix, excluding diagonal self-self scoring): ", n),
             paste("Align score mean: ", x), paste("Align score standard deviation: ", s), paste("Align score median: ", m),
             paste("Align score min: ", rmin), paste("Align score max: ", rmax), paste("Align score 95% confidence interval (upper): ",conf_95_pos),
             paste("Align score 95% confidence interval (lower): ",conf_95_neg), paste("Align score IQR: ", iqr), 
             paste("Align score upper IQR value (Q3): ", iqr_up)), txt_var)
close(txt_var)

####################################
#this is from querying PDB structures pinned at different AA lengths and running through DP 
#(eq = slope of the line of plotting upper IQR (y) against avg. AA len (x) -- see paper for more info)
#T_base <- ((3.657* avg_AA_len) + 260.1)

#this is what users will change if desired -- see figures_3_networks_custom.R and use --networks_custom to use this!
#Tem_custom <- ((3.657* avg_AA_len) + 260.1) + (((3.657* avg_AA_len) + 260.1) * Ts) 
#10-50% scale
t_scale_array <- as.array(c(0.1,0.2,0.3,0.4,0.5))
new_array <- ((3.657*386) + 260.1) + (((3.657* 386) + 260.1) * t_scale_array))

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
     edge.arrow.size=0.5,
     vertex.label.cex=0.5,
     vertex.label.family="Helvetica",
     vertex.label.font=1.5,
     vertex.shape="circle",
     vertex.size=3,
     vertex.color=cluster_colours[membership(community_clustering)],
     vertex.label.color="black",
     edge.width=0.3,
     main="[title]"))
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
     edge.arrow.size=0.5, 
     vertex.label.cex=0.5, 
     vertex.label.family="Helvetica",
     vertex.label.font=1.5,
     vertex.shape="circle", 
     vertex.size=3,
     vertex.color=cluster_colors[membership(community_clustering)], 
     vertex.label.color="black", 
     edge.width=0.4, 
     main= "[title]"))
})

dev.off()

####################################
#create a table of node groups
clu <- components(dv) #TO-TO --- OR TE OR IG ... make this recursive for all 10 permutations, put in a new dir
clu_who <- groups(clu)
lapply(clu_who, function(x) write.table(as.data.frame(x), 'network_cluster_communities.csv', append= T, sep=',' ))

####################################
#need to add in adjust input Final_align.csv mat to have 5 new ones reflecting 10-50% scaling with < cutoff supplanted with 0s in new dir
#useful if people want to use them for downstream/external applications, such as cytoscape for prettier nets or iqtree/figtree/MAFFT&phylo.io for prettier trees
