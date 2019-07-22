#THIS IS A DRAFT, THRESHOLDING METHOD IS ON IT'S WAY

args = commandArgs(trailingOnly=TRUE)
matrix_confirm_temp = args[1]
var_score_to_pass = args[2]
orig_fasta = arg[3]

####################################
#1. calculate orig_fasta average aa length
library(Biostrings)
for_avg_len <- fasta.seqlengths(orig_fasta)
avg_AA_len <- mean(for_avg_len)

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
r <- range(group1[group1!=0])

#2g. 95% confidence interval = x +- (z*(s/sqrt(n)))
# + 
conf_95_pos <- x + (z*(s/sqrt(n)))

####################################
#3. put the conf interval into the equation we derived from PDB querying, see paper methods for details
#This is what people change from wrapper script!
#Ts <-
#default
Ts <- 0.3
Tem <- ((2.824* avg_AA_len) + 504) + (((2.824* avg_AA_len) + 504) * Ts)#this is 30% scaling we found 

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

#EDIT -- circ nets are not included right now
#l <- layout <- layout.reingold.tilford(ig, circular=T) #circ is not useful for connections, but is for name reading / ref for IDs / colour grouping
ll <- layout.fruchterman.reingold(ig, niter=5000) #don't change this to the better layout_with_graphopt below - it's not actually better when the data is contains so many edges like this original will before transformations below
#layout_nicely = same as fruchterman reingold in my tests   

####################################
#5. Apply edgeweight cutoff
te <- delete_edges(ig, E(ig)[weight< Tem])
community_clustering <- fastgreedy.community(te)
cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)
te2 <- layout.fruchterman.reingold(te, niter=5000)


pdf(file = "networks_Rplots3.pdf")

#Now let's finally build those networks
#plot with edge weight applied
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

####################################
#6. delete zero edge nodes of edgeweight-adjusted
#gives more spaced + readable clusters but totherwise same as above
dv <- delete.vertices(te, V(te)[degree(te)==0])

community_clustering <- multilevel.community(dv)
cluster_colors <- rainbow(max(membership(community_clustering)), alpha = 0.5)

lld <- layout_with_graphopt(dv, niter=5000,charge = 0.01) #this is BY FAR the best layout

#plot with edge weight applied and zero edge weight nodes deleted
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

#create a table of node groups
clu <- components(dv) #TO-TO --- OR TE OR IG
clu_who <- groups(clu)
lapply(clu_who, function(x) write.table(as.data.frame(x), 'network_cluster_communities.csv', append= T, sep=',' ))
