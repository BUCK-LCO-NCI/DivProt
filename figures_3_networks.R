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

#EDIT -- circ nets are not included right now
#l <- layout <- layout.reingold.tilford(ig, circular=T) #circ is not useful for connections, but is for name reading / ref for IDs / colour grouping
ll <- layout.fruchterman.reingold(ig, niter=5000) #don't change this to the better layout_with_graphopt below - it's not actually better when the data is contains so many edges like this original will before transformations below
#layout_nicely = same as fruchterman reingold in my tests   

#Apply edgeweight cutoff

#this is the math to compute the threshold expect modified equation value. See the paper for further information.
library(Matrix)
library(stats)

#for upper half of each matrix, excluding the diagonal
group <- matrix_confirm
group[lower.tri(group1,diag=TRUE)] <- 0

#1. n = number of observations (half of mat, excluding diag)
n <- ((((length(group1) /2 ) * (length(group1) /2)) / 2) - (length(group1) /2))

#2. x = mean
x <- mean(group1[group1!=0])

#3. s = standard deviation
s <- sd(group1[group1!=0])

#4. z = z-score constant = 1.960 (95 int) / 1.15 (75 int)
z <- 1.960  

#5. m = median (just for fun)
m <- median(group1[group1!=0])

#6 r = range (just for fun)
r <- range(group1[group1!=0])

# 95% confidence interval = x +- (z*(s/sqrt(n)))
# + 
conf_95_pos <- x + (z*(s/sqrt(n)))

#put the conf interval into the equation we derived from PDB querying, see paper methods for details

#This is what people change from wrapper script!
#Ts <-
#avg_AA_len <- uhhhhhh
#default
Ts <- 0.3
Tem <- ((2.824* avg_AA_len) + 504) + (((2.824* avg_AA_len) + 504) * Ts)#this is 30% scaling we found 


#######
te <- delete_edges(ig, E(ig)[weight< Tem])
community_clustering <- fastgreedy.community(te)
cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)
te2 <- layout.fruchterman.reingold(te, niter=5000)

######
pdf(file = "networks_Rplots3.pdf")

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

#For delete zero edge nodes of edgeweight-adjusted
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
