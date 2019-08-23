#STILL A DRAFT
#nets custom
args = commandArgs(trailingOnly=TRUE)
matrix_confirm_temp = read.csv(args[1], header = TRUE, row.names = 1, check.names = FALSE)
var_score_to_pass = args[2]
orig_fasta = args[3]


dir.create("./network_custom")
dir.create("./network_custom/adjusted_algn_score_csv")
dir.create("./network_custom/community_csv")

#you'll see this is pretty much just a simplified version of the main networks script
Ts <- as.numeric(var_score_to_pass)

#T_base <- ((3.657* avg_AA_len) + 260.1)

#for avg aa len
library(Biostrings)
for_avg_len <- fasta.seqlengths(orig_fasta)
avg_AA_len <- mean(for_avg_len)


Tem_custom <- ((3.657* avg_AA_len) + 260.1) + (((3.657* avg_AA_len) + 260.1) * Ts) 


#1. Network construction
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
ll <- layout.fruchterman.reingold(ig, niter=1000) #don't change this to the better layout_with_graphopt below - it's not actually better when the data is contains so many edges like this original will before transformations below
#layout_nicely = same as fruchterman reingold in my tests   

####################################
#2. Apply edgeweight cutoff
te <- delete_edges(ig, E(ig)[weight< Tem_custom])
community_clustering <- fastgreedy.community(te)
cluster_colours <- rainbow(max(membership(community_clustering)), alpha = 0.5)
te2 <- layout.fruchterman.reingold(te, niter=1000)


pdf(file = "./network_custom/networks_Rplots3_CUSTOM.pdf")

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

lld <- layout_with_graphopt(dv, niter=1000,charge = 0.01) #this is BY FAR the best layout here

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
     main="Custom edge-adjusted network")


dev.off()
##########################
#create a table of node groups
clu <- components(dv) 
clu_who <- groups(clu)
lapply(clu_who, function(x) write.table(as.data.frame(x), './network_custom/community_csv/network_CUSTOM_cluster_communities.csv', append= T, sep=',' ))

##########################
#add make new matrix with 0 supplement at < cutoff
matrix_confirm_temp[matrix_confirm_temp < Tem_custom] <- 0

write.csv(matrix_confirm_temp, file = "./network_custom/adjusted_algn_score_csv/adjusted_network_CUSTOM_Final_3_align_matrix.csv", row.names = TRUE)
