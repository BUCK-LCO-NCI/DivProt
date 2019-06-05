library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

#read in fastq as stringset biostrings object
#fastqish = sysargv from python script
x <- readBStringSet(fastqish, format="fastq",
                      nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE, with.qualities = TRUE) #with.qualities doesn't work, we're going to have to load them in seperately

#check qual is in there - it's not
as.data.frame(x)

#### We made the file in python before this
ww <- read.csv("./just_qual.csv", header = FALSE, row.names = NULL)
www <- data.frame(lapply(ww, trimws), stringsAsFactors = FALSE) #there's leading(?) whitespace in the df that you have to get rid of with this
Tw <- t(www)

quals_tmp <-  BStringSet(Tw)
quals <- PhredQuality(quals_tmp)

#need to read in just qual vals (shoudn't really -- but you do)
full_obj <- QualityScaledBStringSet(x, quals)

length(full_obj)
#PhredQuality = phred vharacter string
qual_mat = qualitySubstitutionMatrices(fuzzyMatch = c(0, 1), alphabetLength = 8L, qualityClass = "PhredQuality", bitScale = 1) #could change bitscle to 0.5 or 2 are common ones

submat <- matrix(qual_mat,8, 8)

my_alphbabet <- str("G","I","E","B","T","S","C","H")

dimnames(submat) <- list(my_alphbabet, my_alphbabet)

phred_mat_fin <- sapply(full_obj, function(i){
  sapply(full_obj, function(j){
    score(pairwiseAlignment(i, j, substitutionMatrix = submat))
  })
})

write.csv(phred_mat_fin, "./phred_align_matrix.csv") #possibly change location

#####################################
#Produce score_mat and prob_mat actual matrics
#####################################
#this is just copied from the top of the older ss_rscript...i suppose that file will turn into just producing the figures

library("tidyr")
library("dplyr")
library("reshape")

#1. for score_mat
##subject file
score_mat_score_table <- read.table("score_mat_score_table.csv", sep = ",", header = FALSE)
colnames(score_mat_score_table) <- c("qseqid", "sseqid", "align_score")
score_mat_score_table$align_score = as.numeric(gsub("Score=", "", score_mat_score_table$align_score)) #gets rid of Score=

####align score
shaped_algn_score <- reshape(score_mat_score_table,idvar = "qseqid", timevar = "sseqid", direction = "wide")
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

write.csv(algnscore_matrix, "score_align_matrix.csv") 


#2. for prob_mat
##subject file
prob_mat_score_table <- read.table("prob_mat_score_table.csv", sep = ",", header = FALSE)
colnames(prob_mat_score_table) <- c("qseqid", "sseqid", "align_score")
prob_mat_score_table$align_score = as.numeric(gsub("Score=", "", prob_mat_score_table$align_score)) #gets rid of Score=

####align score
shaped_algn_score <- reshape(prob_mat_score_table,idvar = "qseqid", timevar = "sseqid", direction = "wide")
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
write.csv(algnscore_matrix, "prob_align_matrix.csv") 


######################################
###SUM ALL THREE MATRICES
######################################
#read in score_mat

score_mat <- read.csv(file = "./score_align_matrix.csv", header = TRUE, row.names = 1)

#read in prob_mat

prob_mat <- read.csv(file = "./prob_align_matrix.csv", header = TRUE, row.names = 1)

#read in phred_mat
#(or ggg var)
phred_mat <- read.csv(file = "./phred_align_matrix.csv"), header = TRUE, row.names = 1)


temp_df <- cbind(score_mat, prob_mat, phred_mat)

fin_df <- sapply(unique(colnames(temp_df)), 
                 function(x) rowSums(temp[, colnames(temp_df) == x, drop = FALSE]))

write.csv(fin_df, file = "final_3_align_matrix.csv")
