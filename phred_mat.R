library(Biostrings)
args = commandArgs(trailingOnly=TRUE)

#read in fastq as stringset biostrings object
#fastq = sysargv from python script
x <- readBStringSet(fastq, format="fastq",
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

write.csv(phred_mat_fin, "./phred_mat_align.cs") #possibly change location

#####################################
#Produce score_mat and prob_mat actual matrics
#####################################






######################################
###SUM ALL THREE MATRICES
######################################
#read in score_mat

score_mat <- read.csv(file = , header = TRUE, row.names = 1)

#read in prob_mat

prob_mat <- read.csv(file = , header = TRUE, row.names = 1)

#read in phred_mat
#(or ggg var)
phred_mat <- read.csv(file = , header = TRUE, row.names = 1)


temp_df <- cbind(score_mat,prob_mat, phred_mat)

fin_df <- sapply(unique(colnames(temp_df)), 
                 function(x) rowSums(temp[, colnames(temp_df) == x, drop = FALSE]))

write.csv(fin_df, file = "final_3_align_matrix.csv")
