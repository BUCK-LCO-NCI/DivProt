#tiny lil' script for aa <-> ss alignment score weighting
#as apparent below, order matters! feed in the aa matrix first
# run as: (sinteractive -> module load R) Rscript post_align_weight.R aa_align_ scores.csv 0.25 ss_align_scores.csv 0.75
#!/usr/bin/env Rscript
#`args<-commandArgs(TRUE)

args <- commandArgs(trailingOnly = TRUE) 

filename1 <- args[1]
filename2 <- args[3]

aa <- read.csv(file = filename1, header = TRUE, row.names = 1)
ss <- read.csv(file = filename2, header = TRUE, row.names = 1)

x = as.numeric(args[2])
y = as.numeric(args[4])

aa_weighted <- (aa * x)
ss_weighted <- (ss * y)

#head(aa_weighted)


temp_df <- cbind(aa_weighted, ss_weighted)
fin_df <- sapply(unique(colnames(temp_df)), 
       function(x) rowSums(temp[, colnames(temp_df) == x, drop = FALSE]))

write.csv(fin_df, file = "aa_ss_weighted_aligments.csv")

