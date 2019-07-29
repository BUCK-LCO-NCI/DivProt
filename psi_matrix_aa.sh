#!/usr/bin/env bash

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>log_aa.out 2>&1

module load blast
module load R

chmod +x psi_matrix.sh

DIRNAME=$2
QUERYFILE=$1
SUBJECTFILE=$1

mkdir -p $DIRNAME
cd $DIRNAME

psiblast -query "../$QUERYFILE" -subject "../$SUBJECTFILE" -outfmt 6 -evalue 100 > PSIBLAST_results.txt
#note the VERY liberal e-value. You can change this if you want, but keep in mind that you will not be able to weight with aa results that are incomplete, i.e. you have no blast results for two seqs (different data sizes cannot be merged) 
#also note though that while e-val 100 is almost comical, in our testing of short (20-55), divergent seqs, even e-val 100 is too strict...as is 10,000. Just keep in mind your data and what can be done with it

module load R #incase you forgot ;)
Rscript ../R_script_aa.R
