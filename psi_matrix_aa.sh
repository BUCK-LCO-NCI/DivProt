#!/usr/bin/env bash

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>log_psi_aa.out 2>&1

module load blast
module load R

chmod +x psi_matrix.sh

DIRNAME=$2
QUERYFILE=$1
SUBJECTFILE=$1

mkdir -p $DIRNAME
cd $DIRNAME

psiblast -query "../$QUERYFILE" -subject "../$SUBJECTFILE" -outfmt 6 > PSIBLAST_results.txt

module load R #incase you forgot ;)
Rscript ../R_script_aa.R
