#!/bin/bash
# example of how to run this script "bash pre_pre_processing.sh Adoma_polyoma_LTandVP1.fasta"

#1: parse input file
mkdir split_out
cd split_out

csplit -s -b "%02d.fasta" -f seq_ -z ../$1 '/>/' {*}
#for fastq: csplit -s -b "%02d.fasta" -f seq_ -z ./testt.fa '/@/' {*} --we should probably integrate this so someone can supply a fastq in additon to fasta

ls seq_* > all_seqs.list #not doing .fasta to avoid input fa getting put in there
#2 put all files as porter5 job submisstions in one .swarm file
ls *.fasta > temp_porter5_submission.swarm

awk '{print "python3 Porter5/Porter5.py -i", $0 , "--cpu 4 --fast"}' temp_porter5_submission.swarm > porter5_submission.swarm

#Porter5 has to be optimised with paths to hhblits, psi-blast, etc...? I'm not sure if it will be a problem?
printf "Now submit a swarm job with the Porter5_swarm_submission.swarm file in ./split_out\n"
