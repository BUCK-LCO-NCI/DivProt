#!/usr/bin/env python3


#Rscript and fqish file need to be in the same dir

from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
import itertools, sys

fastq = sys.argv[1]
#CAN have dictionary created in one step -- on the TODO list
id_list = []
for seq_record in SeqIO.parse(sys.argv[1], "fastq"):
    x = (seq_record.id)
    id_list.append(x)


seq_list = []
for seq_record in SeqIO.parse(fastq, "fastq"):
    x = (seq_record.seq)
    seq_list.append(x)

keys = id_list
values = seq_list
final_dictionary = dict(zip(keys,values))

# TODO:Combine dict create into one
score_table = []
alignment_out = []

for key1 in final_dictionary:
    for key2 in final_dictionary:
            align = pairwise2.align.globalms(final_dictionary[key1],final_dictionary[key2],1,-1,0,0, one_alignment_only=True)
            string = format_alignment(*align[0])
            score_table.append(key1 + "," + key2 + "," + string.splitlines()[3])
            alignment_out.append(key1 + "::" + key2 + "\n" + string)

with open ("struct_score_table.csv", "w") as f:
    for item in score_table:
            f.write("%s\n" % item) #still has 2 char gap ahead of Score but its cool

with open ("alignment_out.txt", "w") as f:
     for item in alignment_out:
             f.write("%s\n" % item) #looks fab

import subprocess
subprocess.check_call(['Rscript', 'R_script_strct.R'], shell=False)
