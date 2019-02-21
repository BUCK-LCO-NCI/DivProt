#If using biowulf, one has to configure an enviornment that contains pandas, as it doesn't seem to be(?) an availible module with python3
#example:
#source /data/belfordak/Buck_lab/conda/etc/profile.d/conda.sh
#conda activate base
#conda activate project_pandas  (just a py env with pandas)
#python

#run as python3 pre_processing.py seq_00.fa.ss3 (except for all file in dir so some sort of loop?) <- TODO
#this script produces the fastqish file from the Porter5 output
#!/usr/bin/env python

import numpy as np
from math import *
import pandas as pd
import sys, re, os

infile = pd.read_csv(sys.argv[1],delimiter = '\t' )

#verified works
maxval= infile.iloc[:,3:-1].max(axis=1) #to work for both ss3 and ss8

phred_calc = (-10)*np.log10(1-(maxval))

infile['phred_calc'] = phred_calc

#for ASCII generation from phred score:
my_dictionary= {
    (0) : '!',
    (1) : '\"',
    (2) : '#',
    (3) : '$',
    (4) : '%',
    (5) : '&',
    (6) : "\'",
    (7) : '(',
    (8) : ')',
    (9) : '*',
    (10) : '+',
    (11) : ',',
    (12) : '-',
    (13) : '.',
    (14) : '/',
    (15) : '0',
    (16) : '1',
    (17) : '2',
    (18) : '3',
    (19) : '4',
    (20) : '5',
    (21) : '6',
    (22) : '7',
    (23) : '8',
    (24) : '9',
    (25) : ':',
    (26) : ';',
    (27) : '<',
    (28) : '=',
    (29) : '>',
    (30) : '?',
    (31) : '@',
    (32) : 'A',
    (33) : 'B',
    (34) : 'C',
    (35) : 'D',
    (36) : 'E',
    (37) : 'F',
    (38) : 'G',
    (39) : 'H',
    (40) : 'I',
    (41) : 'J',
    (42) : 'K'
}

ASCII = infile['phred_calc'].apply(np.floor)
# with var assignment
infile['ASCII']= ASCII.map(my_dictionary) 

#Transposing
infile_transposed_in = infile.T
#print(infile_transposed_in)

fastqish_but_cols = infile_transposed_in.loc[['SS', 'ASCII'],:]

concat_test_temp = pd.Series(fastqish_but_cols.fillna('').values.tolist()).str.join('')
concat_test = np.array(concat_test_temp)

#Part 2: pulling seq names and "+" separator (lines 0 and 2)

infile_name = sys.argv[1]

#1. read in temp_porter5_submission.swarm as list_seq_XX
#note: I'm not sure if this will be the final path (maybe we'll want it back 1 dir?)
#all_seqs.list created in pre_pre_processing script
with open("all_seqs.list") as f: 
    list_seq_XX = f.readlines()
list_seq_XX = [x.strip() for x in list_seq_XX]

y = sorted(list_seq_XX, key = lambda item: int(item.partition('_')[2].partition('.')[0]))
#print(y) above was for sorting so order was ...9,10,11,12 instead of 9,10,100,101, etc. Lib would have been incorrect otherwise

# for getting rid of .fasta. We're going to just use the seq_00 part to search in dict
x = []
for line in y:
    part = line.split(".")
    #x = part[0]
    #print(x)
    x.append(part[0])
#print(x)
list_seq_XX_fin = x

#2. read original.fa as list, but only lines that start with ">"
original_fasta = sys.argv[2]
with open(original_fasta,"r") as f:
    id_temp = []
    for ln in f:
        if ln.startswith(">"):
            x = ln
            id_temp.append(x)
list_seq_fa_orig = [x.strip() for x in id_temp]

#3. dictionary from the two lists
keys = list_seq_XX_fin #id_temp_clean in test
values = list_seq_fa_orig
dict_for_name_appends = dict(zip(keys, values)) 

#3.5 parse file name for key val search
sep = os.path.basename(infile_name).split('.')[0]

#4. if filename = key, print value as line 0 in file
#replace > with @
dict_var = dict_for_name_appends[sep]

temp = list(dict_var)
temp[0] = "@"
dict_var = "".join(temp)
#print(dict_var)

f_out = open((infile_name + '.fastqish'), 'w')
f_out.write(str(dict_var + "\n" + concat_test[0] + "\n" + "+" + "\n" +concat_test[1]))
f_out.close()
