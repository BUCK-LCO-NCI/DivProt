from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
import itertools, sys
from Bio.SubsMat import MatrixInfo
from Bio import Align
import warnings

fastqish = sys.argv[1]

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

###################################################################
####### MATRIX 1: scoring matrix
###################################################################

matrix = {
("G", "C"): -3.0, ("G", "H"): -2.0, ("G", "G"): 5.0, ("G", "I"): -2.0, ("G", "E"): -4.0, ("G", "B"): -3.0, ("G", "T"): -3.0,
("G", "S"): -3.0,("I", "C"): -3.0, ("I", "H"): -2.0, ("I", "G"): -2.0, ("I", "I"): 5.0, ("I", "E"): -4.0, ("I", "B"): -3.0,
("I", "T"): -3.0, ("I", "S"): -3.0,("E", "C"): -3.0, ("E", "H"): -6.0, ("E", "G"): -4.0, ("E", "I"): -4.0, ("E", "E"): 5.0,
("E", "B"): -2.0, ("E", "T"): -3.0, ("E", "S"): -3.0, ("B", "C"): -3.0, ("B", "H"): -3.0, ("B", "G"): -3.0, ("B", "I"): -3.0,
("B", "E"): -2.0, ("B", "B"): 5.0, ("B", "T"): -3.0, ("B", "S"): -3.0, ("T", "C"): -3.0, ("T", "H"): -3.0, ("T", "G"): -3.0,
("T", "I"): -3.0, ("T", "E"): -3.0, ("T", "B"): -3.0, ("T", "T"): 5.0, ("T", "S"): -3.0, ("S", "C"): -3.0, ("S", "H"): -3.0,
("S", "G"): -3.0, ("S", "I"): -3.0, ("S", "E"): -3.0, ("S", "B"): -3.0, ("S", "T"): -3.0, ("S", "S"): 5.0, ("C", "C"): 5.0,
("C", "H"): -3.0, ("C", "G"): -3.0, ("C", "I"): -3.0, ("C", "E"): -3.0, ("C", "B"): -3.0, ("C", "T"): -3.0, ("C", "S"): -3.0
}

from Bio import SubsMat
from Bio.SubsMat import MatrixInfo as matlist

score_table_1 = []
alignment_out_1 = []

def funct_score_mat():
    for key1 in final_dictionary:
        for key2 in final_dictionary:
            if len(final_dictionary[key1]) < 30:
                warnings.warn("You have submitted one or more sequences that contain less than 30 characters. Sequences of this size are typically of low complexity in secondary structure, and thus results regaring them can be less meaningful, amd should be regarded with less confidence")
                #
                align_1 = pairwise2.align.globalds(final_dictionary[key1],final_dictionary[key2],matrix, -10,-1, one_alignment_only=True) #NOTEEEE: The -100 mismatch score is completely arbitrary, but it has to be there as a placeholder, otherwise only matches will be called on from the input matrix, not the mismatches. weird i kno, i mean it's obvioulsy overwritten by the matrix, but weird. oh but gap extend matters
                string_1 = format_alignment(*align_1[0])
                score_table_1.append(key1 + "," + key2 + "," + string_1.splitlines()[3])
                alignment_out_1.append(key1 + "::" + key2 + "\n" + string_1)

#The matrix and list output will be created after matrix 2 (probability martrix) to enable multithreading 
#funct_score_mat()
#with open ("score_mat_score_table.csv", "w") as f:
#    for item in score_table_1:
#            f.write("%s\n" % item)
            
#with open ("score_mat_align_out.txt", "w") as f:
#    for item in alignment_out_1:
#            f.write("%s\n" % item)       
            
#so the above returns an error:
#SystemError: PyEval_EvalFrameEx returned a result with an error set
#The above exception was the direct cause of the following exception:
#but everything is working, the output variables are good, so idk we're good? i'm pretty sure we're good


###################################################################
####### MATRIX 2: substitution matrix (log-odds probability matrix)
###################################################################

##STILL IN PROGRESS##

import os, io
from itertools import combinations, count
from itertools import zip_longest, islice
import uuid
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo


just_fa = []

#make this all a function
with open(sys.argv[1]) as fin:
     paired = zip_longest(*[iter(fin)] * 2, fillvalue='')
     every_other = islice(paired, None, None, 2)
     for lines in every_other:
             line1, line2 = lines
             just_fa.append(list(lines))

for lst in just_fa:
    for j, item in enumerate(lst):
        lst[j] = item.replace('\n', '')
        
just_fa_2 =  list(itertools.combinations(just_fa, 2))

for x in just_fa_2:            
    for n, i in enumerate(x):
        for j,a in enumerate(i):
            if '@' in a:
                x[n][j] = a.replace('@','>')
                
    outname = str(uuid.uuid1()) #+ ".fasta"
    outname_aln = str(outname + ".aln")
    
    with open(outname, "w") as openfile:
        for e in itertools.chain.from_iterable(x):
            openfile.write(e+'\n')
            
    cline = ClustalwCommandline("clustalw2", infile = outname) #doesn't actually need the ".fasta"
    stdout, stderr = cline()
    alignment = AlignIO.read(open(outname_aln), "clustal") #alphabet = alpha)
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    replace_info = summary_align.replacement_dictionary()
    arm = SubsMat.SeqMat(replace_info)
    custom_sub_mat = SubsMat.make_log_odds_matrix(arm)
    
    align_2 = pairwise2.align.globalds(x[0][1],x[1][1], custom_sub_mat, -10,-1, one_alignment_only=True)
    string_2 = format_alignment(*align_2[0]) 
    
    score_f = open("prob_mat_score_table.temp.csv" % i, 'a') #NEED TO CLEAN THIS UP! won't be read into Rscript correctly as is
    score_f.write(x[0][0] + "," + x[1][0] + "," + string_2.splitlines()[3] + "\n")
    score_f.close()
    
    align_f = open("prob_mat_align_out.temp.txt", 'a')
    align_f.write(x[0][0] + "::" + x[1][0] + "\n" + string_2)
    align_f.close()
    
#transformation from negative numbers for score table only -> will do alignment file too

integers = open('test_temp_table.csv', 'r')
third2file = []
smallestInt = float('inf')

for line in integers:
        linesplit = line.strip().split(",")
        third = linesplit[2]
        linesplit2 =third.strip().split("=")
        third2 = (linesplit2[1])
        third2file.append(third2)

fin_t3f = list(map(float, third2file))
min_val = min(fin_t3f)
listlen = len(third2file)

integers = open('prob_mat_score_table.temp.csv', 'r')
for line in integers:
        linesplit = line.strip().split(",")
        third = linesplit[2]
        linesplit2 =third.strip().split("=")
        third2 = (linesplit2[1])
        third2file.append(third2)
        score_f_fin = open("prob_mat_score_table.csv","a")
        l = [float(linesplit2[1]),abs(min_val),float(1)] #this = (original score + abslute val lowest score + 1 ) to make everything positive
        sumvar = sum(l)
        score_f_fin.write(linesplit[0] + "," + linesplit[1] + "," + "  " + "Score=" + str(sumvar) + "\n")
        score_f_fin.close()

        
#delete mat2 unnecessary files now
dir_name = "./" #or subdir if i make things better
dd = os.listdir(dir_name)

for item in dd:
    if item.endswith((".aln", ".dnd", ".temp.csv")): #still need to delete the fasta it makes, that doesn't have a file extension right now 
        os.remove(os.path.join(dir_name, item))            
            
            
            
            
            
            
            
            
#multithreading and running            

#funct_prob_mat()

import threading

if __name__ == "__main__":
    # creating thread
    t1 = threading.Thread(target=funct_score_mat)
    t2 = threading.Thread(target=funct_prob_mat)


    # starting thread 1
    t1.start()
    # starting thread 2
    t2.start()

    # wait until thread 1 is completely executed
    t1.join()
    # wait until thread 2 is completely executed
    t2.join()

    # both threads completely executed
    print("alignmnets 1/3 and 2/3 are done")
 

#######
#Now we'll create the files from both alignments  
######  

#Mat 1
with open ("score_mat_score_table.csv", "w") as f:
    for item in score_table_1:
            f.write("%s\n" % item)
            
with open ("score_mat_align_out.txt", "w") as f:
    for item in alignment_out_1:
            f.write("%s\n" % item)   

            
            
#Mat 2 -- nope, this needs to be defined above with concatination, shouldn't all be saved as a variable
#with open ("prob_mat_score_table.csv", "w") as f:
   # for item in score_table_2:
         #   f.write("%s\n" % item)
            
#with open ("prob_mat_align_out.txt", "w") as f:
   # for item in alignment_out_2:
          #  f.write("%s\n" % item) 
            
###################################################################
####### MATRIX 3: quality matrix (phred from structure prediction) 
###################################################################
#This in in R b/c of Biostrings

#create a just-phred file to read into R
all_lines = []
with open('All_psi_plus_SS8.fastqish') as f:
    x = [v for i, v in enumerate(f, start=1) if i % 4 == 0]
    all_lines.append(x)
    
print(all_lines)


with open('just_qual.csv', mode='w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows(all_lines)

csv_file.close()

#### we're just going to run the R script instead of using RPy2
import subprocess
subprocess.check_call(['Rscript', '--args', '$fastq','phred_mat.R'], stderr=subprocess.STDOUT, shell=False) #pass fastq variable on
#this is fine, EXCEPT FOR reading in the .fastqish file...I need the R script to take the same sysargv[1] as this script...



###################################################################
####### Merge all marticies together --probably just do in the same R script as above
###################################################################

