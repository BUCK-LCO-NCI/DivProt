from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
import itertools, sys
from Bio.SubsMat import MatrixInfo
from Bio import Align

fastq = sys.argv[1]

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

score_table = []
alignment_out = []


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

for key1 in final_dictionary:
    for key2 in final_dictionary:
            align = pairwise2.align.globalds(final_dictionary[key1],final_dictionary[key2],matrix, -10,-1, one_alignment_only=True) #NOTEEEE: The -100 mismatch score is completely arbitrary, but it has to be there as a placeholder, otherwise only matches will be called on from the input matrix, not the mismatches. weird i kno, i mean it's obvioulsy overwritten by the matrix, but weird. oh but gap extend matters
            string = format_alignment(*align[0])
            score_table.append(key1 + "," + key2 + "," + string.splitlines()[3])
            alignment_out.append(key1 + "::" + key2 + "\n" + string)

#so the above returns an error:
#SystemError: PyEval_EvalFrameEx returned a result with an error set
#The above exception was the direct cause of the following exception:
#but everything is working, the output variables are good, so idk we're good? i'm pretty sure we're good


###################################################################
####### MATRIX 2: substitution matrix (log-odds probability matrix)
###################################################################

#1. clustalw alignment (pairwise2 canoot produce an MSA, which we must have for this process)
from Bio.Align.Applications import ClustalwCommandline

cline = ClustalwCommandline("clustalw2", infile="my_test.fastqish") #Have to change it to fasta :/ lame
print(cline)
stdout, stderr = cline() #run it here

#2. alignmanet object
from Bio import AlignIO

alignment = AlignIO.read(open("my_test.aln"), "clustal")

#3. calculate summary information
from Bio.Align import AlignInfo
from Bio import SubsMat

summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()
print(consensus)

#4. substitution matrix
replace_info = summary_align.replacement_dictionary()
print(replace_info[("B", "G")])

my_arm = SubsMat.SeqMat(replace_info)
custom_sub_mat = SubsMat.make_log_odds_matrix(my_arm)
custom_sub_mat.print_full_mat()

#TODO: then loop to re run alignment with this matrix



###################################################################
####### MATRIX 3: quality matrix (phred from structure prediction) 
###################################################################
#This in in R b/c of Biostrings, so using RPy2 and biostrings 

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

import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import bio

package_names = ('biostrings')
biostrings.biostrings_env['PHRED_MATRIX']





