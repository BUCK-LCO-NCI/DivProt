#python3 run_figures.py alignmatrix.csv -heatmap
#python3 run_figures.py alignmatrix.csv -tsne_trees -k 2
#python3 run_figures.py alignmatrix.csv -network original_fasta_name.fa 
#python3 run_figures.py alignmatrix.csv -network_custom original_fasta_name.fa -tm 1.2 

import subprocess
import sys

align_file = sys.argv[1]
script_to_pass = sys.argv[2]
orig_fasta = sys.argv[3]


#this loop if for deining deafult values if -k/-tm are not given 
if script_to_pass == "-tsne_trees":
    var_score_to_pass = sys.argv[4] if len(sys.argv) >= 4 else '1'
elif script_to_pass == "-network_custom":
    var_score_to_pass = sys.argv[4] if len(sys.argv) >= 5 else 'XXXXXX'#'the default cutoff soon-to-come' #NOTEEE TO ME #>=4?


if script_to_pass == "-heatmap":
    subprocess.call(['Rscript', 'figures_1_heatmap.R', align_file], stderr=subprocess.STDOUT, shell=False) #or just Rscript figures_1_heatmap.R
elif script_to_pass == "-tsne_trees":
    subprocess.call(['Rscript', 'figures_2_tsne_trees.R', align_file, var_score_to_pass], stderr=subprocess.STDOUT, shell=False) #pass fastq variable on
elif script_to_pass == "-network":
    subprocess.call(['Rscript', 'figures_3_networks.R', align_file, orig_fasta], stderr=subprocess.STDOUT, shell=False) #pass fastq variable on
elif script_to_pass == "-network_custom":
    subprocess.call(['Rscript', 'figures_4_networks_custom.R', align_file, var_score_to_pass, orig_fasta], stderr=subprocess.STDOUT, shell=False) #pass fastq variable on
