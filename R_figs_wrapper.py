
#python3 run_figures.py alignmatrix.csv -tsne_trees/heatmap/network -k (or -tm for threshold modifier, but actually doesn't matter, just change so people arent confusd) 2
import subprocess
import sys

align_file = sys.argv[1]
script_to_pass = sys.argv[2]

#this loop if for deining deafult values if -k/-tm are not given 
if script_to_pass == "-tsne_trees":
    variable_to_pass = sys.argv[4] if len(sys.argv) >= 4 else '1'
elif scrip t_to_pass == "-network":
    variable_to_pass = sys.argv[4] if len(sys.argv) >= 4 else #'the default cutoff soon-to-come'


if script_to_pass == "-heatmap":
    subprocess.check_call(['Rscript', 'figures_1_heatmap.R'], stderr=subprocess.STDOUT, shell=False) #or just Rscript figures_1_heatmap.R
elif script_to_pass == "-tsne_trees":
    subprocess.check_call(['Rscript', '--args', '$variable_to_pass','figures_2_tsne_trees.R'], stderr=subprocess.STDOUT, shell=False) #pass fastq variable on
elif script_to_pass == "-network":
    subprocess.check_call(['Rscript', '--args', '$variable_to_pass','figures_3_networks.R'], stderr=subprocess.STDOUT, shell=False) #pass fastq variable on
