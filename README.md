# DivProt
Alignment programme for divergent proteins based on amino acid sequence and predicted structure, producing figures that map divergent proteins in proximal phylogenetic space. Alignment 

DivProt works by iteratively aligning sequences, storing the scores (bitscore, e-value, or alignment score for amino acids, and alignment score for secondary structure) in a matrix, and producing figures (heatmap, phylogenetic trees, and networks) to help visualise evolutionary relationships.

For the amino acid method, DivProt uses PSI-BLAST to pull conserved amino acid domains between the input sequences.
For secondary structure, DivProt utilizes Porter5 [https://github.com/mircare/Porter5] to predict structure (using either HHBLITS or both PSI-BLAST and HHBLITS) and custom aligns with a 3-tierd alignemnt strategy that utilises a scoring matrix based on structure grouping, a log-odd probability matrix, and predicted structure probablilty. Amino acid weighting from running your dataset through the PSI-BLAST method can influence your secondary structure scoring (or vise-versa) according to input specifications.

[TODO: figure here]

This is compiled to work on the NIH's Biowulf cluster, and will need adjusting to run locally. When running on biowulf, call an interactive session to avoid any jobs being killed if they take up too much memory + set up a custom enviornment to run python modules not automatically availible (instruction on this are at the top of the pre_pre_proc file).

Anywhere below where you see .ssX, I am referring to the Porter5 output files .ss3 or .ss8. .ssX denotes that you should specify which one you are working with. 


#### Dependencies
1. Python3
  - pandas
  - NumPy
  - Biopython
  - Biotite (if interested in 2D visualisation of structure)
  

2. Porter5
  - hhsuite
  - psiblast
  - uniprot (cp -r /fdb/hhsuite/uniprot20_2016_02 .)
  - uniref90
  
3. R
  - Packages: Biostrings, pheatmap, igraph, dplyr, tidyr,reshape, tsne(?)

4. clustalW

## For running on conserved amino acid domains:

Running the amino acid pipeline is extremely simple. Just upload your protein sequences and run the aligner script for PSIBLAST. Figures will be automatically generated.

Note: the R script currently outs three matrices, the output of your iterative psi-blast alignemnts by evalue, persent identity, and bitscore. However, the figures produced by the script are fed only by the bitscore matrix, as we believe this one to be the most informative of the three pis-blast outputs. You can certainly change this though by going into the R script (R_script_aa) and replacing any instances of "bitscore_matrix" below line 70 with "evalue_matrix" or "pident_matrix". Figures will then be created with that data.

## For running on predicted secondary structure:


### Step 0.
You'll need to set up Porter5 with a single sequence before you can run all in your fasta file. Porter5 requests the paths to its dependencies. I do this by pulling the first sequence in my fasta file, creating a new .fa with that, and submitting that to Porter. This lets me input my paths and set everything up. You only need to do this once.

> Example of how I configure the paths (this will be almost the same for anyone else running on biowulf):

>uniref90 = /fdb/SIFT/uniref90.fa

>hhblits = hhblits [I call the hhsuite module with the job submission, or /usr/local/apps/hhsuite/3.0-beta.3/bin/hhblits without doing that should work too]

>uniprot = /data/belfordak/Buck_lab/Divergent_prots/DivProt/uniprot20_2016_02/uniprot20_2016_02 (cp -r uniprot20_2016_02 db to a dir you have permission to - giving the path to where it's hosted on biowulf will cause an error when you run Porter b/c you don't have permission to access the db in this location)

>psiblast = /usr/local/apps/ncbi-toolkit/21.0.0/bin/psiblast

(If found that your prediction works and you get your .ssX files, you're good to go. If you get an error message, it is likely due to your paths, even though there was probably no error called when you set them up. Also, just fyi, the error messages from Porter don't always obviously direct you to the problem (an error telling you a certain psi-blast related file could not be found doesn't necessarily mean there is anything wrong with your psi-blast directory/module/path).
  
### Step 1. 
Run the pre_pre_processing script on your fasta input file. Both single and multiline fasta are fine.
HOWEVER, you cannot have any non-sequence characters in your file outside of the header. 

~ Note:  You can't have a * at the end of each protein sequence. This will cause Porter5 to fail

Ex.
```
bash pre_pre_processing.sh Adoma_polyoma_LTandVP1.fasta
```

### Step 2. 
Submit a swarm job like the prompt at the end of the pre_pre_processing.sh script says to do with the outfile. Do this from within the DivProt directory (no subdirectory), linking to the location of the file in ./split_out. 

### Step 3.
Run the pre_procesing.py script to take the Porter5 output and produce what I'm calling "fastqish" files.

They're basically fastq files (but they don't contain the line numbers I've included here):

>1. @Arowana_adomavirus_LT | 14901:17525 Reverse
>2. CCCHHHHHHHHHHHHHCCCCCCCCHHHHHHHHHHHHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHHCCCC
>3. "+"
>4. !"#$%'(*++++*'%$"!!!!"##$&'+----,,+++*))&#!!!!!!!#'(+-..--+(%$##$$##"!

You want to run this script on all of the files in the directory (all individual sequences from your input, as Porter5 only takes one at a time)

So in the terminal, run:
```
for f in *.ssX; do python3 ../pre_processing.py $f ../original_fasta.fasta; done
```
Once this is done, you'll want to create one merged file of all the XX.ssX.fastqish files. This file contains the sequence name, secondary structure prediction, and phred score of each position structure.
Create the file with something like:
```
awk 1 *ssX.fastqish > whatever_file_name_you_want (e.x. original_fasta_name.ss3.fastqish)
```
~ Note: we *highly* reccommend using the .ss8 files for meaningful alignments ~

### Step 4
Run the aligner to produce a matrix of alignment scores of all your input structure sequences. 

Run aligner:
```
$ module load R/3.5
$ python3 aligner_finalish.py split_out/original_fasta_name.ss3.fastqish
```

### Step 5
Run the R_figs_wrapper.py script to generate figures after loading an R module. You need to specify which of three figures you want after writing the file in the command line (examples below). Figures were seperated out to accomidate users playing around with papameters without generating lots of unnecessary figures.

The three figures:
1. Heatmap
2. Phylogenetic tree(s) 
3. Networks

2 has the optional varible parameter "-k" for cluster. If your input data is functionally divergent and does not contain a common ancestor (i.e. no relationship between them should reasonably be mapped, and they should not be connected on a tree) then you denote this with the k value when running the script. So denoting -k 3 would produce three phylogenetic trees. Default is 1. Clustering is done on k-means, following TSNE reduction. The cluster plot is included in the  output file.

Example of producing trees:
```
$ python3 run_figures.py input_align.csv -tsne_trees -k 3
```

3 contains an optional variable parameter "-tm" for threshold modifier. The methodolgy to build the networks contains a cutoff value that reduces erroneoud connection between nodes. Without this, given the sensitivity of DivProt, likely every protein would be connected. This cutoff value ..... (see methods section of paper).... . The default value is calculated based on the average length of your input sequence lengths. This value is printed to the termainal, so you can know at what value to start increasing or decreasing, if you wish to view your network with more or less stringency. We encourage the user to read the paper methods section before this, and consider biological meaningfulness when altering such a parameter.

Example of producing networks:
```
$ python3 run_figures.py input_align.csv -networks -tm XXX
```

.... 

#### On the structure:
> XX.ss3: helix (H), strand (E), and coil (C)

> XX.ss8: α-helix (H), π-helix (I), 3<sub>10</sub> helix (G), β-stand (E), β-bridge (B), H bonded turn (T), bend (S), and Coil and others (C)


...That's it for now...

Notes and TODO:
> 1. Add look-up name file for users (seq_07 = input_07_actual_id)
> 2. There currently is not a limit to fasta input size, except figures will only be able to scale up to a certain degree, and will get difficult to read after a certain number of sequences are used
> 3. Add dynamic network edge cutoffs
> 4. Add .log file for ss aligner
> 5. Add multithreading
> 6. Add error messages and --help
> 7. collect Porter output in a new dir (or write something about usr creating and moving into a new dir, then executing pre_processing. will need to change that outfile a little)

