# DivProt
Simple alignment method for divergent proteins based on amino acid sequence and predicted structure, producing figures that map divergent proteins in proximal phyologentic space.

DivProt works by iteritevly aligning sequences, storing the scores (bitscore, e-value, or alignment score for amino acids, and alignment score for secondary structure) in a matrix, and producing figures (eatmap, phylogenetic tress, and networks) to help visualise evolutionary relationships.

For the amino acid method, DivProt uses PSI-BLAST to pull conserved amino acid domains between the input sequences.
For secondary structure, DivProt utilizes Porter5 [https://github.com/mircare/Porter5] to predict structure (using either HHBLITS or both PSI-BLAST and HHBLITS) and custom aligns with [either blomsum62 or custom scoring with account to phred score caluculated from Porter5 output confidence in assignment]

This is compiled to work on the NIH's Biowulf cluster, and will need adjusting to run locally.

User instructions(may change a little)
- will update with more directory info as I think that will be a little confusing

#### Dependencies
1. Python3
  - pandas
  - NumPy
  - Biopython

2. Porter5
  - hhsuite
  - psiblast
  - uniprot (cp -r /fdb/hhsuite/uniprot20_2016_02 .)
  - uniref90
  
3. R

## For running on conserved amino acid domains:

[coming soon]

## For running on predicted secondary structure:


### Step 0.
You'll need to set up Porter5 with a single sequence before you can run all in your fasta file. Porter5 requests the paths to its dependencies. I do this by pulling the first sequence in my fasta file, creating a new .fa with that, and submitting that to Porter. This lets me input my paths and set everything up. You only need to do this once.

> Example of how I configure the paths (this will be almost the same for anyone else running on biowulf):

>uniref90 = /fdb/SIFT/uniref90.fa

>hhblits = hhblits [I call the hhsuite module with the job submission, or /usr/local/apps/hhsuite/3.0-beta.3/bin/hhblits without doing that should work too]

>uniprot = /data/belfordak/Buck_lab/Divergent_prots/DivProt/uniprot20_2016_02/uniprot20_2016_02 (cp -r uniprot20_2016_02 db to a dir you have permission to - giving the path to where it's hosted on biowulf will cause an error when you run Porter b/c you don't have permission to access the db in this location)

>psiblast = /usr/local/apps/ncbi-toolkit/21.0.0/bin/psiblast

(If found that your prediction works and you get your .ssX files, you're good to go. If you get an error message, it is likely due to your paths, even though there was probably no error called when you set them up. Also, just fyi, the error messages from Porter don't always obviously direct you to the problem (an error telling you a certins psi-blast related file could not be found doesn't necessarily mean there is anything wrong with your psi-blast directory/module/path).
  
### Step 1. 
Run the pre_pre_processing script on your fasta input file. Both single and multiline fasta are fine.
Ex.
```
bash pre_pre_processing.sh Adoma_polyoma_LTandVP1.fasta
```

### Step 2. 
Submit a swarm job like the prompt at the end of the pre_pre_processing.sh script says to do with the outfile. Do this from within the Porter5 directory, linking to the location of the file in ~/split_out. 

### Step 3.
Run the pre_procesing.py script to take the Porter5 output and produce what I'm calling "fastqish" files.

They're basically fastq files (but they don't contain the line numbers I've included here):

>1. @Arowana_adomavirus_LT | 14901:17525 Reverse
>2. CCCHHHHHHHHHHHHHCCCCCCCCHHHHHHHHHHHHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHHCCCC
>3. "+"
>4. !"#$%'(*++++*'%$"!!!!"##$&'+----,,+++*))&#!!!!!!!#'(+-..--+(%$##$$##"!

You want to run this script on all of the files in the directory (all individual sequences from your input,as Porter5 only takes on at a time)
So in the terminal, run:
```
for f in *; do python3 pre_procesing.py $f ../original_fasta.fasta; done
```
Once this is done, you'll want to create on meraged file of all the XX.ssX.fastqish files. This file containes the sequence name, secondary structure prediction, and phred score of each position structure.
Create the file with something like:
```
cat *.ssX.fastqish > whatever_file_name_you_want (e.x. original_fasta_name.ss3.fastqish)
```

### Step 4
Run the aligner to produce a metrix of alignment scores of all your input structure sequences. The R file to produce figures will automatically run from the aligner_finalish.py script. They will be outputted into Rplots.pdf.

Run aligner:
```
python3 aligner_finalish.py original_fasta_name.ss3.fastqish
```

#### On the structure:
> XX.ss3: helix (H), strand (E), and coil (C)

> XX.ss8: helix (G), α-helix (H), π-helix (I), β-stand (E), bridge (B), turn (T), bend (S), and others (C)


...That's it for now...

Notes and TODO:
> 1. The aligner currently *does not* take into account phred scores when calculating the alignment score. This may change. If it doesn't, we should change alignment method anyway from custom -> blosum62
> 2. There currently is not a limit to fasta input size, except figures will only be able to scale up to a certain degree, and will get difficult to read after a certain number of sequences are used.
> 3. I'll be benchmarking Porter5 with more threads to see if it makes much of a difference in run time.
> 4. A 3rd method that takes into account both the aa and ss alignment scores for producing a matrix may be be implemented in the future.
> 5. Add .log file for ss aligner
> 6. Add --help
> 7. Add error messages
> 8. collect Porter output in a new dir (or write something about usr creating and moving into a new dir, then executing pre_processing. will need to change that outfile a little)

