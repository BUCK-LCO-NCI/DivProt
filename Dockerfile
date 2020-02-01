FROM ubuntu:latest
LABEL maintainer "anna.k.belford@gmail.com"

#avoid prompts that come from
ARG DEBIAN_FRONTEND=noninteractive

# get the python libs and other software
RUN apt-get update && apt-get install -y \
  git \ 
  python3 \
  python3-setuptools \
  python3-pip \
  python3-numpy \
  python3-requests \
  python3-msgpack \
  python3-matplotlib \
  python3-pandas \
  python3-biopython \
  python-pexpect \
  hhsuite \
  ncbi-blast+ \
  build-essential \
  r-base \
  r-cran-randomforest

#pip seperately 
RUN pip3 install biotite
RUN pip3 install wget

# get DivProt, Porter5
RUN git clone https://github.com/BUCK-LCO-NCI/DivProt

#set DivProt as working directory from now on
WORKDIR /DivProt

#grab porter5
RUN git clone https://github.com/mircare/Porter5/ ./Porter5 

#clean
RUN apt-get autoremove -y && rm -rf /DivProt/readme_figures /DivProt/example_DP_out

#get the R packages
RUN R -e \
  "install.packages('Biostrings',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  "install.packages('pheatmap',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  "install.packages('igraph',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  "install.packages('reshape',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  "install.packages('Matrix',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  "install.packages('Rtsne',dependencies=TRUE, repos='http://cran.rstudio.com/')"

#setting a pathi in the environment --> I do not think this is necessary for hhblits to work for porter, but we'll see
#ENV HHLIB=/hh-suite
#ENV PATH="$HHLIB/bin:$HHLIB/scripts:${PATH}"

#make directory for databases to be deposited in in someone requests uniref and uniprot to be installed to the container (for porter)
RUN mkdir ./dbs

#notes
#check out https://docs.docker.com/get-started/part4/ for info about deploying this container to swarm if you're not on a hpc with the framework already all set up
#i maybe need git clone hhsuite run (?)
