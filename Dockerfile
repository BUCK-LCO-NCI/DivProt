FROM ubuntu:latest
LABEL maintainer "anna.k.belford@gmail.com"

#please don't use this yet...I'm currently working on getting it running with everything we need!

#avoid prompts that come from
ARG DEBIAN_FRONTEND=noninteractive

# get the python libs and other software
RUN apt-get update && apt-get install -y \
  git \ 
  python3 \
  python3-pip \
  python3-setuptools \
  python3-numpy \
  python3-requests \
  python3-msgpack \
  python3-matplotlib \
  python3-pandas \
  python3-biopython \
  hhsuite \
  ncbi-blast+ \
  build-essential \
  r-base \
  r-cran-randomforest

#biotite seperately 
RUN pip3 install biotite

# get DivProt, Porter5
RUN git clone https://github.com/BUCK-LCO-NCI/DivProt

WORKDIR /Divprot
#to put porter inside DP
RUN git clone https://github.com/mircare/Porter5/ 

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


#docker exec = execute a command inside of the container
RUN mkdir '/DivProt/dbs'

#notes
#check out https://docs.docker.com/get-started/part4/ for info about deploying this container to swarm if you're not on a hpc with the framework already all set up
#maybe need git clone hhsuite run
