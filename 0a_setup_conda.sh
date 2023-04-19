#!/bin/bash


conda create -n cdh -y
conda activate cdh

conda install -c conda-forge r-base=4.1.0 -y
conda install -c conda-forge c-compiler -y
conda install -c conda-forge gcc -y
conda install -c anaconda gxx_linux-64 -y
conda install -c anaconda readline -y
conda install -c conda-forge r-biocmanager -y
conda install -c conda-forge python -y
conda install -c conda-forge fs -y
#conda install -c conda-forge r-reticulate -y
#conda install -c conda-forge igraph -y
#conda install -c conda-forge leidenalg -y
#conda install -c conda-forge vtraag -y
#conda install -c conda-forge pandas -y
#conda install -c conda-forge umap -y
#conda install -c conda-forge learn -y
conda install -c conda-forge r-sourcetools -y
conda install -c bioconda freebayes=1.3.4
conda install -c bioconda bcftools=1.12
conda install -c bioconda samtools

pip install vireoSNP