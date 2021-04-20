#!/usr/bin/bash

#---------------------------------------------------------------------
#
# This Shell script downloads the genome variation graphs (VGs) used 
# to test GRAFIMO for the paper 
#   "GRAFIMO: variant and haplotype aware motif scanning via genome
#   variation graphs"
#
# The genome variation graphs have been constructed on the hg38 genome 
# assembly with ~78 millions of genomic variants (SNPs and indels) 
# of 2548 individuals from 1000 Genomes Project phase 3. 
#
#---------------------------------------------------------------------

home_dir=${PWD}

# create a directory containing the VGs to download
mkdir -p vg
cd vg

# download the whole genome variation graph (a VG per chromosome)
# 
# for each VG will be retrieved the XG index and the GBWT index
# the GBWT index allows the tracking of subjects' haplotypes, while
# scanning the VG for motif occurrences
for i in $(seq 1 22; echo "X");
do
    wget http://ncrnadb.scienze.univr.it/vgs/chr${i}.{xg,gbwt}
done

