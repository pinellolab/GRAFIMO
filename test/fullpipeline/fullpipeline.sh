#!/bin/sh

########################################################################
#  prog: fullpipeline.sh
#  author: Manuel Tognon <manu.tognon@gmail.com>
#                        <manuel.tognon@studenti.univr.it>
########################################################################
#
# Test scripts that show the main functionalities of GRAFIMO.
# Available scripts:
#   - fullpipeline.sh
#   - myvgpipeline.sh
#
# In this script will be covered the full pipeline workflow;
#
# We will create the genome graph of the chr22 (hg38 assembly)
# using as variants the corresponding VCF file from the 1000 Genome Project.
#
# The regions we will extract are relative to the best 200 ChIP-Seq peaks
# for CTCF in chromosome 22 (with respect to the associated FDR). The experiment is
# from The Encode Project (ENCFF519CXF).
#
# The motif file used is MA0139.1.jaspar from JASPAR database.
# The file refers to CTCF motif (motif length: 19)
#
# The background file refers to the background probability for nucleotides
# of assembly hg38 (For more details type in your terminal 'cat bg_nt')

chr22=`gunzip chr22.fa.gz; ls chr22.fa`

# get the VCF from The 1000 genome project
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz

# assign mandatory elements to variables
variants=`ls ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz`
bed=`ls 200peaks_ENCFF519CXF_chr22.bed`
CTCF=`ls MA0139.1.jaspar`
bg=`ls bg_nt`


# run the test
grafimo -l $chr22 -c 22 -v $variants --bedfile $bed --motif $CTCF --bgfile $bg --pseudo 0.1 --pvalueT 1e-4

# What are we doing?
#
# We are building the graph genome of the chromosome 22, using as variants those
# present in the downloaded VCF file; remember it is fundamental to
# index your variants (tabix needed); note that if an indexed VCF file is found it
# will be reindexed, in order to avoid possible issues due to different versions of
# the TBI file.
# Using the '-c 22' option we are telling to GRAFIMO to take into account only the
# chromosome 22 in its computations. We can also specify more than one chromosome,
# in this case the genome fasta must be given in input as linear genome. To take
# into account all the chromosomes just don't use the '-c'option.
# We are also declaring that we want to extract the regions defined in the given
# BED file.
# The motif used to score is that defined in MA0139.1.jaspar (both JASPAR and MEME
# format are supported). To build the score matrix for the motif we are also giving
# a background distribution for our alphabet (in our case is the DNA alphabet
# "ACGT") and we are adding a pseudocount, in order to avoid division by zero, while
# processing our motif.
# For further details on the scoring approach used refer to
# xxxx. xxxxx. XXXXX. xxxx
#
# Finally we are declaring that in output we want only sequences that have a score
# with a p-value < 0.0001.
#
# Now you should have in your working directory a new directory named grafimo_out.
# Enter it and you should have three files:
#   - grafimo_out.tsv
#   - grafimo_out.html
#   - grafimo_out.gff
#
# These files are a summarize the results obtained by GRAFIMO pipeline.

echo "Content of the directory grafimo_out:"
tree grafimo_out_*/


