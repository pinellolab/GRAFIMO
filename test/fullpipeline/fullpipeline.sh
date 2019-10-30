#!/bin/sh

########################################################################
#  prog: fullpipeline.sh
#  author: Manuel Tognon <manu.tognon@gmail.com>
#
#
# Test scripts that show the main functionalities of GRAFIMO.
# Available scripts:
#   - fullpipeline.sh
#   - myvgpipeline.sh
#
# In this script will be covered the full pipeline;
#
# We will create the genome graph of the chr22 (hg19 assembly)
# using as variants the VCF file from the 1000 Genome Project
#
# The regions we will extract are relative to the best 200 ChIP-Seq peaks
# for CTCF in chromosome 22. The experiment is from The Encode Project
# (ENCFF519CXF).
#
# The motif file used is MA0139.1.jaspar from JASPAR database.
# The file refers to CTCF motif (length of motif: 19)
#
# The background file refers to the background probability for nucleotides
# of assembly hg19 (For more details type in your console 'cat bg_nt')

chr22=`gunzip chr22.fa.gz; ls chr22.fa`

# get the VCF from The 1000 genome project
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz

# assign it to a variable
variants=`ls ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz`
bed=`ls 200peaks_ENCFF519CXF_chr22.bed`
CTCF=`ls MA0139.1.jaspar`
bg=`ls bg_nt`


# run the test
grafimo --linear_genome $chr22 --chroms 22 --vcf $variants --bedfile $bed --motif $CTCF --bgfile $bg --pseudo 0.1 --pvalueT 1e-3

# What we did?
#
# We are building the graph genome of the chromosome 22, using as variants to
# represent those present in the given VCF file; remember it is fundamental to
# index your variants (tabix needed); note that you can also index your VCF
# before giving the non-indexed VCF to grafimo (the indexing step will be
# skipped).
# We are also declaring that we want to extract the regions defined in the given
# BED file.
# The motif used to score is that defined in MA0139.1.jaspar (for now only JASPAR
# format is supported). To build the score for each position we are also giving the
# background distribution for our alphabet (in our case is the DNA alphabet "ACGT")
# and a pseudocount, which will be added to each entry of the motif PWM.
# For further details on the scoring formula used refer to
# xxxx. xxxxx. Bioinformatics. xxxx
#
# Finally we are declaring that in output we want only sequences that have a score
# with a p-value <= 0.001.
#
# Now you should have in your working directory a new directory named grafimo_out.
# Enter it and you should have three files:
#   - grafimo_out.tsv
#   - grafimo_out.html
#   - grafimo_out.gff
#
# These files are a summarize the results obtained by GRAFIMO pipeline.

echo "Content of the directory grafimo_out:"
tree grafimo_out/


