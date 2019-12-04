#!/bin/sh

########################################################################
#  prog: myvgpipeline.sh
#  author: Manuel Tognon <manu.tognon@gmail.com>
#
#
# Test scripts that show the main functionalities of GRAFIMO.
# Available scripts:
#   - fullpipeline.sh
#   - myvgpipeline.sh
#
# In this script will be covered the pipeline where user can give to grafimo
# its own genome graphs (VG or XG format);
#
# We give you the chr22 (hg19 assembly) genome graph, created using hg19 as
# reference and using as variants a VCF file from the 1000 Genome Project.
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

chr22=`tar xvzf chr22.xg.zip; ls chr22.xg`
bed=`ls 200peaks_ENCFF519CXF_chr22.bed`
CTCF=`ls MA0139.1.jaspar`
bg=`ls bg_nt`

# run the example
grafimo -d ./ --chroms 22 --bedfile $bed --motif $CTCF --bgfile $bg --pseudo 0.1 

# What did we do?
#
# We gave to grafimo a directory containing a genome graph (the current directory
# './' )that we had available. In alternative we colud have given in input a
# path to a whole genome graph. Then we extract from chr22.xg the
# regions defined in the BED file 200peaks_ENCFF519CXF_chr22.bed (collection of
# 200 regions for CTCF ChIP-Seq peaks, assembly hg19).
#
# We also declared that the hits have to be scored using data available in the
# motif MA0139.1.jaspar (for now only JASPAR file format is supported by GRAFIMO)
# from JASPAR database.
# With the --bgfile option we are giving a file that tells GRAFIMO which
# is the background distribution to use in the scoring step; moreover the --pseudo
# option declares that to each entry of the motif PWM will be added the specified
# pseudocount
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
