#!/bin/sh

########################################################################
#  prog: myvgpipeline.sh
#  author: Manuel Tognon <manu.tognon@gmail.com>
#                        <manuel.tognon@studenti.univr.it>
########################################################################
#
# Test scripts that show the main functionalities of GRAFIMO.
# Available scripts:
#   - fullpipeline.sh
#   - myvgpipeline.sh
#
# In this script will be covered the pipeline without the genome graph
# creation. In this pipline users can give to grafimo both a set of genome
# graphs in which the motif will be searched and the whole genome graph
# (both XG and VG format are allowed). Note that is suggested to use the first
# approach, rather than giving the whole genome graph in input, because the
# search will be much more demanding and could cause memory issues on the less
# powerful systems.
#
# We will download the chr22 (hg38 assembly) genome graph, created using hg38 as
# reference and using as variants a VCF file from the 1000 Genome Project (for
# further details read fullpipeline.sh).
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

# assign mandatory elements to variables
chr22=`tar xvzf chr22.xg.zip; ls chr22.xg`
bed=`ls 200peaks_ENCFF519CXF_chr22.bed`
CTCF=`ls MA0139.1.jaspar`
bg=`ls bg_nt`

# run the example
grafimo -d ./ --chroms 22 --bedfile $bed --motif $CTCF --bgfile $bg --pseudo 0.1 --pvalueT 1e-4

# What are we doing?
#
# We are giving to grafimo a directory containing the chr22 genome graph (the
# current directory './' ), previously downloaded. In alternative we can give in
# input the path to a whole genome graph. Then we are extracting from chr22.xg the
# regions defined in the BED file 200peaks_ENCFF519CXF_chr22.bed (collection of
# 200 regions for CTCF ChIP-Seq peaks ordered by FDR, assembly hg38).
#
# We are also declaring that the hits have to be scored using data available in the
# motif MA0139.1.jaspar (both JASPAR and MEME file format are supported by GRAFIMO)
# from JASPAR database.
# With the --bgfile option we are giving a file that tells GRAFIMO which
# is the background distribution to use in the motif preprocessing step; moreover
# the --pseudo option declares that to each entry of the motif PWM will be added
# the specified pseudocount, in order to avoid division by zero during the motif
# preprocessing.
#
# Finally we are telling to GRAFIMO to store only hits which have a p-value
# associated to their score which is < 0.0001.
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

