#!/bin/bash

#-----------------------------------------------------------------------
#
# This shell script shows how to search for motif occurrences on genome 
# variation graphs (VGs) by using GRAFIMO.
#
# This example searches a toy motif called 'example.meme' (see 
# data/example.png for the motif logo) on a genome variation graph. The 
# considered genome has two chromosomes 'x' and 'y'. The corresponding 
# VGs (XG and GBWT indexes are in 'data/mygenome/' directory). 
#
# The motif will be searched on the regions defined in the
# 'data/regions.bed' BED file. The BED file contains 86 regions where 
# the motif will be searched.
#
# Note that all the functionalities showed here apply also for motif 
# search on single whole genome variation graphs. You need to use `-g` 
# option instead of `-d`, followed by the path to your graph.
#
#----------------------------------------------------------

# let us assign the motif and the BED file to variables
motif="data/example.meme"
bedfile="data/regions.bed"

# Now we can search the motif on our genome variation graph.
#
# Note that with `-d` option we specify the directory containing
# the VGs of chromosomes (XG and GBWT indexes).
#
# In the VGs are encoded two haplotypes.
grafimo findmotif -d data/mygenome/ -m $motif -b $bedfile

# The results are stored in the directory 'grafimo_out_PID_MOTIFID'.
# 
# The directory contains a TSV report, an HTML report and a GFF3 report. 
# The GFF3 can be loaded on UCSC genome browser asnd visualized as custom 
# track.

# GRAFIMO computes a log-likelihood score, a P-value and a q-value
# for each retrieved potential motif occurrence. Such measures are 
# weighted by a background probability distribution. By default, GRAFIMO 
# assumes a uniform background distribution. The user can provide
# a custom background distribution by passing it to GRAFIMO as a text 
# file like 'bg_nt' file in 'data/' directory. Let us assume we want 
# that the background distribution is the one described in 'bg_nt' file 
# and we want store the results in 'grafimo_out_bg'. We use `-k` and `-o` 
# options when calling GRAFIMO from command-line, in order to do this.
grafimo findmotif -d data/mygenome -m $motif -b $bedfile -k data/bg_nt -o grafimo_out_bg

# By default GRAFIMO retrieves only those potential motif occurrences
# whose P-value is smaller than 1e-4. Let us assume we want to
# retrieve those motif occurrences with a P-value smaller than 0.05
# rather than 1e-4 and store the results in 'grafimo_out_05'. This 
# can be done using `-t` and `-o` options.
grafimo findmotif -d data/mygenome -m $motif -b $bedfile -t 0.05 -o grafimo_out_05

# By default, GRAFIMO apply the threshold on the P-values, but it is
# also possible to apply it on q-values. Let us assume we want to
# retrieve only those potential motif occurrences with q-value smaller 
# than 1e-4 and we want store the results in 'grafimo_out_qval'. This 
# can be done using `--qvalueT`, `-t` and `-o` options, when calling 
# GRAFIMO from command-line.
grafimo findmotif -d data/mygenome -m $motif -b $bedfile --qvalueT -t 1e-4 -o grafimo_out_qval

# GRAFIMO returns only those potential motif occurrences which are 
# observed in at least one of the haplotypes embedded in the VGs. 
# GRAFIMO allows the user to study the recombinats that can be obtained 
# from the genomic variants used for building the VG. To do so, it is 
# sufficient to add the `--recomb` option when calling GRAFIMO from 
# command-line. Let us assume we are interested in studying all the 
# possible recombinants that can be obtained from the VCF variants, and 
# we want store the results in 'grafimo_out_recomb'.
grafimo findmotif -d data/mygenome -m $motif -b $bedfile --recomb -o grafimo_out_recomb

# Let us now assume we are interested in retrieving only those motif
# occurrences which can be found on chromosome 'x' and we want store 
# the results in 'grafimo_out_chrx' directory. This can be done by using 
# `--chroms-find` and `-o` options
grafimo findmotif -d data/mygenome -m $motif -b $bedfile --chroms-find x -o grafimo_out_chrx


# See https://github.com/pinellolab/GRAFIMO/wiki for additional details

