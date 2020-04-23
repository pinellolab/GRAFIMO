#!/bin/sh

########################################################################
#  prog: vg_scan.sh
#  author: Manuel Tognon <manu.tognon@gmail.com>
#                        <manuel.tognon@studenti.univr.it>
########################################################################
#
# Test scripts to show the main functionalities of GRAFIMO (GRAph-based Finding
# of Individual Motif Occurrences).
#
# Available scripts:
#   - vg_creation.sh
#   - vg_scan.sh
#
# In this script will be bruefly explained how to scan VG(s) using GRAFIMO.
#
# We will scan the genome variation graphs created during the previous tutorial
# (vg_creation.sh).
#
# First cpy all the VGs in the current directory (if you prefer you can keep them in
# their original location)
cp ../vg_creation/chr*.xg .

# In the current directory there are a file named MA0139.1.meme, which represents
# the PWM for CTCF, downloaded from JASPAR database (entry MA0139.1), a file named
# ENCFF633TEW_2000.bed.zip (you need to uncompress it before running GRAFIMO
# analysis), which represents 2000 genomic coordinates, where CTCF motif
# will be searched, and an another file named `bg_nt`, which represents
# a background distribution for the genome (to define your own background file,
# refer to http://meme-suite.org/doc/bfile-format.html).
#
# Note that genomic coordinates defined in ENCFF633TEW_2000.bed.gz are ChIP-seq
# narrow peaks regions, available on the ENCODE project database; these coordinates
# were mapped on hg38 genome assembly and found on cell line K562.

# uncompress the BED file
gunzip ENCFF633TEW_2000.bed.zip

# Scan the VGs
grafimo findmotif -d ./ -b ENCFF633TEW_2000.bed -m MA0139.1.meme -k bg_nt

# Note that is possible to constrain the scan to a subset of the chromosome, using
# the -c option. Also you can define the threshold for the sequences (-t option)
# and decide if apply it on the P-values or on the q-values (--qvalueT option).
# To avoid problems during log-odds matrix computation is added a pseudocount to
# PWM's entry (default 0.1). You can set your pseudocount value by using -p option.
#
# For example calling GRAFIMO on the same data used before, but using a threshold
# value of 1e-3 on the q-values, constraining the scan to chromosome 20 and 21, and
# using a pseudocount of 0.001:
#   grafimo findmotif -d ./ -b ENCFF633TEW_2000.bed.gz -m MA0139.1.meme -k bg_nt -c 20 21 -p 0.001 -t 1e-3 --qvalueT
#
# The output of grafimo is a directory, named grafimo_out_JOBID_MOTIFID (default),
# containing three files:
#   - TSV results summary
#   - HTML results summary
#   - GFF (suitable for input to UCSC genome browser)
#
# Note, also, that is possible to scan a whole genome variation graph (-g option),
# but this will result in a slower running time of GRAFIMO.
# For example:
#   grafimo findmotif -g /path/to/my/vg.vg -b ENCFF633TEW_2000.bed -m MA0139.1.meme -k bg_nt

