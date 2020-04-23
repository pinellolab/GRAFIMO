#!/bin/sh

########################################################################
#  prog: vg_creation.sh
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
# In this script will be briefly explained how to create your VG from your data.
#
# We will create the genome variation graph (VG) for the whole human genome,
# creating a VG for each chromosome, using the genome reference (hg38 assembly) and
# using simple variants (SNPs and indel), from all the subject in the 1000 Genomes
# project.
#
# Instead of creating a single whole genome variation graphs, is created a VG per
# chromosome. This allows a execution for GRAFIMO, and allows the user to run VG
# also on its own laptop.
#
# Using this approach there is no drawback, since the results are equal to the ones
# obtained using a whole genome variation graph.

# first step:
#   download the reference genome from the UCSC server (hg38 assembly)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz  # uncompress the gz

# second step:
#   download the variants (indels + SNPs) from the 1000 GP
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz

# third step:
#   build the genome variation graph (VG)
grafimo buildvg -l hg38.fa -v ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz

# Note that if the user wants to build the VG only for a subset of chromosome
# it is possible to define this subset by using the -c option, in the call to
# GRAFIMO. For example if we want to build the VG only for chromosome 2, 18 and X,
# we should type:
#   grafimo -l hg38.fa -v ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz -c 2 18 X
#

