#!/bin/bash

#-----------------------------------------------------------------------
#
# This shell script shows how to build a genome variation graph (VG) 
# with GRAFIMO, starting from a reference genome and a VCF file.
#
# During this tutorial you will create a VG using reference FASTA file
# data/xy.fa and the VCF data/xy2.vcf.gz.
#
# Note that these files are simple toy examples.
#
# xy.fa contains the sequences of two toy chromosomes named 'x' and 'y'.
#
# The VCF xy2.vcf.gz contains genetic variants (SNPs and indels) 
# affecting both 'x' and 'y' sequences. For each chromosome are encode 
# two haplotypes.
#
#-----------------------------------------------------------------------

# assign the reference genome FASTA and the VCF to variables
reference="data/xy.fa"
vcf="data/xy2.vcf.gz"

# GRAFIMO will create a genome variation graph for each chromosome 
# available in the reference FASTA file, in our example 'x' and 'y'
#
# The two resulting VGs will be indexed, and the XG index and the GBWT 
# index are obtained.
#
# The latter index allows to recover haplotypes when scanning the VG.

# Now, we can call GRAFIMO to build our genome VG 
grafimo buildvg -l $reference -v $vcf

# As you can see, in the current directory four new files 'x.xg', 
# 'x.gbwt', 'y.xg' and 'y.gbwt' are created.
#
# These are the VGs of chromosomes 'x' and 'y', enriched with genomic 
# variants contained in the given VCF.
#
# We can use `vg view` command to inspect the structure of the VGs.
vg view -dp x.xg - | dot -Tpng -o x.png 
vg view -dp y.xg - | dot -Tpng -o y.png 

# Now, let us assume we are interested in creating the VG
# only for chromosome 'x' from our reference FASTA
#
# Moreover, let us assume we want to store the resulting 
# VG in the 'chrx_vg' directory
#
# This can be easily by using `--chroms-build`  and `-o` options 
# when calling GRAFIMO from command-line
mkdir -p chrx_vg  # create the output directory
grafimo buildvg -l $reference -v $vcf --chroms-build x -o chrx_vg

# You can notice that in the 'chrx_vg', only the VG
# for chromosome 'x' has been created.

# Let us assume we are constructing the VGs for a group of control subjects.
# To identify such VGs we'd like to add the prefix 'control' to VGs names.
#
# This can be easily done by using the '--chroms-prefix-build' option:
mkdir -p control_vg  # suppose we store the resulting VGs in control_vg directory
grafimo buildvg -l data/xy.fa -v data/xy2.vcf.gz -o control_vg --chroms-prefix-build control

# But if the user prefer he can define a space separated file to rename
# each chromosome VG as he/she wish. On the first column should be the 
# chromosomes names as used in reference FASTA, while on the second 
# column there should be the new names to assign to the VGs.
# Refer to data/namemap.txt for an example.
#
# Thus let's rename 'x' VG to awesomeVG and 'y' to fancyVG:
mkdir -p namemap_example  # suppose we store the results in namemap_example
grafimo buildvg -l data/xy.fa -v data/xy2.vcf.gz -o namemap_example --chroms-namemap-build data/namemap.txt

# Now, let us assume we would like to build the genome VG (a variation 
# graph for each chromosome) but by using a fresh VCF index (TBI).
#
# This can be done by adding `--reindex` option to GRAFIMO's command-line 
# call. Let us store the result in a new directory called 'vg_fresh_index'.
mkdir -p vg_fresh_index  # create the output directory
grafimo buildvg -l $reference -v $vcf -o vg_fresh_index --reindex

# As you can notice, in 'vg_fresh_index' the VGs for both chromosomes
# 'x' and 'y' (XG and GBWT indexes) are contained built using a fresh 
# VCF index. 
#
# Their structure is identical to the original VGs.

# See https://github.com/pinellolab/GRAFIMO/wiki for additional details

