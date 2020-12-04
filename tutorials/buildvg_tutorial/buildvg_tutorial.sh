#!/bin/bash

#----------------------------------------------------------
#
# This shell script shows how to build a genome variation
# graph (VG) with GRAFIMO, starting from a reference genome 
# and a VCF file
#
# This example creates a reference VG file by using as the FASTA
# file data/xy.fa and the VCF file xy2.vcf.gz. These files are 
# toy examples.
#
# data/xy.fa contains the sequences of two toy chromosomes, 
# named 'x' and 'y'. 
#
# The VCF file xy2.vcf.gz contains variants (SNPs and 
# indels) affecting both 'x' and 'y' chromosomes. 
# Two haplotypes are available for each chromosome.
#
#----------------------------------------------------------

# assign the reference genome FASTA and the VCF file to
# variables
reference="data/xy.fa"
vcf="data/xy2.vcf.gz"



# GRAFIMO will create a genome variation graph for each
# chromosome available in the reference FASTA file, in
# our example 'x' and 'y'
#
# The two resulting VGs will be indexed, and the XG
# index and the GBWT index are obtained.
#
# The latter index allows to retrieve the haplotypes when
# scanning the VG

# Now, we can call GRAFIMO to build our genome VG 
grafimo buildvg -l $reference -v $vcf

# As you can see, in the current directory four new files 
# 'chrx.xg', 'chrx.gbwt', 'chry.xg' and 'chry.gbwt' are 
# created
#
# These are the VGs of chromosomes 'x' and 'y', enriched with 
# the genomic variants contained in the given VCF file
#
# We can use `vg view` command to inspect the structure 
# of the VGs
vg view -dp chrx.xg - | dot -Tpng -o chrx.png 
vg view -dp chry.xg - | dot -Tpng -o chry.png 



# Now, let us assume we are interested in creating the VG
# only for chromosome 'x' from our reference FASTA
#
# Moreover, let us assume we want to store the resulting 
# VG in the 'chrx_vg' directory
#
# This can be easily by using `-c`  and `-o` options 
# when calling GRAFIMO from command-line
mkdir -p chrx_vg  # create the output directory

grafimo buildvg -l $reference -v $vcf -c x -o chrx_vg

# You can notice that in the 'chrx_vg', only the VG
# for chromosome 'x' has been created.



# Now, let us assume we would like to build the 
# genome VG (a variation graph for each chromosome)
# but by using a fresher VCF index (TBI) file.
#
# This can be done by adding `--reindex` option to 
# GRAFIMO's command-line call. Let us store the result
# in a new directory called 'vg_fresh_index'
mkdir -p vg_fresh_index  # create the output directory

grafimo buildvg -l $reference -v $vcf -o vg_fresh_index --reindex

# As you can see, in 'vg_fresh_index' the VGs for both chromosomes
# 'x' and 'y' (XG and GBWT indexes) are contained built
# using a fresher VCF index. 
#
# Their structure is identical to the original VGs



# See https://github.com/InfOmics/GRAFIMO/wiki for additional details

