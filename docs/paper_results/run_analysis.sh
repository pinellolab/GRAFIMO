#!/usr/bin/bash

#---------------------------------------------------------------------
#
# This shell script reproduce the analysis presented in 
#   "GRAFIMO: variant and haplotype aware motif scanning via genome
#   variation graphs"
#
# Are saerched the potential occurrences of three transcription factor
# binding sites (TFBS) motif:
#   - CTCF (JASPAR ID MA0139.1)
#   - ATF3 (JASPAR ID MA0605.2)
#   - GATA1 (JASPAR ID MA0035.4)
#
# The three motifs are searched on a genome variation graph (hg38 
# genome assembly) constructed with genomic variants (SNPs and indels) 
# of 2548 individuals from the 1000 Genomes Project phase 3.
# 
# The TFBS motifs are searched on the top 3000 genomic regions 
# corresponding to ChIP-seq peak regions for each motif repectively. 
# The peak regions have been obtained from the ENCODE Project data 
# portal.  
# The ChIP-seq peaks were obtained from different cell lines. 
#
# In the analysis were kept only those potential motif occurrences
# whose P-value was < 1e-4.
#
# For further details refer to the main paper and supplementary 
# materials
#---------------------------------------------------------------------

cd tf_motifs
myhome=${PWD}

## run the analysis for CTCF
cd CTCF
ctcf_home=${PWD}

# search on data from Cell line A549
cd Cell_line_A549

grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all -t 1

cd $ctcf_home

# search on data from Cell line GM12878
cd Cell_line_GM12878

grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all -t 1

cd $ctcf_home

# search on data from Cell line HepG2
cd Cell_line_HepG2

grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all -t 1

cd $ctcf_home

# search on data from Cell line K562
cd Cell_line_K562

grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all -t 1

cd $ctcf_home

# search on data from Cell line MCF-7
cd Cell_line_MCF-7

grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all -t 1

cd $ctcf_home
cd $myhome

## run the analysis for ATF3
cd ATF3
atf3_home=${PWD}

# search on data from Cell line H1
cd Cell_line_H1

grafimo findmotif -d ../../../vg -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all -t 1

cd $atf3_home

# search on data from Cell line HepG2
cd Cell_line_HepG2

grafimo findmotif -d ../../../vg -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all -t 1

cd $atf3_home

# search on data from Cell line K562
cd Cell_line_K562

grafimo findmotif -d ../../../vg -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all -t 1

#cd $atf3_home
cd $myhome

## run the analysis for GATA1
cd GATA1
gata1_home=${PWD}

# search on data from Cell line K562
cd Cell_line_K562

grafimo findmotif -d ../../../vg -m ../MA0035.4.meme -b ENCFF811YFQ_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_ENCFF811YFQ -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0035.4.meme -b ENCFF811YFQ_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all_ENCFF811YFQ -t 1

grafimo findmotif -d ../../../vg -m ../MA0035.4.meme -b ENCFF939ODZ_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_ENCFF939ODZ -t 1e-4
grafimo findmotif -d ../../../vg -m ../MA0035.4.meme -b ENCFF939ODZ_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_all_ENCFF939ODZ -t 1

cd $myhome
