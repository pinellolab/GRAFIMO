#!/bin/bash

########################################################################
#
# run_analysis_pop.sh
#
# recover potential motif occurrences and samples where they can be found

cd tf_motifs
myhome=${PWD}

cd CTCF
ctcf_home=${PWD}

cd Cell_line_A549
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop -t 1e-4 --track-samples -j 16

cd $ctcf_home

cd Cell_line_GM12878
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop -t 1e-4 --track-samples -j 16

cd $ctcf_home

cd Cell_line_HepG2
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop -t 1e-4 --track-samples -j 16

cd $ctcf_home

cd Cell_line_K562
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop -t 1e-4 --track-samples -j 16

cd $ctcf_home

cd Cell_line_MCF-7
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg -m ../MA0139.1.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop -t 1e-4 --track-samples -j 16

cd $ctcf_home
cd $myhome

cd ATF3
atf3_home=${PWD}

cd Cell_line_H1
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg  -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop -t 1e-4 --track-samples -j 16

cd $atf3_home

cd Cell_line_HepG2
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg  -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop -t 1e-4 --track-samples -j 16

cd $atf3_home

cd Cell_line_K562
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg  -m ../MA0605.2.meme -b *_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop -t 1e-4 --track-samples -j 16

cd $atf3_home
cd $myhome

cd GATA1
gata1_home=${PWD}

cd Cell_line_K562
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg  -m ../MA0035.4.meme -b ENCFF811YFQ_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop_ENCFF811YFQ -t 1e-4 --track-samples -j 16
grafimo findmotif -d /media/data/users/mtognon/Desktop/xg  -m ../MA0035.4.meme -b ENCFF939ODZ_3000.bed -k ../../bg_nt --chroms-prefix-find chr -o grafimo_out_pop_ENCFF939ODZ -t 1e-4 --track-samples -j 16

cd $gata1_home
cd $myhome
