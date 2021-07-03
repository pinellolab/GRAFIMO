# How to reproduce results presented in "GRAFIMO: variant and haplotype aware motif scanning on pangenome graph"

The current directory contains the data and shell scripts to reproduce the analysis presented in the paper "GRAFIMO: variant and haplotype aware motif scanning on pangenome graph".

## Prerequisites

- Unix shell
- GRAFIMO and the required dependencies should be installed (for further information see https://github.com/pinellolab/GRAFIMO or https://github.com/InfOmics/GRAFIMO)
- Jupyter Lab or Jupyter Notebook (with Python3 interpreter)

## How to reproduce the analysis

**1.** The first step is to download the pangenome variation graphs (VGs) we used during the GRAFIMO analysis. We created the VGs by enriching the hg38 reference genome with genetic variants (SNPs and indels) of 2548 individuals from the 1000 Genomes Project (1000GP).

**NB** Our analysis were carried out with VG v1.31.0

To download the VGs, type
```
bash get_vgs.sh
```

The VGs will be stored in a newly created directory, called ```vg```.

**2.** Then, we can search the potential motif occurrences on the VG (it can take a while)
```
bash run_analysis.sh
```

For simplicity, we already provided the top 3000 ChIP-seq optimal IDR thresholded peak regions for each transcription factor binding site motif (```./tf_motifs/TFNAME/Cell_line_NAME/EXPERIMENTCODE_3000.bed).

The ChIP-seq optimal IDR thresholded peak regions were originally stored in BigBED file format. To obtain a suitable input for GRAFIMO (requires that genomic regions are given in BED format), we first need to convert them in the corresponding BED files.

Given that ```bigBedToBed``` (http://hgdownload.cse.ucsc.edu/admin/exe/) tool is available on your Unix PATH, to convert the BigBED to the corresponding BED, type
```
bigBedToBed /path/to/BIGBED.bigbed /path/to/BED.bed
```

In our analysis we considered only those features beloning to the canonical chromosomes (chr1, chr2, chrX, etc.). To filter the other features we provided a Python script
```
python3 process_bed.py /path/to/BED.bed /path/to/BED_filt.bed
```

To keep only the top 3000 regions (sorted by *q*-value):
```
cat /path/to/BED_filt.bed | sort -rk9 | head -n3000 > /path/to/BED_3000.bed
```

**3.** Once the analysis is completed, we can reproduce the results presented in the paper by running ```grafimo_results_analysis.ipynb``` IPython Notebook.

## Reproducing population specific binding site analysis

To reproduce the analysis on the fraction of population specific binding sites recovered only on individual haplotypes are required: 

- GRAFIMO v1.1.5 (still in alpha) available [here](https://github.com/pinellolab/GRAFIMO/tree/track-samples)

- VG v1.32.1 or later

Once the requirements have been satisfied, the search can start (it will take a while):
```
bash run_analysis_pop.sh
```

Once the search has been completed, the analysis can be reproduced by running the IPython notebook ```population_specific_bs.ipynb```.
