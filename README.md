# GRAFIMO
GRAph-based Finding of Individual Motif Occurrences

![image](https://user-images.githubusercontent.com/51021763/80205317-b1535b80-862a-11ea-98bb-7e9680da7881.png)

## Introduction

Transcription factors (TFs) are proteins that promote or reduce the expression of genes by binding nearby short DNA sequences known as transcription factor binding sites (TFBS). While several tools have been developed to report potential occurrences of TFBS in linear DNA sequences, no tool exists to find them in genome variation graphs (VGs). VGs are sequence-labeled graphs that can represent many genomes in a single compact data structure. Scanning for binding sites in VGs allows us to capture how genomic variation affects the binding landscape of a factor in a population of individuals.

GRAFIMO is a command-line tool for the scanning of known TF DNA motifs represented as Position Weight Matrices (PWMs) in VGs. GRAFIMO implements features commonly available in motif scanning tools for linear genomes. Essentially, it extends the standard PWM scanning procedure to consider variation and alternative haplotypes encoded in a VG. Using GRAFIMO, we recover several binding sites that are missed when using only the reference genome, and which could constitute individual-specific binding events.


## GRAFIMO installation
**dependencies needed**
- vg (v1.21.0 or later) (https://github.com/vgteam/vg). Note that vg must be in your PATH.
- tabix (https://github.com/samtools/tabix)
- graphviz (https://www.graphviz.org/)

**Building and installation**

To build and install GRAFIMO follow the following steps:
```
git clone https://github.com/pinellolab/GRAFIMO.git
cd GRAFIMO
python3 setup.py install --user
grafimo --help
```

If the help is correctly shown, you have successfully built GRAFIMO.

We successfully built GRAFIMO also using a Windows Subsystem for Linux.

**Docker image**

Note that docker must be installed in your machine:
- https://docs.docker.com/docker-for-windows/install/ (Windows)
- https://docs.docker.com/docker-for-mac/install/ (MacOS)

There are two available ways to get GRAFIMO running on your machine using docker:
- run the image from the Docker Hub
- build the image by your own and run it
 
First option:

 make the pull of the image from the docker hub:
 ```
 docker pull pinellolab/grafimo
 docker run -it pinellolab/grafimo grafimo --help
 ```
If the help is correctly printed, then the pull has been successfully completed and you should be able to run the image on your machine.
 
Second option:
 
clone GRAFIMO repository, enter it, build the docker image and run it by typing the following commands:
 ```
 git clone https://github.com/pinellolab/GRAFIMO.git
 cd GRAFIMO
 docker build -t grafimo .
 docker run -it grafimo grafimo --help
 ```
 If the help is correctly printed, then the image has been successfully built.

## Usage

GRAFIMO performs motif scanning on the VG genome variation-graph data-structure.  
The given motifs are searched in the genomic coordinates specified in the input BED file.

GRAFIMO scans a precomputed genome variation graph or a set of pre-computed genome graphs, e.g. one for each chromosome (suggested approach), for the occurrences of the given motif. It also allows the user to compute the VGs for each chromosome from scratch, given the reference FASTA file and a set of genomic variants in a VCF file.

The results are stored in three different file formats: TSV (can be viewed using Microsoft Excel or MacOS Numbers), HTML (can be viewed by any browser) and GFF (data can be loaded on UCSC genome browser).

With GRAFIMO is also possible to compute the VG from your data (given a reference genome in FASTA format and a VCF file containing genomic variants). Will be created a set of genome variation graphs (one for each chromosome). This approach is equivalent to create a single whole genome variation graph, but allows a faster GRAFIMO analysis, in terms of running time. 

**Input**

The following three tables summarize the available input arguments for GRAFIMO.
The first one shows the arguments to build the genome variation graph from your data. 
The second table shows the arguments to scan a precomputed VG or a set of genome variation-graphs.

*GRAFIMO building of genome variation graphs*

To build your genome variation graph with grafimo, as a set of VGs (one for each chromosome), type ```grafimo buildvg```, followed by the following arguments:

Options                   | Parameter      | Description                  | Default behavior 
--------------------------- | -------------- | ---------------------------- | ------------
`-l` | `LINEAR-GENOME` | Path to the linear genome (the FA or FASTA formats are required) [mandatory] | No default 
`-v` | `VCF`| Path to the VCF file that contains the variants the user wants to consider in the variation graph. The VCF must be compressed (e.g. VCF.vcf.gz) [mandatory] | No default
`-c`| `1 2 X ...` | List of chromosomes for which the variation graph will be created | Is created the variation graph for all the chromosomes
`--cores` | `NCORES` | Number of cores to use during VG construction | Takes all the available cores

*GRAFIMO scanning of pre-computed VG(s)*

To scan a pre-computed VG (set of VGs) type ```grafimo findmotif```, followed by the following arguments:

Options                   | Parameter      | Description                  | Default behavior 
--------------------------- | -------------- | ---------------------------- | ------------
`-g` | `GRAPH-GENOME` | Path to a pre-computed genome variation graph. It must be in XG or VG format. The former is preferred over the latter (**NB** this argumennt is an alternative to the `-d` option) [mandatory] | No default 
`-d` | `GRAPH-GENOME-DIR`| Path to the directory that contains a set of VGs, on which the motif will be searched. Is assumed that has been constructed a VG for each chromosome (**NB** this is an alternative to `-g`option) [mandatory] | No default
`-b` | `BEDFILE` | Path to the BED file containing a set genomic coordinates where the motif will be searched [mandatory] | No default
`-m` | `MOTIF` | Path to the motif file motif to search on VG(s). Both MEME and JASPAR formats are allowed [mandatory] | No default
`-h`|  | Print the help and exit | No default
`--version`|  | Print GRAFIMO current version and exit | No default
`--cores` | `NCORES` | Number of cores to use during GRAFIMO analysis | Takes all the available cores
`-c`| `1 2 X ...` | List of chromosomes where the motif will be searched | The motif is searched on all the chromosomes
`-k`| `BACKGROUND` | Path to a background file, which specifies the source of a 0-order background model to convert a probability matrix to a log-odds score matrix and to use in P-values estimation for corresponding scores. The guidelines to specify this file can be found at http://meme-suite.org/doc/bfile-format.html | Assume a uniform distribution for the background
`-p`| `PSEUDOCOUNT` | Pseudocount that will be added to each motif count (probability), to avoid problems during log-odds matrix computation due to possible divisions by zero | 0.1
`-t` | `THRESHOLD` | Hits with a P-value (by default) or a q-value higher than the defined threshold won't be returned in the results | 1e-4
`-q`|  | If used the q-values will not be computed | q-values are computed
`-r`|  | Only sequences belonging to the forward strand will be returned in the results | Are returned sequences belonging to both forward and reverse strands 
`-f`|  | Print results on terminal, without creating the three results files | The results are summarized in three files stored in a newly created directory
`-o`| `OUTDIR` | Name of the directory where the results will be stored | Is created a directory named `grafimo_out_JOBID_MOTIFNAME`
`--qvalueT` |  | If used the threshold is apllied on q-values rather than on P-values | The threshold is applied on P-values
`--top-graphs`| `GRAPHS_NUM` | For the first `GRAPHS_NUM` regions (ordered by P-value) containing an hit, will be returned their PNG image
`--verbose`|  | Prints additional informations about GRAFIMO execution 

**Output**

GRAFIMO outputs the results of its analysis in three different file formats: TSV, HTML and GFF3. The three files are stored in the specified directory (`-o` option) or, by default, in a directory named `grafimo_out_JOBID_MOTIFNAME`.

The three files will be named `grafimo_out.*`:
- grafimo_out.tsv, TSV file containing the scored sequences; this file format can be easily viewed and managed with Microsoft Excel or MacOS Numbers
- grafimo_out.html, HTML version of the TSV file; this file can be viewed using the browser you like most
- grafimo_out.gff, a GFF3 file that allow the user to load the results on the UCSC genome browser as custom tracks

## Advanced examples

To take a deeper look into GRAFIMO functionalities enter the `test` directory by typing:
```
cd GRAFIMO/test
```

Let's start with an example showing how GRAFIMO creates a set of genome variation-graphs (one for each chromosome) from your data.

Enter the directory `vg_creation`, by typing
```
cd vg_creation
```

First, we need a reference genome and a set of variants to create a genome variation-graph for each chromosome. Then, we need to download the genome FASTA file from the UCSC server by typing
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz  # uncompress the gz
```
and a VCF file containing SNPs and indel from the 1000 Genomes
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz
```

To build your first genome variation graph (VG for each chromosome approach) with GRAFIMO, type
```
grafimo buildvg -l hg38.fa \
-v ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz 
```
GRAFIMO builds a variation graph for each chromosome enriched with the variants defined in the VCF given in input. If your analysis should be limited to a subset of chromosomes, just use the ```-c``` option. By default will be used all the cores available on your machine. To specify a custom number of cores, use the ```--cores```.
All the previously described steps are summarized in the script ```vg_creation.sh```, available in `vg_creation` folder.


If you have a precomputed VG genome variation-graph or a set of VGs (e.g. one for each chromosome), you can follow the following steps.

To test this functionality of GRAFIMO go back to the ```test``` directory and enter ```vg_scan``` directory by typing 
```
cd ../vg_scan
```

First we copy the variation graphs obtained in the first part of the tutorial, in the current directory by typing
```
cp ../vg_creation/chr*.xg .
```
Now, the variation graphs should be available in the current directory. 
In ```vg_scan``` directory you'll find a file named `MA0139.1.meme`, which represents the PWM for CTCF, downloaded from JASPAR database (entry MA0139.1) a file named `ENCFF633TEW_2000.bed.zip` (you need to uncompress it before running GRAFIMO analysis), which represents 2000 genomic coordinates where CTCF motif will be searched, and an another file named `bg_nt`, which represents a background distribution for the genome (to define your own background file, refer to http://meme-suite.org/doc/bfile-format.html). Note that genomic coordinates defined in `ENCFF633TEW_2000.bed.zip` are ChIP-seq narrow peaks regions, available on the ENCODE project database; these coordinates were mapped on hg38 genome assembly and found on cell line K562.

To run GRAFIMO analysis type in your terminal
```
grafimo findmotif -d ./ -b ENCFF633TEW_2000.bed -m MA0139.1.meme -k bg_nt
```

GRAFIMO scans the variation-graphs in the specified directory (`-d` option) for CTCF motif occurrences in the regions defined in the given BED file. From these regions are obtained sequences of length 19, where 19 is the length of CTCF motif, with a sliding window approach, and each sequence is scored using the motif PWM given to GRAFIMO. Note that the search can be constrained to a subset of chromosomes, by using the `-c` option.
The resulting scores are then filtered by P-value < 1e-4 (you can define your threshold and if it works on P-values or q-values). 

The output is a directory (named ```grafimo_out_JOBID_MA0139.1```) , containing the results, summarized in the files:
 - grafimo_out.tsv
 - grafimo_out.html
 - grafimo_out.gff

All the previously described steps are summarized in the script ```vg_scan.sh```, available in `vg_scan` folder.

