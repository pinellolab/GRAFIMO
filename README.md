# GRAFIMO
Graph-based Find Individual Motif Occurrences

![workflow](https://user-images.githubusercontent.com/51021763/66830153-1a0f9c00-ef55-11e9-9702-fa6ba0f2cf61.png)

## GRAFIMO installation and usage
**dependencies needed**
- vg (https://github.com/vgteam/vg)
- tabix (https://github.com/samtools/tabix)
- graphviz (https://www.graphviz.org/)

**pip installation**

To build GRAFIMO on your machine, follow these steps:
```
git clone https://github.com/pinellolab/GRAFIMO.git
cd GRAFIMO
python3 setup.py install
grafimo --help
```

If the help is shown, you have succesfully built GRAFIMO

**docker image**

This option is suggested for MacOS and Windows users.

Note that docker must be installed in your machine:
- https://docs.docker.com/docker-for-windows/install/ (Windows)
- https://docs.docker.com/docker-for-mac/install/ (MacOS)

There are two available ways to get grafimo running on your machine using docker:
- build the image by your own and run it
- run the image from the Docker Hub

First option:
 
 clone GRAFIMO repository, enter it, build the docker image and run it by typing the following commands:
 ```
 git clone https://github.com/pinellolab/GRAFIMO.git
 cd GRAFIMO
 docker build -t grafimo .
 docker run -it grafimo
 grafimo --help
 ```
 If the help is correctly printed, then the image has been succesfully built.
 
Second option:

 make the pull of the image from the docker hub:
 ```
 docker pull pinellolab/grafimo
 docker run -it pinellolab/grafimo
 grafimo --help
 ```
 If the help is correctly printed, then the pull has been succesfully completed and you should be able to run the image on      your machine.

## Usage

Here is a guide to the input parameters for GRAFIMO. The first two tables show the parameters that are mandatory for the 
corresponding pipelines.
The last table shows the optional parameters that can be given to GRAFIMO.
Then a brief explanation of the output will be provided.

Finally are given some advanced examples to try GRAFIMO.

**Input**

*With variation graph creation*

Options                   | Parameter      | Description                  | Default behavior 
--------------------------- | -------------- | ---------------------------- | ------------
`-l` | `LINEAR-GENOME` | Path to the linear genome (the FA or FASTA formats are required) | No default 
`-v` | `VCF`| Path to the VCF file that contains the variants the user wants to consider in the variation graph. The VCF must be compressed (e.g. VCF.vcf.gz) | No default
`-b` | `BEDFILE` | Path to the bedfile containing the ChIP-seq narrow peaks | No default
`-m` | `MOTIF` | Path to the motif file to use to score the possible binding sites for a certain transcription factor. Motif files in both the JASPAR and MEME format can be given in input | No default

*Without variation graph creation*

Options                   | Parameter      | Description                  | Default behavior 
--------------------------- | -------------- | ---------------------------- | ------------
`-g` | `GRAPH-GENOME` | Path to a genome graph that is already available. It must be in XG or VG format. The former is preferred over the latter (**NB** this is an alternative to the `-d` option) | No default 
`-d` | `GRAPH-GENOME-DIR`| Path to the directory that contains different genome graphs from which the ChIP-seq peaks will be extracted. Here we are assuming that a graph genome has been constructed for each chromosome (**NB** this is an alternative to `-g`option) | No default
`-b` | `BEDFILE` | Path to the bedfile containing the ChIP-seq narrow peaks | No default
`-m` | `MOTIF` | Path to the motif file to use to score the possible binding sites for a certain transcription factor. Motif files in both the JASPAR and MEME format can be given in input | No default

*Options*

Options                   | Parameter      | Description                  | Default behavior 
--------------------------- | -------------- | ---------------------------- | ------------
`-h`|  | Print the help and exit | No default
`--version`|  | Print GRAFIMO's version and exit | No deafault
`--cores` | `NCORES` | Defines the number of cores of your machine to use running GRAFIMO | Takes all the available cores
`-c`| `1 2 X ...` | List of chromosomes for which the variation graph will be created | Is created the variation graph for all the chromosomes
`-k`| `BACKGROUND` | 	Path to the background file specifying the source of a 0-order background model to convert a probability matrix to a log-odds score matrix and to use in estimating the p-values of match scores. The background model normalizes for biased distribution of individual letters in the sequences. If both strands are being scored, the background model is modified by averaging the frequencies of letters and their reverse complements. The guidelines to specify this file can be found at http://meme-suite.org/doc/bfile-format.html | Assume a uniform distribution for the background
`-p`| `PSEUDOCOUNT` | Pseudocount that will be added to each count of the motif, to avoid problems in the computation of the log odds score matrix, due to possible divisions by zero | Is used 0.1
`-t` | `THRESHOLD` | Hits with a p-value higher than the defined threshold won't be inserted in the final result summary | A threshold of 1e-4 is used
`-q`|  | Specify if the q-value will not be computed for the scored sequences | The q-value is computed
`-r`|  | Specify if only sequences form the forward strand have to be considered | Sequences from both the forward and reverse strand will be scored
`-o`| `OUTDIR` | Specify the name of the directory where the results will be stored | A directory named `grafimo_out` will be created

**Output**

- grafimo_out.tsv, TSV file containing the scored sequences; that format can be easily viewed and managed with Excel
- grafimo_out.html, HTML version of the TSV file
- grafimo_out.gff, a GFF3 file that allow the user to see the results in the UCSC genome browser

**Advanced examples**

If you want to try advanced examples to take a look into GRAFIMO functionlities
enter the ```test``` directory by typing:
```
cd test
```

Let's start with an example showing the entire pipeline with the creation of the variation graph.
Enter the directory ```fullpipeline```.
Run the script ```fullpipeline.sh```, by typing on your terminal:
```
chmod +x fullpipeline.sh
./fullpipeline
```

In the directory you'll find a file ```chr22.fa```, that is the fasta file for the chromosome 22 
(hg38 assembly), a narrowPeak BED file ```200peaks_ENCFF519CXF_chr22.bed```, in this file
there are the first 200 best peaks (ordered by q-value) for the transcription factor
CTCF from ChIP-seq data by The ENCODE Project (ENCFF519CXF.bed) related to chr22, a 
motif file ```MA0139.1.jaspar``` representing the CTCF motif from JASPAR database 
(http://jaspar.genereg.net/) and ```bg_nt```. The
latter file represents a background distribution for the alphabet.
 
By running ```fullpipeline.sh``` you'll also get a VCF file ```ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz.vcf.gz```
from the 1000 Genome Project, used to build the variation graph of chr22.

*What is GRAFIMO doing?*

GRAFIMO build a variation graph from the reference and the VCF given in input. Then
it will use the indexed version of the graph, in order to speed-up the following queries.
The next step uses the BED barrowPeak file to query regions from the variation graph.
Finally, from those regions are obtained sequences of length L, where L is the length of 
the motif, with a sliding window approach, and each sequence is scored using the motif
given to GRAFIMO.
The scores are finally filtered by p-value: since we have set a p-value threshold of 1e-3,
all the sequences with a p-value <= 1e-3 will be kept in the result.

The output is a directory (named ```grafimo_out``` by default) that contains the results
made of:
 - grafimo_out.tsv
 - grafimo_out.html
 - grafimo_out.gff


If you have your own variation graph available in VG or XG format you can skip
the variation graph creation step.

To test this pipeline go back to the ```test``` directory and enter ```myvgpipeline```
with:
```
cd myvgpipeline
```
In the current directory there is a toy example of the variation graph of the chromosome 22,
named ```chr22.xg```. Remember that also the VG format is allowed, but its
processing will take more time since we need to index it for the following queries;
you'll also find a BED narrowPeak file ```200peaks_ENCFF519CXF_chr22.bed```, in this file
there are the first 200 best peaks for CTCF from ChIP-seq data by The ENCODE
Project (ENCFF519CXF.bed) related to chr22, a motif file ```MA0139.1.jaspar``` 
representing the CTCF motif from JASPAR database (http://jaspar.genereg.net/) and 
```bg_nt```. The atter file represents a background distribution for the alphabet.

Now we are ready to try this pipeline. Type on the terminal:
```
chmod +x myvgpipeline.sh
./myvgpipeline
```

*What is GRAFIMO doing?*

With this pipeline we skip the variation graph creation, instead we query the 
regions defined in the BED narrowPeak file on the given variation graph.
In our example we are using the option ```--graph_genome_dir``` this means that 
we are giving to GRAFIMO the path to the directory which contains the 
graph genomes we want use. If you want to score only certain genome graphs 
(as in the example) is better to use the option ```--chroms``` that specifies 
which chromosome you want to query. If you want to query all the chromosome just
don't add this argument.
Then, from the obtained subgraphs are extracted as done before the sequences of 
length L, where L is the length of the motif. Everyone of these sequences are scored and 
filtered using a p-value thrshold, like in the other pipeline.

The output is a directory (named ```grafimo_out``` by default) hat contains the results
made of:
 - grafimo_out.tsv
 - grafimo_out.html
 - grafimo_out.gff
