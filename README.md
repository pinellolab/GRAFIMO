# GRAFIMO
Graph-based Finding of Individual Motif Occurrences

![workflow](https://user-images.githubusercontent.com/51021763/66830153-1a0f9c00-ef55-11e9-9702-fa6ba0f2cf61.png)

## GRAFIMO installation
**dependencies needed**
- vg (v1.21.0 or later) (https://github.com/vgteam/vg)
- tabix (https://github.com/samtools/tabix)
- graphviz (https://www.graphviz.org/)

**Building and installation**

To build and install GRAFIMO follow the following steps:
```
git clone https://github.com/pinellolab/GRAFIMO.git
cd GRAFIMO
python3 setup.py install
grafimo --help
```

If the help is correctly shown, you have succesfully built GRAFIMO.

**Docker image**

This option is suggested for MacOS and Windows users.

Note that docker must be installed in your machine:
- https://docs.docker.com/docker-for-windows/install/ (Windows)
- https://docs.docker.com/docker-for-mac/install/ (MacOS)

There are two available ways to get grafimo running on your machine using docker:
- run the image from the Docker Hub
- build the image by your own and run it
 
First option:

 make the pull of the image from the docker hub:
 ```
 docker pull pinellolab/grafimo
 docker run -it pinellolab/grafimo
 grafimo --help
 ```
 If the help is correctly printed, then the pull has been succesfully completed and you should be able to run the image on      your machine.
 
Second option:
 
 clone GRAFIMO repository, enter it, build the docker image and run it by typing the following commands:
 ```
 git clone https://github.com/pinellolab/GRAFIMO.git
 cd GRAFIMO
 docker build -t grafimo .
 docker run -it grafimo
 grafimo --help
 ```
 If the help is correctly printed, then the image has been succesfully built.

## Usage

GRAFIMO performs motif scanning on the VG genome-graph structure.  
The given motifs are saerched in the regions specified in the input BED files.

GRAFIMO presents many additional features to the motif scanning, beside the score, it gives a p-value, a q-value and other important informations, like the chromosome where the hit was found or its starting and ending positions.

GRAFIMO allows the user to follow two possible pipelines. The first one creates the genome-graphs for each chromosome from the given genome and VCF file (note that it is possible to build also a single or a variable number of genome-graphs, it is not mandatory to build the whole genome) and then search the given motifs.
The second pipeline performs the motif scanning step on a given set of genome-graphs or on a whole genome-graph, searching the motif in the regions defined in the given BED file.

The results are stored in three different file formats: TSV (can be viewed using Microsoft Excel or MacOS Numbers), HTML (can be viewed by any browser) and GFF (data can be loaded on UCSC genome browser).

Here there is a guide to the input parameters for GRAFIMO. The first two tables show the parameters that are mandatory for the 
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
`-m` | `MOTIF` | Path to the motif file to use to score the possible binding sites for a certain transcription factor. Motif files in both the MEME and JASPAR format can be given in input | No default

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
`--cores` | `NCORES` | Number of cores to use running GRAFIMO | Takes all the available cores
`-c`| `1 2 X ...` | List of chromosomes for which the variation graph will be created (or on which the motif will be searched)| Is created the variation graph for all the chromosomes (the motif is searched on all the chromosomes)
`-k`| `BACKGROUND` | 	Path to the background file specifying the source of a 0-order background model to convert a probability matrix to a log-odds score matrix and to use in estimating the p-values of match scores. The background model normalizes for biased distribution of individual letters in the sequences. The guidelines to specify this file can be found at http://meme-suite.org/doc/bfile-format.html | Assume a uniform distribution for the background
`-p`| `PSEUDOCOUNT` | Pseudocount that will be added to each count of the motif, to avoid problems in the computation of the log-odds score matrix, due to possible divisions by zero | Is used 0.1
`-t` | `THRESHOLD` | Hits with a P-value (by default) or a q-value higher than the defined threshold won't be returned in the results | A threshold of 1e-4 is used
`-q`|  | Specify if the q-value will not be computed for the scored sequences | The q-value is computed
`-r`|  | Specify if only sequences from the forward strand have to be scanned | Sequences from both the forward and reverse strand are scanned 
`-f`|  | Print results directly on the terminal without storing them | The results are stored in a directory and not printed on the terminal
`-o`| `OUTDIR` | Specify the name of the directory where the results are stored | A directory named `grafimo_out_JOBID_MOTIFNAME` is created
`--qvalueT` |  | If set to true, the threshold is applied on the q-values | The threshold is applied on the P-values
`--top-graphs`| `GRAPHS_NUM` | Are created the PNG images of the first `GRAPHS_NUM` regions (in the results order)
`--verbose`|  | Prints a lot of additional informations about the current run 

**Output**

- grafimo_out.tsv, TSV file containing the scored sequences; that format can be easily viewed and managed with Microsoft Excel or MacOS Numbers
- grafimo_out.html, HTML version of the TSV file; this file can be viewed using the browser you like most
- grafimo_out.gff, a GFF3 file that allow the user to load and view the results on the UCSC genome browser

**Advanced examples**

If you want to try advanced examples to take a deeper look onto GRAFIMO functionalities
enter the ```test``` directory by typing:
```
cd GRAFIMO/test
```

Let's start with an example showing the entire pipeline with the creation of the variation graph and the motif scanning.
Enter the directory ```fullpipeline```.
Run the script ```fullpipeline.sh```, by typing on your terminal:
```
chmod +x fullpipeline.sh 
./fullpipeline
```

In the ```fullpipeline``` directory you'll find a file ```chr22.fa```, that is the fasta file for the chromosome 22 
(hg38 assembly), a narrowPeak BED file ```200peaks_ENCFF519CXF_chr22.bed```, in this file
there are the first 200 best peaks (ordered by q-value) for the transcription factor
CTCF from ChIP-seq data by The ENCODE Project (ENCFF519CXF.bed) related to chr22, a 
motif file ```MA0139.1.jaspar``` representing the CTCF motif from JASPAR database 
(http://jaspar.genereg.net/) and ```bg_nt```. The
latter file represents a background distribution for the alphabet.
 
By running ```fullpipeline.sh``` you'll also get a VCF file ```ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz.vcf.gz```
from The 1000 Genome Project, used to build the variation graph of chr22.

*What is GRAFIMO doing?*

GRAFIMO build a variation graph from the reference and the VCF given in input. Then
it will use the indexed version of the graph to query the regions defined in the BED file.
Finally, from these regions are obtained sequences of length L, where L is the length of 
the motif, with a sliding window approach, and each sequence is scored using the motif
given to GRAFIMO.
The scores are finally filtered by p-value: since we have set a p-value threshold of 1e-4,
all the sequences with a p-value <= 1e-4 will be kept in the results.

The output is a directory (named ```grafimo_out_JOBID_MA0139.1```) that contains the results
made of:
 - grafimo_out.tsv
 - grafimo_out.html
 - grafimo_out.gff


If you have your own whole genome-graph or a set of genome-graphs, in VG or XG format, you can skip
the variation graph creation step.

To test this pipeline go back to the ```test``` directory and enter ```myvgpipeline```
with:
```
cd myvgpipeline
```

First we copy the variation graph of chromomosome 22 that we created before by typing:
```
cp ../fullpipeline/chr22.xg .
```
Now the variation graph of the chromosome 22 is available in the current directory. 
Remember that also the VG (e.g. ```chr22.vg```) format is allowed, but it will be 
indexed to perform the following queries.
You'll also find, as before, a BED narrowPeak file ```200peaks_ENCFF519CXF_chr22.bed```, in this file
there are the first 200 best peaks for CTCF from ChIP-seq data by The ENCODE
Project (ENCFF519CXF.bed) related to chr22, a motif file ```MA0139.1.jaspar``` 
representing the CTCF motif from JASPAR database (http://jaspar.genereg.net/) and 
```bg_nt```. Again, the latter file represents a background distribution for the alphabet.

Now we are ready to try this pipeline. Type on the terminal:
```
chmod +x myvgpipeline.sh
./myvgpipeline
```

*What is GRAFIMO doing?*

With this pipeline we skip the genome-graph creation, instead we directly query it, to get the 
regions defined in the BED narrowPeak file.
In our example we are using the option ```--graph_genome_dir``` this means that 
we are giving to GRAFIMO the path to the directory which contains the 
genome-graphs we want to query. If you want to score only certain genome graphs 
(as in the example) you can use the option ```--chroms```, which tells GRAFIMO 
which chromosome you want to query. If you want to query all the chromosome just
don't use this argument.
Then, from the obtained subgraphs are extracted, as before, the sequences of 
length L, where L is the length of the motif. Everyone of these sequences are scored and 
filtered using a p-value threshold, like in the previous example.

The output is a directory (named ```grafimo_out``` by default) that contains the results
made of:
 - grafimo_out.tsv
 - grafimo_out.html
 - grafimo_out.gff
