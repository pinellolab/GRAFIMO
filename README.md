# GRAFIMO
Graph-based Find Individual Motif Occurrences

![workflow](https://user-images.githubusercontent.com/51021763/66830153-1a0f9c00-ef55-11e9-9702-fa6ba0f2cf61.png)

## GRAFIMO installation and usage
**dependencies needed**
- vg (https://github.com/vgteam/vg)
- tabix (https://github.com/samtools/tabix)

**pip installation**

**docker image**

This option is suggested for MacOS and Windows users.

Note that docker must be installed in your machine:
- https://docs.docker.com/docker-for-windows/install/ (Windows)
- https://docs.docker.com/docker-for-mac/install/ (MacOS)

There are two available ways to get grafimo running on your machine using docker:
- build the image by your own and run it
- run the image from the Docker Hub

First option:
 clone GRAFIMO repository, enter it, build the docker image and run it following these commands:
 ```
 git clone https://github.com/InfOmics/GRAFIMO.git
 cd GRAFIMO
 docker build -t grafimo .
 docker run -it grafimo
 grafimo --help
 ```
 If the help is correctly printed, then the image has been succesfully built.
 
 Second option:
  

## Installation test
 If everything went right type:
 ```
 grafimo --help
 ```
 if the help is correctly printed, then the tool is installed.

 To test the functionalities of GRAFIMO download the directory test.
 Enter it, by typing:
 ```
 cd test
 ```

 Then type:
 ```
 grafimo_test
 ```
If at the end appears the message:
"All tests passed! Enjoy GRAFIMO for your research!"
All the tests were passed.

Remember to exit from the test directory with:
```
cd ..
```

## Usage

Here is a brief guide to help you using GRAFIMO.

We provide the user with two possible workflows:
- with the variation graph creation
- using the user variation graph

**with variation graph creation**

Input:
- linear reference genome in .fa format (e.g. hg38.fa)
- VCF file compressed (e.g. my_vcf.vcf.gz)
- BED file containing the regions to extract from the genome
- JASPAR motif file (available at http://jaspar.genereg.net/)
- [OPTIONS]

Output:
- TSV file containing the sequences statistically significant with their scores
- HTML file containing the sequences statistically significant with their scores
- GFF file containing the hits found by GRAFIMO

Example:

```
grafimo --linear_genome hg.fa --vcf vcf.vcf.gz --bedfile bed.bed --motif m.jaspar [OPTIONS]
```

**without variation graph creation**

Input:
- your variation graph both in VG and XG format (better if in XG)
- BED file containg the regions to extract from the genome
- JASPAR motif file (available at http://jaspar.genereg.net/)
- [OPTIONS]

Output:
- TSV file containing the sequences statistically significant with their scores
- HTML file containing the sequences statistically significant with their scores
- GFF file containing the hits found by GRAFIMO

Example:
```
grafimo --graph_genome my_vg.xg --bedfile my_bed.bed --motif my_motif.jaspar
```

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

What is GRAFIMO doing?

GRAFIMO build a variation graph from the reference and the VCF given in input. Then
it will use the indexed version of the graph, in order to speed-up the following queries.
The next step uses the BED barrowPeak file to query regions from the variation graph.
Finally, from those regions are obtained sequences of length L, where L is the length of 
the motif, with a sliding window approach, and each sequence is scored using the motif
given to GRAFIMO.
The scores are finally filtered by p-value: since we have set a p-value threshold of 1e-3,
all the sequences with a p-value <= 1e-3 will be kept in the result.

The output is a directory (named ```grafimo_out``` by default) that contains the results
obtained in three formats:
 - TSV
 - HTML
 - GFF
 
The latter format allows the user to load the results obtained with GRAFIMO to a genome
browser.


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

What is GRAFIMO doing?

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

The output is a directory (named ```grafimo_out``` by default) that contains the results
obtained in three formats:
 - TSV
 - HTML
 - GFF
 
The latter format allows the user to load the results obtained with GRAFIMO to a genome
browser.
