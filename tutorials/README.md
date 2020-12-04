# GRAFIMO tutorials

In this directory are contained two hands-on tutorials on **how to search motif occurrences in VGs** and on **how to build a genome variation graph (VG) with GRAFIMO**.

## Searching for motif occurrences in genome variation graphs

The directory `findmotif_tutorial` contains a Shell script with an hands-on tutorial on how to search potential motif occurrences in
VGs with GRAFIMO. 

The content of `findmotif_tutorial` is:
- **`findmotif_tutorial.sh`**: Shell script describing how to search motif occurrences on genome variation graphs. In the script are shown 
a simple motif search with GRAFIMO and some advanced options. Given that GRAFIMO, VG, Tabix and GraphViz are installed in your system and are 
reachable via Unix `$PATH`, you can run the script to inspect the output.
- **`data` directory**: directory containing the data used during the tutorial steps

Please note that data used for the tutorial are simple toy examples.

## Building genome variation graphs with GRAFIMO

The directory `buildvg_tutorial` contains a Shell script with an hands-on tutorial on how to build genome variation graphs with GRAFIMO,
by using a reference genome FASTA file and a VCF file.

The content of `buildvg_tutorial` is:
- **`buildvg_tutorial.sh`**: Shell script describing how to build genome variation graphs with GRAFIMO. In the script are shown 
a simple motif search with GRAFIMO and some advanced options. Given that GRAFIMO, VG, Tabix and GraphViz are installed in your 
system and are reachable via Unix `$PATH`, you can run the script to inspect the output.
- **`data` directory**: directory containing data used during the tutorial steps

Please note that data used for the tutorial are simple toy examples.
