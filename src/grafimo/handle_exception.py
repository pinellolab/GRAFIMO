"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it

Script that defines funxctions to handle exceptions that can arise from
an uncorrect usage of the tool.

"""

import sys

def throw_arg_errors():
    
    sys.stderr.write('ERROR: -- incorrect argument given. Type fimovg --help for help\n')
    return 1

def throw_args_num_error():
    
    sys.stderr.write('ERROR: -- incorrect number of arguments given. Type fimovg --help for help\n')
    return 1

def throw_wget_genome_failure():
    
    sys.stderr.write('ERROR: -- unable to download the reference genome from UCSC. ')
    sys.stderr.write('Try to download it manually at ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz\n')
    return 1
    
def decompress_genome_error():
    
    sys.stderr.write("ERROR: -- unable to decompress genome. Make sure you have command 'gunzip'\n")
    return 1

def throw_wget_vcf_failure():
    
    sys.stderr.write('ERROR: -- unable to download the VCF from The 1000 genome project. ')
    sys.stderr.write('Try to download it manually at http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz\n')
    return 1

def throw_empty_motif_error():
    
    sys.stderr.write('ERROR: -- the motif matrix is empty\n')
    return 1

def throw_motif_matrix_not_pd_dataframe():
    
    sys.stderr.write('ERROR: -- the motif matrix is not an instance of pandas Dataframe\n')
    return 1

def throw_incorrect_width_error():
    
    sys.stderr.write('ERROR: -- invalid value for motif width assigned\n')
    return 1

def throw_not_bool_error():
    
    sys.stderr.write('ERROR: -- the value used is not an instance of bool\n')
    return 1
    
def throw_not_str_error():
    
    sys.stderr.write('ERROR: -- the value used is not an instance of str\n')
    return 1 

def throw_empty_TF_name_error():
    
    sys.stderr.write('ERROR: -- no name given for the transcription factor motif\n')
    return 1

def throw_not_jaspar_error():
    
    sys.stderr.write('ERROR: -- the motif file given is not in jaspar format. You can find them at http://jaspar.genereg.net\n')
    return 1

def throw_motif_file_not_read_error():
    
    sys.stderr.write('ERROR: -- not possible to read the given motif file\n')
    return 1

def throw_rev_matrix_not_requested_error():
    
    sys.stderr.write('ERROR: -- attempt to initialize a reverse matrix even if not requested\n')
    return 1

def throw_not_Motif_instance():
    
    sys.stderr.write('ERROR: -- the given value is not an instance of Motif\n')
    return 1

def throw_scaled_matrix_not_found():
    
    sys.stderr.write("ERROR: -- the Motif object doesn't have a scaled matrix associated\n")
    return 1

def throw_no_reverse_error():
    
    sys.stderr.write('ERROR: -- not allowed attempt to compute the reverse complement motif matrix\n')
    return 1

def throw_not_dna_alphabet_error():
    
    sys.stderr.write('ERROR: -- the alphabet is not recognized as the DNA alphabet [ A C G T]\n')
    return 1

def throw_not_able_to_read_bg_file():
    
    sys.stderr.write('ERROR: -- not able to read the backgroud file\n')
    return 1

def throw_subprocess_error():
    
    sys.stderr.write('ERROR: -- error while executing a shell command\n')
    return 1

def throw_bedfile_reading_error():
    
    sys.stderr.write('ERROR: -- not possible to read the bedfile\n')
    return 1

def throw_not_xg_error():
    
    sys.stderr.write('ERROR: -- is required the graph genome in xg format\n')
    return 1

def throw_wrong_pvalueT_err():
    
    sys.stderr.write('ERROR: -- pvalue must be a float value\n ')
    return 1

def throw_keys_diff_err():
    
    sys.stderr.write('ERROR: -- the dictionary keys are different\n')
    return 1

def throw_not_list_err():
    
    sys.stderr.write('ERROR: -- the object given is not an instance of list\n')
    return 1

def throw_objs_list_length():
    
    sys.stderr.write('ERROR: -- the length of the object list is null\n')
    return 1

def throw_no_xg_and_vg_err():
    
    sys.stderr.write('ERROR: -- the given file is neither a vg or an xg\n')
    return 1

def throw_not_wrong_input_err():
    
    sys.stderr.write('ERROR: -- incorrect input, please check it\n')
    return 1


    

    

    
    
    

