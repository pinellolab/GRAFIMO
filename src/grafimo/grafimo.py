"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it


Here are defined the steps of the two available pipelines:
    
    - with vg creation
    - only scoring (indexing of the vg if required [vg => xg])

The final results can be visualized directly on the terminal if 
text_only == True (only in TSV format).
By default the results are written to the repository given by the user.

"""

from . import subgraphs_extraction as sge
from . import motif as mtf
from . import motif_set as mtf_set
from . import paths_scoring as ps
from . import vgCreation as vgc
from . import objs_writer as ow
import time

# version of GRAFIMO
__version__ = '0.10'


# pipeline that creates the graphs of the chromosomes beside the motif discovery
def with_vg_pipeline(cores, linear_genome, vcf, chroms, bedfile, motifs, bgfile, 
                         pseudo, pvalueT, no_reverse, qvalue, text_only, dest, top_graphs, 
                         pipeline, verbose = False):
    
    gplus = False # prevent unexpected behaviors

    printWelcomeMsg("WITH_VG_CREATION")

    ### create the genome graphs ###
    if verbose:
        print("\nCreating the variation graph for chromosomes" + str(chroms))
    start_vgc = time.time()
    
    vg_loc = vgc.create_vg(chroms, linear_genome, vcf) # create the vg

    end_vgc = time.time()
    if verbose:
        print("Elapsed time for the creation of the genome graphs:", str(end_vgc - start_vgc))

    motifSet = mtf_set.MotifSet()
    for motif in motifs:
        if verbose:
            print("\nBuilding the motif " + motif)
        start_mc = time.time()
        
        m = mtf.get_motif_pwm(motif, bgfile, pseudo, no_reverse) # create the motif
        motifSet.addMotif(m)

        end_mc = time.time()
        if verbose:
            print("Elapsed time for the processing of motif " + motif + ":", str(end_mc - start_mc))

    # from the motif set get the list of the processed motifs
    motif_lst = motifSet.getMotifsList()

    # perform mpotif scanning for each processed motif
    for m in motif_lst:
        
        start_sge = time.time()

        data = sge.get_data(vg_loc, bedfile, m.getWidth(), pipeline, gplus, 
                                chroms, cores, verbose) # extract the region peaks

        end_sge = time.time()
        if verbose:
            print("Elapsed time to extract regions from " + bedfile + "for motif " + motif + 
                    ":", str(end_sge - start_sge))

        start_scan = time.time()
        
        df = ps.scoreGraphsPaths(data, m, pvalueT, cores, no_reverse, qvalue) # motif scanning

        end_scan = time.time()
        if verbose:
            print("Elapsed time to extract regions from " + bedfile + "for motif " + motif + 
                    ":", str(end_scan - start_scan))

        if text_only:
            # print directly on the terminal the results in a TSV-like manner
            print(df)

        else:
            objs_towrite=[df] # initialize the list of objects to save
            ow.writeresults(objs_towrite, dest, m.getMotifID(), top_graphs, vg_loc) # write results
    

# pipeline that performs the motif discovery on the given genome graph(s)
def without_vg_pipeline(cores, graph_genome, bedfile, motifs, bgfile, pseudo, 
                            pvalueT, no_reverse, qvalue, text_only, dest, top_graphs,
                            pipeline, gplus=False, chroms=[], verbose = False):
    
    printWelcomeMsg("WITHOUT_VG_CREATION")
    
    motifSet = mtf_set.MotifSet()
    for motif in motifs:
        
        if verbose:
            print('\nBuilding the motif ' + motif)

        start_mc = time.time()
        
        m = mtf.get_motif_pwm(motif, bgfile, pseudo, no_reverse) # create the motif
        motifSet.addMotif(m) # add the motif

        end_mc = time.time()
        if verbose:
            print("Elapsed time for the processing of motif " + motif + ":", str(end_mc - start_mc))

    motif_lst = motifSet.getMotifsList()
    for m in motif_lst:

        start_sge = time.time()

        data = sge.get_data(graph_genome, bedfile, m.getWidth(), pipeline, gplus, 
                                chroms, cores) # extract the region peaks
        
        end_sge = time.time()
        if verbose:
            print("Elapsed time to extract regions from " + bedfile + "for motif " + motif + 
                    ":", str(end_sge - start_sge))
        
        start_scan = time.time()

        df = ps.scoreGraphsPaths(data, m, pvalueT, cores, no_reverse, qvalue) # motif scanning
        
        end_scan = time.time()
        if verbose:
            print("Elapsed time to extract regions from " + bedfile + "for motif " + motif + 
                    ":", str(end_scan - start_scan))

        if text_only:
            # print directly on the terminal the results in a TSV-like manner
            print(df)

        else:
            objs_towrite=[df] # initialize the list of objects to save
            ow.writeresults(objs_towrite, dest, m.getMotifID(), top_graphs, graph_genome)


def printWelcomeMsg(pipeline):
    """
        Prints the intro message for GRAFIMO
        ----
        Params:
            pipeline (str) : chosen pipeline
        ----
        Returns:
             None
    """
    for _ in range(50):
        print('*', end='')
    print()
    print("\tWELCOME TO GRAFIMO v", __version__, sep='')
    print()
    print("Beginning the " + pipeline + " pipeline")
    print()

    for _ in range(50):
        print('*', end='')
    print('\n')

