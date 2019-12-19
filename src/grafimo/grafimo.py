"""

@author: Manuel Tognon

@email: manu.tognon@gmail.com
@email: manuel.tognon@studenti.univr.it


Script where are defined the software version working and the two different
pipelines available:
    
    - with vg creation
    - only scoring (indexing of the vg if required [vg => xg])

"""

from . import subgraphs_extraction as sge
from . import motif as mtf
from . import motif_set as mtf_set
from . import paths_scoring as ps
from . import vgCreation as vgc
from . import objs_writer as ow
from grafimo.utils import isMEME_ff
import sys
import os

__version__ = '0.8'

def with_vg_pipeline(cores, linear_genome, vcf, chroms, bedfile, motifs, bgfile, 
                         pseudo, pvalueT, no_reverse, qvalue, text_only, dest, top_graphs, 
                         pipeline, verbose = False):
    
    gplus = False # prevent unexpected behaviors

    printWelcomeMsg("WITH_VG_CREATION")

    if verbose:
        print("\nCreating the variation graph for chromosomes" + str(chroms))
    vg_loc = vgc.create_vg(chroms, linear_genome, vcf) # create the vg

    motifSet = mtf_set.MotifSet()
    for motif in motifs:
        if verbose:
            print("\nBuilding the motif " + motif)
        
        m = mtf.get_motif_pwm(motif, bgfile, pseudo, no_reverse) # create the motif
        motifSet.addMotif(m)

    motif_lst = motifSet.getMotifsList()
    for m in motif_lst:
        data = sge.get_data(vg_loc, bedfile, m.getWidth(), pipeline, gplus, 
                                chroms, cores, verbose) # extract the region peaks
        df = ps.scoreGraphsPaths(data, m, pvalueT, cores, no_reverse, qvalue) # scoring

        if text_only:
            # print directly on the terminal the results in a TSV-like manner
            print(df)

        else:
            objs_towrite=[df] # initialize the list of objects to save
            ow.writeresults(objs_towrite, dest, m.getMotifID(), top_graphs, vg_loc) # write results
    
def without_vg_pipeline(cores, graph_genome, bedfile, motifs, bgfile, pseudo, 
                            pvalueT, no_reverse, qvalue, text_only, dest, top_graphs,
                            pipeline, gplus=False, chroms=[], verbose = False):
    
    printWelcomeMsg("WITHOUT_VG_CREATION")
    
    motifSet = mtf_set.MotifSet()
    for motif in motifs:
        if verbose:
            print('\nBuilding the motif ' + motif)
        
        m = mtf.get_motif_pwm(motif, bgfile, pseudo, no_reverse) # create the motif
        motifSet.addMotif(m) # add the motif

    motif_lst = motifSet.getMotifsList()
    for m in motif_lst:
        data = sge.get_data(graph_genome, bedfile, m.getWidth(), pipeline, gplus, 
                                chroms, cores) # extract the region peaks
        df = ps.scoreGraphsPaths(data, m, pvalueT, cores, no_reverse, qvalue)

        if text_only:
            # print directly on the terminal the results in a TSV-like manner
            print(df)

        else:
            objs_towrite=[df] # initialize the list of objects to save
            ow.writeresults(objs_towrite, dest, m.getMotifID(), top_graphs, graph_genome)

def printWelcomeMsg(pipeline):
    """
        Prints the initial message for GRAFIMO
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