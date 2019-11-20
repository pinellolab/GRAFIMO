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
from . import paths_scoring as ps
from . import vgCreation as vgc
from . import objs_writer as ow

__version__='0.8'

def with_vg_pipeline(cores, linear_genome, vcf, chroms, bedfile, motif, bgfile, 
                         pseudo, pvalueT, no_reverse, qvalue, dest, pipeline):
    
    gplus=False # prevent unexpected behaviors

    printWelcomeMsg("WITH_VG_CREATION")
    
    ##TO DO:q-value computation

    vg_loc=vgc.create_vg(chroms, linear_genome, vcf) # create the vg
    m=mtf.get_motif_pwm(motif, bgfile, pseudo, no_reverse) # create the motif
    data=sge.get_data(vg_loc, bedfile, m.getWidth(), pipeline, gplus, 
                          chroms, cores) # extract the region peaks
    df=ps.scoreGraphsPaths(data, m, pvalueT, cores, no_reverse, qvalue) # scoring
    
    objs_towrite=[df] # initialize the list of objects to save
    ##TO DO: matplotlib or seaborn plots 
    ##TO DO: add plots objects to objs_toWrite() 
    ow.writeresults(objs_towrite, dest) # write results
    
def without_vg_pipeline(cores, graph_genome, bedfile, motif, bgfile, pseudo, 
                            pvalueT, no_reverse, qvalue, dest, pipeline, gplus=False, chroms=[]):
    
    printWelcomeMsg("WITHOUT_VG_CREATION")
    
    ##TO DO:q-value computation
    
    m=mtf.get_motif_pwm(motif, bgfile, pseudo, no_reverse)
    data=sge.get_data(graph_genome, bedfile, m.getWidth(), pipeline, gplus, 
                        chroms, cores)
    df=ps.scoreGraphsPaths(data, m, pvalueT, cores, no_reverse, qvalue)
    
    objs_towrite=[df] # initialize the list of objects to save
    ##TO DO: matplotlib or seaborn plots 
    ##TO DO: add plots objects to objs_toWrite() 
    ow.writeresults(objs_towrite, dest)

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
    for _ in range(35):
        print('*', end='')
    print()
    print("\tWELCOME TO GRAFIMO v", __version__, sep='')
    print()
    print("Beginning the " + pipeline + " pipeline")
    print()

    for _ in range(35):
        print('*', end='')
    print('\n')