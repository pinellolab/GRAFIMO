#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 12:07:48 2019

@author: manuel
"""

import motif as mtf
import sys

if __name__=='__main__':
    
    try:
        motif_file=sys.argv[1]
        bg=sys.argv[2]
        pseudo=float(sys.argv[3])
    
        m=mtf.build_motif_pwm_jaspar(motif_file, bg, pseudo)
   
    
        sm=mtf.scale_pwm(m.getMotif_matrix())
        m.setMotif_matrix_scaled(sm)
        
        m=mtf.comp_pval_mat(m)
        
        print(m.getMotifID())
        print(m.getMotifName())
        
        print()
        
        print(m.getMotif_matrix())
        print(m.getMotif_matrix_scaled())
        print(m.getMotif_pval_mat(), len(m.getMotif_pval_mat()), sum(m.getMotif_pval_mat()))
        
    except:
        
        sys.stderr.write('motif.py: test not passed\n')
        
    else:
        
        sys.stderr.write('motif.py: test passed\n')
    
    
    