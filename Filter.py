import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc



def Check(x,dc,H):
    
    #Settings
    Timers=True
    tic = time.perf_counter()

    
    Hsum=H.sum(axis=1).tolist()
    
    new_dc=H.multiply(x).multiply(dc).sum(axis=0).transpose()/(Hsum*x)
    
    
    if Timers:
        toc = time.perf_counter()  
        print('Check Time: '+str(toc-tic))
    return new_dc






























