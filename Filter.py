import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc



def Check(x,dc,H):
    
    #Settings
    Timers=True
    tic = time.perf_counter()
    
    new_dc=np.transpose((np.sum(np.transpose(x)*H[:,:]*np.transpose(dc), axis=1))/(np.transpose(x[:])*np.sum(H[:,:], axis=1)))

    

    if Timers:
        toc = time.perf_counter()  
        print('Check Time: '+str(toc-tic))
    return new_dc































