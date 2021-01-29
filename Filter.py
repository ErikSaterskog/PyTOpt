import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc



def Check(x,dc,H):
    
    #Settings
    Timers=True
    tic = time.perf_counter()
    

    #Initialize 
    nElem=np.size(x,0)
    new_dc=np.zeros([np.size(dc),1])
    
    #breakpoint()
    for elem in range(0,nElem):
        #breakpoint()
        new_dc[elem][0]=(np.sum(np.transpose(x)*H[elem,:]*np.transpose(dc)))/(x[elem]*np.sum(H[elem,:]))


    if Timers:
        toc = time.perf_counter()  
        print('Check Time: '+str(toc-tic))
    return new_dc
































