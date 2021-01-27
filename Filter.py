import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc



def Check(eDof,coords,rmin,x,dc):
    
    #Settings
    Timers=True
    tic = time.perf_counter()
    
    #Initialize 
    nDof=np.max(eDof)
    nElem=np.size(eDof,0)
    ex=coord[0]
    ey=coord[1]
    
    new_dc = np.zeros([nElem,1])
    
    #Find centers of elements
    
    
    
    
    
    
    np.ix_(eDof[elem,:],eDof[elem,:])
    
    
    
    
    #Loop through all elements
    
    #Find elements within rmin of current element
    
    
    for i in range(0,nelx):
        for j in range(0,nely):
            _sum = 0
            
            
            for k in range(int(max([i-np.floor(rmin),0])),int(min([i+np.floor(rmin)+1,nelx]))):
                for l in range(int(max([j-np.floor(rmin),0])),int(min([j+np.floor(rmin)+1,nely]))):
                    
                    
                    
                    
                    fac = rmin-np.sqrt((i-k)**2+(j-l)**2)
                    _sum = _sum + max(fac,0)
                    new_dc[j,i] = new_dc[j,i] + max(fac,0)*x[l,k]*dc[l,k]
                    
                    
            new_dc[j,i] = new_dc[j,i]/(x[j,i]*_sum)
            
    if Timers:
        toc = time.perf_counter()  
        print('Check Time:')
        print(toc-tic)        
    return new_dc



































