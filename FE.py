import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc



def FE(x,SIMP_penal,eDof,coord,fixDofs,F):

    #Settings
    E=210*1e9
    v=0.3
    ptype=2         #ptype=1 => plane stress, ptype=2 => plane strain
    ep=[ptype,1]    #ep[ptype, thickness]  
    Timers=True     #Print Timers 
    
    #Check sizes
    nDof=np.max(eDof)
    nElem=np.size(eDof,0)
    ex=coord[0]
    ey=coord[1]
    
    
    #Initialize Vecotors and Matrices
    K = np.zeros([nDof,])
    #F = np.zeros([nDof,1])
    U = np.zeros([nDof,1])
    
    
    #Check element type
    if len(eDof[0,:])==6:   #Triangular Element
        Tri=True
    elif len(eDof[0,:])==8:
        Tri=False           #Use Quad Instead
    else:
        raise Exception('Unrecognized Element Shape, Check eDof Matrix')
    
    
    #Start timer
    tic1 = time.perf_counter()
    
    #Sparse Matrices?
    #K=csr_matrix(K)
    #F=csr_matrix(F)
    
    #Linear Elastic Constitutive Matrix
    D=cfc.hooke(ptype, E, v)
    
    #ASSEMBLE, should be done using coo_matrix() if possible
    if Tri:  #Tri Elements
        for elem in range(0,nElem):            
            edofIndex=np.ix_(eDof[elem,:],eDof[elem,:])                        #Finding the indexes from eDof
            Ke=cfc.plante(ex[elem,:],ey[elem,:],ep,D)                          #Element Stiffness Matrix for Triangular Element
            K[edofIndex] = K[edofIndex] + x[elem]**SIMP_penal*Ke
    else:    #Quad Elements
        for elem in range(0,nElem):            
            edofIndex=np.ix_(eDof[elem,:],eDof[elem,:])                        #Finding the indexes from eDof
            Ke=cfc.plani4e(ex[elem,:],ey[elem,:],ep,D)                         #Element Stiffness Matrix for Quad Element
            K[edofIndex] = K[edofIndex] + x[elem]**SIMP_penal*Ke
            

    toc1 = time.perf_counter()

    
    #F[1,0]=-1e6
    
    #Add Boundary Conditions
        
    #fixDofs = np.array(fixDofs)        
            
    allDofs = [i for i in range(0,nElem)]
    freeDofs = np.setdiff1d(allDofs, fixDofs)
    
    tic2 = time.perf_counter()
    Ue = spsolve(K[np.ix_(freeDofs,freeDofs)],F[np.ix_(freeDofs)])
    toc2 = time.perf_counter()
    
    j = -1
    for i in Ue:
        j = j + 1
        U[freeDofs[j]] = i 


    if Timers==True:
        print('FE, ASSEM:')
        print(toc1-tic1)
        print('FE, SOLVE:')
        print(toc2-tic2)

    return U


















