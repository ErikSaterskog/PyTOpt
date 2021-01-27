import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc



def FE(x,SIMP_penal,eDof,coord):

    
    #Check sizes
    nDof=np.max(eDof)
    nElem=np.size(eDof,0)
    ex=coord[0]
    ey=coord[1]
    
    
    #Initialize Vecotors and Matrices
    K = np.zeros([nDof,])
    F = np.zeros([nDof,1])
    U = np.zeros([nDof,1])
    
    fixdofs = []
    
    #Check element type
    if len(eDof[0,:])==4:   #Triangular Element
        Tri=True
    elif len(eDof[0,:])==5:
        Tri=False           #Use Quad Instead
    else:
        raise Exception('Unrecognized Element Shape, Check eDof Matrix')
    
    
    #Start timer
    tic = time.perf_counter()
    
    #K=csr_matrix(K)
    #F=csr_matrix(F)
    
    #ASSEMBLE, should be done using coo_matrix() if possible
    for elem in range(0,nElem):            
            edofIndex=np.ix_(eDof[elem,:],eDof[elem,:])      #Finding the indexes from eDof
            if Tri:
                Ke=cfc.plante(ex[elem,:],ey[elem,:],ep,D)    #Element Stiffness Matrix for Triangular Element
            else:
                Ke=cfc.plani4e(ex[elem,:],ey[elem,:],ep,D)   #Element Stiffness Matrix for Quad Element
            K[edofIndex] = K[edofIndex] + x[elem]**SIMP_penal*Ke
            

    toc = time.perf_counter()
    print('FE, ASSEM:')
    print(toc-tic)
    F[1,0]=-1e6
    

    
    tic = time.perf_counter()
        
    for i in range(0,2*(nely+1),2):
        fixdofs.append(i)
        
    fixdofs.append(2*(nelx+1)*(nely+1)-1)
    fixdofs = np.array(fixdofs)        
            
    alldofs = [i for i in range(0,2*(nely+1)*(nelx+1))]
    freedofs = np.setdiff1d(alldofs, fixdofs)
    

    Ue = spsolve(K[np.ix_(freedofs,freedofs)],F[np.ix_(freedofs)])

    
    j = -1
    for i in Ue:
        j = j + 1
        U[freedofs[j]] = i 

    toc = time.perf_counter()
    print('FE, SOLVE:')
    print(toc-tic)

    return U


















