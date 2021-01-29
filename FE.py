import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc



def _FE(x,SIMP_penal,eDof,coords,fixDofs,F,ep,mp):

    #Settings
    E=mp[0]
    v=mp[1]
    ptype=ep[0]
    Timers=False     #Print Timers 
    
    #Check sizes
    nDof=np.max(eDof)
    nElem=np.size(eDof,0)
    nx=coords[:,0]
    ny=coords[:,1]
    
    #Initialize Vecotors and Matrices
    K = np.zeros([nDof,nDof])
    #F = np.zeros([nDof,1])
    U = np.zeros([nDof,1])
    
    
    #Check element type
    if len(eDof[0,:])==6:   #Triangular Element
        Tri=True
        elemX=np.zeros([nElem,3])
        elemY=np.zeros([nElem,3])
    elif len(eDof[0,:])==8:
        Tri=False           #Use Quad Instead
        elemX=np.zeros([nElem,4])
        elemY=np.zeros([nElem,4])
    else:
        raise Exception('Unrecognized Element Shape, Check eDof Matrix')


    #Find The coordinates for each element's nodes
    for elem in range(0,nElem):
        
        nNode=np.ceil(np.multiply(eDof[elem,:],0.5))-1
        nNode=nNode.astype(int)
        
        elemX[elem,:]=nx[nNode[0:8:2]]
        elemY[elem,:]=ny[nNode[0:8:2]]

    
    #Start timer
    tic1 = time.perf_counter()
    
    #Sparse Matrices??
    #K=csr_matrix(K)
    #F=csr_matrix(F)
    
    #Linear Elastic Constitutive Matrix
    D=cfc.hooke(ptype, E, v)
    
    #ASSEMBLE, should be done using coo_matrix() if possible
    if Tri:  #Tri Elements
        for elem in range(0,nElem):  
            #breakpoint()
            edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)                        #Finding the indexes from eDof
            Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D)                    #Element Stiffness Matrix for Triangular Element
            K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke
    else:    #Quad Elements
        for elem in range(0,nElem):  
            edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)                        #Finding the indexes from eDof
            Ke=cfc.plani4e(elemX[elem,:],elemY[elem,:],ep,D)                   #Element Stiffness Matrix for Quad Element
            #breakpoint()
            K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke[0]
            

    toc1 = time.perf_counter()

    
    #F[1,0]=-1e6
    
    #Add Boundary Conditions
        
    #fixDofs = np.array(fixDofs)        
    

    allDofs = [i for i in range(0,nDof)]        
    freeDofs = np.setdiff1d(allDofs, fixDofs)
    breakpoint()


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


















