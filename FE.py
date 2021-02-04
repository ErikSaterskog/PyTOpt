import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc
import Mod_Hook as mh



def _FE(x,SIMP_penal,eDof,coords,fixDofs,F,ep,mp):

    #Settings
    E=mp[0]
    v=mp[1]
    ptype=ep[0]
    Timers=True     #Print Timers 
    
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

    tic2 = time.perf_counter()
    Ue = spsolve(K[np.ix_(freeDofs,freeDofs)],F[np.ix_(freeDofs)])
    toc2 = time.perf_counter()
    
    j = -1
    for i in Ue:
        j = j + 1
        U[freeDofs[j]] = i 


    if Timers==True:
        print('FE, Assem.: '+str(toc1-tic1))
        print('FE, Solve:  '+str(toc2-tic2))
        

    return U








def _FE_NL(x,SIMP_penal,eDof,coords,fixDofs,F,ep,mp):
    
    #Settings
    E=mp[0]
    v=mp[1]
    ptype=ep[0]
    Timers=True     #Print Timers 
    err=1e9
    TOL=1e6
    
    #Check sizes
    nDof=np.max(eDof)
    nElem=np.size(eDof,0)
    nx=coords[:,0]
    ny=coords[:,1]
    
    #Initialize Vecotors and Matrices
    K = np.zeros([nDof,nDof])
    #F = np.zeros([nDof,1])
    U = np.zeros([nDof,1])
    
    D=cfc.hooke(ptype, E, v)
    
    
    allDofs = [i for i in range(0,nDof)]        
    freeDofs = np.setdiff1d(allDofs, fixDofs)
    
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

    
    
    #check element shape
    if Tri==True:
    
        
        for elem in range(0,nElem):    #Gissning med linjärt fall
            Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D)
            K=cfc.assem(eDof,K,Ke)

    
        Ue = spsolve(K[np.ix_(freeDofs,freeDofs)],F[np.ix_(freeDofs)])
            
        j = -1
        for i in Ue:
            j = j + 1
            U[freeDofs[j]] = i
            
        while err>TOL:
    
            ed=cfc.extractEldisp(eDof,U)

            for elem in range(0,nElem):    
                eps = np.zeros([6,])
                discard,eps_2D=cfc.plants(elemX[elem,:],elemY[elem,:],ep,D,ed[elem,:]) 
                eps[0:4] = eps_2D
                sig,D_new = mh._mod_hook(eps,mp)
                Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D_new[np.ix_([0,1,3],[0,1,3])])
                K=cfc.assem(eDof,K,Ke)
            g= np.matmul(K[np.ix_(freeDofs,freeDofs)],U[np.ix_(freeDofs)])-F[np.ix_(freeDofs)]
            err = np.linalg.norm(g)
            print(err)
            
            U[np.ix_(freeDofs)] = U[np.ix_(freeDofs)] - spsolve(K[np.ix_(freeDofs,freeDofs)],g).reshape(48,1)
            
            
    else:
        raise Exception ('QUADS NOT YET IMPLEMENTED')
        
    return U,D

#
#
#
#
#
#
#
#












