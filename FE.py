import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc
import Mod_Hook as mh
import Plani4s



def FE(x,SIMP_penal,eDof,coords,fixDofs,F,ep,mp):

    #Settings
    E=mp[0]
    v=mp[1]
    ptype=ep[0]
    Timers=True     #Print Timers 
    
    
    nDof,nElem,K,U,Tri,elemX,elemY=init(eDof,coords)
    
    
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
            edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)                        #Finding the indexes from eDof
            Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D)                    #Element Stiffness Matrix for Triangular Element
            K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke
    else:    #Quad Elements
        for elem in range(0,nElem):  
            edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)                        #Finding the indexes from eDof
            Ke=cfc.plani4e(elemX[elem,:],elemY[elem,:],ep,D)                   #Element Stiffness Matrix for Quad Element
            K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke[0]
            

    toc1 = time.perf_counter()


    allDofs = range(nDof)        
    freeDofs = np.setdiff1d(allDofs, fixDofs)

    tic2 = time.perf_counter()
    U[np.ix_(freeDofs)] = spsolve(K[np.ix_(freeDofs,freeDofs)],F[np.ix_(freeDofs)]).reshape(len(freeDofs),1)
    toc2 = time.perf_counter()
    

    if Timers==True:
        print('FE, Assem.: '+str(toc1-tic1))
        print('FE, Solve:  '+str(toc2-tic2))
        

    return U








def _FE_NL(x,SIMP_penal,eDof,coords,fixDofs,F,ep,mp):
    
    #Settings
    E=mp[0]
    v=mp[1]
    ptype=ep[0]
    err=1e9
    TOL=1e-6
    Timers=True     #Print Timers 
    
    nDof,nElem,K,U,Tri,elemX,elemY=init(eDof,coords)
    
    D=cfc.hooke(ptype, E, v)
    
    allDofs = range(nDof)       
    freeDofs = np.setdiff1d(allDofs, fixDofs)

        
    for elem in range(nElem):    #Gissning med linjärt fall
        edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1) 
        if Tri:
            Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D)
        else:
            Ke=cfc.plani4e(elemX[elem,:],elemY[elem,:],ep,D)[0]
        K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke

    
    U[np.ix_(freeDofs)] = spsolve(K[np.ix_(freeDofs,freeDofs)],F[np.ix_(freeDofs)]).reshape(len(freeDofs),1)
          
    while err>TOL:
        K = np.zeros(np.shape(K))
        ed=cfc.extractEldisp(eDof,U) 

        for elem in range(nElem):
            if Tri:
                edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)
                eps = np.zeros([6,])
                eps_2D=cfc.plants(elemX[elem,:],elemY[elem,:],ep,D,ed[elem,:])[1] 
                eps[0:4] = eps_2D
                D_new = mh._mod_hook(eps,mp)[1]
                Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D_new[np.ix_([0,1,2,3],[0,1,2,3])])
                K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke
            else:
                edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)
                eps=Plani4s.plani4s(elemX[elem,:],elemY[elem,:],ep,ed[elem,:]) 
                D_new = mh._mod_hook(eps,mp)[1]
                Ke=cfc.plani4e(elemX[elem,:],elemY[elem,:],ep,D_new[np.ix_([0,1,2,3],[0,1,2,3])])[0]
                K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke
                
        R= np.matmul(K[np.ix_(freeDofs,freeDofs)],U[np.ix_(freeDofs)])-F[np.ix_(freeDofs)]
        err = np.linalg.norm(R)
        print(err)
            
        U[np.ix_(freeDofs)] = U[np.ix_(freeDofs)] - spsolve(K[np.ix_(freeDofs,freeDofs)],R).reshape(len(freeDofs),1)
            
        
    return U



def init(eDof,coords):
    
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
    for elem in range(nElem):
        
        nNode=np.ceil(np.multiply(eDof[elem,:],0.5))-1
        nNode=nNode.astype(int)
        
        elemX[elem,:]=nx[nNode[0:8:2]]
        elemY[elem,:]=ny[nNode[0:8:2]]
        
    return nDof,nElem,K,U,Tri,elemX,elemY











