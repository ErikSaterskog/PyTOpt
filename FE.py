
import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc
import Mod_Hook as mh
import Plani4s
import elem3n



def FE(x,SIMP_penal,eDof,coords,fixDofs,F,ep,mp):

    #Settings
    E=mp[0]
    v=mp[1]
    ptype=ep[0]
    Timers=True     #Print Timers 
    
    nDof,nElem,K,U,Tri,elemX,elemY=init(eDof,coords)
    
    tic1 = time.perf_counter()       #Start timer
    
    allDofs = range(nDof)        
    freeDofs = np.setdiff1d(allDofs, fixDofs)
    D=cfc.hooke(ptype, E, v)         #Linear Elastic Constitutive Matrix
    
    #ASSEMBLE, should be done using coo_matrix() if possible
    if Tri:  #Tri Elements
        for elem in range(nElem):  
            edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)                    #Finding the indexes from eDof
            Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D)               #Element Stiffness Matrix for Triangular Element
            K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke
    else:    #Quad Elements
        for elem in range(nElem):  
            edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)                    #Finding the indexes from eDof
            Ke=cfc.plani4e(elemX[elem,:],elemY[elem,:],ep,D)                   #Element Stiffness Matrix for Quad Element
            K[edofIndex] = K[edofIndex] + x[elem][0]**SIMP_penal*Ke[0]
            

    toc1 = time.perf_counter()
    tic2 = time.perf_counter()
    
    U[np.ix_(freeDofs)] = spsolve(K[np.ix_(freeDofs,freeDofs)],F[np.ix_(freeDofs)]).reshape(len(freeDofs),1)
    toc2 = time.perf_counter()
    

    if Timers==True:
        print('FE, Assem.: '+str(toc1-tic1))
        print('FE, Solve:  '+str(toc2-tic2))
        

    return U


def _FE_NL(x,SIMP_penal,eDof,coords,fixDofs,F,ep,mp):
    """
    INPUT:
        x          - element densities, design variables
        SIMP_penal - Penaltyfactor preventing x to be between 0-1
        eDof       - Element degrees of freedom
        coord      - Node coordinates both in x- and y-direction
        fixDofs    - Degrees of freedom prescribed by boundary condition
        F          - Forcevector
        ep         - element parameters
                        - ptype(if plane strain or plain stress)
                        - t thickness
                        - ir intregration rule
        mp         - Material parameters consisting of Young's modulus and 
                     Poisson's ratio.
                     
                     
    OUTPUT:
        U          - Displacement vector
    
      
    """
    #Settings
    E=mp[0]
    v=mp[1]
    ptype=ep[0]
    err=1e9     # Setting an error, arbritary big.
    TOL=1e-6    # Setting a resonable low tolerance. 
    
    ep=[2,1,2,2]   #TEMPORARY TESTING
    
     # Collecting necessary variables as the K stiffness matrix for example.
    nDof,nElem,K,U,Tri,elemX,elemY=init(eDof,coords)
    
    U = FE(x,SIMP_penal,eDof,coords,fixDofs,F,ep,mp)
    
    
    # Define all degrees of freedom and free degrees of freedom.
    allDofs = range(nDof)       
    freeDofs = np.setdiff1d(allDofs, fixDofs)

    #Newton iteration loop until convergens.
        # Starting by calculating the strain and then define the constitutive 
        # relation when having a nonlinear material model. Reassemble the nonlinear
        # stiffness matrix K. Checking the the residual.
    while err>TOL:
        
        K=np.zeros(np.shape(K))
        R = np.zeros(np.shape(F))
        ed=cfc.extractEldisp(eDof,U)
        
        for elem in range(nElem):
            
            edofIndex=np.ix_(eDof[elem,:]-1,eDof[elem,:]-1)

            Ke, fint, fext, stress, epsilon=elem3n.elem3n((ed[elem,:]), elemX[elem,:], elemY[elem,:], ep, mp) #här kna man skicka in en materiafunktion istället för att definera den i elem3n
            
            K[edofIndex]=K[edofIndex]+Ke
        
            R[np.ix_(eDof[elem,:]-1)]=R[np.ix_(eDof[elem,:]-1)]+fint-fext-F[np.ix_(eDof[elem,:]-1)]
            

        err = np.linalg.norm(R[freeDofs])
        #print(err)
           
            
        U[np.ix_(freeDofs)] = U[np.ix_(freeDofs)] - spsolve(K[np.ix_(freeDofs,freeDofs)],R[freeDofs]).reshape(len(freeDofs),1)
                
    
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























