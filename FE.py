# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:52:55 2021

@author: Daniel
"""



import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc
import Mod_Hook as mh
import Plani4s
import elem3n
import elem4n


class _FE():
    
    
    def __init__(self,eDof,coords,mp,fixDofs):
        self.mp = mp
        self.E=mp[0]
        self.v=mp[1]
        
        self.nDof=np.max(eDof)
        self.nElem=np.size(eDof,0)
        nx=coords[:,0]
        ny=coords[:,1]
        
        #Initialize Vecotors and Matrices
        #K = np.zeros([nDof,nDof])
        self.K = coo_matrix((self.nDof,self.nDof))
        #F = np.zeros([nDof,1])
        self.U = np.zeros([self.nDof,1])
        
        
        #Check element type
        if len(eDof[0,:])==6:   #Triangular Element
            self.Tri=True
            self.elemX=np.zeros([self.nElem,3])
            self.elemY=np.zeros([self.nElem,3])
        elif len(eDof[0,:])==8:
            self.Tri=False           #Use Quad Instead
            self.elemX=np.zeros([self.nElem,4])
            self.elemY=np.zeros([self.nElem,4])
        else:
            raise Exception('Unrecognized Element Shape, Check eDof Matrix')
    
    
        #Find The coordinates for each element's nodes
        for elem in range(self.nElem):
            
            nNode=np.ceil(np.multiply(eDof[elem,:],0.5))-1
            nNode=nNode.astype(int)
            
            self.elemX[elem,:]=nx[nNode[0:8:2]]
            self.elemY[elem,:]=ny[nNode[0:8:2]]
        
        allDofs = range(self.nDof)        
        self.freeDofs = np.setdiff1d(allDofs, fixDofs)
        self.eDof = eDof
        
    
        
    def fe(self,x,SIMP_penal,F,ep):
        #Settings
        
        ptype=ep[0]
        Timers=True     #Print Timers 
        
        
        tic1 = time.perf_counter()       #Start timer
        
        
        D=cfc.hooke(ptype, self.E, self.v)         #Linear Elastic Constitutive Matrix
        
        row=[]
        col=[]
        data=[]
        
        #ASSEMBLE, should be done using coo_matrix() if possible
        
        for elem in range(self.nElem):  
                
            edofIndex=(self.eDof[elem,:]-1).tolist() 
            if self.Tri:
                Ke=cfc.plante(self.elemX[elem,:],self.elemY[elem,:],ep[0:2],D)               #Element Stiffness Matrix for Triangular Element
            else:
                Ke=cfc.plani4e(self.elemX[elem,:],self.elemY[elem,:],ep,D)[0] 
            row.extend(edofIndex*len(edofIndex))
            col.extend(np.repeat(edofIndex,len(edofIndex)))
            data.extend(np.reshape(Ke*x[elem][0]**SIMP_penal,np.size(Ke)).tolist()[0])
    
        K=coo_matrix((data,(row,col)),shape=(self.nDof,self.nDof))
        K=K.tocsc()
        
        toc1 = time.perf_counter()
        tic2 = time.perf_counter()
        
        self.U[np.ix_(self.freeDofs)] = spsolve(K[np.ix_(self.freeDofs,self.freeDofs)],F[np.ix_(self.freeDofs)]).reshape(len(self.freeDofs),1)
        toc2 = time.perf_counter()
        
    
        if Timers==True:
            print('FE, Assem.: '+str(toc1-tic1))
            print('FE, Solve:  '+str(toc2-tic2))
            
    
        return self.U
    
    
    def fe_nl(self,x,SIMP_penal,F,ep):
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
        err=1e9     # Setting an error, arbritary big.
        TOL=1e-6    # Setting a resonable low tolerance. 
        
        ep=[2,1,2,2]   #TEMPORARY TESTING
        
        
        U = FE.fe(self,x,SIMP_penal,F,ep)
        
        
        #Newton iteration loop until convergens.
            # Starting by calculating the strain and then define the constitutive 
            # relation when having a nonlinear material model. Reassemble the nonlinear
            # stiffness matrix K. Checking the the residual.
        while err>TOL:
            
            K=np.zeros(np.shape(self.K))
            R = np.zeros(np.shape(F))
            fextGlobal = np.zeros(np.shape(F))
            fintGlobal = np.zeros(np.shape(F)) 
            #ed=cfc.extractEldisp(self.eDof,U)
            dr = np.zeros([np.size(self.eDof,1),self.nElem])
            
            
            for elem in range(self.nElem):
                edofIndex2D=np.ix_(self.eDof[elem,:]-1,self.eDof[elem,:]-1)
                edofIndex1D=np.ix_(self.eDof[elem,:]-1)
                Ue = U[edofIndex1D]
                if self.Tri:
                    Ke, fint, fext, stress, epsilon=elem3n.elem3n(Ue.reshape(np.size(self.eDof,1),), self.elemX[elem,:], self.elemY[elem,:], ep, self.mp) #här kna man skicka in en materiafunktion istället för att definera den i elem3n
                else:
                    Ke, fint, fext, stress, epsilon=elem4n.elem4n(Ue.reshape(np.size(self.eDof,1),), self.elemX[elem,:], self.elemY[elem,:], ep, self.mp) #här kna man skicka in en materiafunktion istället för att definera den i elem3n
                
                
                fext+=F[edofIndex1D]
                
                K[edofIndex2D]=K[edofIndex2D]+Ke*x[elem][0]**SIMP_penal
                
                R[edofIndex1D]=R[edofIndex1D]+fint*x[elem][0]**SIMP_penal-fext
                
                dr[:,elem] = (SIMP_penal*x[elem][0]**(SIMP_penal-1)*fint).reshape(np.size(self.eDof,1),)
                
                fextGlobal[edofIndex1D]+=fext
                fintGlobal[edofIndex1D]+=fint
                
            err = np.linalg.norm(R[self.freeDofs])
            print(err)
               
            U[np.ix_(self.freeDofs)] = U[np.ix_(self.freeDofs)] - spsolve(K[np.ix_(self.freeDofs,self.freeDofs)],R[self.freeDofs]).reshape(len(self.freeDofs),1)
                    
        
        return U,K,dr
    
    

















