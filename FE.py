

import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc
import Mod_Hook as mh
import Plani4s
import elem3n
import elem4n
import MaterialModelSelection as MMS


class FE():
    
    
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
            
            self.elemX=np.zeros([self.nElem,3])
            self.elemY=np.zeros([self.nElem,3])
        elif len(eDof[0,:])==8: #Quad Element
            
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
        
    
        
    def fe(self,x,SIMP_penal,F,ep,fun):
        #Settings
        epLin=ep.copy()
        epLin[3]=1
        Timers=True     #Print Timers 
        
        tic1 = time.perf_counter()       #Start timer
        
        row=[]
        col=[]
        data=[]
        
        #ASSEMBLE, should be done using coo_matrix() if possible
        
        for elem in range(self.nElem):  
                
            #part 1
            #ticedof1 = time.perf_counter()
            edofIndex=(self.eDof[elem,:]-1).tolist() 
            edofIndex2D=np.ix_(self.eDof[elem,:]-1,self.eDof[elem,:]-1)
            #ticedof2 = time.perf_counter()
            
            #part 2
            #ticedof = time.perf_counter()
            Ke, fint, fext, stress, epsilon=fun(np.zeros(np.size(self.eDof,1),), self.elemX[elem,:], self.elemY[elem,:], epLin, self.mp) #här kna man skicka in en materiafunktion istället för att definera den i elem3n
            #ticedof = time.perf_counter()
            
            #part 3
            3ticedof = time.perf_counter()
            Ke=np.matrix(Ke)               
            #ticedof = time.perf_counter()
            
            #part 4
            #ticedof = time.perf_counter()
            row.extend(edofIndex*len(edofIndex))
            col.extend(np.repeat(edofIndex,len(edofIndex)))
            data.extend(np.reshape(Ke*x[elem][0]**SIMP_penal,np.size(Ke)).tolist()[0])
            #ticedof = time.perf_counter()
        
        #part 5
        #ticedof = time.perf_counter()
        K=coo_matrix((data,(row,col)),shape=(self.nDof,self.nDof))
        K=K.tocsc()
        #ticedof = time.perf_counter()
        
        
        
        toc1 = time.perf_counter()
        tic2 = time.perf_counter()
        
        self.U[np.ix_(self.freeDofs)] = spsolve(K[np.ix_(self.freeDofs,self.freeDofs)],F[np.ix_(self.freeDofs)]).reshape(len(self.freeDofs),1)
        toc2 = time.perf_counter()
        
    
        if Timers==True:
            print('FE, Assem.: '+str(toc1-tic1))
            print('FE, Solve:  '+str(toc2-tic2))
            
    
        return self.U
    
    
    def fe_nl(self,x,SIMP_penal,F,ep,fun):
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
        err=1e9                  # Setting an error, arbritary big.
        TOL=1e-11*max(abs(F))    # Setting a resonable low tolerance. 
        
        U = FE.fe(self,x,SIMP_penal,F,ep,fun)
        

        lambdaF = U.copy()
        
        index1D=np.ix_(self.freeDofs)
        index2D=np.ix_(self.freeDofs,self.freeDofs)
        
        newtonIt = 0

        if ep[3]==1:  #Check if linear.
            return U,[],[]
            
        
        #Newton iteration loop until convergens.
        while err>TOL:
            row=[]
            col=[]
            data=[]
            newtonIt +=1
            R = np.zeros(np.shape(F)) 
            dR = np.zeros([self.nElem,np.size(self.eDof,1)])
            fextGlobal = F.copy()
            fintGlobal = np.zeros(U.shape)
            
            for elem in range(self.nElem):
                
                edofIndexList=(self.eDof[elem,:]-1).tolist()
                edofIndex1D=np.ix_(self.eDof[elem,:]-1)
                
                Ue = U[edofIndex1D]
                
                Ke, fint, fext, stress, epsilon=fun(Ue.reshape(np.size(self.eDof,1),), self.elemX[elem,:], self.elemY[elem,:], ep, self.mp) 
                Ke=np.matrix(Ke)
                
                fextGlobal[edofIndex1D]+=fext
                fintGlobal[edofIndex1D]+=fint*x[elem][0]**SIMP_penal
                dR[elem,:] = (SIMP_penal*x[elem][0]**(SIMP_penal-1)*fint).reshape(np.size(self.eDof,1),)
                
                row.extend(edofIndexList*len(edofIndexList))
                col.extend(np.repeat(edofIndexList,len(edofIndexList)))
                data.extend(np.reshape(Ke*x[elem][0]**SIMP_penal,np.size(Ke)).tolist()[0])
    
            
            K=coo_matrix((data,(row,col)),shape=(self.nDof,self.nDof))
            K=K.tocsc()
              
            R=fintGlobal-fextGlobal
            err = np.linalg.norm(R[self.freeDofs])

            
            U[index1D] = U[index1D] - spsolve(K[index2D],R[self.freeDofs]).reshape(len(self.freeDofs),1)
                    

        lambdaF[index1D] = -spsolve(K[index2D],F[self.freeDofs]).reshape(len(self.freeDofs),1)
              
        
        print('N.iters:    ' + str(newtonIt))
        print('Final error:' + str(err))
        return U,dR,lambdaF
    
    

















