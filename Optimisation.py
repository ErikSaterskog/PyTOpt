
import numpy as np
import calfem.core as cfc 
import Filter
from MMA_DTU_funct import mmasub,kktcheck


class Optimisation:
    
    def __init__(self):
        
        pass
    
    def OC(self,nel,x,volfrac,dc):
        l1,l2,step,damping=0,1e5,0.2,0.5
        while (l2-l1) > 1e-4:
            lmid = 0.5*(l2+l1)
            
            FirstMin = np.minimum.reduce([(x+step),np.array(-np.sign(dc))*x*np.array(np.power((abs(dc)/lmid),damping))])
            SecondMin = np.minimum.reduce([np.ones([nel,1]), FirstMin])
            FirstMax = np.maximum.reduce([(x-step),SecondMin])
            xnew = np.maximum.reduce([0.001*np.ones([nel,1]),FirstMax])
            
            if (sum(sum(xnew)) - volfrac*nel) > 0:
                l1 = lmid
            else:
                l2 = lmid
    
        return xnew

    def mma(self,bc,MP,f,Edof,elemx,elemy,x,SIMP_penal,ep,elementType,materialFun,FEM,el_type,D,eq,weightMatrix,volFrac,ObjectFun):
        # -*- coding: utf-8 -*-
        """
        This function is uses Method of moving asymptotes along with KKt check to 
        optimize the objective function.
        
           Input:
            bc   : Boundary condition               
            MP   : Material parameter [E,A,Sigma_y] 
            f    : Loads                            
            Edof : Element degrees of freedom       
            ex   : Elements position in x           
            ey   : Elements position in y           
            Ndof : Number of degrees of freedom     
            
          Output:
             xval: vector with the optimized areas
             ed  : displacements
    
        """


       
        
        #algorithmic parameters and initial guess
        nelem = np.size(Edof,0)
        
        #initialization
        eeen = np.ones((nelem,1))
        eeem = np.ones((1,1))
        zerom = np.zeros((1,1))
        xmin = 1.e-2*eeen   #lower bound on x
        xmax = 1*eeen       #upper bound on x
        move = 1.0
        c = 1000*eeem
        d = zerom.copy()
        a0 = 0
        a=np.zeros((1,1))
        xold1 = x.copy()
        xold2 = x.copy()
        low = xmin.copy()
        upp = xmax.copy()
        
        maxoutit =250
        kkttol = 0	
        dc = x.copy()
        
        # Calculate function values and gradients of the objective and constraints functions
    
        
        U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x,SIMP_penal,f,ep,elementType,materialFun)
        
          
        f0val, dc = ObjectFun(nelem, ep, el_type, elemx, elemy, D, eq, U, Edof, fext_tilde, fextGlobal, SIMP_penal, x, dc, dR, freedofs, K)
        
        dc = Filter.Filter(x,dc,weightMatrix)
        
        
    
        dfdx = eeen.copy()
        dfdx = dfdx.reshape(1,nelem)
        fval = sum(x)-volFrac*len(x)
        
        # The iterations starts
        kktnorm = kkttol+10
        outit = 0
        
    #-----------------------------------------------------------------------------
        
        while (kktnorm > kkttol) and (outit < maxoutit):
            outit += 1
            
            
            dc = np.array(dc)
            x,y,z,lam,xsi,eta,mu,zet,s,low,upp =  \
                mmasub(1,nelem,outit,x,xmin,xmax,xold1,xold2,f0val,dc,fval,dfdx,low,upp,a0,a,c,d,move)
            dc = np.matrix(dc)                                                                
            # update
            xold2 = xold1.copy()
            xold1 = x.copy()
            
    #------------------------------------------------------------------------------        
            # Re-calculate function values and gradients of the objective and constraints functions
            U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x,SIMP_penal,f,ep,elementType,materialFun)
            f0val = np.matmul(f.T,U)  
              
            f0val, dc = ObjectFun(nelem, ep, el_type, elemx, elemy, D, eq, U, Edof, fext_tilde, fextGlobal, SIMP_penal, x, dc, dR, freedofs, K)
            
            dc = Filter.Filter(x,dc,weightMatrix)
            
        
            dfdx = eeen.copy()
            dfdx = dfdx.reshape(1,nelem)
            fval = sum(x)-volFrac*len(x)
    #----------------------------------------------------------------------------- 
       
            # The residual vector of the KKT conditions is calculated
            dc = np.array(dc)
            residu,kktnorm,residumax = \
                kktcheck(1,nelem,x,y,z,lam,xsi,eta,mu,zet,s,xmin,xmax,dc,fval,dfdx,a0,a,c,d)
            dc = np.matrix(dc)
            print(f0val)
        return f0val,x,eps_h

