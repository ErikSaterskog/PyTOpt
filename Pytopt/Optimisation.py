
import numpy as np
from Pytopt import Filter
from Pytopt.MMA_fun import mmasub,kktcheck


def OC(x, volfrac, G0, dG0, loop):
    nel=len(x)
    l1,l2,step,damping=0,1e5,0.02,0.5
    while (l2-l1) > 1e-4:
        lmid = 0.5*(l2+l1)
        
        FirstMin = np.minimum.reduce([(x+step),np.array(-np.sign(dG0))*x*np.array(np.power((abs(dG0)/lmid),damping))])
        SecondMin = np.minimum.reduce([np.ones([nel,1]), FirstMin])
        FirstMax = np.maximum.reduce([(x-step),SecondMin])
        xnew = np.maximum.reduce([0.001*np.ones([nel,1]),FirstMax])
        
        if (sum(sum(xnew)) - volfrac*nel) > 0:
            l1 = lmid
        else:
            l2 = lmid
            
    return xnew

def MMA(x, volFrac, G0, dG0, loop):
    """
    This function is uses Method of moving asymptotes along with KKT check to 
    optimize the objective function.
    
       Input:
       
        f       : Loads                            
        edof    : Element degrees of freedom       
        elemx   : Elements position in x           
        elemy   : Elements position in y           
        
        
      Output:
         f0val  : The object function value 
         x      : Vector with the design variables
         eps_h  : Vector with the element hydrostatic strain 

    """
    
    #Initialization ##############################
    
    nelem = np.size(x,0)     
    eeen = np.ones((nelem,1))
    xmin = 1.e-2*eeen           #Lower bound on x
    xmax = eeen                 #Upper bound on x
    move = 0.01
    c = 1000 
    d = np.zeros((1,1))
    a0 = 1
    a=np.zeros((1,1))
    xold1 = x.copy()
    xold2 = x.copy()
    low = xmin.copy()
    upp = xmax.copy()
    
    
    ###########################################

    dfdx = eeen.copy()
    dfdx = dfdx.reshape(1,nelem)
    fval = sum(x)-volFrac*len(x)
    outit=0                        
    maxoutit=10
    kkttol = 0	
    kktnorm = 10

    while (kktnorm > kkttol) and (outit < maxoutit):
        outit += 1

        dG0 = np.array(dG0)
        x,y,z,lam,xsi,eta,mu,zet,s,low,upp = mmasub(1, nelem, outit, x, xmin, xmax, xold1, xold2, G0, dG0, fval, dfdx, low, upp, a0, a, c, d, move)
    
        # update
        xold2 = xold1.copy()
        xold1 = x.copy()
        
        dfdx = eeen.copy()
        dfdx = dfdx.reshape(1,nelem)
        fval = sum(x)-volFrac*len(x)
#----------------------------------------------------------------------------- 
   
        # The residual vector of the KKT conditions is calculated
        residu,kktnorm,residumax = kktcheck(1,nelem,x,y,z,lam,xsi,eta,mu,zet,s,xmin,xmax,dG0,fval,dfdx,a0,a,c,d)

    return x



