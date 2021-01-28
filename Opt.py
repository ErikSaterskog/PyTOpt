import numpy as np

class Optimisation:
    
    def __init__(self):
        
        pass
    
    def OC(self,nelx,nely,x,volfrac,dc):
        l1,l2,move=0,1e5,0.2
        while (l2-l1) > 1e-4:
            lmid = 0.5*(l2+l1)
            xnew = np.maximum.reduce([0.001*np.ones([nely,nelx]),np.maximum.reduce([np.maximum.reduce([(x-move),np.minimum.reduce([np.ones([nely,nelx]),np.minimum.reduce([(x+move),x*(np.sqrt(-dc/lmid))]) ]) ]) ])])
    
            if (sum(sum(xnew)) - volfrac*nelx*nely) > 0:
                l1 = lmid
            else:
                l2 = lmid
    
        return xnew
