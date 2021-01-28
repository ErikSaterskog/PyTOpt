import numpy as np

class Optimisation:
    
    def __init__(self):
        
        pass
    
    def OC(self,nel,x,volfrac,dc):
        l1,l2,move=0,1e5,0.2
        while (l2-l1) > 1e-4:
            lmid = 0.5*(l2+l1)
            FirstMin = np.minimum.reduce([(x+move),np.transpose(np.transpose(x)*(np.sqrt(-dc/lmid)))])
            SecondMin = np.minimum.reduce([np.ones([nel,1]), FirstMin])
            FirstMax = np.maximum.reduce([(x-move),SecondMin])
            xnew = np.maximum.reduce([0.001*np.ones([nel,1]),FirstMax])
            if (sum(sum(xnew)) - volfrac*nel) > 0:
                l1 = lmid
            else:
                l2 = lmid
    
        return xnew
