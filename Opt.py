import numpy as np

import nlopt 
 


class Optimisation:
    
    def __init__(self):
        
        pass
    
    def OC(self,nel,x,volfrac,dc):
        l1,l2,move=0,1e5,0.2
        while (l2-l1) > 1e-4:
            lmid = 0.5*(l2+l1)
            FirstMin = np.minimum.reduce([(x+move),x*(np.sqrt(-dc/lmid))])
            SecondMin = np.minimum.reduce([np.ones([nel,1]), FirstMin])
            FirstMax = np.maximum.reduce([(x-move),SecondMin])
            xnew = np.maximum.reduce([0.001*np.ones([nel,1]),FirstMax])
            if (sum(sum(xnew)) - volfrac*nel) > 0:
                l1 = lmid
            else:
                l2 = lmid
    
        return xnew

    def mma(self,nel,x,volfrac,dc):
        
        opt = nlopt.opt("NLOPT_LD_MMA", 10) #Choosing optmisation algorithm
        def f(x, grad):                     #f is the objective function
            if grad.size > 0:               #grad is kind of the number of elements
                pass
        return f
    
        opt.set_min_objective(f)
        lb = 1                              #Placeholders for now
        ub = 100                            #Placeholders for now
        opt.set_lower_bounds(lb)
        opt.set_upper_bounds(ub)
        
        
        # Both fc and h is in the same form as f. 
        opt.add_inequality_constraint(fc, tol=0)
        opt.add_equality_constraint(h, tol=0)
        
        def c(result, x, grad):
            if grad.size > 0:
                pass
                #...set grad to gradient, in-place...
                #result[0] = ...value of c_0(x)...
                #result[1] = ...value of c_1(x)...
        return c
        
        opt.add_inequality_mconstraint(c, tol)
        opt.add_equality_mconstraint(c, tol)
        
        #Later here should stopping criterion be implemented
        
        xnew = opt.optimize(x)
        
        return xnew