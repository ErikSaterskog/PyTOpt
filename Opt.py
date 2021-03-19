
import numpy as np
import nlopt 
import calfem.core as cfc 
import Filter


class Optimisation:
    
    def __init__(self):
        
        pass
    
    def OC(self,nel,x,volfrac,dc):
        l1,l2,step,damping=0,1e5,0.2,0.5
        while (l2-l1) > 1e-4:
            lmid = 0.5*(l2+l1)
            
            FirstMin = np.minimum.reduce([(x+step),x*np.array(np.power((-dc/lmid),damping))])
            SecondMin = np.minimum.reduce([np.ones([nel,1]), FirstMin])
            FirstMax = np.maximum.reduce([(x-step),SecondMin])
            xnew = np.maximum.reduce([0.001*np.ones([nel,1]),FirstMax])
            
            if (sum(sum(xnew)) - volfrac*nel) > 0:
                l1 = lmid
            else:
                l2 = lmid
    
        return xnew

    def mma(self,nelem,SIMP_penal,edof,f,ep,elemx,elemy,D,weightMatrix,volFrac,x,elementType,el_type,FEM,materialFun):
        
        def Objfun(x,grad):
            x = x.reshape(nelem,1)
            grad = grad.reshape(nelem,1)
            global U
            U,dR,lambdaF,sig_VM = FEM.fe_nl(x,SIMP_penal,f,ep,elementType,materialFun)
            c = np.matmul(f.T,U)            
            for elem in range(nelem):
                if ep[3]==1:
                    
                    if el_type==2:
                        Ke=cfc.plante(elemx[elem,:],elemy[elem,:],ep[0:2],D) 
                    else:
                        Ke=cfc.plani4e(elemx[elem,:],elemy[elem,:],ep,D)[0]   
                        
                    Ue = U[np.ix_(edof[elem,:]-1)]
                    
                    grad[elem] = -SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    
                else:
                    lambdaFe = lambdaF[np.ix_(edof[elem,:]-1)]
                    grad[elem] = np.matmul(lambdaFe.T,dR[elem,:].reshape(np.size(edof,1),1))
            
            
            grad[:] = Filter.Check(x,grad,weightMatrix)
            grad = grad.reshape(len(x),)
            x = x.reshape(len(x),)
            
            ce = np.array(c[0]).reshape([1,])
            print('G0:'+str(ce[0]))
            print('------------------------------')
            return ce[0]
        
        def VolCon(x,grad):
            grad[:] = 4**x.copy()
            x = x.reshape(nelem,1)
            grad = grad.reshape(nelem,1)
            grad[:] = Filter.Check(x,grad,weightMatrix)
            grad = grad.reshape(len(x),)
            x = x.reshape(len(x),)
            return sum(x)-volFrac*len(x)
        
        opt = nlopt.opt(nlopt.LD_MMA, len(x))           #Choosing optmisation algorithm
        opt.set_min_objective(Objfun)
        lb = np.zeros(np.size(x))+0.001                    
        ub = np.ones(np.size(x))                           
        opt.set_lower_bounds(lb)
        opt.set_upper_bounds(ub)
        opt.add_inequality_constraint(VolCon,0)
        opt.set_xtol_rel(1e-2)
        opt.set_maxeval(150)
        x = opt.optimize(x.reshape(len(x),))
        
        return x