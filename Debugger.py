

import FE
import numpy as np



"""Numerical Sensitivity Analysis"""

def num_Sens_Anal(x,SIMP_penal,edof,coords,bc,f,ep,mp,numElem):
    
    eps=1e-7
    dc = x.copy() 
    
    for elem in range(0,numElem):
        
        x[elem]=x[elem]-eps
        U1 = FE._FE(x,SIMP_penal,edof,coords,bc,f,ep,mp)
    
        x[elem]=x[elem]+(eps*2)
        U2 = FE._FE(x,SIMP_penal,edof,coords,bc,f,ep,mp)
    

        G01=np.matmul(np.transpose(f),U1)
        G02=np.matmul(np.transpose(f),U2)
        
        dc[elem]=(G02-G01)/(eps*2)
    
    return dc














