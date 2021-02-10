

import FE



"""Numerical Sensitivity Analysis"""

def num_Sens_Anal(x,SIMP_penal,edof,coords,bc,f,ep,mp,numElem):
    
    eps=1e-6
    dc = x.copy() 
    
    for elem in range(0,numElem):
        
        x[elem]=x[elem]-eps
        U1 = FE._FE(x,SIMP_penal,edof,coords,bc,f,ep,mp)
    
        x[elem]=x[elem]+eps
        U2 = FE._FE(x,SIMP_penal,edof,coords,bc,f,ep,mp)
    
        G01=f*U1
        G02=f*U2
        
        dc[elem]=(G01-G02)/eps
    
    return dc














