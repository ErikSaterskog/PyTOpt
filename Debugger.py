

import FE_test
import numpy as np



"""Numerical Sensitivity Analysis"""

def num_Sens_Anal(x,SIMP_penal,edof,coords,bc,f,ep,mp,numElem):
    
    eps=1e-7
    dc = x.copy() 
    FEM = FE_test.FE(edof,coords,mp,bc)    
    
    for elem in range(0,numElem):
        
        x[elem]=x[elem]-eps
        U1 = FEM.fe(x,SIMP_penal,f,ep).copy()

        x[elem]=x[elem]+(eps*2)
        U2 = FEM.fe(x,SIMP_penal,f,ep).copy()

        x[elem]=x[elem]-(eps)

        G01=np.matmul(np.transpose(f),U1)
        G02=np.matmul(np.transpose(f),U2)
        
        dc[elem]=(G02-G01)/(eps*2)
    return dc

def num_Sens_Anal_NL(x,SIMP_penal,edof,coords,bc,f,ep,mp,numElem):
    
    eps=1e-7
    dc = x.copy() 
    FEM = FE_test.FE(edof,coords,mp,bc)   
    
    for elem in range(0,numElem):
        
        x[elem]=x[elem]-eps       
        U,K,dr = FEM.fe_nl(x,SIMP_penal,f,ep)
        U1=U.copy()
        
        x[elem]=x[elem]+(eps*2)       
        U,K,dr = FEM.fe_nl(x,SIMP_penal,f,ep)
        U2=U.copy()
        
        x[elem]=x[elem]-(eps)

        G01=np.matmul(np.transpose(f),U1)
        G02=np.matmul(np.transpose(f),U2)
        
        dc[elem]=(G02-G01)/(eps*2)
    
    return dc














