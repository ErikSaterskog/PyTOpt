

from Pytopt import FE
import numpy as np



"""Numerical Sensitivity Analysis"""


def num_Sens(x, SIMP_penal, edof, coords, bc, f, ep, mp, numElem, elementFun, materialFun, eq, eps):

    dc = x.copy() 
    FEM = FE.FE(edof,coords,mp,bc)  
    
    for elem in range(0,numElem):
        
        x[elem]=x[elem]-eps       
        U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x, SIMP_penal, f, ep, elementFun, materialFun, eq)
        U1=U.copy()
        
        x[elem]=x[elem]+(eps*2)       
        U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x, SIMP_penal, f, ep, elementFun, materialFun, eq)
        U2=U.copy()
        
        x[elem]=x[elem]-(eps)

        G01=np.matmul(np.transpose(f),U1)
        G02=np.matmul(np.transpose(f),U2)
        
        dc[elem]=(G02-G01)/(eps*2)
    
    return dc














