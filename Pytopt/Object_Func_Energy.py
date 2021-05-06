"""
Calculates the objective function energy and its derivate.

Inputs:
    nelem       -number of elements
    ep          -element parameters
    el_type     -element type
    elemx       -x-coordinates for element corners
    elemy       -y-coordinates for element corners
    D           -constitutive matrix
    eq          -Body forces
    U           -Displacement vector
    edof        -Element degrees of freedom
    fext_tilde  -External force vector with only body forces
    fextGlobal  -External force vector
    SIMP_const  -Solid isotropic material with penalisation constant
    x           -Design variables
    dG0         -Derivative of objective function
    dR          -Derivative of Residual from FEM
    freedofs    -Free degree of freedoms
    K           -Stiffness matrix
    
Output:
    G0          -Objective Function
    dG0         -Derivative of objective function
    
    

Written 2021-05
Made By: Daniel Pettersson & Erik Säterskog
"""


import numpy as np
import calfem.core as cfc
from scipy.sparse.linalg import spsolve


def Energy(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_const, x, dG0, dR, freedofs, K):
    
    
    
    index1D=np.ix_(freedofs)
    index2D=np.ix_(freedofs,freedofs)
    if ep[3]==False:
        lambdaF = np.zeros(U.shape)
        lambdaF[index1D] = -spsolve(K[index2D],fextGlobal[freedofs]).reshape(len(freedofs),1)
    
    for elem in range(nelem):
        if ep[3]:
                    
            if el_type==2:
                Ke=cfc.plante(elemx[elem,:],elemy[elem,:],ep[0:2],D) #Element Stiffness Matrix for Triangular Element
                if eq is not None:
                    eqe=[eq,eq,eq]
                    eqe=np.array(eqe).reshape(1,6)
                else:
                        eqe=np.zeros([1,6])
            else:
                Ke=cfc.plani4e(elemx[elem,:],elemy[elem,:],ep,D)[0]  #Element Stiffness Matrix for Quad Element
                if eq is not None:
                    eqe=[eq,eq,eq,eq]
                    eqe=np.array(eqe).reshape(1,8)
                else:
                    eqe=np.zeros([1,8])

            Ue = U[np.ix_(edof[elem,:]-1)]
            fext_tildee = fext_tilde[np.ix_(edof[elem,:]-1)]
            dG0[elem] = np.matmul(fext_tildee.T,Ue)  -  SIMP_const*x[elem][0]**(SIMP_const-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    
                    
        else:
            if eq is not None:
                eqe=[eq,eq,eq]
                eqe=np.array(eqe).reshape(1,6)
            else:
                eqe=np.zeros([1,6])

            
            lambdaFe = lambdaF[np.ix_(edof[elem,:]-1)]
            
            Ue = U[np.ix_(edof[elem,:]-1)]
            fext_tildee = fext_tilde[np.ix_(edof[elem,:]-1)]
            dG0[elem] = np.matmul(fext_tildee.T,Ue) + np.matmul(lambdaFe.T,dR[elem,:].reshape(np.size(edof,1),1))
                    
                
        G0=np.matmul(fextGlobal.T,U)
            
    return G0, dG0










