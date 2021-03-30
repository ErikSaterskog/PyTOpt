


import numpy as np
from scipy.sparse.linalg import spsolve


def Displacement(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K):
    
    
    
    index1D=np.ix_(freedofs)
    index2D=np.ix_(freedofs,freedofs)
    
    lambdaF = np.zeros(U.shape)
    lambdaF[index1D] = -2*spsolve(K[index2D],U[index1D]).reshape(len(freedofs),1)

    for elem in range(nelem):  
        if eq is not None:
            eqe=[eq,eq,eq]
            eqe=np.array(eqe).reshape(1,6)
        else:
            eqe=np.zeros([1,6])

        lambdaFe = lambdaF[np.ix_(edof[elem,:]-1)]
        

        dG0[elem] = np.matmul(lambdaFe.T,dR[elem,:].reshape(np.size(edof,1),1))
                    
        if dG0[elem] >0:
            print(str(elem) + ':' +str(dG0[elem]))
            dG0[elem] = 0
                
        G0=np.matmul(fextGlobal.T,U)

    return G0, dG0










