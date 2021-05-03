


import numpy as np
import calfem.core as cfc
from scipy.sparse.linalg import spsolve


def Displacement(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K):
    
    
    
    index1D=np.ix_(freedofs)
    index2D=np.ix_(freedofs,freedofs)
    
    lambdaF = np.zeros(U.shape)
    lambdaF[index1D] = -2*spsolve(K[index2D],U[index1D]).reshape(len(freedofs),1)

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
            lambdaFe = lambdaF[np.ix_(edof[elem,:]-1)]
            dG0[elem] =  SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(lambdaFe.T, np.matmul(Ke,Ue))
                
                    
        else:
                        
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
                    
            
    G0=np.matmul(U.T,U)
    return G0, dG0










