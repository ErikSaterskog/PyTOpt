


import numpy as np
import calfem.core as cfc



def Energy(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, SIMP_penal, x, dG0, lambdaF, dR):
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
            dG0[elem] = np.matmul(fext_tildee.T,Ue)  -  SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    
            if dG0[elem] >0:
                print(str(elem) + ':' +str(dG0[elem]))
                dG0[elem] = 0
                    
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
                    
            if dG0[elem] >0:
                print(str(elem) + ':' +str(dG0[elem]))
                dG0[elem] = 0
                
            
    return dG0










