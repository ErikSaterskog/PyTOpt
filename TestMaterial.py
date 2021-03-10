
import numpy as np

def head(eps,mp):
    
    sig = mat(eps,mp)
    
    D = numD(eps,sig,mp)
    
    return sig.reshape(6,1),D

def mat(eps,mp):
    
   
    E1,E3        = mp[0]
    nu1,nu3      = mp[1]
    eps_y = mp[2]
    G1 = E1/(2*(1+nu1))
    K1 = E1/(3*(1-2*nu1))
    G3 = E3/(2*(1+nu3))
    K3 = E3/(3*(1-2*nu3))
    
    E2 = (E1+E3)/2
    nu2 = (nu1+nu3)/2
    G2 = E2/(2*(1+nu2))
    K2 = E2/(3*(1-2*nu2))
    
    eps_h = sum(eps[:2])/3
    
    I_v = np.array([1,1,1,0,0,0])
    I_vT = np.array([[1],[1],[1],[0],[0],[0]])
    I_s = np.diag([1,1,1,0.5,0.5,0.5])
    I_sdev = (I_s - 1/3*I_v*I_vT)
    
       
    if eps_h >= 0.5*eps_y:
        k = max([min([1-0.5*eps_y/eps_h,1]),0])
        if k > 0.5:
            k2 = max([min([1-eps_y/eps_h,1]),0])
            eps1 = eps*(1-k)*(1-k2)
            eps2 = eps*k*(1-k2)
            eps3 = eps*k2
            eps_dev1 = np.matmul(I_sdev,eps1)
            eps_dev2 = np.matmul(I_sdev,eps2)
            eps_dev3 = np.matmul(I_sdev,eps3)
        else:
            k2 =0
            eps1 = eps*(1-k)*(1-k2)
            eps2 = eps*k*(1-k2)
            eps3 = eps*k2
            eps_dev1 = np.matmul(I_sdev,eps1)
            eps_dev2 = np.matmul(I_sdev,eps2)
            eps_dev3 = np.matmul(I_sdev,eps3)
    else:
        k,k2=0,0
        eps1 = eps*(1-k)*(1-k2)
        eps2 = eps*k*(1-k2)
        eps3 = eps*k2
        eps_dev1 = np.matmul(I_sdev,eps1)
        eps_dev2 = np.matmul(I_sdev,eps2)
        eps_dev3 = np.matmul(I_sdev,eps3)
        
    sigma = 2*G1*eps_dev1 + K1*np.matmul(I_v*I_vT,eps1) + 2*G2*eps_dev2 + K2*np.matmul(I_v*I_vT,eps2) + 2*G3*eps_dev3 + K3*np.matmul(I_v*I_vT,eps3)
    
    return sigma





def Dfun(E,nu):
    return E/((1+nu)*(1-2*nu))*np.array([[1-nu,nu,nu,0],[nu,1-nu,nu,0],[nu,nu,1-nu,0],[0,0,0,(1-2*nu)/2]])

def numD(eps,sig,mp):
    delta = 1e-7
    D = np.zeros([6,6])
    for i in range(0,6):
        eps2 = eps.copy()
        eps2[i] = eps[i] + delta
        sig2 = mat(eps2,mp)
        
        D[i,:] = (sig2-sig)/delta
    return D
    


