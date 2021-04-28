
import numpy as np

def mod_hooke(eps,mp):
    if max(abs(eps)) > 0.3:
        raise Exception('ERROR ERROR Part is breaking')
    else:
        pass
    E = mp['E']
    nu = mp['nu']
    G = np.zeros([3,1])
    G[0] = E/(2*(1+nu))
    G[1] = -G[0].copy()*10
    G[2] = G[0].copy()
    K = E/(3*(1-2*nu))
    
    I_v = np.array([1,1,1,0,0,0])
    I_vT = np.array([[1],[1],[1],[0],[0],[0]])
    I_s = np.diag([1,1,1,0.5,0.5,0.5])
    I_sdev = (I_s - 1/3*I_v*I_vT)
    
    eps_dev = np.matmul(I_sdev,eps)
    
    
    sigma = np.zeros([6,])
    D = np.zeros([6,6])
    Eps_sum = sum(eps_dev*eps_dev)
    
    for i in range(0,3):
        sigma = sigma + 2*G[i]*eps_dev*(Eps_sum)**i
        D = D + 2*G[i]*I_sdev*(Eps_sum)**i
    sigma = sigma + K*np.matmul(I_v*I_vT,eps)
    D  = D + K*I_v*I_vT
    
    sigma.shape = (6,1)
    return sigma, D
