"""
Bilinear material model.
Inputs:
    eps     -epsilon, the strain
    mp[E,nu,eps_y]
             E          - Young's modulus
             nu         - Poission's ratio
             eps_y      - Yielding strain for bilinear material model

Outputs:
    sigma   -Stress
    D       -Constitutive matrix


Written 2021-05
Made By: Daniel Pettersson & Erik Säterskog
"""
import numpy as np

def Bilinear(eps,mp):
    
    E1       = mp['E']
    nu1      = mp['nu']
    eps_y    = mp['eps_y']
    
    E2 = E1*0.05
   
    G = E1/(2*(1+nu1))
    K1 = E1*G/(3*(3*G-E1))
    K2 = E2*G/(3*(3*G-E2))
    
    eps_h = sum(eps[:2])/3
    
    I_v = np.array([[1,1,1,0,0,0]])
    I_vT = np.array([[1],[1],[1],[0],[0],[0]])
    I_s = np.diag([1,1,1,0.5,0.5,0.5])
    I_sdev = (I_s - 1/3*I_v*I_vT)
    
    
    if eps_h > eps_y:
        k = 1-eps_y/eps_h
        dkdeps = 3*eps_y/(3*eps_h)**2
        
    else:
        k = 0
        dkdeps = 0
        
    eps1 = eps*(1-k)
    eps2 = eps*k
    eps_dev1 = np.matmul(I_sdev,eps1)
    eps_dev2 = np.matmul(I_sdev,eps2)
    
        
    sigma = 2*G*eps_dev1 + K1*np.matmul(I_v*I_vT,eps1) + 2*G*eps_dev2 + K2*np.matmul(I_v*I_vT,eps2)

    E1_v = (2*G*I_sdev + K1*I_v*I_vT)
    E2_v = (2*G*I_sdev + K2*I_v*I_vT)
    
    D = (E1_v*(1-k)+E2_v*k
         -np.matmul(E1_v,(eps.reshape(6,1)*(I_v*dkdeps)))
         +np.matmul(E2_v,(eps.reshape(6,1)*(I_v*dkdeps))))
    
    return sigma.reshape(6,1),D


    


