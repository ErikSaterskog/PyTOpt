# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 10:25:11 2021

@author: Daniel
"""



import numpy as np

def elastic(eps,mp):
    E = mp[0]
    nu = mp[1]
    G = E/(2*(1+nu))
    K = E/(3*(1-2*nu))
    
    I_v = np.array([1,1,1,0,0,0])
    I_vT = np.array([[1],[1],[1],[0],[0],[0]])
    I_s = np.diag([1,1,1,0.5,0.5,0.5])
    I_sdev = (I_s - 1/3*I_v*I_vT)
        
    sigma = np.zeros([6,])
    D = np.zeros([6,6])
        
    D  =K*I_v*I_vT+ 2*G*I_sdev
    sigma = np.matmul(D,eps)
    sigma.shape = (6,1)
    return sigma, D


def mod_hooke(eps,mp):
    if max(abs(eps)) > 0.3:
        raise Exception('ERROR ERROR Part is breaking')
    else:
        pass
    E = mp[0]
    nu = mp[1]
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

def Bilinear(eps,mp):
    
   
    E1       = mp[0]
    nu1      = mp[1]
    eps_y    = mp[2]
    
    E2 = E1*0.1
    
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

def numD(eps,sig,mp,mat):
    delta = 1e-7
    D = np.zeros([6,6])
    for i in range(0,6):
        eps2 = eps.copy()
        eps2[i] = eps[i] + delta
        sig2 = mat(eps2,mp)[0]
        
        D[i,:] = (sig2-sig)/delta
    return D
    

