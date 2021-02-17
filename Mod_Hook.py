# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:04:30 2021

@author: Daniel
"""
import numpy as np

def _mod_hook(eps,mp):
    E = mp[0]
    nu = mp[1]
    G = np.zeros([2,1])
    G[0] = E/(2*(1+nu))
    G[1] = -G[0].copy()*0.1
    K = E/(3*(1-2*nu))
    
    I_v = np.array([1,1,1,0,0,0])
    I_vT = np.array([[1],[1],[1],[0],[0],[0]])
    I_s = np.diag([1,1,1,0.5,0.5,0.5])
    I_sdev = (I_s - 1/3*I_v*I_vT)
    
    eps_dev = np.matmul(I_sdev,eps)
    
    
    sigma = np.zeros([6,])
    D = np.zeros([6,6])
    Eps_sum = sum(eps_dev*eps_dev)
    
    for i in range(0,2):
        sigma = sigma + 2*G[i]*eps_dev*(1e4*Eps_sum)**i
        D = D + 2*G[i]*I_sdev*(1e4*Eps_sum)**i
    sigma = sigma + K*np.matmul(I_v*I_vT,eps)
    D  = D + K*I_v*I_vT
    
    sigma.shape = (6,1)
    return sigma, D
