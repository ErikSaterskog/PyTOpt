# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:23:12 2021

@author: Daniel
"""

import numpy as np
import calfem.core as cfc

def mat(eps,mp):
    
    E_ten      = mp[0]
    E_com       = mp[0]
    nu_ten      = 0.3
    nu_com      = 0.3
    low_lim=0.0005
    high_lim=0.0025
    sig_y = 220e6
    G = E_com/(2*(1+nu_com))
    K_com = E_com/(3*(1-2*nu_com))
    K_ten = E_ten/(3*(1-2*nu_ten))
    
    eps_h = sum(eps[:2])/3
    
    I_v = np.array([1,1,1,0,0,0])
    I_vT = np.array([[1],[1],[1],[0],[0],[0]])
    I_s = np.diag([1,1,1,0.5,0.5,0.5])
    I_sdev = (I_s - 1/3*I_v*I_vT)
    
    eps_dev = np.matmul(I_sdev,eps)
    eps_vM = np.sqrt(2/3)*np.sqrt(np.matmul(eps.T,np.matmul(I_sdev,eps)))
    
    eps_h = sum(eps[:2])/3

        
    if eps_h >= sig_y/(3*G): 
        print('Comp')
        Gstar = sig_y/(3*eps_h)
        #eps_hy = np.zeros(np.shape(eps))
        #eps_hy[np.where(eps>0)] = sig_y/(3*G)
        sigma_ten = 2*Gstar*eps_dev + K_ten*np.matmul(I_v*I_vT,eps)*(sig_y/(3*G*eps_h))
        sigma_com = eps*0
        D = 2*Gstar*I_sdev - 4*Gstar/(3*eps_h**2)*np.matmul(eps_dev.reshape(6,1),eps_dev.reshape(1,6)) + K_ten*I_v*I_vT*(sig_y/(3*G*eps_h))
    else:
        print('ten')
        ten_quote = 0
        sigma_com = 2*G*eps_dev + K_ten*np.matmul(I_v*I_vT,eps)
        sigma_ten = eps*0
        D = Dfun(E_com,nu_com)


    sigma = sigma_com + sigma_ten
    
    return sigma.reshape(6,1), D





def Dfun(E,nu):
    return E/((1+nu)*(1-2*nu))*np.array([[1-nu,nu,nu,0],[nu,1-nu,nu,0],[nu,nu,1-nu,0],[0,0,0,(1-2*nu)/2]])




