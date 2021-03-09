# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:23:12 2021

@author: Daniel
"""

import numpy as np
import calfem.core as cfc

def mat(eps,mp):
    
   
    E       = mp[0]
    nu      = 0.3
    low_lim=0.0005
    high_lim=0.0025
    sig_y = 220e6
    G = E/(2*(1+nu))
    K = E/(3*(1-2*nu))
    
    eps_h = sum(eps[:2])/3
    
    I_v = np.array([1,1,1,0,0,0])
    I_vT = np.array([[1],[1],[1],[0],[0],[0]])
    I_s = np.diag([1,1,1,0.5,0.5,0.5])
    I_sdev = (I_s - 1/3*I_v*I_vT)
    
    eps_dev = np.matmul(I_sdev,eps)
    eps_vM = np.sqrt(2/3)*np.sqrt(np.matmul(eps.T,np.matmul(I_sdev,eps)))
    
        
    if eps_h >= sig_y/(3*G): 
        Gstar = sig_y/(3*eps_h)
        sigma = 2*Gstar*eps_dev + K*np.matmul(I_v*I_vT,eps)
        Estar = 2*Gstar*(1+nu)
        D = Dfun(Estar,nu)
        #D = 2*Gstar*I_sdev - 4*sig_y/(9*eps_h**3)*np.matmul(eps_dev.reshape(6,1),eps_dev.reshape(1,6)) + K*I_v*I_vT
    else:
        sigma = 2*G*eps_dev + K*np.matmul(I_v*I_vT,eps)
        D = Dfun(E,nu)


    
    
    return sigma.reshape(6,1), D





def Dfun(E,nu):
    return E/((1+nu)*(1-2*nu))*np.array([[1-nu,nu,nu,0],[nu,1-nu,nu,0],[nu,nu,1-nu,0],[0,0,0,(1-2*nu)/2]])




