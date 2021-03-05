# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:23:12 2021

@author: Daniel
"""

import numpy as np
#import calfem.core as cfc

def sick(eps,mp):
    
    E_ten      = mp[0]
    E_com       = E_ten*7
    nu_ten      = 0.3
    nu_com      = 0.3
    # lame_max    = E_max*nu_max/((1+nu_max)*(1-2*nu_max))
    # lame_min    = E_min*nu_min/((1+nu_min)*(1-2*nu_min))
    # my_max      = E_max/(2*(1+nu_max)) 
    # my_min      = E_min/(2*(1+nu_min))
    # sig=np.zeros([3,])
    
    epsilon = np.matrix([[eps[0],eps[3]/2,eps[5]/2],[eps[3]/2,eps[1],eps[4]/2],[eps[5]/2,eps[4]/2,eps[2]]])
    
    lamd,nor  = np.linalg.eig(epsilon)
    eps_h = sum(eps[:2])/3
    if eps_h>0.0005: 
        print('ten')
        D = E_ten/((1+nu_ten)*(1-2*nu_ten))*np.array([[1-nu_ten,nu_ten,nu_ten,0],[nu_ten,1-nu_ten,nu_ten,0],[nu_ten,nu_ten,1-nu_ten,0],[0,0,0,(1-2*nu_ten)/2]])
        
    else:
        print('Comp')
        D = E_com/((1+nu_com)*(1-2*nu_com))*np.array([[1-nu_com,nu_com,nu_com,0],[nu_com,1-nu_com,nu_com,0],[nu_com,nu_com,1-nu_com,0],[0,0,0,(1-2*nu_com)/2]])
    
    #sig = D_1*epskvot + D_2*epsrest Detta ska implementeras kanske
        
    sigma = np.zeros(np.shape(eps))
    sigma[:4] = np.matmul(D,eps[:4])
    return sigma.reshape(6,1), D


