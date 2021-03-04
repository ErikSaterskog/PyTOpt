# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:23:12 2021

@author: Daniel
"""

import numpy as np
import calfem.core as cfc

def sick(eps,mp):
    
    E_max       = mp[0]
    E_min       = E_max*7/6
    nu_max      = 0.3
    nu_min      = 0.3
    lame_max    = E_max*nu_max/((1+nu_max)*(1-2*nu_max))
    lame_min    = E_min*nu_min/((1+nu_min)*(1-2*nu_min))
    my_max      = E_max/(2*(1+nu_max)) 
    my_min      = E_min/(2*(1+nu_min))
    sig=np.zeros([3,])
    
    epsilon = np.matrix([[eps[0],eps[3]/2,eps[5]/2],[eps[3]/2,eps[1],eps[4]/2],[eps[5]/2,eps[4]/2,eps[2]]])
    
    lamd,nor  = np.linalg.eig(epsilon)
    
    if lamd[0]+nu_max*lamd[1]>0 and lamd[1]+nu_max*lamd[0]>0:
        
        D = E_max/((1+nu_max)*(1-2*nu_max))*np.array([[1-nu_max,nu_max,nu_max,0],[nu_max,1-nu_max,nu_max,0],[nu_max,nu_max,1-nu_max,0],[0,0,0,(1-2*nu_max)/2]])
        
    elif lamd[0]+nu_min*lamd[1]<=0 and lamd[1]+nu_min*lamd[0]<=0:
        
        D = E_min/((1+nu_min)*(1-2*nu_min))*np.array([[1-nu_min,nu_min,nu_min,0],[nu_min,1-nu_min,nu_min,0],[nu_min,nu_min,1-nu_min,0],[0,0,0,(1-2*nu_min)/2]])
        
    elif lamd[0]+nu_min*lamd[1]>0 and lamd[1]+nu_max*lamd[0]<=0:
        
        D = E_max/((1+nu_max)*(1-2*nu_max))*np.array([[1-nu_min,nu_min,nu_min,0],[nu_min,1-nu_min,nu_min,0],[nu_min,nu_min,1-nu_min,0],[0,0,0,(1-2*nu_min)/2]])
        
    elif lamd[0]+nu_max*lamd[1]<=0 and lamd[1]+nu_min*lamd[0]>0:
        
        D = E_min/((1+nu_min)*(1-2*nu_min))*np.array([[1-nu_max,nu_max,nu_max,0],[nu_max,1-nu_max,nu_max,0],[nu_max,nu_max,1-nu_max,0],[0,0,0,(1-2*nu_max)/2]])
        
    
    for i in range(0,3):
        lam = np.ma.array(lamd, mask=False)
        lam.mask[i] = True
        if lamd[i]>0:
                
            sig[i] = lamd[i]*(2*my_max+lame_max)+lame_max*lam.sum()
        else:
            sig[i] = lamd[i]*(2*my_min+lame_min)+lame_min*lam.sum()
        
    sigma_mat = sig[0]*nor[:,0]*nor[:,0].T +sig[1]*nor[:,1]*nor[:,1].T+sig[2]*nor[:,2]*nor[:,2].T
    sigma = np.array([sigma_mat[0,0],sigma_mat[1,1],sigma_mat[2,2],sigma_mat[0,1],sigma_mat[1,2],sigma_mat[0,2]])        
    
    
    return sigma.reshape(6,1), D


