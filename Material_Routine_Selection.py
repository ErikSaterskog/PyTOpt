# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 10:25:11 2021

@author: Daniel
"""



import numpy as np
import Material_Elastic
import Material_Bilinear
import Material_ModifiedHooke

def Elastic(eps,mp):
    return Material_Elastic.elastic(eps,mp) 

def ModifiedHooke(eps,mp):
    return Material_ModifiedHooke.mod_hooke(eps,mp)

def Bilinear(eps,mp):
    return Material_Bilinear.Bilinear(eps,mp)

def numD(eps,sig,mp,mat):
    delta = 1e-7
    D = np.zeros([6,6])
    for i in range(0,6):
        eps2 = eps.copy()
        eps2[i] = eps[i] + delta
        sig2 = mat(eps2,mp)[0]
        
        D[i,:] = (sig2-sig)/delta
    return D
    

