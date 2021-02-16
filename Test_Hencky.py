# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:29:40 2021

@author: Daniel
"""

import numpy as np
import matplotlib.pyplot as plt
import Mod_Hook as MH

E = 210e9
nu= 0.3

mp = [E,nu,0]


eps11 = np.linspace(0,0.01,100)

eps = np.zeros([6,100])
sigma = eps.copy()
eps[1,:] = eps11
for i in range(1,100):
    sigma[:,i],D = MH._mod_hook(eps[:,i],mp);

plt.figure()
plt.plot(eps[1,:],sigma[1,:])

plt.plot([eps[1,1],eps[1,-1]],[sigma[1,1],sigma[1,-1]])
