# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:29:40 2021

This script shows the different between the modified nonlinear materialmodel
and ordinary linear hookes law.

Just run. You could ajust the strain(eps11) by increasing or decreasing the
second input to np.linespace.


"""

import numpy as np
import matplotlib.pyplot as plt
import Mod_Hook as MH
import calfem.core as cfc

E = 210e9
nu= 0.3

mp = [E,nu,0]

D_lin = cfc.hooke(2,E,nu)

eps11 = np.linspace(0,0.01,1000)

eps = np.zeros([6,1000])
sigma = eps.copy()
sig_Lin = eps.copy()
eps[0,:] = eps11
for i in range(1,1000):
    sigma[:,i],D = MH._mod_hook(eps[:,i],mp)
    sig_Lin[0,i] = D_lin[0,0] * eps[0,i]

plt.figure()
plt.plot(eps[0,:],sigma[0,:])

plt.plot(eps[0,:],sig_Lin[0,:])
