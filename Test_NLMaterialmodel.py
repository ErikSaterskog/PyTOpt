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
import Material_Bilinear as MB

E = 210e9
nu= 0.3

mp = [E,nu,0.002]
numEval=1000

eps11 = np.linspace(-0.002,0.005,numEval)

eps = np.zeros([6,numEval])
sigma = eps.copy()
eps[0,:] = eps11
for i in range(1,numEval):
    
    sigma_temp,D = MB.head(eps[:,i],mp)
    sigma[:,i]=sigma_temp.reshape(6)

plt.figure()
plt.plot(eps[0,1:],sigma[0,1:])

