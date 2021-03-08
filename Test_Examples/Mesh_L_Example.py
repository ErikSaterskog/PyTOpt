# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:52:07 2021

@author: Daniel
"""


#Really simple test
import numpy as np
import calfem.geometry as cfg
import calfem.vis as cfv
import Main
import elastic as el


g = cfg.Geometry()

g.point([0,0])                  #0
g.point([1,0])                  #1
g.point([1,0.5],marker=9)       #2
g.point([0.5,0.5])              #3
g.point([0.5,1])                #4
g.point([0,1])                  #5

g.spline([0, 1])
g.spline([1, 2])
g.spline([2, 3])
g.spline([3, 4])
g.spline([4, 5],marker=4)
g.spline([5, 0])

g.surface([0, 1, 2, 3,4,5])

force = [-1e7,9,2] #First magnitude, second marker,third direction
bmarker = 4

E = 210e9 # Young's modulus
nu = 0.3 #Poisson's ratio

mp = [E,nu]

volFrac = 0.3 # Constraint on 50% volume
meshSize=0.04 # The average length of one element. 
rMin = meshSize*np.sqrt(2)*0.5 # How aggressive the filter should be. Smaller -> less aggressive
changeLimit=0.01 # How small change between two optmisation we allow before stopping.
el_type = 2   #2-Tri,  3-Quad


ep=[2,1,2,2]    #ep[ptype, thickness, integration rule(only used for QUAD),linear(1)/nonlinear(2)]  
SIMP_penal = 3
method='OC'
Debug=False


settings = [volFrac,meshSize, rMin, changeLimit, SIMP_penal, method, Debug]


#sick.sick(epsilon,mp)
#mh.mod_hook(epsilon, mp)
materialFun=el.elastic

Main.Main(g, el_type, force, bmarker, settings, mp, ep, materialFun)







