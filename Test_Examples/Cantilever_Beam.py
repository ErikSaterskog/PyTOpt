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

g = cfg.Geometry()

g.point([0,0])                 #0
g.point([1,0])        #1
g.point([1,0.5],marker=9)               #2               #2
g.point([1,1])               #2
g.point([0,1])               #3
g.point([0,0.6])               #3
g.point([0,0.4])               #3


g.spline([0, 1],marker=0)
g.spline([1, 2],marker=1)
g.spline([2, 3],marker=2)
g.spline([3, 4],marker=3)
g.spline([4, 5],marker=4)
g.spline([5, 6],marker=5)
g.spline([6, 0],marker=6)


g.surface([0, 1, 2, 3, 4,5,6])

force = [-1e8,9,2] #First magnitude, second marker, third direction
bmarker = 5


E = 210e9 # Young's modulus
nu = 0.3 #Poisson's ratio

mp = [E,nu]

volFrac = 0.2 # Constraint on 50% volume
meshSize=0.04 # How fine mesh we want. 1 is only one element and 0 is infinity.
rMin = meshSize*np.sqrt(2)*0.5 # How aggressive the filter should be. Smaller -> less aggressive
changeLimit=0.005 # How small change between two optmisation we allow before stopping.
el_type = 3   #2-Tri,  3-Quad



settings = [volFrac,meshSize,rMin,changeLimit]


Main._Main(g,el_type,force,bmarker,settings,mp)










