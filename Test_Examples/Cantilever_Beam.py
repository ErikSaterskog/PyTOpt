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

el_type = 2   #2-Tri,  3-Quad
g = cfg.Geometry()

g.point([0,0])                  #0
g.point([1,0])                  #1
g.point([1,0.25],marker=9)      #2
g.point([0,0.25])               #3

g.spline([0, 1],marker=0)
g.spline([1, 2],marker=1)
g.spline([2, 3],marker=2)
g.spline([3, 0],marker=3)

g.surface([0, 1, 2, 3])

force = [-1e6,9,2] #First magnutude, second marker,third direction
bmarker = 3

Main._Main(g,el_type,force,bmarker)










