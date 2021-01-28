# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:52:07 2021

@author: Daniel
"""


#Really simple test
import numpy as np
import calfem.geometry as cfg
import calfem.vis as cfv
import Mesh



el_type = 3
g = cfg.Geometry()

g.point([0,0])          #0
g.point([1,0])          #1
g.point([1,1],marker=9)       #2
g.point([0,1])    #3

g.spline([0, 1],marker=0)
g.spline([1, 2],marker=1)
g.spline([2, 3],marker=2)          
g.spline([3, 0],marker=3)



g.surface([0, 1, 2, 3])

force = [-1000,9,2] #First magnutude, second marker,third direction
bmarker = 3

_mesh = Mesh.Mesh(g,0.25)

coords, edof, dofs, bdofs = _mesh.quad()



cfv.drawMesh(coords, edof, 2, el_type)
cfv.showAndWait()
#Main.mainer(g,el_type,force,bmarker)










