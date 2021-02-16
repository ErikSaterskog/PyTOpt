# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:52:07 2021

The example start by import the necessary libaries. After that is the geometry constructed.
Do not change these. 

The marker indicates an ID for the point or line. This ID is later used to 
define boundary conditions and prescribed forces. Change these if you are sure
about how they works and when you want a diffrerent problem. 

The force vector is consisting of three values. The first value is the applied
force magnitude and can be altered by the user. However, the later two should not 
be touched. As the first of them inicated where the load is applied and the last
in what direction the force is pointing. A value of 2 indicates in y-direction.

b-marker inidicates what line should be prescirbed.

E and nu are material parameters and can be altered after what material of interest.

VolFrac    - How many procent of the volume should the final solution have contra
             the original.
meshsize   - How "large" each element should be. This is a scale value between 0
             and 1 where 1 means few and large elements and 0 means infinite many
             elements.
rMin       - For the filtering, how large radius should the filter take into 
             account. With a low value the filter will only notice the closest 
             neighbour for each element. A large value will make the filter take 
             a large elements into account when filtering each element.
changeLimit- For OC as optimisation method, what tolerance for the change 
             between iteration is sufficiant.
el_type    - 2 means triangular elements and 3 means quad elements.
Linear     - True or False


Then we call on the Main module to start the optimisation.


"""


#Really simple test
import numpy as np
import calfem.geometry as cfg
import calfem.vis as cfv
import Main

g = cfg.Geometry()

g.point([0,0])                 #0
g.point([1,0])        #1
g.point([1,0.4],marker=9)               #2               #2
g.point([1,0.8])               #2
g.point([0,0.8])               #3
g.point([0,0.5])               #3
g.point([0,0.3])               #3


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

volFrac = 0.4 # Constraint on 50% volume
meshSize=0.02 # The average length of one element. 
rMin = meshSize*np.sqrt(2)*0.5 # How aggressive the filter should be. Smaller -> less aggressive
changeLimit=0.005 # How small change between two optmisation we allow before stopping.
el_type = 2   #2-Tri,  3-Quad
Linear = False


settings = [volFrac,meshSize,rMin,changeLimit,Linear]


Main._Main(g,el_type,force,bmarker,settings,mp)










