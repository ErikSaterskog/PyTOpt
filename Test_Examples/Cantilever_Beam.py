# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:52:07 2021

The example starts by importing the necessary libaries. After that is the geometry constructed.
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
ep         -
             t         - thickness 
             linear    - True if linear, else false
             el_type    - 2 means triangular elements and 3 means quad elements.
method     - 
             OC - Optimal Criterion methon
             MMA - Method of moving asymptotes
debug      - True/False if the sensitivity should be checked numerically.
materialFun- Determine which material model that should be used. The user can 
             add ones own material model as long as the input is a strain
             and a material parametervector 

Then we call on the Main module to start the optimisation.
"""

import calfem.geometry as cfg
import PyTOpt
import Material_Routine_Selection as mrs
import Object_Func_Selection as ofs
g = cfg.Geometry()

g.point([0,0])                 #0
g.point([1,0])                 #1
g.point([1,0.4],marker=9)      #2               
g.point([1,0.8])               #3
g.point([0,0.8])               #4
g.point([0,0.5])               #5
g.point([0,0.3])               #6


g.line([0, 1],marker=0)
g.line([1, 2],marker=1)
g.line([2, 3],marker=2)
g.line([3, 4],marker=3)
g.line([4, 5],marker=4)
g.line([5, 6],marker=5)
g.line([6, 0],marker=6)


g.surface([0, 1, 2, 3, 4,5,6])

force = [-1e6,9,2]      #First magnitude, second marker, third direction
bmarker = 5
eq=[0,0]


E = 210e9               # Young's modulus
nu = 0.3                #Poisson's ratio
eps_y = 0


mp = [E,nu,eps_y]

volFrac = 0.3           # Constraint on volume
meshSize=0.1            # The average length of one element. 
rMin = meshSize*0.7     # How aggressive the filter should be. Smaller -> less aggressive
changeLimit=0.0        # How small change between two optmisation we allow before stopping.
ep=[1,True,2]          #ep[thickness, linear(True)/nonlinear(False),2-Tri,  3-Quad]  
SIMP_penal = 3
method='OC'
Debug=False

settings = [volFrac,meshSize, rMin, changeLimit, SIMP_penal, method, Debug]

materialFun = mrs.Bilinear
ObjectFun = ofs.Energy

PyTOpt.Main(g, force, bmarker, settings, mp, ep, materialFun, ObjectFun, eq)














