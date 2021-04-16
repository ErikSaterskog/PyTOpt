# Importing Modules
import calfem.geometry as cfg
import Pytopt.PyTOpt as PyTOpt
from Pytopt import Material_Routine_Selection as mrs
from Pytopt import Object_Func_Selection as ofs
#####################

g = cfg.Geometry()

g.point([0,0])   
g.point([1,0.5],marker=1)
g.point([0,1])


g.line([0, 1])
g.line([1, 2])
g.line([2, 0],marker=2)


g.surface([0 ,1, 2])

force = [-1e5, 1, 2]
bmarker = 2
eq = [0,0]

volFrac = 0.3 
meshSize=0.1
rMin = meshSize*0.7 
changeLimit=0.01 
SIMP_penal = 3
method='OC'
Debug=False
E = 210e9 
nu = 0.3 
mp = [E,nu]

ep=[1,True,2]    #ep[thickness, linear(True)/nonlinear(False),2-Tri,  3-Quad]  

settings = [volFrac,meshSize, rMin, changeLimit, SIMP_penal, method, Debug]

# Material model and Objective function
materialFun = mrs.Elastic
ObjectFun = ofs.Energy
#######################

# Calling the optimisation
PyTOpt.Main(g, force, bmarker, settings, mp, ep, materialFun, ObjectFun, eq)
#######################











