# Importing Modules
import calfem.geometry as cfg
import Pytopt.PyTOpt as PyTOpt
from Pytopt import Material_Routine_Selection as mrs
from Pytopt import Object_Func_Selection as ofs
from Pytopt import Optimisation as opt
#####################

# Creating geometry
g = cfg.Geometry()

g.point([0,0])   
g.point([1,0.49])
g.point([1,0.5],marker=1)
g.point([0,1])

g.line([0, 1])
g.line([1, 2])
g.line([2, 3])
g.line([3, 0],marker=2)

g.surface([0 ,1, 2, 3])
#####################
# Forces and boundary conditions
force = [-1e5, 1, 2]
bmarker = 2
eq = [0,0]
#####################
# Material parameters
E = 210e9               # Young's modulus
nu = 0.3                # Poisson's ratio
eps_y = 0               # Strain border for Bilinear material model
mp = {'E':E,'nu':nu,'eps_y':eps_y}
materialFun = mrs.Elastic # Material model
#####################
# Settings, Objective function and Optimisation routine
ep=[1,True,3]
settings = {'volFrac':0.3,'meshSize':0.08,'rmin':0.08*0.7,'changeLimit': 0.01,'SIMP_penal':3}
ObjectFun = ofs.Energy
OptFun = opt.MMA
#####################
# Calling the optimisation
PyTOpt.Main(g, force, bmarker, mp, ep, materialFun, ObjectFun, OptFun,settings,eq)
#######################









