"""
This is an example of a bridge topology optimisation. The example utilises the program
PyTOpt written by Daniel Pettersson and Erik Säterskog

Settings explanation:
force      - [Magnitude, Force Marker, Direction]
bmarker    - Boundary marker
eq         - Body forces,  [0, -9.81*7750] for steel on earth
mp[E,nu,eps_y]
             E          - Young's modulus
             nu         - Poission's ratio
             eps_y      - Yielding strain for bilinear material model
ep[thickness, linear, el_type]
             thickness  - thickness of the 2D material
             linear     - True-linear, False-nonlinear
             el_type    - 2 indicates triangular elements and 3 indicates
             quad elements.
Settings[volfrac, meshSize, rmin, changelimit, SIMP_const
             VolFrac    - Fraction of max material allowed in optimisation.
             meshsize   - Avarage element size. [m]
             rMin       - Filter radius. Large values gives smudged solutions. [m]
             changeLimit- Determines at what change the optimisation breaks
             SIMP_const - Solid isotropic material penalisation constant
materialFun- Determines which material model that should be used.
ObjectivFun- Determines which objective function that should be used.
Optfun     - Determines which optimisation algorithm that should be used.


Written 2021-05
Made By: Daniel Pettersson & Erik Säterskog
"""

# Importing Modules
import calfem.geometry as cfg
import Pytopt.PyTOpt as PyTOpt
from Pytopt import Material_Routine_Selection as mrs
from Pytopt import Object_Func_Selection as ofs
from Pytopt import Optimisation as opt

# Creating geometry
g = cfg.Geometry()

g.point([0,0])                 
g.point([30,0])
g.point([35,0])
g.point([65,0])
g.point([70,0])
g.point([100,0])
g.point([100,25])                 
g.point([0,25])               

g.line([0, 1],marker=0)
g.line([1, 2],marker=1)
g.line([2, 3],marker=2)
g.line([3, 4],marker=3)
g.line([4, 5],marker=4)
g.line([5, 6],marker=5)
g.line([6, 7],marker=6)
g.line([7, 0],marker=7)

g.surface([0, 1, 2, 3, 4, 5, 6, 7])

# Forces and boundary conditions
force = [7e4,6,2]                               
bmarker = [1,3,5,7]                             
eq = [0, 0]

# Material parameters
mp = {'E':210e9,'nu':0.3,'eps_y':0 }            
materialFun = mrs.Bilinear                      

# Settings, Objective function and Optimisation routine
ep=[1,True,2]
settings = {'volFrac':0.3,'meshSize':2,'rmin':2*0.7,'changeLimit': 0.01,'SIMP_const':3}
ObjectFun = ofs.Energy
OptFun = opt.OC

# Calling the optimisation
PyTOpt.Main(g, force, bmarker, mp, ep, materialFun, ObjectFun, OptFun,settings,eq)

