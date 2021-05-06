"""
This is an example of a Clamped Beam topology optimisation. The geometry is
rectangular with a clamped boundary on the left and right side and a downwards force
in the middle on the top side. The example utilises the program PyTOpt written
by Daniel Pettersson and Erik Säterskog

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
             Bilinear/Elastic/ModifiedHooke
ObjectivFun- Determines which objective function that should be used.
             Compliance/Displacement
Optfun     - Determines which optimisation algorithm that should be used.
             OC/MMA

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
g.point([1,0])                 
g.point([1,0.4])                             
g.point([0.5,0.4],marker=9)    
g.point([0,0.4])               

g.line([0, 1],marker=0)
g.line([1, 2],marker=1)
g.line([2, 3],marker=2)
g.line([3, 4],marker=3)
g.line([4, 0],marker=4)

g.surface([0, 1, 2, 3, 4])

# Forces and boundary conditions
force = [-4e7,9,2]                                  
bmarker = [1,4]                                    
eq = [0,0]                                           

# Material parameters
mp = {'E':210e9,'nu':0.3,'eps_y':0 }                
materialFun = mrs.Bilinear                         

# Settings, Objective function and Optimisation routine
ep=[1,True,2]
settings = {'volFrac':0.3,'meshSize':0.08,'rmin':0.08*0.7,'changeLimit': 0.01,'SIMP_const':3}
ObjectFun = ofs.Compliance
OptFun = opt.OC

# Calling the optimisation
PyTOpt.Main(g, force, bmarker, mp, ep, materialFun, ObjectFun, OptFun, settings, eq)

