"""
The example starts by importing the necessary libaries. After that is the 
geometry constructed. The marker indicates an ID for the point or line. 
This ID is later used to define boundary conditions and prescribed forces. 
Change these if you are sure about how they works. 

The force vector is consisting of three values. The first value is the applied
force magnitude and can be altered by the user. The second inicated where the load 
is applied and the last in what direction the force is pointing. A value of 2 
indicates in y-direction. b-marker inidicates what line should be prescirbed.
E, nu and eps_y are material parameters and can be altered after what material 
of interest.

VolFrac    - How many procent of the volume should the final solution have contra
             the original.
meshsize   - How large the average element should be. 
rMin       - For the filtering, how large radius should the filter take into 
             account. With a low value the filter will only notice the closest 
             neighbour for each element. A large value will make the filter take 
             a large elements into account when filtering each element.
changeLimit- For OC as optimisation method, what tolerance for the change 
             between iteration is sufficiant.

ep         -
             t          - thickness 
             linear     - True-linear, False-nonlinear
             el_type    - 2 means triangular elements and 3 means quad elements.
method     - 
             OC - Optimal Criterion methon
             MMA - Method of moving asymptotes
debug      - True/False if the sensitivity should be checked numerically.
materialFun- Determine which material model that should be used. The user can 
             add ones own material model as long as the input is a strain
             and a material parametervector 
ObjectivFun- Determine which objective function that should be used. The user can 
             add ones own objective function as long as the input is the same
             as for the already exisitng objective funtions.

Then we call on the Main module to start the optimisation.
"""


# Importing Modules
import calfem.geometry as cfg
import PyTOpt
import Material_Routine_Selection as mrs
import Object_Func_Selection as ofs
#################

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
g.point([100,50])
g.point([0,50])

g.line([0, 1],marker=0)
g.line([1, 2],marker=1)
g.line([2, 3],marker=2)
g.line([3, 4],marker=3)
g.line([4, 5],marker=4)
g.line([5, 6],marker=5)
g.line([6,8],marker=6)
g.line([8,9],marker=7)
g.line([9,7],marker=8)
g.line([7, 0],marker=9)
g.line([6,7],marker=10)

g.surface([0,1,2,3,4,5,10,9])
g.surface([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
#######################

# Forces and boundary conditions
force = [-7e4,10,2] #First magnitude, second marker, third direction
bmarker = [5,9,7]
eq = [0,-9.81*7750]
#######################

# Material parameters
E = 210e9       # Young's modulus
nu = 0.3        # Poisson's ratio
eps_y = 1e-6    # Strain border for Bilinear material model
mp = [E,nu,eps_y]
######################

# Settings
volFrac = 0.3       
meshSize=2       
rMin = meshSize*0.7 
changeLimit=0.01    
ep=[1,True,2]       
SIMP_penal = 3
method='OC'
Debug=False
settings = [volFrac,meshSize, rMin, changeLimit, SIMP_penal, method, Debug]
######################

# Material model and Objective function
materialFun = mrs.Elastic
ObjectFun = ofs.Energy
######################

# Calling the optimisation
PyTOpt.Main(g, force, bmarker, settings, mp, ep, materialFun, ObjectFun, eq)
######################