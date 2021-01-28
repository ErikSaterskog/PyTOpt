"""
This file will call an example, mesh it and solve the optimisation
"""

import numpy as np
import calfem.utils as cfu
import Mesh
import calfem.vis as cfv
import FE
import Opt
import Filter
import matplotlib.pyplot as plt

def _Main(g,el_type,force,bmarker):
     
    """MESHING"""
    _mesh = Mesh.Mesh(g,0.05)
    
    if el_type == 2:
        coords, edof, dofs, bdofs = _mesh.tri()
    elif el_type ==3:
        coords, edof, dofs, bdofs = _mesh.quad()
    else:
        print("Wrong el_type!")
    
    cfv.drawMesh(coords, edof, 2, el_type)
    plt.pause(2)
    
    """ Denote forces and boundary conditions """
    nDofs = np.max(edof)
    f = np.zeros([nDofs,1])
    
    bc = np.array([],'i')
    bcVal = np.array([],'f')
    
    bc, bcVal = cfu.applybc(bdofs, bc, bcVal, bmarker, value=0.0, dimension=0)
    
    
    cfu.applyforce(bdofs, f, force[1], force[0], force[2])
    
    
    """ Optimisation """
    try:
        x = Opt._Opt()
    except:
        print("Optimisation is not yet implemented")
    SIMP_penal = 3
    nElem=np.size(edof,0)
    x =np.zeros([nElem,1])+0.5
    try: 
        U =FE._FE(x,SIMP_penal,edof,coords,bdofs,f)
    except:
        print("FE is not done yet!")
    
    cfv.draw_element_values(x, coords, edof, 2, el_type, 
                      draw_elements=True, draw_undisplaced_mesh=True, 
                      title="Density", magnfac=25.0)

    
