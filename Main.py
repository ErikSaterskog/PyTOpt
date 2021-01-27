"""
This file will call an example, mesh it and solve the optimisation
"""

import numpy as np
import calfem.utils as cfu
import Mesh
import calfem.vis as cfv
import FE


def mainer(g,el_type,force,bmarker):
       
    _mesh = Mesh.Mesh(g,0.05)
    
    if el_type == 2:
        coords, edof, dofs, bdofs = _mesh.tri()
    elif el_type ==3:
        coords, edof, dofs, bdofs = _mesh.quad()
    else:
        print("Wrong el_type!")
    
    nDofs = np.max(edof)
    f = np.zeros([nDofs,1])
    
    bc = np.array([],'i')
    bcVal = np.array([],'f')
    
    bc, bcVal = cfu.applybc(bdofs, bc, bcVal, bmarker, value=0.0, dimension=0)
    
    
    cfu.applyforce(bdofs, f, force[1], force[0], force[2])
    
    SIMP_penal = 3
    nElem=np.size(edof,0)
    x =np.zeros([nElem,1])
    U =FE._FE(x,SIMP_penal,edof,coords,bdofs,f)
        
    cfv.drawMesh(coords, edof, 2, el_type)
    #cfv.figure()
