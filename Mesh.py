

import numpy as np
import calfem.mesh as cfm

class Mesh:
    
    def __init__(self,g,elsize):
        self.mesh = cfm.GmshMesh(g)
        self.mesh.dofs_per_node = 2  # Degrees of freedom per node.
        self.mesh.el_size_factor = elsize  # Factor that changes element sizes.
        
    
    def quad(self):
        self.mesh.el_type = 3
        
        # mesh creation 
        coords, edof, dofs, bdofs, elementmarkers = self.mesh.create()
        
        return coords, edof, dofs, bdofs
    
        
    def tri(self):
        self.mesh.el_type = 2
        
        coords, edof, dofs, bdofs, elementmarkers = self.mesh.create()
        return coords, edof, dofs, bdofs
        
    
    
