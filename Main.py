

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
    
    
    """ Denote forces and boundary conditions """
    nDofs = np.max(edof)
    f = np.zeros([nDofs,1])
    
    bc = np.array([],'i')
    bcVal = np.array([],'f')
    
    bc, bcVal = cfu.applybc(bdofs, bc, bcVal, bmarker, value=0.0, dimension=0)
    
    
    cfu.applyforce(bdofs, f, force[1], force[0], force[2])
    
    
    """ Optimisation """
    try:
        change = 2
        loop = 0
        SIMP_penal = 3
        nElem=np.size(edof,0)
        x =np.zeros([nElem,1])+0.5
        x[3] = 0.9 # Tillagd för att Lindemann dividerar med 0. Kommer tas bort när optimeringen är inkluderad.
        while change > 0.01:
            loop = loop + 1
            xold = x.copy()
            
            U = FE._FE(x,SIMP_penal,edof,coords,bdofs,f)
            
            print('FE:')
            print(toc-tic)
            
            
            tic = time.perf_counter()
            Ke = lk()
            c = 0
            dc = np.zeros([nely,nelx])
            for ely in range(0,nely):
                for elx in range(0,nelx):
                    n1 = (nely+1)*(elx)+ely
                    n2 = (nely+1)*(elx+1)+ely
                    Ue = U[[2*n1,2*n1+1, 2*n2,  2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2,  2*n1+3]]
                    #c = c + x[ely,elx]**penal*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    #c = c + (x[ely,elx]**penal*np.transpose(Ue).dot(Ke.dot(Ue)))
                    dc[ely,elx] = -penal*x[ely,elx]**(penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    #dc[ely,elx] = -penal*x[ely,elx]**(penal-1)*(np.transpose(Ue).dot(Ke.dot(Ue)))
            
            toc = time.perf_counter()
            print('C+Sens:')
            print(toc-tic)
        
        
        
            tic = time.perf_counter()
            dc = Check(nelx,nely,rmin,x,dc)
            toc = time.perf_counter()
            print('Check:')
            print(toc-tic)
            
            tic = time.perf_counter()
            x = OC(nelx,nely,x,volfrac,dc)
            toc = time.perf_counter()
            print('OC:')
            print(toc-tic)
            
            change =np.max(np.max(abs(x-xold)))
            print('------------')
        
        
    except:
        print("Optimisation is not yet implemented")
    
    try: 
        U =FE._FE(x,SIMP_penal,edof,coords,bdofs,f)
    except:
        print("FE is not done yet!")
    
    """ Visualisation """
    
    c = cfv.draw_element_values(U, coords, edof, 2, el_type, 
                      draw_elements=True, draw_undisplaced_mesh=False, 
                      title="Density", magnfac=25.0)
    
    cfv.showAndWait()

    
