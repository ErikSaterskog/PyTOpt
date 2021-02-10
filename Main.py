

import numpy as np
import calfem.utils as cfu
import Mesh
import calfem.vis as cfv
import FE
import Opt
import Filter
import time
import matplotlib.pyplot as plt
import calfem.core as cfc
from scipy.sparse import csr_matrix
import Mod_Hook as mh
import Plani4s
import Debugger

def _Main(g,el_type,force,bmarker,settings,mp):
    
    #Settings
    E=mp[0]#210*1e9
    v=mp[1]#0.3
    Linear = True
    ptype=2         #ptype=1 => plane stress, ptype=2 => plane strain
    ep=[ptype,1,2]    #ep[ptype, thickness, integration rule(only used for QUAD)]  
    change = 2
    loop = 0
    SIMP_penal = 3
    volFrac,meshSize, rMin, changeLimit = settings
    
    
    
    """ Meshing """
    
    _mesh = Mesh.Mesh(g,meshSize)
    
    if el_type == 2:
        coords, edof, dofs, bdofs = _mesh.tri()
    elif el_type ==3:
        coords, edof, dofs, bdofs = _mesh.quad()
    else:
        print("Wrong Element Type!")
    
    

    """ Denote Forces and Boundary Conditions """
    
    nDofs = np.max(edof)
    f = np.zeros([nDofs,1])
    bc = np.array([],'i')
    bcVal = np.array([],'f')
    
    bc, bcVal = cfu.applybc(bdofs, bc, bcVal, bmarker, value=0.0, dimension=0)
    bc=bc-1   #Fix a calfem bug
    
    cfu.applyforce(bdofs, f, force[1], force[0], force[2])
    
    
    
    """ Initialization Cont."""

    nElem=np.size(edof,0)
    x =np.zeros([nElem,1])+volFrac
        
    #Check sizes, Initialize
    nElem=np.size(edof,0)
    nx=coords[:,0]
    ny=coords[:,1]
    elemCenterX=np.zeros([nElem,1])
    elemCenterY=np.zeros([nElem,1])
    print('Initializing...')
    print('Number of elements: '+str(nElem))
    
    
    #Check element type
    if len(edof[0,:])==6:   #Triangular Element
        Tri=True
        elemX=np.zeros([nElem,3])
        elemY=np.zeros([nElem,3])
    elif len(edof[0,:])==8:
        Tri=False           #Use Quad Instead
        elemX=np.zeros([nElem,4])
        elemY=np.zeros([nElem,4])
    else:
        raise Exception('Unrecognized Element Shape, Check eDof Matrix')
    
    #Find The coordinates for each element's nodes
    for elem in range(0,nElem):
            
        nNode=np.ceil(np.multiply(edof[elem,:],0.5))-1
        nNode=nNode.astype(int)
            
        elemX[elem,:]=nx[nNode[0:8:2]]
        elemY[elem,:]=ny[nNode[0:8:2]]
        elemCenterX[elem]=np.mean(elemX[elem])
        elemCenterY[elem]=np.mean(elemY[elem])

    #Linear Elastic Constitutive Matrix
    D=cfc.hooke(ptype, E, v)


    #Create weighting matrix for Filter
    weightMatrix=np.zeros([nElem,nElem])

    for elem in range(0,nElem):
        xDist=elemCenterX-elemCenterX[elem]
        yDist=elemCenterY-elemCenterY[elem]
        dist=np.sqrt(xDist**2+yDist**2)            #Calculates the distance from the current element to all others
        #for elemOther in range(0,nElem):           #Checks which are inside the radius rMin
        #breakpoint()
        weightMatrix[:,elem]=np.maximum(rMin-dist,np.zeros([nElem,1]))[:,0]
            
    #weightMatrix=csr_matrix(weightMatrix)




        """ MAIN LOOP """

    while change > changeLimit:
        
        loop = loop + 1
        xold = x.copy()
        
        
        #"""LINEAR"""
        if Linear == True:
            U = FE._FE(x,SIMP_penal,edof,coords,bc,f,ep,mp)  #FEA
        
        
            dc = xold.copy() 
            
            tic=time.perf_counter()
            
            if Tri:  #Tri Elements
                for elem in range(0,nElem):  
                    Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D)   #!THIS COULD BE PLACED OUTSIDE OF LOOP!               #Element Stiffness Matrix for Triangular Element
                    Ue = U[np.ix_(edof[elem,:]-1)]
                    dc[elem] = -SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    
            else:    #Quad Elements
                for elem in range(0,nElem):            
                    Ke=cfc.plani4e(elemX[elem,:],elemY[elem,:],ep,D)    #!THIS COULD BE PLACED OUTSIDE OF LOOP!           #Element Stiffness Matrix for Quad Element
                    Ue = U[np.ix_(edof[elem,:]-1)]
                    dc[elem] = -SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke[0],Ue))
            
            
            
            if True and loop==1:
                dc_Num=Debugger.num_Sens_Anal(x,SIMP_penal,edof,coords,bc,f,ep,mp,nElem)
            
                plt.plot(range(0,nElem),dc_Num-dc)
                plt.xlabel('Element')
                plt.ylabel('dc difference')
            
            
        #"""NON LINEAR"""
        else:
            U = FE._FE_NL(x,SIMP_penal,edof,coords,bc,f,ep,mp)  #FEA
            dc = xold.copy() 
            
            tic=time.perf_counter()
            
            ed = cfc.extractEldisp(edof, U)
        
            if Tri:  #Tri Elements
                for elem in range(0,nElem):  
                    eps = np.zeros([6,])
                    eps_2D=cfc.plants(elemX[elem,:],elemY[elem,:],ep,D,ed[elem,:])[1] 
                    eps[0:4] = eps_2D
                    D_new = mh._mod_hook(eps,mp)[1]
                    Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D_new[np.ix_([0,1,2,3],[0,1,2,3])])
                    Ue = U[np.ix_(edof[elem,:]-1)]
                    dc[elem] = -SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    
            else:    #Quad Elements
                for elem in range(0,nElem):            
                    eps=Plani4s.plani4s(elemX[elem,:],elemY[elem,:],ep,ed[elem,:]) 
                    D_new = mh._mod_hook(eps,mp)[1]
                    Ke=cfc.plani4e(elemX[elem,:],elemY[elem,:],ep,D_new[np.ix_([0,1,2,3],[0,1,2,3])])[0]
                    Ue = U[np.ix_(edof[elem,:]-1)]
                    dc[elem] = -SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    
                    
                
        
        toc=time.perf_counter()

        dc = Filter.Check(x,dc,weightMatrix)
        
        try:
            x = Opt.Optimisation().OC(nElem,x,volFrac,dc)
        except:
            print("Optimisation is not yet implemented")
        

        change =np.max(np.max(abs(x-xold)))
        print('Sens. Anal.:'+str(toc-tic))
        print('Change:     '+str(change))
        print('Iteration:  '+str(loop))
        print('---------------------------')
        if loop == 100:                                                          # If alternating
            break
        
        
        
        
        
    """ Visualisation """
    
    cfv.draw_element_values(x, coords, edof, 2, el_type,displacements=U,
                      draw_elements=True, draw_undisplaced_mesh=False, 
                      title="Density", magnfac=1.0)
    
    cfv.showAndWait()

    
