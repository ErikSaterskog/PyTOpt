
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
from scipy.sparse import csc_matrix, csr_matrix, lil_matrix
import Mod_Hook as mh
import Plani4s
import Debugger
import elem3n
from scipy.sparse.linalg import spsolve
import elem4n
import MaterialModelSelection as MMS

def Main(g,el_type,force,bmarker,settings,mp,ep, materialFun):
    
    #Initiating
    change = 2
    loop = 0
    
    
    #Settings
    E=mp[0]
    v=mp[1]
    try:
        volFrac,meshSize, rMin, changeLimit,SIMP_penal,method,Debug = settings
    except:
        try:
            volFrac,meshSize, rMin, changeLimit,SIMP_penal,method = settings
            Debug = False
        except:
            try:
                volFrac,meshSize, rMin, changeLimit,SIMP_penal = settings
                Debug = False
                method ='OC'
            except:
                raise Exception('Too few inputet settings. Requires a least 5 inputs.')
    
    """ Meshing """
    
    mesh = Mesh.Mesh(g,meshSize)
    
    if el_type == 2:
        coords, edof, dofs, bdofs = mesh.tri()
        if ep[3]==1:
            elementFun = MMS.LinTri
        else:
            elementFun = MMS.Tri    
    elif el_type ==3:
        coords, edof, dofs, bdofs = mesh.quad()
        if ep[3]==1:
            elementFun = MMS.LinQuad
        else:
            elementFun = MMS.Quad
    else:
        print("Wrong Element Type!")
    
    

    """ Denote Forces and Boundary Conditions """
    
    nDofs = np.max(edof)
    
    f = np.zeros([nDofs,1])
    
    
    bc = np.array([],'i')
    bcVal = np.array([],'f')
    
    try:
        for bcmarker in bmarker:
            bc, bcVal = cfu.applybc(bdofs, bc, bcVal, bcmarker, value=0.0, dimension=0)
    except:
        bc, bcVal = cfu.applybc(bdofs, bc, bcVal, bmarker, value=0.0, dimension=0)
    
    bc=bc-1   #Fix a calfem bug
    
    try:
        for fcmarker in force[1]:
            cfu.applyforce(bdofs, f, fcmarker, force[0], force[2])
    except:
        cfu.applyforce(bdofs, f, force[1], force[0], force[2])
    
    
    """ Initialization Cont."""

    nElem=np.size(edof,0)
    
    x =np.zeros([nElem,1])+1#volFrac
        
    #Check sizes, Initialize
    nx=coords[:,0]
    ny=coords[:,1]
    elemCenterX=np.zeros([nElem,1])
    elemCenterY=np.zeros([nElem,1])
    print('Initializing...')
    print('Number of elements: '+str(nElem))
    
    
    #Check element type
    if len(edof[0,:])==6:   #Triangular Element
        elemX=np.zeros([nElem,3])
        elemY=np.zeros([nElem,3])
    elif len(edof[0,:])==8: #Quad Element
        elemX=np.zeros([nElem,4])
        elemY=np.zeros([nElem,4])
    else:
        raise Exception('Unrecognized Element Shape, Check eDof Matrix')
    
    #Find The coordinates for each element's nodes
    for elem in range(nElem):
            
        nNode=np.ceil(np.multiply(edof[elem,:],0.5))-1
        nNode=nNode.astype(int)
            
        elemX[elem,:]=nx[nNode[0:8:2]]
        elemY[elem,:]=ny[nNode[0:8:2]]
        elemCenterX[elem]=np.mean(elemX[elem])
        elemCenterY[elem]=np.mean(elemY[elem])

    #Linear Elastic Constitutive Matrix
    D=cfc.hooke(ep[0], E, v)


    #Create weighting matrix for Filter
    weightMatrix=lil_matrix((nElem,nElem))

    ticH=time.perf_counter()
    for elem in range(0,nElem):
        xDist=elemCenterX-elemCenterX[elem]
        yDist=elemCenterY-elemCenterY[elem]
        dist=np.sqrt(xDist**2+yDist**2)            #Calculates the distance from the current element to all others
        
        weightMatrix[:,elem]=np.maximum(rMin-dist,np.zeros([nElem,1]))
            
    tocH=time.perf_counter()
    print('H:'+str(tocH-ticH))


    FEM = FE.FE(edof,coords,mp,bc)

    """ MAIN LOOP """
    if method == 'OC':
        while change > changeLimit:
            
            loop = loop + 1
            xold = x.copy()
            dc = xold.copy() 
            
        
            U,dR,lambdaF,sig_VM = FEM.fe_nl(x, SIMP_penal, f, ep, elementFun, materialFun)
                        
            tic=time.perf_counter()
            
            for elem in range(nElem):
                if ep[3]==1:
                    
                    if el_type==2:
                        Ke=cfc.plante(elemX[elem,:],elemY[elem,:],ep[0:2],D)   #!THIS COULD BE PLACED OUTSIDE OF LOOP!               #Element Stiffness Matrix for Triangular Element
                    else:
                        Ke=cfc.plani4e(elemX[elem,:],elemY[elem,:],ep,D)[0]    #!THIS COULD BE PLACED OUTSIDE OF LOOP!           #Element Stiffness Matrix for Quad Element
                        
                    Ue = U[np.ix_(edof[elem,:]-1)]
                    dc[elem] = -SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    
                else:
                    lambdaFe = lambdaF[np.ix_(edof[elem,:]-1)]
                    dc[elem] = np.matmul(lambdaFe.T,dR[elem,:].reshape(np.size(edof,1),1))
       
            if Debug and loop==1:
                dc_Num=Debugger.num_Sens_Anal(x,SIMP_penal,edof,coords,bc,f,ep,mp,nElem,elementFun)

                plt.plot(range(0,nElem),(1-dc_Num/dc)*100)
                plt.xlabel('Element')
                plt.ylabel('Percent Error')    
            
            
            toc=time.perf_counter()

            dc = Filter.Check(x,dc,weightMatrix)

            ticOpt=time.perf_counter()
            x = Opt.Optimisation().OC(nElem,x,volFrac,dc)
            tocOpt=time.perf_counter()
    
            change =np.max(np.max(abs(x-xold)))
            
            print('Sens. Anal.:'+str(toc-tic))
            print('Opt:        '+str(tocOpt-ticOpt))
            print('Change:     '+str(change))
            print('Iteration:  '+str(loop))
            print('---------------------------')
            if loop == 50:                                                          # If alternating
                break
            
        
    else:
        
        x = Opt.Optimisation().mma(nElem,SIMP_penal,edof,f,ep,elemX,elemY,D,weightMatrix,volFrac,x,elementFun,el_type,FEM,materialFun)
        x = x.reshape(nElem,1)
        U,dR,lambdaF,sig_VM = FEM.fe_nl(x,SIMP_penal,f,ep,elementFun,materialFun)
        #raise Exception('Not implemented yet!')
            
            
    """ Visualisation """
    
    cfv.draw_element_values(x, coords, edof, 2, el_type,displacements=None,
                      draw_elements=True, draw_undisplaced_mesh=False, 
                      title="Density", magnfac=1.0,clim=(0,1))
    
    cfv.showAndWait()

    
