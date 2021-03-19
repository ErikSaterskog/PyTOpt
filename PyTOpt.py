"""
This the main-module and is called by the user with the inputs
g           - Geometry object
el_type     - 2 if triangles, 3 if quad elements
force       - A vector with the magnitude and position of the prescribed forces.
bmarker     - A vector with the boundaries and prescribed displacements.
settings    - Includes settings for restrictions on the optimisation and the filtering.
                volfrac     - How much material is allowed
                meshSize    - How big the average element should be.
                rMin        - In how large radius the element should accknowledge other elements in the filter
                changeLimit - A tolerance for the optimisation
                SIMP_penal  - A penalty factor forcing elementdensity towards zero or one
                method      - Which optimisation method should be used.
                Debug       - A option to check the sensitivity analysis against numerical sensitivity.
mp          - Material parameters
ep          - Element parameters
                ptype       - 2 if plain strain
                t           - Thickness
                ir          - Integration rule
                non/lin     - 1 if linear analysis, 2 if nonlinear analysis
materialFun - A material model with strain and mp as input and stress and consitutive matrix as output.

The module ends by plotting the result of 

"""
import numpy as np
import calfem.utils as cfu
import Mesh
import FE
import Optimisation as Opt
import Filter
import time
import matplotlib.pyplot as plt
import calfem.core as cfc
from scipy.sparse import lil_matrix
import Debugger
import Element_Routine_Selection as ERS
import json


def Main(g,force,bmarker,settings,mp,ep, materialFun, eq=None):
    
    #Initiating
    change = 2
    loop = 0
    ticGlobal=time.perf_counter()
    save=False
    
    
    #Settings
    E = mp[0]
    nu = mp[1]
    el_type=ep[2]
    ptype=2
    intRule=2
    ep=[ptype, ep[0], intRule, ep[1], ep[2]]
    plt.rcParams["figure.dpi"] = 200
    
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
        if ep[3]:
            elementFun = ERS.LinTri
        else:
            elementFun = ERS.Tri    
    elif el_type == 3:
        coords, edof, dofs, bdofs = mesh.quad()
        if ep[3]:
            elementFun = ERS.LinQuad
        else:
            elementFun = ERS.Quad
    else:
        print("Wrong Element Type!")
    
    
    """ Denote Forces and Boundary Conditions """
    
    f = np.zeros([np.max(edof),1])
        
    bc = np.array([],'i')
    bcVal = np.array([],'f')
    
    try:
        for bcmarker in bmarker:
            bc, bcVal = cfu.applybc(bdofs, bc, bcVal, bcmarker, value=0.0, dimension=0)
    except:
        bc, bcVal = cfu.applybc(bdofs, bc, bcVal, bmarker, value=0.0, dimension=0)
    
    bc=bc-1   
    
    try:
        for fcmarker in force[1]:
            cfu.applyforce(bdofs, f, fcmarker, force[0], force[2])
    except:
        cfu.applyforce(bdofs, f, force[1], force[0], force[2])
    
    
    """ Initialization Cont."""

    nelem = np.size(edof,0)
    
    x = np.zeros([nelem,1])+volFrac
        
    #Check sizes, Initialize
    nx=coords[:,0]
    ny=coords[:,1]
    elemCenterx=np.zeros([nelem,1])
    elemCentery=np.zeros([nelem,1])
    print('Initializing...')
    print('Number of elements: '+str(nelem))
    
    
    #Check element type
    if len(edof[0,:])==6:   #Triangular Element
        elemx=np.zeros([nelem,3])
        elemy=np.zeros([nelem,3])
    elif len(edof[0,:])==8: #Quad Element
        elemx=np.zeros([nelem,4])
        elemy=np.zeros([nelem,4])
    else:
        raise Exception('Unrecognized Element Shape, Check edof Matrix')
    
    #Find The coordinates for each element's nodes
    for elem in range(nelem):
            
        nnode=np.ceil(np.multiply(edof[elem,:],0.5))-1
        nnode=nnode.astype(int)
            
        elemx[elem,:]=nx[nnode[0:8:2]]
        elemy[elem,:]=ny[nnode[0:8:2]]
        elemCenterx[elem]=np.mean(elemx[elem])
        elemCentery[elem]=np.mean(elemy[elem])

    #Linear Elastic Constitutive Matrix
    D=cfc.hooke(ep[0], E, nu)


    #Create weighting matrix for Filter
    weightMatrix=lil_matrix((nelem,nelem))

    ticH=time.perf_counter()
    for elem in range(0,nelem):
        xdist=elemCenterx-elemCenterx[elem]
        ydist=elemCentery-elemCentery[elem]
        dist=np.sqrt(xdist**2+ydist**2)            
        
        weightMatrix[:,elem]=np.maximum(rMin-dist,np.zeros([nelem,1]))
            
    tocH=time.perf_counter()
    print('H:'+str(tocH-ticH))

    #Initiate the FEM
    FEM = FE.FE(edof,coords,mp,bc)

    """ MAIN LOOP """
    if method == 'OC':
        while change > changeLimit:
            
            loop = loop + 1
            xold = x.copy()
            dc = xold.copy() 
            
        
            U,dR,lambdaF,sig_VM,fext_tilde = FEM.fe_nl(x, SIMP_penal, f, ep, elementFun, materialFun, eq)
                        
            tic=time.perf_counter()
            
            for elem in range(nelem):
                if ep[3]:
                    
                    if el_type==2:
                        Ke=cfc.plante(elemx[elem,:],elemy[elem,:],ep[0:2],D) #Element Stiffness Matrix for Triangular Element
                        if eq is not None:
                            eqe=[eq,eq,eq]
                            eqe=np.array(eqe).reshape(1,6)
                        else:
                            eqe=np.zeros([1,6])
                    else:
                        Ke=cfc.plani4e(elemx[elem,:],elemy[elem,:],ep,D)[0]  #Element Stiffness Matrix for Quad Element
                        if eq is not None:
                            eqe=[eq,eq,eq,eq]
                            eqe=np.array(eqe).reshape(1,8)
                        else:
                            eqe=np.zeros([1,8])

                    Ue = U[np.ix_(edof[elem,:]-1)]
                    fext_tildee = fext_tilde[np.ix_(edof[elem,:]-1)]
                    dc[elem] = np.matmul(fext_tildee.T,Ue)  -  SIMP_penal*x[elem][0]**(SIMP_penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                    
                    if dc[elem] >0:
                        print(str(elem) + ':' +str(dc[elem]))
                        dc[elem] = 0
                    
                else:
                    if eq is not None:
                            eqe=[eq,eq,eq]
                            eqe=np.array(eqe).reshape(1,6)
                    else:
                            eqe=np.zeros([1,6])
                    lambdaFe = lambdaF[np.ix_(edof[elem,:]-1)]
                    Ue = U[np.ix_(edof[elem,:]-1)]
                    dc[elem] = np.matmul(eqe,Ue) + np.matmul(lambdaFe.T,dR[elem,:].reshape(np.size(edof,1),1))
                    
                    if dc[elem] >0:
                        print(str(elem) + ':' +str(dc[elem]))
                        dc[elem] = 0
                        
            if Debug and loop==1:
                dc_Num=Debugger.num_Sens_Anal(x,SIMP_penal,edof,coords,bc,f,ep,mp,nelem,elementFun)

                plt.plot(range(0,nelem),(1-dc_Num/dc)*100)
                plt.xlabel('Element')
                plt.ylabel('Percent Error')    
            
            
            toc=time.perf_counter()

            dc = Filter.Filter(x,dc,weightMatrix)

            ticOpt=time.perf_counter()
            x = Opt.Optimisation().OC(nelem,x,volFrac,dc)
            tocOpt=time.perf_counter()
    
            change = np.max(np.max(abs(x-xold)))
            
            print('Sens. Anal.:'+str(toc-tic))
            print('Opt:        '+str(tocOpt-ticOpt))
            print('Change:     '+str(change))
            print('Iteration:  '+str(loop))
            print('---------------------------')
            
            if loop == 50:
                break
            
            if loop % 5 == 0: 
                fig, ax = plt.subplots()
                for j in range(0,nelem):
                    ax.fill(elemx[j,:], elemy[j,:], color = [1,1,1]*(1-x[j]))
                ax.axis('equal')
                plt.show() 
                
        
    elif method =='MMA':
        
        x = Opt.Optimisation().mma(nelem,SIMP_penal,edof,f,ep,elemx,elemy,D,weightMatrix,volFrac,x,elementFun,el_type,FEM,materialFun,eq)
        x = x.reshape(nelem,1)
        #U,dR,lambdaF,sig_VM = FEM.fe_nl(x,SIMP_penal,f,ep,elementFun,materialFun,eq)
     
    else:
        raise Exception('No Optimisation method of that names')
            
            
    """ Visualisation """
    
    tocGlobal=time.perf_counter()
    timeMin=(tocGlobal-ticGlobal)/60
    print('Total computation time: '+str(int(timeMin))+'m '+str(round(np.mod(timeMin,1)*60,1))+'s')
    
    if save:
        data = {'x': x.tolist(),'coords': coords.tolist(),'edof': edof.tolist(), 'el_type': el_type}
        with open('Saved_Results/myfile.json', 'w') as outfile:
            json.dump(data, outfile, indent=4)
    
    fig, ax = plt.subplots()
    for j in range(0,nelem):
        ax.fill(elemx[j,:], elemy[j,:], color = [1,1,1]*(1-x[j]))
        
    ax.axis('equal')
    plt.show()    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
