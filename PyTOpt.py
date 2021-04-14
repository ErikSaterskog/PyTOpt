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

## IMPORTING MODULES ############################
import numpy as np
import calfem.utils as cfu
import Mesh
import FE
import Optimisation as Opt
import Filter
import time
import matplotlib.pyplot as plt
import calfem.core as cfc
from scipy.sparse import coo_matrix
import Debugger
import Element_Routine_Selection as ERS
import json
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import Material_Routine_Selection as mrs

################################################


def Main(g, force, bmarker, settings, mp, ep, materialFun, ObjectFun, eq=None):
    
    
    """ Settings """
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
    
    """------------------------------------------"""
    
    
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
    
    """---------------------------------------------"""
    
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
    
    """---------------------------------------------"""
    
    """ Initialization"""
    
    change = 2
    loop = 0
    ticGlobal=time.perf_counter()
    save=False
    
    nelem = np.size(edof,0)
    x = np.zeros([nelem,1])+volFrac # Startguess
    
    #Initiate the FEM and Linear Elastic Constitutive Matrix
    FEM = FE.FE(edof,coords,mp,bc)
    D=cfc.hooke(ep[0], E, nu)
    
    #Check sizes, Initialize
    nx=coords[:,0]
    ny=coords[:,1]
    elemCenterx=np.zeros([nelem,1])
    elemCentery=np.zeros([nelem,1])
    print('Initializing...')
    print('Number of elements: '+str(nelem))
    
    #Find The coordinates for each element's nodes
    if len(edof[0,:])==6:   #Triangular Element
        elemx=np.zeros([nelem,3])
        elemy=np.zeros([nelem,3])
    elif len(edof[0,:])==8: #Quad Element
        elemx=np.zeros([nelem,4])
        elemy=np.zeros([nelem,4])
    else:
        raise Exception('Unrecognized Element Shape, Check edof Matrix')
        
    for elem in range(nelem):
            
        nnode=np.ceil(np.multiply(edof[elem,:],0.5))-1
        nnode=nnode.astype(int)
            
        elemx[elem,:]=nx[nnode[0:8:2]]
        elemy[elem,:]=ny[nnode[0:8:2]]
        elemCenterx[elem]=np.mean(elemx[elem])
        elemCentery[elem]=np.mean(elemy[elem])

    """---------------------------------"""

    """ Creating weighting matrix for Filter"""
    data = []
    row = []
    col = []
    
    ticH=time.perf_counter()
    for elem in range(0,nelem):
        xdist=elemCenterx-elemCenterx[elem]
        ydist=elemCentery-elemCentery[elem]
        dist=np.sqrt(xdist**2+ydist**2)            
        weightdata = np.maximum(rMin-dist,np.zeros([nelem,1]))
        data.extend(weightdata[np.where(weightdata>0)])
        row.extend(np.where(weightdata>0)[0])
        col.extend(np.repeat(elem,len(np.where(weightdata>0)[0])))
        
    
    weightMatrix=coo_matrix((data,(row,col)),shape=(nelem,nelem))
       
    tocH=time.perf_counter()
    print('H:'+str(tocH-ticH))
    """-------------------------------------"""
    

    """ MAIN LOOP """
    if method == 'OC':
        while change > changeLimit:
            
            loop = loop + 1
            xold = x.copy()
            dG0 = xold.copy() 
            
            #FE Calculation
            U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x, SIMP_penal, f, ep, elementFun, materialFun, eq)
                        
            #Object Function Calculation
            G0, dG0 = ObjectFun(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K)
            
            if Debug and loop==1:
                dG0_Num=Debugger.num_Sens_Anal(x,SIMP_penal,edof,coords,bc,f,ep,mp,nelem,elementFun)
                plt.plot(range(0,nelem),(1-dG0_Num/dG0)*100)
                plt.xlabel('Element')
                plt.ylabel('Percent Error')

            #Apply Filter
            dG0 = Filter.Filter(x,dG0,weightMatrix)

            x = Opt.Optimisation().OC(nelem,x,volFrac,dG0)

            change = np.max(np.max(abs(x-xold)))
            
            print('G0:         '+str(float(G0)))
            print('Change:     '+str(change))
            print('Iteration:  '+str(loop))
            print('---------------------------')
            
            if loop == 50:
                break
            
            if loop % 5 == 0: 
                patches = []
                fig, ax = plt.subplots()
                for j in range(0,nelem):
                    polygon = Polygon(np.transpose([elemx[j,:], elemy[j,:]]))
                    patches.append(polygon)

                p = PatchCollection(patches, cmap=matplotlib.cm.Greys)
                p.set_array(np.transpose(x)[0])
                ax.add_collection(p)
                ax.axis('equal')
                plt.show() 
                
        
    elif method =='MMA':
        G0,x,eps_h = Opt.Optimisation().mma(f,edof,elemx,elemy,x,SIMP_penal,ep,elementFun,materialFun,FEM,el_type,D,eq,weightMatrix,volFrac,ObjectFun)
    else:
        raise Exception('No Optimisation method of that name. Try OC or MMA.')
    
    """--------------------------------------------------"""        
            
    """ Visualisation """
    
    # Time and objective function value
    tocGlobal=time.perf_counter()
    timeMin=(tocGlobal-ticGlobal)/60
    print('Total computation time: '+str(int(timeMin))+'m '+str(round(np.mod(timeMin,1)*60,1))+'s')
    print('Final G0: '+str(float(G0)))
    
    # Possibility to save settings
    if save:
        data = {'x': x.tolist(),'coords': coords.tolist(),'edof': edof.tolist(), 'el_type': el_type}
        with open('Saved_Results/myfile.json', 'w') as outfile:
            json.dump(data, outfile, indent=4)
    
    # Plot the final optimised structure
    patches = []
    fig, ax = plt.subplots()
    for j in range(0,nelem):
        polygon = Polygon(np.transpose([elemx[j,:], elemy[j,:]]))
        patches.append(polygon)

    p = PatchCollection(patches, cmap=matplotlib.cm.Greys)
    p.set_array(np.transpose(x)[0])
    ax.add_collection(p)
    ax.axis('equal')
    plt.show()  
    
    # If linear, find the hydrostatic strain
    if ep[3]:
        ep[3] = False
        materialFun = mrs.Elastic 
        if el_type == 2:
            elementFun = ERS.Tri    
        elif el_type == 3:
            elementFun = ERS.Quad
        U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x, SIMP_penal, f, ep, elementFun, materialFun, eq)

    eps_h[np.where(x<0.1)]=0
    
    # Plot the hydrostatic strain over the optimised structure.
    patches = []
    fig, ax = plt.subplots()
    for j in range(0,nelem):
        polygon = Polygon(np.transpose([elemx[j,:], elemy[j,:]]))
        patches.append(polygon)

    p = PatchCollection(patches, cmap=matplotlib.cm.RdBu)
    p.set_array(np.transpose(eps_h)[0])
    plt.colorbar(p)
    p.set_clim(-max(abs(eps_h))/2, max(abs(eps_h))/2)
    ax.add_collection(p)
    ax.axis('equal')
    plt.title('Hydrostatic Strain')
    plt.show()  

    """------------------------------------------------------"""
    
    
    
    
    
    
    
    
    
    
    
    
    
    
