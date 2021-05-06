"""
This the main-module and is called by the user.
Inputs:
    g           - Geometry object
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
    mp[E,nu,eps_y]
             E          - Young's modulus
             nu         - Poission's ratio
             eps_y      - Yielding strain for bilinear material model
    ep[thickness, linear, el_type]
             thickness  - thickness of the 2D material
             linear     - True-linear, False-nonlinear
             el_type    - 2 indicates triangular elements and 3 indicates
             quad elements.
    materialFun - A material model with strain and mp as input and stress and consitutive matrix as output.
    ObjectFun   - An objective function with multiply inputs.

Output:
    Plots showing the density of the optimised result.

Written 2021-05
Made By: Daniel Pettersson & Erik Säterskog
"""

## IMPORTING MODULES ############################
import time
import numpy as np
import calfem.utils as cfu
import calfem.core as cfc
from scipy.sparse import coo_matrix
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from Pytopt import Element_Routine_Selection as ERS
from Pytopt import Material_Routine_Selection as mrs
from Pytopt import Mesh, FE, Filter, Debugger
################################################


_defaultSettings = {'volFrac':0.5,'meshSize':0.1,'rmin':0.1*0.7,'changeLimit': 0.01,'SIMP_const':3,'Debug':False}

def Main(g, force, bmarker, mp, ep, materialFun, ObjectFun, OptFun, settingsdict={}, eq=None, maxiter=30):
    
    
    """ Settings """
    E = mp['E']
    nu = mp['nu']
    el_type=ep[2]
    ptype=2
    intRule=2
    ep=[ptype, ep[0], intRule, ep[1], ep[2]]
    plt.rcParams["figure.dpi"] = 200
    settings=[]
    for setting in _defaultSettings.keys():
        if settingsdict.get(setting) != None:
            settings.append(settingsdict.get(setting))
        else:
            settings.append(_defaultSettings.get(setting))
    
    
    volFrac,meshSize, rMin, changeLimit,SIMP_const,Debug = settings
    
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

    
    nelem = np.size(edof,0)
    x = np.zeros([nelem,1])+volFrac # Startguess
    
    #Initiate the FEM and Linear Elastic Constitutive Matrix
    FEM = FE.FE(edof,coords,mp,bc)
    D=cfc.hooke(ep[0], E, nu)
    G0List=[]
    Areae = []
    
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
        
        PQ = np.sqrt((elemx[elem,1]-elemx[elem,0])**2 + (elemy[elem,1]-elemy[elem,0])**2)
        PR = np.sqrt((elemx[elem,2]-elemx[elem,0])**2 + (elemy[elem,2]-elemy[elem,0])**2)
        QR = np.sqrt((elemx[elem,2]-elemx[elem,1])**2 + (elemy[elem,2]-elemy[elem,1])**2)
        if len(edof[0,:])==6:
            semiper = (PR+PQ+QR)/2
            Areae.append(np.sqrt(semiper*(semiper-PR)*(semiper-PQ)*(semiper-QR)))
        else:
            
            RS = np.sqrt((elemx[elem,3]-elemx[elem,2])**2 + (elemy[elem,3]-elemy[elem,2])**2)
            PS = np.sqrt((elemx[elem,3]-elemx[elem,0])**2 + (elemy[elem,3]-elemy[elem,0])**2)

            semiper1 = (PS+RS+PR)/2
            semiper2 = (PR+PQ+QR)/2

            A1 = np.sqrt(semiper1*(semiper1-PR)*(semiper1-RS)*(semiper1-PS))
            A2 = np.sqrt(semiper2*(semiper2-PR)*(semiper2-PQ)*(semiper2-QR))

            Areae.append(A1+A2)
    """---------------------------------"""
    Areae = (Areae/sum(Areae)).reshape(nelem,1)  
    
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
    while change > changeLimit:
            
        loop = loop + 1
        xold = x.copy()
        dG0 = xold.copy() 
            
        #FE Calculation
        U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x, SIMP_const, f, ep, elementFun, materialFun, eq)
        
        #Object Function Calculation
        G0, dG0 = ObjectFun(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_const, x, dG0, dR, freedofs, K)
        G0List.append(G0[0][0])

        
        if Debug and loop==1:
            error_Std=[]
            error_Mean=[]
            epsList=np.logspace(-7,-2.9, num=70)
            for eps in epsList:
                dG0_Num=Debugger.num_Sens(x, SIMP_const, edof, coords, bc, f, ep, mp, nelem, elementFun, materialFun, eq, eps)
                error=(1-dG0_Num/dG0)
                error_Std.append(np.std(error))
                error_Mean.append(np.mean(error))
                print(eps)
            
            plt.figure()

            plt.errorbar(np.log10(epsList), error_Mean, yerr=error_Std)
            plt.xlabel('Log10 of Perturbation')
            plt.ylabel('Error')

            
        #Apply Filter
        dG0 = Filter.Filter(x, dG0, weightMatrix)

        x = OptFun(x, volFrac, G0, dG0, Areae)

        change = np.max(np.max(abs(x-xold)))
    
            
        print('G0:         '+str(float(G0)))
        print('Change:     '+str(change))
        print('Opt.Iters:  '+str(loop))
        print('---------------------------------')
            
        if loop == maxiter:
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
            
            

        if loop % 5 == -1:
            rMin=np.maximum(rMin/1.2,0.001)
            
            data = []
            row = []
            col = []
            for elem in range(0,nelem):
                xdist=elemCenterx-elemCenterx[elem]
                ydist=elemCentery-elemCentery[elem]
                dist=np.sqrt(xdist**2+ydist**2)            
                weightdata = np.maximum(rMin-dist,np.zeros([nelem,1]))
                data.extend(weightdata[np.where(weightdata>0)])
                row.extend(np.where(weightdata>0)[0])
                col.extend(np.repeat(elem,len(np.where(weightdata>0)[0])))
            
            weightMatrix=[]
            weightMatrix=coo_matrix((data,(row,col)),shape=(nelem,nelem))
     
    """--------------------------------------------------"""        
            
    """ Visualisation """
    
    # Time and objective function value
    tocGlobal=time.perf_counter()
    timeMin=(tocGlobal-ticGlobal)/60
    print('Total computation time: '+str(int(timeMin))+'m '+str(round(np.mod(timeMin,1)*60,1))+'s')
    print('Final G0: '+str(float(G0)))
    print('---------------------------------')
    print(G0List)

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
        U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x, SIMP_const, f, ep, elementFun, materialFun, eq)
    
    print('Max eps_h:  '+str(np.max(eps_h)))
    print('Min eps_h:  '+str(np.min(eps_h)))
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
