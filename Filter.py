import numpy as np
import time
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
import calfem.core as cfc



def Check(eDof,coords,dofs,rMin,x,dc):
    
    #Settings
    Timers=True
    tic = time.perf_counter()
    
    #Check element type
    if len(eDof[0,:])==6:   #Triangular Element
        Tri=True
    elif len(eDof[0,:])==8:
        Tri=False           #Use Quad Instead
    else:
        raise Exception('Unrecognized Element Shape, Check eDof Matrix')
    

    #Initialize 
    nDof=np.max(eDof)
    nElem=np.size(eDof,0)
    nx=coords[:,0]
    ny=coords[:,1]
    new_dc=np.zeros([np.size(dc),1])
    
    #Find Elements Coordinates
    if Tri:
        elemX=np.zeros([nElem,3])
        elemY=np.zeros([nElem,3])
    else:
        elemX=np.zeros([nElem,4])
        elemY=np.zeros([nElem,4])
    
    
    elemCenterX=np.zeros([nElem,1])
    elemCenterY=np.zeros([nElem,1])
    
    
    #Find Centers
    for elem in range(0,nElem):
        
        nNode=np.ceil(np.multiply(eDof[elem,:],0.5))-1
        nNode=nNode.astype(int)
        
        elemX[elem,:]=nx[nNode[0:8:2]]
        elemY[elem,:]=ny[nNode[0:8:2]]
        elemCenterX[elem]=np.mean(elemX[elem])
        elemCenterY[elem]=np.mean(elemY[elem])
    

    for elem in range(0,nElem):
        _sum=0
        sumHxf=0
        xDist=elemCenterX-elemCenterX[elem]
        yDist=elemCenterY-elemCenterY[elem]
        dist=np.sqrt(xDist**2+yDist**2)     #Calculates the distance from the current element to all others
        for elemOther in range(0,nElem):    #Checks which are inside the radius rMin
            if dist[elemOther]<rMin:
                H=dist[elemOther]
                _sum=_sum+H
                sumHxf = sumHxf + H*x[elemOther]*dc[elemOther]
        new_dc[elem][0]=(sumHxf)/(x[elem]*_sum)


    if Timers:
        toc = time.perf_counter()  
        print('Check Time:')
        print(toc-tic)        
    return new_dc
































