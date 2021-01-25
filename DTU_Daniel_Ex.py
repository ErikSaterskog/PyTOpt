# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 15:08:32 2021

@author: Daniel
"""
import numpy as np
import calfem.core as cfc
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.sparse import csc_matrix, linalg


def top(nelx,nely,volfrac,penal,rmin):
    x =np.zeros([nely,nelx])
    x[0:nely,0:nelx] = volfrac
    loop = 0
    change = 1
    
    while change > 0.01:
        loop = loop + 1
        xold = x.copy()
        U = FE(nelx,nely,x,penal)
        
        Ke = lk()
        c = 0
        dc = np.zeros([nely,nelx])
        for ely in range(0,nely):
            for elx in range(0,nelx):
                n1 = (nely+1)*(elx)+ely
                n2 = (nely+1)*(elx+1)+ely
                Ue = U[[2*n1,2*n1+1, 2*n2,  2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2,  2*n1+3]]
                c = c + x[ely,elx]**penal*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
                dc[ely,elx] = -penal*x[ely,elx]**(penal-1)*np.matmul(np.transpose(Ue), np.matmul(Ke,Ue))
    
        dc = Check(nelx,nely,rmin,x,dc)
        
        x = OC(nelx,nely,x,volfrac,dc)
        change =np.max(np.max(abs(x-xold)))
        print(str(change) + '\n')
        
        viridis = cm.get_cmap('Greys', 29)
        
        
        if loop == 1:
            fig, ax = plt.subplots(1, 1, figsize=(6,3),constrained_layout=True)
        
        psm = ax.pcolormesh(np.flip(x,0), cmap=viridis, rasterized=True, vmin=0, vmax=1)
        fig.colorbar(psm, ax=ax)
        plt.pause(1e-6)
        #plt.colormaps(gray)
        #plt.imagesc(-x)
        #plt.axis(equal)
        #plt.axis(tight)
        #plt.axis(off)
        
        plt.show()
        
        if loop > 100:
            return x
        
    return x

def OC(nelx,nely,x,volfrac,dc):
    l1,l2,move=0,1e5,0.2
    while (l2-l1) > 1e-4:
        lmid = 0.5*(l2+l1)
        xnew = np.maximum.reduce([0.001*np.ones([nely,nelx]),np.maximum.reduce([np.maximum.reduce([(x-move),np.minimum.reduce([np.ones([nely,nelx]),np.minimum.reduce([(x+move),x*(np.sqrt(-dc/lmid))]) ]) ]) ])])

        if (sum(sum(xnew)) - volfrac*nelx*nely) > 0:
            l1 = lmid
        else:
            l2 = lmid

    return xnew
def Check(nelx,nely,rmin,x,dc):
    new_dc = np.zeros([nely,nelx])
    for i in range(0,nelx):
        for j in range(0,nely):
            _sum = 0
            for k in range(int(max([i-np.floor(rmin),0])),int(min([i+np.floor(rmin)+1,nelx]))):
                for l in range(int(max([j-np.floor(rmin),0])),int(min([j+np.floor(rmin)+1,nely]))):
                    fac = rmin-np.sqrt((i-k)**2+(j-l)**2)
                    _sum = _sum + max(fac,0)
                    new_dc[j,i] = new_dc[j,i] + max(fac,0)*x[l,k]*dc[l,k]
            new_dc[j,i] = new_dc[j,i]/(x[j,i]*_sum)
    return new_dc

def FE(nelx,nely,x,penal):
    Ke = lk()
    #K = csc_matrix((2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1))).toarray()
    #F = csc_matrix((2*(nely+1)*(nelx+1),1)).toarray()
    K = np.array(np.zeros([2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1)]))
    F = np.zeros([2*(nely+1)*(nelx+1),1])
    U = F.copy()
    fixdofs = []
    for elx in range(0,nelx):
        for ely in range(0,nely):
            n1 =(nely+1)*(elx)+ely
            n2 = (nely+1)*(elx+1)+ely
            edof = [2*n1,2*n1+1, 2*n2,  2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2,  2*n1+3]
            K[np.ix_(edof,edof)] = K[np.ix_(edof,edof)] + x[ely,elx]**penal*Ke

    F[1,0]=-1e8
    
    for i in range(0,2*(nely+1),2):
        fixdofs.append(i)
        
    fixdofs.append(2*(nelx+1)*(nely+1)-1)
            
            
    alldofs = [i for i in range(0,2*(nely+1)*(nelx+1))]
    freedofs = np.setdiff1d(alldofs, fixdofs)
    
    Inv_K = np.linalg.inv(K[np.ix_(freedofs,freedofs)])
    
    U[np.ix_(freedofs)] = np.matmul(Inv_K,F[freedofs])    
    return U


def lk():
    E = 210*10**9
    nu = 0.3
    ex = [0,1,1,0,0.5,1,0.5,0]
    ey = [0,0,1,1,0,0.5,1,0.5]
    ep = [1,2]
    #D  = [[(1-nu)*E/((1+nu)*(1-2*nu)), (1-2*nu)/2*E/((1+nu)*(1-2*nu))],[(1-2*nu)/2*E/((1+nu)*(1-2*nu)),(1-nu)*E/((1+nu)*(1-2*nu))]]
    #Ke = cfc.flw2i8e(ex,ey,ep,D)
    
    k = np.array([1/2-nu/6, 1/8+nu/8,
                  -1/4-nu/12, -1/8+3*nu/8,
                  -1/4+nu/12, -1/8-nu/8,
                  nu/6, 1/8-3*nu/8])
    Ke = E/(1-nu**2)*np.mat([[k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]],
                             [k[1],k[0],k[7],k[6],k[5],k[4],k[3],k[2]],
                             [k[2],k[7],k[0],k[5],k[6],k[3],k[4],k[1]],
                             [k[3],k[6],k[5],k[0],k[7],k[2],k[1],k[4]],
                             [k[4],k[5],k[6],k[7],k[0],k[1],k[2],k[3]],
                             [k[5],k[4],k[3],k[2],k[1],k[0],k[7],k[6]],
                             [k[6],k[3],k[4],k[1],k[2],k[7],k[0],k[5]],
                             [k[7],k[2],k[1],k[4],k[3],k[6],k[5],k[0]]])
    #Edof = np.zeros([(nely+1)*nelx+1)])
    
    
    
    return Ke
    
    
    
x = top(20,10,0.5,5.0,0.5)
print(str(x))
