# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 15:08:32 2021

@author: Daniel
"""
import numpy as np
import calfem.core as cfc
import matplotlib.pyplot as plt
import time
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.sparse import csc_matrix, linalg, csr_matrix, lil_matrix, coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib.image import NonUniformImage



def top(nelx,nely,volfrac,penal,rmin):
    x =np.zeros([nely,nelx])
    x[0:nely,0:nelx] = volfrac
    loop = 0
    change = 1
    
    while change > 0.01:
        loop = loop + 1
        xold = x.copy()
        
        
        tic = time.perf_counter()
        U = FE(nelx,nely,x,penal)
        toc = time.perf_counter()
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
        #print(str(change) + '\n')
        
        
    #    if loop > 50:
    #        fig, ax = plt.subplots(1, 1, figsize=(6,3),constrained_layout=True)
    #        viridis = cm.get_cmap('Greys', 29)
    #        psm = ax.pcolormesh(np.flip(x,0), cmap=viridis, rasterized=True, vmin=0, vmax=1)
    #        fig.colorbar(psm, ax=ax)
    #        plt.pause(1e-6)
    #        plt.show()
    #        return x
    #fig, ax = plt.subplots(1, 1, figsize=(6,3),constrained_layout=True)
    #viridis = cm.get_cmap('Greys', 29)
    #psm = ax.pcolormesh(np.flip(x,0), cmap=viridis, rasterized=True, vmin=0, vmax=1)
    #fig.colorbar(psm, ax=ax)
    #plt.pause(1e-6)
    #plt.show()    
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
    K = np.zeros([2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1)])
    F = np.zeros([2*(nely+1)*(nelx+1),1])
    U = F.copy()
    fixdofs = []
    edof=[]
    tic = time.perf_counter()
    for elx in range(0,nelx):
        for ely in range(0,nely):
            n1 =(nely+1)*(elx)+ely
            n2 = (nely+1)*(elx+1)+ely
            edof=(2*n1,2*n1+1, 2*n2,  2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2,  2*n1+3)
            #edof.extend([2*n1,2*n1+1, 2*n2,  2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2,  2*n1+3])
            edofIndex=np.ix_(edof,edof)
            K[edofIndex] = K[edofIndex] + x[ely,elx]**penal*Ke
            

    toc = time.perf_counter()
    print('FE, ASSEM:')
    print(toc-tic)
    F[1,0]=-1e6
    
    #K=csr_matrix(K)
    #F=csr_matrix(F)
    
    tic = time.perf_counter()
        
    for i in range(0,2*(nely+1),2):
        fixdofs.append(i)
        
    fixdofs.append(2*(nelx+1)*(nely+1)-1)
    
    fixdofs = np.array(fixdofs)        
            
    alldofs = [i for i in range(0,2*(nely+1)*(nelx+1))]
    freedofs = np.setdiff1d(alldofs, fixdofs)
    

    Ue = spsolve(K[np.ix_(freedofs,freedofs)],F[np.ix_(freedofs)])

    
    j = -1
    for i in Ue:
        j = j + 1
        U[freedofs[j]] = i 
    return U

    toc = time.perf_counter()
    print('FE, SOLVE:')
    print(toc-tic)

def lk():
    E = 210*1e9
    nu = 0.3
    ex = [0,1,1,0,0.5,1,0.5,0]
    ey = [0,0,1,1,0,0.5,1,0.5]
    ep = [1,2]
    
    
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
    
    
    
    
    return Ke
    

tic = time.perf_counter()
x = top(40,30,0.5,3.0,1.5)
toc = time.perf_counter()
print(toc-tic)


plt.matshow(x)
plt.show()










