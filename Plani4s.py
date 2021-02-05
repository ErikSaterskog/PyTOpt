# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 09:16:26 2021

@author: Daniel


Addition to the missing plani4s in CALFEM-Python

"""

import numpy as np
import numpy.linalg

def gauss_quad(ir):
    if ir ==1:
        GP = np.array([[0] , [0]])
        W = np.array([[2], [2]])
    elif ir ==2:
        g,w = 0.577,1
        GP = np.array([[-g,g,-g,g],[-g,-g,g,g]])
        W = np.array([[w,w,w,w],[w,w,w,w]]) 
    elif ir ==3:
        raise Exception('Not yet implemented')
    else:
        raise Exception('Not correct value of ir')


    wp =W[0,:]*W[1,:]
    xsi = GP[0,:]
    eta = GP[1,:]
    return xsi,eta,wp

def shape_functions(xsi,eta,ir):
    NGP = ir*ir
    N = np.zeros([len(xsi),4])
    dNr = np.zeros([NGP*2,4])
    
    N[:,0] = (1-xsi)*(1-eta)/4
    N[:,1] = (1+xsi)*(1-eta)/4
    N[:,2] = (1+xsi)*(1+eta)/4
    N[:,3] = (1-xsi)*(1+eta)/4

    dNr[0:NGP*2:2,0] = -(1-eta)/4
    dNr[0:NGP*2:2,1] =  (1-eta)/4
    dNr[0:NGP*2:2,2] =  (1+eta)/4
    dNr[0:NGP*2:2,3] = -(1+eta)/4
    
    dNr[1:NGP*2+1:2,0] = -(1-xsi)/4
    dNr[1:NGP*2+1:2,1] = -(1+xsi)/4
    dNr[1:NGP*2+1:2,2] =  (1+xsi)/4
    dNr[1:NGP*2+1:2,3] =  (1-xsi)/4

    
    return N,dNr   
    

def plani4s(ex,ey,ep,ed):


    ptype,t,ir = ep
    
    if ptype == 1:
        raise Exception("Not implemented plane stress yet")
    elif ptype ==2:
        eps_2d = np.zeros([3,])
        xsi,eta,wp = gauss_quad(ir)
        N,dNr = shape_functions(xsi,eta,ir)
        JacTran = np.matmul(dNr,np.array([ex,ey]).reshape(4,2))
        ngp = ir*ir
        for i in range(0,ngp):
            ind = [2*i,2*i+1]
            try:
                detJ = np.linalg.det(JacTran[ind,:])
            except:
                raise Exception("Determinant too small!")
            invJac = np.linalg.inv(JacTran[ind,:])
            dNx  = np.matmul(invJac,dNr[ind,:])
            
            B = np.zeros([3,ngp*2])
            
            B[0,0:ngp*2:2]= dNx[0,:]
            B[1,1:ngp*2:2]= dNx[1,:]
            B[2,0:ngp*2:2]= dNx[1,:]
            B[2,1:ngp*2:2]= dNx[0,:]
            
            eps_2d = eps_2d + np.matmul(B,ed)
        
        eps = np.zeros([6,])
        eps[[0,1,3],] = eps_2d
        return eps
    else:
        raise Exception("No such ptype..")
 



ir =2
ex = [0,10,10,0]
ey = [0,0,10,10]
ep = [2,1,2]
ed = [0.1,0.1,0,0.3,0,0,0.7,0.2]

eps = plani4s(ex,ey,ep,ed)


