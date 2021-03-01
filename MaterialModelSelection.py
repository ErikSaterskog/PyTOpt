
import numpy as np
import elem3n
import elem4n
import calfem.core as cfc
import Plani4s as P4

def Tri(ue,ex,ey,ep,mp,eq=None):
    return elem3n.elem3n(ue, ex, ey, ep, mp, eq=None)

def Quad(ue,ex,ey,ep,mp,eq=None):
    return elem4n.elem4n(ue, ex, ey, ep, mp, eq=None)

def LinTri(ue,ex,ey,ep,mp,eq=None):
    D = cfc.hooke(ep[0],mp[0],mp[1])
    Ke,fe = cfc.plante(ex, ey, ep[:2], D, eq=[0,0])
    #sig,eps = cfc.plants(ex, ey, ep, D, ue) 
    #fint = np.matmul(Ke,ue)
    return Ke,[],fe,[],[]

def LinQuad(ue,ex,ey,ep,mp,eq=None):
    D = cfc.hooke(ep[0],mp[0],mp[1])
    Ke,fe = cfc.plani4e(ex, ey, ep, D, eq=[0,0])
    #sig,eps = P4.plani4s(ex,ey,ep,ue)
    #fint = np.matmul(Ke,ue)
    return Ke,[],fe,[],[]
