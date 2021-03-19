
import Element_Tri_Routine
import Element_Quad_Routine
import calfem.core as cfc


def Tri(ue,ex,ey,ep,mp, materialFun, eq=None):
    if eq==None:
        eq=[0,0] 
    return Element_Tri_Routine.Element_Tri_Routine(ue, ex, ey, ep, mp, materialFun, eq)

def Quad(ue,ex,ey,ep,mp, materialFun, eq=None):
    if eq==None:
        eq=[0,0] 
    return Element_Quad_Routine.Element_Quad_Routine(ue, ex, ey, ep, mp, materialFun, eq)

def LinTri(ue,ex,ey,ep,mp, materialFun, eq=None):
    if eq==None:
        eq=[0,0] 
    D = cfc.hooke(ep[0],mp[0],mp[1])
    Ke,fe = cfc.plante(ex, ey, ep[:2], D, eq)
    return Ke,[],fe.reshape(6,1),[],[]

def LinQuad(ue,ex,ey,ep,mp, materialFun, eq=None):
    if eq==None:
        eq=[0,0] 
    D = cfc.hooke(ep[0],mp[0],mp[1])
    Ke,fe = cfc.plani4e(ex, ey, ep, D, eq)
    return Ke,[],fe.reshape(6,1),[],[]
