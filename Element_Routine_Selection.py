
import Element_Tri_Routine
import Element_Quad_Routine
import calfem.core as cfc


def Tri(ue,ex,ey,ep,mp, materialFun, eq=None):
    return Element_Tri_Routine.Element_Tri_Routine(ue, ex, ey, ep, mp, materialFun, eq=None)

def Quad(ue,ex,ey,ep,mp, materialFun, eq=None):
    return Element_Quad_Routine.Element_Quad_Routine(ue, ex, ey, ep, mp, materialFun, eq=None)

def LinTri(ue,ex,ey,ep,mp,eq=None):
    D = cfc.hooke(ep[0],mp[0],mp[1])
    Ke,fe = cfc.plante(ex, ey, ep[:2], D, eq=[0,0])
    return Ke,[],fe,[],[]

def LinQuad(ue,ex,ey,ep,mp,eq=None):
    D = cfc.hooke(ep[0],mp[0],mp[1])
    Ke,fe = cfc.plani4e(ex, ey, ep, D, eq=[0,0])
    return Ke,[],fe,[],[]
