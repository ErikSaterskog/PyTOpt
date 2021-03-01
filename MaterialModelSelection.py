

import elem3n
import elem4n
import calfem.core as cfc

def Tri(ue,ex,ey,ep,mp,eq=None):
    return elem3n.elem3n(ue, ex, ey, ep, mp, eq=None)

def Quad(ue,ex,ey,ep,mp,eq=None):
    return elem4n.elem4n(ue, ex, ey, ep, mp, eq=None)

def LinTri(ex,ey,ep,D,eq=None):
    return cfc.plante(ex, ey, ep, D, eq=None)

def LinQuad(ex, ey, ep, D, eq=None):
    return cfc.plani4e(ex, ey, ep, D, eq=None)