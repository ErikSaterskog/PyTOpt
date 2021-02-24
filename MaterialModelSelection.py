# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 12:10:37 2021

@author: Daniel
"""
import elem3n
import elem4n

def Tri(ue,ex,ey,ep,mp,eq=None):
    return elem3n.elem3n(ue, ex, ey, ep, mp, eq=None)

def Quad(ue,ex,ey,ep,mp,eq=None):
    return elem4n.elem4n(ue, ex, ey, ep, mp, eq=None)

