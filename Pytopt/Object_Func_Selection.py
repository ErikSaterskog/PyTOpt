"""
Selection script for the different objective functions


Written 2021-05
Made By: Daniel Pettersson & Erik Säterskog
"""

from Pytopt import Object_Func_Compliance, Object_Func_Displacement


def Compliance(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K):
    return Object_Func_Compliance.Compliance(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K)

def Displacement(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K):
    return Object_Func_Displacement.Displacement(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K)


