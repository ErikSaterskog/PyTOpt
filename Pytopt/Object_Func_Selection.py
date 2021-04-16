from Pytopt import Object_Func_Energy, Object_Func_Displacement



def Energy(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K):
    return Object_Func_Energy.Energy(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K)

def Displacement(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K):
    return Object_Func_Displacement.Displacement(nelem, ep, el_type, elemx, elemy, D, eq, U, edof, fext_tilde, fextGlobal, SIMP_penal, x, dG0, dR, freedofs, K)


