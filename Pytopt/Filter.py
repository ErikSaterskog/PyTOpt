"""
Low-pass filter.
Inputs:
    x   -Design variables
    dc  -derivative of objective function
    H   -Weight matrix
Outputs:
    new_dc  -Modified derivative of objective function


Written 2021-05
Made By: Daniel Pettersson & Erik Säterskog
"""
def Filter(x,dG0,H):
    
    Hsum=H.sum(axis=1).tolist()
    new_dc=H.multiply(x).multiply(dG0).sum(axis=0).transpose()/(Hsum*x)

    return new_dc






























