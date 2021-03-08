
import numpy as np
#import calfem.core as cfc

def mat(eps,mp):
    
    sigma = np.zeros(np.shape(eps))
    
    E_ten      = mp[0]
    E_com       = E_ten*3
    nu_ten      = 0.3
    nu_com      = 0.3
    
    low_lim=0.002
    
    
    eps_h = sum(eps[:2])/3
    
    D_ten=Dfun(E_ten,nu_ten)
    D_com=Dfun(E_com,nu_com) 

    x=1200*(eps_h-low_lim)
    
    ten_quote=1-np.exp(-x)/(1+np.exp(-x))
    comp_quote=1-ten_quote
    
    D=D_com*comp_quote + D_ten*ten_quote 
    sigma[:4] = np.matmul(D,eps[:4])
    
    return sigma.reshape(6,1), D


def Dfun(E,nu):
    return E/((1+nu)*(1-2*nu))*np.array([[1-nu,nu,nu,0],[nu,1-nu,nu,0],[nu,nu,1-nu,0],[0,0,0,(1-2*nu)/2]])















