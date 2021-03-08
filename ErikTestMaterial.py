
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

def mathencky(eps,mp):
    import calfem.core as cfc
    
    sigma = np.zeros(np.shape(eps))
    eps=eps.reshape(6,1)
    
    E      = mp[0]
    nu     = 0.3
    sy0    = 220e6
    
    delta=np.array([[1],[1],[1],[0],[0],[0]])
    G=E/(2*(1+nu))
    K=E/(3*(1-2*nu))*0.1
    
    I_s=np.diag([1,1,1,0.5,0.5,0.5])
    I_sdev=I_s-(1/3)*delta*delta.transpose()
    
    #eps_VM=np.sqrt(2/3)*np.sqrt(eps*I_sdev*eps.transpose())
    
    if sum(eps[:2])/3>0:
        sign=1
    else:
        sign=-1
    eps_VM=sign*np.sqrt(2/3)*np.sqrt(np.matmul(np.matmul(eps.transpose(),I_sdev),eps))
    
    
    if eps_VM<=(sy0/(3*G)):
        G_star=G
        dGde=0
    else:
        G_star=sy0/(3*eps_VM)
        dGde=-sy0/(3*(eps_VM**2))
        
    eps_dev=I_sdev*eps  
    sigma=(2*G_star*I_sdev+K*(delta*delta.transpose()))*eps
    
    if eps_VM < sy0/(3*G):
        D=cfc.hooke(2,E,nu)
    else:
        D=(2*G_star*I_sdev+K*(delta*delta.transpose())+dGde*(4/3)*(eps_dev*eps_dev.transpose())/eps_VM)
    
    sigma=sigma[0,:]

    
    return sigma.reshape(6,1), D















