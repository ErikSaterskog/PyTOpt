import numpy as np
import numpy.linalg
import Mod_Hook as mh
from scipy.sparse.linalg import spsolve
import elastic as el

def gauss_quad(ir):
    if ir ==1:
        GP = np.array([[0] , [0]])
        W = np.array([[2], [2]])
    elif ir ==2:
        g,w = 0.577,1
        GP = np.array([[-g,g,-g,g],[-g,-g,g,g]])
        W = np.array([[w,w,w,w],[w,w,w,w]]) 
    elif ir ==3:
        raise Exception('Not yet implemented')
    else:
        raise Exception('Not correct value of ir')


    wp =W[0,:]*W[1,:]
    xsi = GP[0,:]
    eta = GP[1,:]
    return xsi,eta,wp

def shape_functions(xsi,eta,ir):
    NGP = ir*ir
    N = np.zeros([len(xsi),4])
    dNr = np.zeros([NGP*2,4])
    
    N[:,0] = (1-xsi)*(1-eta)/4
    N[:,1] = (1+xsi)*(1-eta)/4
    N[:,2] = (1+xsi)*(1+eta)/4
    N[:,3] = (1-xsi)*(1+eta)/4

    dNr[0:NGP*2:2,0] = -(1-eta)/4
    dNr[0:NGP*2:2,1] =  (1-eta)/4
    dNr[0:NGP*2:2,2] =  (1+eta)/4
    dNr[0:NGP*2:2,3] = -(1+eta)/4
    
    dNr[1:NGP*2+1:2,0] = -(1-xsi)/4
    dNr[1:NGP*2+1:2,1] = -(1+xsi)/4
    dNr[1:NGP*2+1:2,2] =  (1+xsi)/4
    dNr[1:NGP*2+1:2,3] =  (1-xsi)/4

    
    return N,dNr   
    

def elem4n(ue, ex, ey, ep, mp, materialFun, eq=None):
    
    ptype   = ep[0]             # Which analysis type?
    t       = ep[1]             # Element thickness
    ir      = ep[2]  
    ngp= ir*ir                  # Integration rule and number of gauss points
    matmod  = ep[3]
    xsi,eta,wp = gauss_quad(ir)
    N,dNr = shape_functions(xsi,eta,ir)
    JacTran = np.matmul(dNr,np.array([ex,ey]).T)
    
    b=np.zeros([2,1])
    
    if not eq is None:
        b[0] = eq[0]
        b[1] = eq[1]
    
    
    Ke      = np.zeros([8,8]);   #Preallocate Ke
    fint    = np.zeros([8,1]);   #Preallocate fint
    fext    = np.zeros([8,1]);   #Preallocate fext
    stress  = np.zeros([6, ngp]);#Preallocate stress
    
    
    if ptype == 1:
        raise Exception("Not implemented plane stress yet")
    elif ptype ==2:
       
        for i in range(0,ngp):
            ind = [2*i,2*i+1]
            try:
                detJ = np.linalg.det(JacTran[ind,:])
            except:
                raise Exception("Determinant too small!")

            dNx=spsolve(JacTran[ind,:],dNr[ind,:])

            
            B = np.zeros([3,ngp*2])
            
            B[0,0:ngp*2:2]= dNx[0,:]
            B[1,1:ngp*2:2]= dNx[1,:]
            B[2,0:ngp*2:2]= dNx[1,:]
            B[2,1:ngp*2:2]= dNx[0,:]
            
            N2=np.zeros([2,8])

            N2[0,np.ix_([0,2,4,6])] = N[i,:]
            N2[1,np.ix_([1,3,5,7])]  = N[i,:]
            
            epsilon = np.zeros([6,])
            epsilon[np.ix_([0,1,3])] = np.matmul(B,ue)
            
            #Calculate material response at current gauss point
            [sigma, dsde] = materialFun(epsilon, mp)

                
            stress[:, i] = sigma.reshape(6,)                                   #Save stress for current gauss point
        
           #Calculate the gauss point's contribution to element stiffness and forces
            Dm=dsde[np.ix_([0, 1, 3],[0, 1, 3])]                               # Components for plane strain
            Ke=Ke+np.matmul(np.matmul(B.T,Dm),B)*detJ*wp[i]*t                  # Stiffness contribution
            fint=fint+np.matmul(B.T,sigma[np.ix_([0,1,3])])*wp[i]*detJ*t       # Internal force vector 
            fext=fext+np.matmul(N2.T,b)*detJ*wp[i]*t                           # External force vector
    
    else:
        raise Exception('Only plane strain ep(1)=ptype=2 allowed (unless ep(2)=2, then ep(1)=3 is allowed)');
    
    
    return Ke, fint, fext, stress[:,0],epsilon
            
