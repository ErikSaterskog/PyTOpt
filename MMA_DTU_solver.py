def mma_solver(bc,MP,f,Edof,elemx,elemy,x,SIMP_penal,ep,elementType,materialFun,FEM,el_type,D,eq,weightMatrix,volFrac):
    # -*- coding: utf-8 -*-
    """
    This function is uses Method of moving asymptotes along with KKt check to 
    optimize the objective function.
    
       Input:
        bc   : Boundary condition               
        MP   : Material parameter [E,A,Sigma_y] 
        f    : Loads                            
        Edof : Element degrees of freedom       
        ex   : Elements position in x           
        ey   : Elements position in y           
        Ndof : Number of degrees of freedom     
        
      Output:
         xval: vector with the optimized areas
         ed  : displacements

    """
    import numpy as np
    from MMA_DTU_funct import mmasub,kktcheck
    import Object_Func_Selection as ofs
    import Filter
    ObjectFun = ofs.Energy
    
    m=1
    
    #algorithmic parameters and initial guess
    nelem = np.size(Edof,0)
    
    #initialization
    eeen = np.ones((nelem,1))
    eeem = np.ones((m,1))
    zerom = np.zeros((m,1))
    xmin = 1.e-1*eeen   #lower bound on x
    xmax = 1*eeen       #upper bound on x
    move = 1.0
    c = 1000*eeem
    d = zerom.copy()
    a0 = 0
    a=np.zeros((m,1))
    xold1 = x.copy()
    xold2 = x.copy()
    low = xmin.copy()
    upp = xmax.copy()
    
    maxoutit =150
    kkttol = 0	
    dc = x.copy()
    
    # Calculate function values and gradients of the objective and constraints functions

    
    U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x,SIMP_penal,f,ep,elementType,materialFun)
    
      
    f0val, dc = ObjectFun(nelem, ep, el_type, elemx, elemy, D, eq, U, Edof, fext_tilde, fextGlobal, SIMP_penal, x, dc, dR, freedofs, K)
    
    dc = Filter.Filter(x,dc,weightMatrix)
    
    

    dfdx = eeen.copy()
    dfdx = dfdx.reshape(m,nelem)
    fval = sum(x)-volFrac*len(x)
    
    # The iterations starts
    kktnorm = kkttol+10
    kktnorm = kkttol+10
    outit = 0
    
#-----------------------------------------------------------------------------
    
    while (kktnorm > kkttol) and (outit < maxoutit):
        outit += 1
        
        
        dc = np.array(dc)
        xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp =  \
            mmasub(m,nelem,outit,x,xmin,xmax,xold1,xold2,f0val,dc,fval,dfdx,low,upp,a0,a,c,d,move)
        dc = np.matrix(dc)                                                                
        # update
        xold2 = xold1.copy()
        xold1 = x.copy()
        x = xmma.copy()
        
#------------------------------------------------------------------------------        
        # Re-calculate function values and gradients of the objective and constraints functions
        U, dR, sig_VM, fext_tilde, fextGlobal, eps_h, freedofs, K = FEM.fe_nl(x,SIMP_penal,f,ep,elementType,materialFun)
        f0val = np.matmul(f.T,U)  
          
        f0val, dc = ObjectFun(nelem, ep, el_type, elemx, elemy, D, eq, U, Edof, fext_tilde, fextGlobal, SIMP_penal, x, dc, dR, freedofs, K)
        
        dc = Filter.Filter(x,dc,weightMatrix)
        
        
    
        dfdx = eeen.copy()
        dfdx = dfdx.reshape(m,nelem)
        fval = sum(x)-volFrac*len(x)
#----------------------------------------------------------------------------- 
   
        # The residual vector of the KKT conditions is calculated
        dc = np.array(dc)
        residu,kktnorm,residumax = \
            kktcheck(m,nelem,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,dc,fval,dfdx,a0,a,c,d)
        dc = np.matrix(dc)
    
    return f0val,x,eps_h
