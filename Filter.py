




def Check(nelx,nely,rmin,x,dc):
    
    
    
    new_dc = np.zeros([nely,nelx])
    for i in range(0,nelx):
        for j in range(0,nely):
            _sum = 0
            for k in range(int(max([i-np.floor(rmin),0])),int(min([i+np.floor(rmin)+1,nelx]))):
                for l in range(int(max([j-np.floor(rmin),0])),int(min([j+np.floor(rmin)+1,nely]))):
                    fac = rmin-np.sqrt((i-k)**2+(j-l)**2)
                    _sum = _sum + max(fac,0)
                    new_dc[j,i] = new_dc[j,i] + max(fac,0)*x[l,k]*dc[l,k]
            new_dc[j,i] = new_dc[j,i]/(x[j,i]*_sum)
    return new_dc



































