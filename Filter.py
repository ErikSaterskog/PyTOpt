
def Filter(x,dc,H):
    
    Hsum=H.sum(axis=1).tolist()
    new_dc=H.multiply(x).multiply(dc).sum(axis=0).transpose()/(Hsum*x)

    return new_dc






























