import time


def Filter(x,dc,H):
    
    #Settings
    Timers=True
    tic = time.perf_counter()

    Hsum=H.sum(axis=1).tolist()
    new_dc=H.multiply(x).multiply(dc).sum(axis=0).transpose()/(Hsum*x)
    
    if Timers:
        toc = time.perf_counter()  
        print('Check Time: '+str(toc-tic))
    return new_dc






























