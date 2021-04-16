import calfem.vis as cfv
import json
import numpy as np


with open('Test_Examples/Saved_Results/myfile.json','r') as f:
    file=json.load(f)

x=np.array(file['x'])
coords=np.array(file['coords'])
edof=np.array(file['edof'])
el_type=file['el_type']

cfv.draw_element_values(x, coords, edof, 2, el_type,displacements=None,draw_elements=True, draw_undisplaced_mesh=False,title="Density", magnfac=1.0,clim=(0,1))
cfv.showAndWait()