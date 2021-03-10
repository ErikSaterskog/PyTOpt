
import calfem.geometry as cfg
import elastic as el
import PyTOpt

g = cfg.Geometry()

g.point([0,0])   
g.point([1,0.5],marker=1)
g.point([0,1])


g.line([0, 1])
g.line([1, 2])
g.line([2, 0],marker=2)


g.surface([0 ,1, 2])

force = [-1e5, 1, 2]
bmarker = 2

volFrac = 0.3 
meshSize=0.1
rMin = meshSize*0.7 
changeLimit=0.01 
SIMP_penal = 3
method='OC'
Debug=False

E = 210e9 
nu = 0.3 
mp = [E,nu]

ep=[2,1,2,True,2]

settings = [volFrac,meshSize, rMin, changeLimit, SIMP_penal, method, Debug]


materialFun=el.elastic

PyTOpt.Main(g, force, bmarker, settings, mp, ep, materialFun)














