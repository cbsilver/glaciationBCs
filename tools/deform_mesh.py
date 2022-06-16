#import glacierclass as glc
#import crustclass as crc
from glaciationBCs import glacierclass as glc	#glacial objects
from glaciationBCs import crustclass as crc 	#crustal objects

import numpy as np
import meshio


name = 'basinmesh_domain'
name = 'basinmesh_boundary_top'
name = 'simTH_Bense_ts_170_t_20017500.000000'
name = 'simTH_Bense_ts_160_t_20015000.000000'
mesh=meshio.read(name + '.vtu')

mesh.points

#(Set=="TH"): # units: kg, m, a, K
L_dom = 1150000 #m
L_max = 575000 #m
H_max = 3200 #m
x_0 = 0.0 #m
t_0 = 20000000 #a
t_1 = t_0 + 12500 #a
t_2 = t_1 +  5000 #a
t_3 = t_2 +  5000 #a
t_4 = t_3 + 10000 #a

glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)

# time for last LGM
t = t_0+17500
t = t_0+15000
y_scale = 20

for point in mesh.points:
	x = point[0]
	u_y = y_scale * glacier.local_deflection_heuristic(x,t)
	print(u_y)
	point += [0,u_y,0]

mesh.write(name + '_def.vtu')

# for all meshes
