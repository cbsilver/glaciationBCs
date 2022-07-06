from glaciationBCs import glacierclass as glc	#glacial objects
from glaciationBCs import crustclass as crc 	#crustal objects
from glaciationBCs import airclass as air		# aerial objects
import numpy as np
import importlib
from importlib import reload

s_a = 365.25*24*3600 #=31557600 seconds per year

Set = "T"

T_N = 266.15 #K
T_S = 276.15 #K
T_C = 8 #K

if (Set=="M"): # units: kg, m, s, K
	L_dom = 120000 #m
	L_max = 0.7*L_dom
	H_max = 200 #m
	x_0 = -0.5*L_dom
	t_0 = 0.00 #s
	t_1 = 1.0000 #s
	t_2 = 2
	t_3 = 3
	t_4 = 4

if (Set=="T"): # units: kg, m, a, K
	L_dom = 1150000 #m
	L_max = 575000 #m
	H_max = 3200 #m
	x_0 = 0.0 #m
	t_0 = 17500 #a
	t_1 = t_0 + 12500 #a
	t_2 = t_1 +  5000 #a
	t_3 = t_2 +  5000 #a
	t_4 = t_3 + 10000 #a

if (Set=="TH"): # units: kg, m, a, K
	L_dom = 1150000 #m
	L_max = 575000 #m
	H_max = 3200 #m
	x_0 = 0.0 #m
	t_0 = 20000000 #a
	t_1 = t_0 + 12500 #a
	t_2 = t_1 +  5000 #a
	t_3 = t_2 +  5000 #a
	t_4 = t_3 + 10000 #a

if (Set=="HM"): # units: kg, m, s, K
	L_dom = 1150000 #m
	L_max = 575000 #m
	H_max = 3200 #m
	x_0 = 0.0 #m
	t_0 = 0.0 * 32500 * s_a #a
	t_1 = t_0 + 12500 * s_a #a
	t_2 = t_1 +  5000 * s_a #a
	t_3 = t_2 +  5000 * s_a #a
	t_4 = t_3 + 10000 * s_a #a
	
glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
#importlib.reload(glaciationBCs)
importlib.reload(glc)
glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)

glacier.plot_deflection_rate()
"""

glacier.read_data_GIA()

# find index for y=-7500m
i=3

# stage: glacial advance
t0=0.0
tE=12.5
tRange = np.linspace(t0,tE,25+1)
glacier.lineplot_evolution_uxuy(i,tRange)

glacier.interpolateX_data_uxuy(8001,0,11.8)
glacier.interpolateX_data_uxuy(1150000,-50,32.5)

importlib.reload(crc)
crust = crc.crust()
path2data = '/home/cbs/Forschung/Simulations/OpenGeoSys/SedimentaryBasinBense/HM/dataGIA/'
crust = crc.crust(path2data)

crust.ylineplot_evolution_uxuy(i,tRange)
crust.xlineplot_evolution_uxuy(1150000,tRange)

crust.interpolateX_data_uxuy(8001,0,11.8)
crust.interpolateY_data_uxuy(1150000,-6100,12.25)  

crust.check_interpolationY_uxuy(32.5)
crust.check_interpolationX_uxuy(32.5)
crust.check_interpolationY_uxuy(12.5)

crust.ux_data.loc[230] # get row with index 230
crust.ux_data.loc[[229,230,231]]

a[start:stop]  # items start through stop-1
a[start:]      # items start through the rest of the array
a[:stop]       # items from the beginning through stop-1
a[:]           # a copy of the whole array

xi = self.xvalues[0:231] # from element 0 to 231-1=230

uxt1 = crust.ux_data[str(2.5)].to_numpy(dtype=float)
uxt1yi = uxt1[0:10]
a = uxt1yi
a+=2
uxt1yi
uxt1[0:10]

import copy
uxt1yi = copy.copy(uxt1)
uxt1yi = uxt1.copy()

"""
