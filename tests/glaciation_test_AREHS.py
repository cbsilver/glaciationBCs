import importlib
import pyvista as pv
from glaciationBCs import glacierclass_AREHS as glc	# glacier
from glaciationBCs import crustclass_AREHS as crc 	# earth crust
from glaciationBCs import repoclass_AREHS as dgr	# repository
from glaciationBCs import airclass_AREHS as air		# atmosphere
from glaciationBCs.constants_AREHS import *
from glaciationBCs.pythonBCsOGS_AREHS import uvw_bounds, L_max

"""
dimension = 2
# Ton-Nord 2D
u_min = 9000		 #m
u_max = 20950		 #m
v_min =-2216.03		 #m
v_max = 67.0103		 #m
H_max = 700			 #m
repo_size = 3000  	 #m
L_dom = u_max - u_min
u_0 = u_max
L_max = 0.8 * L_dom
"""

importlib.reload(dgr)
importlib.reload(crc)
importlib.reload(glc.tcr)

dim = 3 # or 2
# xmin, xmax, ymin, ymax, zmin, zmax
bounds = (9000, 20950, -2216.03, 67.0103, 0., 0.)
repo_size = 3000

repo = dgr.repo(BE_Q, BE_z, BE_f, HA_Q, HA_z, HA_f, BE_vol, HA_vol, 
	repo_size, t_inter_BE, t_inter_HA, t_filled, dim)
repo.plot_evolution()
repo.print_max_load()

u_max = uvw_bounds(dim, bounds, 0)[0]
glacier = glc.glacier(L_max(dim, bounds), H_max, u_max, t_)
glacier.plot_evolution()
glacier.plot_evolving_shape()

air = air.air(T_ini, T_min, t_)
air.plot_evolution()

v_min, v_max = uvw_bounds(dim, bounds, 1)
crust = crc.crust(q_geo, v_min, v_max, T_ini, T_bot)
crust.plot_profile(T_ini)
crust.plot_profile_evolution()
#crust.plot_lithostatic_stress()
