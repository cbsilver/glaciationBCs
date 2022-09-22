import pyvista as pv
import numpy as np
import pandas as pd
import importlib

from glaciationBCs import glacierclass_AREHS as glc	# glacier
from glaciationBCs import crustclass_AREHS as crc 	# earth crust
from glaciationBCs import repoclass_AREHS as dgr	# repository
from glaciationBCs import airclass_AREHS as air		# atmosphere
from glaciationBCs.constants_AREHS import *

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
import model_properties_dummy as props

u_min = props.coords_min[0]
u_max = props.coords_max[0]
v_min = props.coords_min[1]
v_max = props.coords_max[1]
L_max = 0.8 * (u_max - u_min)


importlib.reload(dgr)
importlib.reload(crc)
importlib.reload(glc.tcr)

repo = dgr.repo(BE_Q, BE_z, BE_f, HA_Q, HA_z, HA_f, BE_vol, HA_vol, 
	props.repo_size, t_inter_BE, t_inter_HA, t_filled, props.dimension)
repo.plot_evolution()
repo.print_max_load()

glacier = glc.glacier(L_max, H_max, u_max, t_)
glacier.plot_evolution()
glacier.plot_evolving_shape()

air = air.air(T_ini, T_min, t_)
air.plot_evolution()

crust = crc.crust(q_geo, v_min, v_max, T_ini, T_bot)
crust.plot_profile(T_ini, props)
crust.plot_profile_evolution(props)
crust.plot_lithostatic_stress(props)
