from glaciationBCs import glacierclass_AREHS as glc	# glacier
from glaciationBCs import crustclass_AREHS as crc 	# earth crust
from glaciationBCs import repoclass_AREHS as dgr	# repository
from glaciationBCs import airclass_AREHS as air		# atmosphere

from glaciationBCs.constants_AREHS import *

import importlib

importlib.reload(dgr)
importlib.reload(glc.tcr)

repo = dgr.repo(BE_Q, BE_z, BE_f, HA_Q, HA_z, HA_f, BE_vol, HA_vol, lrepo, t_inter_BE, t_inter_HA, t_filled)
repo.plot_evolution()
repo.print_max_load()

glacier = glc.glacier(L_dom, L_max, H_max, u_0, t_)
glacier.plot_evolution()
glacier.plot_evolving_shape()

air = air.air(T_ini, T_min, t_)
air.plot_evolution()

crust = crc.crust(q_geo, v_min, v_max, T_ini, T_bot)
crust.plot_profile(T_ini)
crust.plot_profile_evolution()
