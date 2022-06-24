import glacierclass_AREHS as glc	# glacier
import crustclass_AREHS as crc 	# earth crust
import repoclass_AREHS as dgr	# repository
import airclass_AREHS as air		# atmosphere

import importlib

from constants_AREHS import *

importlib.reload(dgr)

repo = dgr.repo(BE_Q, BE_z, BE_f, HA_Q, HA_z, HA_f, BE_vol, HA_vol, dgr_area, t_inter_BE, t_inter_HA, t_filled)
repo.plot_evolution()
repo.print_max_load()

glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_)
glacier.plot_evolution()

air = air.air(T_ini, T_min, t_)
air.plot_evolution()
