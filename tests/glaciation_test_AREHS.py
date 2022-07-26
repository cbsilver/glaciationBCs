import pyvista as pv
import numpy as np
import pandas as pd

# the following is automatically done by the toy.py when calling the workflow
# mesh = pv.read("~/glaciationBCs/tests/south.vtu")
#
# bounds = []
# for m_id in np.unique(mesh.cell_data["MaterialIDs"]):
#     # get ymin and ymax
#     bounds += mesh.threshold((m_id, m_id), "MaterialIDs").bounds[2:4]
# bounds = list(np.sort(np.unique(bounds))[::-1])
#
# f = open("south_layer_bounds.py", "w")
# f.write(f"import numpy as np\nsouth_layer_bounds = np.array({bounds})")
# f.close()
#
# toys_file = '~/dgr/input/parameters/toys.csv'
# toy_id = 9
# dfs = pd.read_csv(toys_file)
# dfs = dfs.set_index('set_id')
# df_toy = dfs.xs(toy_id)
#
# parameters_file = '~/dgr/input/parameters/parameterTable_rev03.csv'
# dfm = pd.read_csv(parameters_file, delimiter=",")
# df_layers = pd.merge(df_toy, dfm, how='left', left_on=['model_id', 'layer_id'], right_on=['model_id', 'layer_id'])
# rho_array = df_layers["solid_mass_density_[kg/m^3]_mean"].values.tolist()[1:]
# poro_array = df_layers["Porosity_[-]_mean"].values.tolist()[1:]
# f = open("layer_props.py", "w")
# f.write("import numpy as np\n" +
#          f"rho_array = np.array({rho_array})\n" +
#          f"poro_array = np.array({poro_array})")
# f.close()


from glaciationBCs import glacierclass_AREHS as glc	# glacier
from glaciationBCs import crustclass_AREHS as crc 	# earth crust
from glaciationBCs import repoclass_AREHS as dgr	# repository
from glaciationBCs import airclass_AREHS as air		# atmosphere

from glaciationBCs.constants_AREHS import *

import importlib

importlib.reload(dgr)
importlib.reload(crc)
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
# crust.plot_lithostatic_stress()s
