# List of all relevant constants
# Physical units: kg, m, s, K

import numpy as np
from dgr import model_interface as model
from glaciationBCs import coord_control_AREHS as uvw# coordinates

# Physical constants in units: kg, m, s, K
gravity = 9.81		 #m/s²
rho_wat = 1000		 #kg/m³
rho_ice = 900		 #kg/m³
c_p_wat = 4280		 #J/kg/K

# Numerical constants
s_a = 365.25*24*3600 #=31557600 seconds per year
eps = 1.e-8

# settings / set up
dimension = model.dimension
plotinput = False

# Geomodel-specific parameters
coord_ctrl = uvw.coord_control(model.dimension)
coords_min = (model.xmin, model.ymin, model.zmin)
coords_max = (model.xmax, model.ymax, model.zmax)
u_min, v_min, w_min = coord_ctrl.assign_coordinates(coords_min)
u_max, v_max, w_max = coord_ctrl.assign_coordinates(coords_max)

# Target edge-length of repos, resulting repo will have a shorter edge-length
# actual length/area is determined by the resulting mesh
salz_kissen_2d_repo = 3000 #m
salz_kissen_3d_repo = 3000 #m²
ton_nord_2d_repo = 3000 #m
ton_nord_3d_repo = 3000 #m²
repos = {(2, 1): ton_nord_2d_repo, (2, 3): salz_kissen_2d_repo,
         (3, 1): ton_nord_3d_repo, (3, 3): salz_kissen_3d_repo}

drepo_crop = repos[(dimension, model.model_id)]
drepo = 1.
# the actual drepo value is only written after execution of crop_repo.py
# but drepo here has to exist beforehand for the workflow to work
if hasattr(model, 'drepo'):
    drepo = model.drepo

# these define the vertical bounds for cropping the repository
# in the end only the upper bound is important, since we extract the top 
# boundary, the bottom bound only has to be large enough to catch at least 
# one element worth of thickness
ton_nord_repo_v_bounds = (-890, -2200)
salz_repo_v_bounds = (-600, -800)
repo_v_bounds = {1: ton_nord_repo_v_bounds, 3: salz_repo_v_bounds}
repo_vmin, repo_vmax = repo_v_bounds[model.model_id]

# parameters for glacier
L_dom = u_max - u_min#m
L_max = 0.8 * L_dom	 #m
H_max = 700			 #m
u_0 = u_min 		 #m

# Key points for time control
t_0 = 0.0 * 0.000 * s_a #s
t_1 = t_0 +  5000 * s_a #s
t_2 = t_1 + 10000 * s_a #s
t_3 = t_2 + 20000 * s_a #s
t_4 = t_3 + 20000 * s_a #s
t_5 = t_4 + 30000 * s_a #s
t_6 = t_5 + 10000 * s_a #s
t_ = [t_0, t_1, t_2, t_3, t_4, t_5, t_6]

# Thermal Parametrization
q_geo = 0.065		 #W/m²   set to 0 to see the effect clearly
T_ini = 273.15 + 8.5 #K
T_min = 273.15 - 1.5 #K
T_bot = 345.783      #K

# Radioactive waste (Jobmann et al., 2017 - Projekt Ansicht):
# Parameters RK-BE
BE_Q = np.array([842.65, 1269.66, 3895.17, 8308.18, 42363.74])		# W/m³
BE_z = np.array([3.58e-11, 2.95e-10, 9.12e-10, 1.04e-8, 2.62e-8])	# 1/s
BE_f = (2/3)*0.05

# Parameters RK-HA:
HA_Q = np.array([7583.16, 2412.91, 2458.56, 2546.25, 7231.62])		# W/m³
HA_z = np.array([8.28e-9, 6.12e-10, 7.15e-10, 7.15e-10, 8.28e-9])	# 1/s
HA_f = (2*0.182/0.72)*0.05

# Total volume:
BE_vol = 7632.0 #m³
HA_vol = 1342.8 #m³

# Times:
t_zwischen_BE = 23*s_a
t_zwischen_HA = 30*s_a
t_voll = 80*s_a
# => end of stepwise filling = heat source starts in the model

t_inter_BE = 23*s_a
t_inter_HA = 30*s_a
t_filled = 80*s_a