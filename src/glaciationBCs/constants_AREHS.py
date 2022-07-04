# Physical units: kg, m, s, K

# Physical constants in units: kg, m, s, K
gravity = 9.81 #m/s²

# Numerical constants
s_a = 365.25*24*3600 #=31557600 seconds per year
eps = 1.e-8

# Parametrization
q_geo = 0.1			 #W/m²
T_ini = 273.15 + 8.5 #K
T_min = 273.15 - 1.5 #K

# Geomodel-specific
x_min = 9000		 #m
x_max = 20950		 #m
y_min =-1000		 #m
y_max = 0.0			 #m

lrepo = 3000		 #m / m²
xrmin = (x_max+x_min)/2 - lrepo/2
xrmax = (x_max+x_min)/2 + lrepo/2
yrmin = -1270 - 320/2.
yrmax = -1250 + 320/2.

L_dom = x_max - x_min#m
L_max = 0.8 * L_dom	 #m
H_max = 700			 #m
x_0 = x_min 		 #m

# Key points for time control
t_0 = 0.0 * 0.000 * s_a #s
t_1 = t_0 +  5000 * s_a #s
t_2 = t_1 + 10000 * s_a #s
t_3 = t_2 + 20000 * s_a #s
t_4 = t_3 + 20000 * s_a #s
t_5 = t_4 + 30000 * s_a #s
t_6 = t_5 + 10000 * s_a #s
t_ = [t_0, t_1, t_2, t_3, t_4, t_5, t_6]

# Radioactive waste (Jobmann et al., 2017 - Projekt Ansicht):
# Parameters RK-BE
BE_Q = [842.65, 1269.66, 3895.17, 8308.18, 42363.74]	# W/m³
BE_z = [3.58e-11, 2.95e-10, 9.12e-10, 1.04e-8, 2.62e-8]	# 1/s
BE_f = (2/3)*0.05

# Parameters RK-HA:
HA_Q = [7583.16, 2412.91, 2458.56, 2546.25, 7231.62]	# W/m³
HA_z = [8.28e-9, 6.12e-10, 7.15e-10, 7.15e-10, 8.28e-9]	# 1/s
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

# settings / set up
# flag for 2D / 3D
plotinput = False
