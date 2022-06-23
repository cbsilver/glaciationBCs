# Physical units: kg, m, s, K

# hier nur zeitabhängige RBn
# voneinander unabhängige Methoden in eigene Klasse/eigene Datei stecken
# bei Baumstruktur: gemeinsame (Basis-)Klasse oben drüber
# bei zirkulären Abhängigkeiten: Methoden/Objekte in eine Klasse
# Bedeutung von x,y,z durch georef. Modell

# Physical constants in units: kg, m, s, K
gravity = 9.81 #m/s²
fricnum = 0.2

# Numerical constants
s_a = 365.25*24*3600 #=31557600 seconds per year
eps = 1.e-8

#TODO
"""
rho_ice = 900 #kg/m³
rho_wat =1000 #kg/m³
T_under = 273.15 + 0.5 #K
"""
q_geo = 0.1 #W/m²

# Choose parametrization
T_ini = 273.15 + 8.5 #K
T_min = 273.15 - 1.5 #K

coord_min = 9000 	#m
coord_max = 20950	#m
L_dom = coord_max - coord_min #m
L_max = 0.8* L_dom	#m
H_max = 700		#m
x_0 = coord_min #m
t_0 = 0.0 * 0.000 * s_a #s
t_1 = t_0 +  5000 * s_a #s
t_2 = t_1 + 10000 * s_a #s
t_3 = t_2 + 20000 * s_a #s
t_4 = t_3 + 20000 * s_a #s
t_5 = t_4 + 30000 * s_a #s
t_6 = t_5 + 10000 * s_a #s

#Jobmann et al., 2017 - Projekt Ansicht:
#Parameters RK-BE
BE_Q1 = 842.65
BE_z1 = 3.58e-11

BE_Q2 = 1269.66
BE_z2 = 2.95e-10

BE_Q3 = 3895.17
BE_z3 = 9.12e-10

BE_Q4 = 8308.18
BE_z4 = 1.04e-8

BE_Q5 = 42363.74
BE_z5 = 2.62e-8

BE_Q = [842.65, 1269.66, 3895.17, 8308.18, 42363.74]	# W/m³
BE_z = [3.58e-11, 2.95e-10, 9.12e-10, 1.04e-8, 2.62e-8]	# 1/s

#Parameters RK-HA:
HA_Q1 = 7583.16
HA_z1 = 8.28e-9

HA_Q2 = 2412.91
HA_z2 = 6.12e-10

HA_Q3 = 2458.56
HA_z3 = 7.15e-10

HA_Q4 = 2546.25
HA_z4 = 7.15e-10

HA_Q5 = 7231.62
HA_z5 = 8.28e-9

HA_Q = [7583.16, 2412.91, 2458.56, 2546.25, 7231.62]	# W/m³
HA_z = [8.28e-9, 6.12e-10, 7.15e-10, 7.15e-10, 8.28e-9]	# 1/s

#Total volume:
BE_vol = 7632.0 #m³
HA_vol = 1342.8 #m³

#Total surface:
dgr_area = 3000.0 #m / m²

#Times:
t_zwischen_BE = 23*s_a
t_zwischen_HA = 30*s_a
t_voll = 80*s_a 
# => end of stepwise filling = heat source starts in the model

t_inter_BE = 23*s_a
t_inter_HA = 30*s_a
t_filled = 80*s_a

# settings / set up
# flag for 2D / 3D
plotinput = 0
