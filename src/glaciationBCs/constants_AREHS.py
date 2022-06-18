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

rho_ice = 900 #kg/m³
rho_wat =1000 #kg/m³
T_under = 273.15 + 0.5 #K

# Choose parametrization
T_N = 266.15 #K
T_S = 276.15 #K
T_C = 8 #K
Set = "TH"

L_dom = 120000 #m
L_max = 0.7*L_dom
H_max = 200 #m
x_0 = -0.5*L_dom
t_0 = 0.00 #s
t_1 = 1.0000 #s
t_2 = 2
t_3 = 3
t_4 = 4

# settings / set up
# flag for 2D / 3D
