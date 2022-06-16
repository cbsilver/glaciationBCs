# Physical units: kg, m, s, K

# voneinander unabhängige Methoden in eigene Klasse/eigene Datei stecken
# bei Baumstruktur: gemeinsame (Basis-)Klasse oben drüber
# bei zirkulären Abhängigkeiten: Methoden/Objekte in eine Klasse

# Physical constants in units: kg, m, s, K
gravity = 9.81 #m/s²
fricnum = 0.2
b_sub = 0.300 # dim-less proportionality constant glacier height - subsidence
b_reb = 0.225 # dim-less proportionality constant glacier height - rebounding

# Numerical constants
s_a = 365.25*24*3600 #=31557600 seconds per year
eps = 1.e-8

rho_ice = 900 #kg/m³
rho_wat =1000 #kg/m³
T_under = 273.15 + 0.5 #K
