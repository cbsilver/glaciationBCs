# Data model of the evolving crust (thermal and mechanical farfield)
# Physical units: kg, m, s, K

import numpy as np

from glaciationBCs.constants_AREHS import gravity
from glaciationBCs.constants_AREHS import rho_wat
from glaciationBCs.constants_AREHS import c_p_wat

class crust():
	# class variables:
	v_fluid = 1e-11	#m/s

	# constructor
	def __init__(self, q_geo, y_min, y_max, T_ini, T_bot):
		# instance variables: owned by instances of the class, can be different for each instance
		self.q_geo = q_geo
		self.y_min = y_min
		self.y_max = y_max
		self.T_bot = T_bot
		self.T_ini = T_ini
	
	def geothermal_heatflux(self):
		return [0.0, self.q_geo, 0.0]

	def displacement_below(self):
		return [0.0, 0.0, 0.0]
		
	def displacement_aside(self):
		return [0.0, 0.0, 0.0]

	def geothermal_temperature(self, y, T_atm):
		# linear profile according to geothermal heatflux
		DT = self.T_bot - T_atm
		Dy = (self.y_min - self.y_max)
		return DT/Dy * (y - self.y_max) + T_atm

	def lateral_heatflux(self, y, T_atm):
		# linear profile according to ???
		q_max = self.v_fluid * (T_atm - self.T_ini) * rho_wat * c_p_wat
		Dy = (self.y_min - self.y_max)
		return q_max/Dy * (y - self.y_max)

	def hydrostatic_pressure(self, y):
		# linear profile according to gravity
		p_pore = rho_wat * gravity * (self.y_max - y)
		return p_pore
