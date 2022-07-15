# Data model of the evolving crust (thermal and mechanical farfield)
# Physical units: kg, m, s, K

import numpy as np

from glaciationBCs.constants_AREHS import gravity
from glaciationBCs.constants_AREHS import rho_wat

class crust():
	# class variables:
	# -

	# constructor
	def __init__(self, q_geo, y_min, y_max, T_bot):
		# instance variables: owned by instances of the class, can be different for each instance
		self.q_geo = q_geo
		self.y_min = y_min
		self.y_max = y_max
		self.T_bot = T_bot
	
	def geothermal_heatflux(self):
		return [0.0, self.q_geo, 0.0]

	def displacement_below(self,x,y,z,t):
		return [0.0, 0.0, 0.0]
		
	def displacement_aside(self,x,y,z,t):
		return [0.0, 0.0, 0.0]

	def geothermal_temperature(self,x,y,z,t,T_atm):
		# linear profile according to geothermal heatflux
		DT = self.T_bot - T_atm
		Dy = (self.y_min - self.y_max)
		return DT/Dy * (y - self.y_max) + T_atm
		
		T_crust = 0
		return T_crust

	def hydrostatic_pressure(self,x,y,z,t):
		# linear profile according to gravity
		p_pore = rho_wat * gravity * (self.y_max - y)
		return p_pore
