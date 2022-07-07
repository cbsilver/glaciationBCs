# Data model of the evolving crust (thermal and mechanical farfield)
# Physical units: kg, m, s, K

import numpy as np

from glaciationBCs.constants_AREHS import gravity
from glaciationBCs.constants_AREHS import y_max # TODO
from glaciationBCs.constants_AREHS import y_min # TODO
from glaciationBCs.constants_AREHS import T_bot # TODO
from glaciationBCs.constants_AREHS import T_ini # TODO


class crust():
	# class variables:
	rho_wat =1000 #kg/m³ TODO

	# constructor
	def __init__(self, q_geo):
		# instance variables: owned by instances of the class, can be different for each instance
		self.q_geo = q_geo
	
	def geothermal_heatflux(self):
		return [0.0, self.q_geo, 0.0]

	def displacement_below(self,x,y,z,t):
		#TODO
		return [0.0, 0.0, 0.0]
		
	def displacement_aside(self,x,y,z,t):
		#TODO
		return [0.0, 0.0, 0.0]

	def geothermal_temperature(self,x,y,z,t,T_top):
		# linear profile according to geothermal heatflux
		DT = T_bot - T_top
		Dy = (y_min - y_max)
		return DT/Dy * (y - y_max) + T_ini
		
		T_crust = 0
		return T_crust

	def hydrostatic_pressure(self,x,y,z,t):
		# linear profile according to gravity
		p_pore = self.rho_wat * gravity * (y_max - y)
		return p_pore
