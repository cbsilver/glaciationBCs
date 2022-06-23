# Climate model of the evolving atmosphere
# parameterized analytical function for (cyclic) temperature evolution
# Physical units: kg, m, s, K

import numpy as np

from constants_AREHS import *

class air():
	# class variables: owned by the class itself, shared by all instances of the class
	pressure  = 0.e3 #Pa
	
	def __init__(self, T_ini, T_min):
		# instance variables
		self.T_ini = T_ini
		self.T_min = T_min

	def temperature(self):
		return self.T_median
	
	# linear temperature profile from north to south
	def temperature_profile(self, x, t):

		return self.T_ini * time_factor_air(t)
