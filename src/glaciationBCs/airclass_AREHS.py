# Climate model of the evolving atmosphere
# parameterized analytical function for (cyclic) temperature evolution
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt

from glaciationBCs import time_control_AREHS as tcr
from glaciationBCs.constants_AREHS import s_a

class air():
	# class variables: owned by the class itself, shared by all instances of the class
	pressure  = 0.e3 #Pa
	
	def __init__(self, T_ini, T_min, t_):
		# instance variables
		self.T_ini = T_ini
		self.T_min = T_min
		
		self.t_ = t_		
		T_ = [T_ini, T_ini, T_min, T_min, T_min, T_min, T_ini]
		self.tcr = tcr.time_control(t_, T_)
	
	# linear temperature profile from north to south
	def temperature(self, t):
		return self.tcr.function_value(t)
		
	def plot_evolution(self):
		self.tcr.plot_evolution()
