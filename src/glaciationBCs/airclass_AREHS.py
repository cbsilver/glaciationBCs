# Climate model of the evolving atmosphere
# parameterized analytical function for (cyclic) temperature evolution
# Physical units: kg, m, s, K

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from math import pi, sin, cos, sinh, cosh, sqrt, exp

class air():
	# class variables: owned by the class itself, shared by all instances of the class
	pressure  = 0.e3 #Pa
	hydrohead = 0.e3 #m
	
	def __init__(self, L_dom, T_north0, T_south0, T_rise, t_0, t_1, t_2=0, t_3=0, t_4=0):
		# instance variables
		self.L_dom = L_dom
		self.T_north0 = 266.15 #K = -7°C in the north (left)
		self.T_south0 = 276.15 #K = +3°C in the south (right)
		self.T_median = 0.5*(T_north0 + T_south0)
		self.T_rise = 8 #K
		self.T_drop = 0 #K
		self.t_0 = t_0
		self.t_1 = t_1
		self.t_2 = t_2
		self.t_3 = t_3
		self.t_4 = t_4

	def temperature(self):
		return self.T_median
	
	# linear temperature profile from north to south
	def temperature_profile(self, x, t):
		if (     0.0 < t <= self.t_0):# pre-glacial surface temperature decrease
			DT = self.T_drop * (1 - t / self.t_0)
			T_north = self.T_north0 + DT
			T_south = self.T_south0 + DT
		if (self.t_0 < t <= self.t_2):# glacial surface temperature distribution
			T_north = self.T_north0
			T_south = self.T_south0
		if (self.t_2 < t <= self.t_3):	# steadily rising surface temperature
			DT = self.T_rise * (t-self.t_2) / (self.t_3-self.t_2)
			T_north = self.T_north0 + DT
			T_south = self.T_south0 + DT
		if (self.t_3 < t <= self.t_4):	# interglacial surface temperature distribution
			T_north = self.T_north0 + self.T_rise
			T_south = self.T_south0 + self.T_rise

		return T_north + (T_south-T_north) * (x/self.L_dom)
