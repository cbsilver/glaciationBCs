# Data model of the evolving crust (thermal and mechanical farfield)
# Physical units: kg, m, s, K

import copy as cp
import numpy as np
import matplotlib.pyplot as plt
import functools

class crust():
	# class variables:
		
	# constructor
	def __init__(self):
		# instance variables: owned by instances of the class, can be different for each instance
		self.q_geo = 0.05
	
	def geothermal_heatflux(self):
		#TODO
		return 0
		
	def displacement_bottom(self,x,y,z,t):
		#TODO
		return 0.0
