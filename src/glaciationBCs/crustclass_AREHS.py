# Data model of the evolving crust (thermal and mechanical farfield)
# Physical units: kg, m, s, K

import numpy as np

#from glaciationBCs.constants_AREHS import gravity

class crust():
	# class variables:
		
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
