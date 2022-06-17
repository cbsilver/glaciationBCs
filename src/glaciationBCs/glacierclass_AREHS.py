# Model of the evolving glacier extensions (length and height)
# parameterized analytical function for glacier geometry
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from glaciationBCs import constants_AREHS

from math import pi, sin, cos, sinh, cosh, sqrt, exp

class glacier():
	# class variables: owned by the class itself, static, shared by all class instances
	rho_ice = 900 #kg/m³
	rho_wat =1000 #kg/m³
	T_under = 273.15 + 0.5 #K
	
	# constructor
	def __init__(self, L_dom, L_max, H_max, x_0):
		# instance variables
		self.L_dom = L_dom
		self.L_max = L_max
		self.H_max = H_max
		self.x_0 = x_0

	def normalstress(self, x, t):
		return -self.rho_ice * gravity * self.local_height(x,t)
		
	def tangentialstress(self, x, t):
		return fricnum * self.normalstress(x, t)	

	def pressure(self, x, t):
		return -self.normalstress(x,t) 
    
	def temperature(self, x, t):
		return self.T_under
    
	# analytical function for the glacier's shape
	def local_height(self,x,y,z,t):
		# TODO coords = swap(coords) # y->x, z->y
	
		l = self.length(t)
		if l==0:
			return 0
		else:
			xi = (x-self.x_0) / l
			if xi<=1: 
				return self.height(t) * ((1 - (xi**2.5)**1.5))
			else: 
				#print("Warning: local coordinate must not be greater 1, but is ", xi)
				return 0

	# piecewise linear laws for the evolution of the glacier's dimensions
	def height(self, t):

		return self.H_max * time_factor_glacier(t)

	def length(self, t):

		return self.L_max * time_factor_glacier(t)

	# analytical function for the glacier meltwater production
	def local_meltwater(self,x,t):
		# constant flux at a temperate glacier base
		q = 6e-3 * 1 / s_a # = 6mm/a
		# constant flux at a frozen glacier base
		q = 0.0
		
		return q

	# auxiliary functions
	def print_max_load(self):
		print("Maximal normal stress due to glacier load: ")
		print(self.normalstress(0,self.t_1)/1e6, "MPa")
		
	def plot_evolution(self):
		tRange = np.linspace(self.t_1,self.t_0,11)
		fig,ax = plt.subplots()
		ax.set_title('Glacier evolution') #'Gletschervorschub'
		for t in tRange:
			xRange = np.linspace(self.x_0, self.x_0 + self.length(t),110)
			yRange = np.empty(shape=[0])
			for x in xRange:
				y = self.local_height(x,t)	
				yRange = np.append(yRange,y)
			ax.plot(xRange,yRange,label='t=$%.2f $ ' %t)
			ax.fill_between(xRange, 0, yRange)
		ax.set_xlabel('$x$ / m')
		ax.set_ylabel('height / m')
		ax.grid()
		fig.legend()
		plt.show()
	

