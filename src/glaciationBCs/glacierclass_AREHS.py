# Model of the evolving glacier extensions (length and height)
# parameterized analytical function for glacier geometry
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt
from glaciationBCs import time_control_AREHS as tcr

from glaciationBCs.constants_AREHS import gravity
from glaciationBCs.constants_AREHS import s_a

class glacier():
	# class variables: owned by the class itself, static, shared by all class instances
	rho_ice = 900 #kg/m³
	rho_wat =1000 #kg/m³
	T_under = 273.15 + 0.5 #K
	fricnum = 0.2
	qf_melt = 6e-3 * 1 / s_a # = 6mm/a
	
	# constructor
	def __init__(self, L_dom, L_max, H_max, x_0, t_):
		# instance variables
		self.L_dom = L_dom
		self.L_max = L_max
		self.H_max = H_max
		self.x_0 = x_0
		self.t_ = t_
		
		H_ = [0.0, 0.0, 0.0, 0.0, H_max, H_max, 0.0]
		L_ = [0.0, 0.0, 0.0, 0.0, L_max, L_max, 0.0]
		self.tcr_h = tcr.time_control(t_, H_)
		self.tcr_l = tcr.time_control(t_, L_)
		
	def normalstress(self, x, t):
		return -self.rho_ice * gravity * self.local_height(x,t)
		
	def tangentialstress(self, x, t):
		return self.fricnum * self.normalstress(x, t)	

	def pressure(self, x, t):
		return -self.normalstress(x,t) 
    
	def temperature(self, x, t):
		return self.T_under
    
	# analytical function for the glacier's shape
	def local_height(self,x,t):
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
		return self.tcr_h.function_value(t)

	def length(self, t):
		return self.tcr_l.function_value(t)

	# analytical function for the glacier meltwater production
	def local_meltwater(self,x,t):
		# constant flux at a temperate glacier base
		q = qf_melt
		# constant flux at a frozen glacier base
		q = 0.0
		
		return q

	# auxiliary functions
	def print_max_load(self):
		print("Maximal normal stress due to glacier load: ")
		print(self.normalstress(self.x_0,self.t_[5])/1e6, "MPa")
		
	def plot_evolving_shape(self):
		tRange = np.linspace(self.t_[6],self.t_[0],11)
		fig,ax = plt.subplots()
		ax.set_title('Glacier evolution') #'Gletschervorschub'
		for t in tRange:
			xRange = np.linspace(self.x_0, self.x_0 + self.length(t),110)
			yRange = np.empty(shape=[0])
			for x in xRange:
				y = self.local_height(x,t)	
				yRange = np.append(yRange,y)
			ax.plot(xRange,yRange,label='t=$%.2f $ ' %(t/s_a))
			#ax.fill_between(xRange, 0, yRange)
		ax.set_xlabel('$x$ / m')
		ax.set_ylabel('height / m')
		ax.grid()
		fig.legend()
		# fig.savefig("glacier_test.png")
		plt.show()
	
	def plot_evolution(self):
		self.tcr_h.plot_evolution()
		self.tcr_l.plot_evolution()

