# Model of the evolving glacier extensions (length and height)
# parameterized analytical function for glacier geometry
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from math import pi, sin, cos, sinh, cosh, sqrt, exp

# Physical constants in units: kg, m, s, K
gravity = 9.81 #m/s²
fricnum = 0.2
b_sub = 0.300 # dim-less proportionality constant glacier height - subsidence
b_reb = 0.225 # dim-less proportionality constant glacier height - rebounding

# Numerical constants
s_a = 365.25*24*3600 #=31557600 seconds per year
eps = 1.e-8

stages = {
	0 : "pre-glaciation period", # in the glacial (cold) period
    1 : "glacier advance",
    2 : "glacier dormancy",
    3 : "glacier retreat",
    4 : "post-glaciation period", # interglacial (warm) period
    }

# Optional: define deflection mode
deflection_modes = {
	0 : "heuristic approximation from glacier height (Bense)",
    1 : "analytical smooth elastic bending line (Silbermann)",
    }
deflection_mode = 0


class glacier():
	# class variables: owned by the class itself, static, shared by all class instances
	rho_ice = 900 #kg/m³
	rho_wat =1000 #kg/m³
	T_under = 273.15 + 0.5 #K
	
	# constructor
	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2=0, t_3=0, t_4=0):
		# instance variables
		self.L_dom = L_dom
		self.L_max = L_max
		self.H_max = H_max
		self.x_0 = x_0
		self.t_0 = t_0; #print("t0 = ", t_0)
		self.t_1 = t_1; #print("t1 = ", t_1)
		self.t_2 = t_2; #print("t2 = ", t_2)
		self.t_3 = t_3; #print("t3 = ", t_3)
		self.t_4 = t_4; #print("t4 = ", t_4)

	def stagecontrol(self, t):
		print("t = ", t)
		if (     0.0 < t <= self.t_0):
			return stages[0]
		if (self.t_0 < t <= self.t_1):
			return stages[1]
		if (self.t_1 < t <= self.t_2):
			return stages[2]
		if (self.t_2 < t <= self.t_3):
			return stages[3]
		if (self.t_3 < t <= self.t_4):
			return stages[4]
		return "undefined glacier stage"

	def normalstress(self, x, t):
		return -self.rho_ice * gravity * self.local_height(x,t)
		
	def tangentialstress(self, x, t):
		return fricnum * self.normalstress(x, t)	
    
	def hydrohead(self, x, t):
		# equals self.pressure(x,t) / (self.rho_wat * gravity)
		return self.rho_ice/self.rho_wat * self.local_height(x,t)

	def pressure(self, x, t):
		return -self.normalstress(x,t) 
    
	def temperature(self, x, t):
		return self.T_under
    
	# analytical function for the glacier's shape
	def local_height(self,x,t):
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
	
	def local_height_rate(self,x,t):
		h = self.height(t)
		l = self.length(t)
		if l==0:
			return 0
		else:
			xi = (x-self.x_0) / l
			if xi<=1:
				doth = self.height_rate(t)
				dotl = self.length_rate(t)
				part1 = doth / h * ((1 - (xi**2.5)**1.5))
				part2 = dotl / l * 15/4.0 * ((1 - (xi**2.5)**0.5)) * (xi**2.5)
				return h * (part1 + part2)
			else:
				#print("Warning: local coordinate must not be greater 1, but is ", xi)
				return 0

	# piecewise linear laws for the evolution of the glacier's dimensions
	def height(self, t):
		if (     0.0 < t <= self.t_0):
			return 0.0
		if (self.t_0 < t <= self.t_1):
			return self.H_max * (t-self.t_0) / (self.t_1-self.t_0)
		if (self.t_1 < t <= self.t_2):
			return self.H_max
		if (self.t_2 < t <= self.t_3):
			return self.H_max * (1 - (t-self.t_2) / (self.t_3-self.t_2))
		if (self.t_3 < t <= self.t_4):
			return 0.0
		return 0.0

	def height_rate(self, t):
		if (     0.0 < t <= self.t_0):
			return 0.0
		if (self.t_0 < t <= self.t_1):
			return self.H_max / (self.t_1-self.t_0)
		if (self.t_1 < t <= self.t_2):
			return 0.0
		if (self.t_2 < t <= self.t_3):
			return self.H_max / (-(self.t_3-self.t_2))
		if (self.t_3 < t <= self.t_4):
			return 0.0
		return 0.0

	def length(self, t):
		if (     0.0 < t <= self.t_0):
			return 0.0
		if (self.t_0 < t <= self.t_1):
			return self.L_max * (t-self.t_0) / (self.t_1-self.t_0)
		if (self.t_1 < t <= self.t_2):
			return self.L_max
		if (self.t_2 < t <= self.t_3):
			return self.L_max * (1 - (t-self.t_2) / (self.t_3-self.t_2))
		if (self.t_3 < t <= self.t_4):
			return 0.0
		return 0.0

	def length_rate(self, t):
		if (     0.0 < t <= self.t_0):
			return 0.0
		if (self.t_0 < t <= self.t_1):
			return self.L_max / (self.t_1-self.t_0)
		if (self.t_1 < t <= self.t_2):
			return 0.0
		if (self.t_2 < t <= self.t_3):
			return self.L_max / (-(self.t_3-self.t_2))
		if (self.t_3 < t <= self.t_4):
			return 0.0
		return 0.0

	# analytical function for the lithosphere deflection due to glacier load
	def local_deflection(self,x,t):
		if (deflection_mode == 0):
			deflection = self.local_deflection_heuristic(x,t)
		if (deflection_mode == 1):
			deflection = self.local_deflection_elastic(x,t)
		return deflection
		
	# heuristic approximation from glacier height (Bense)
	def local_deflection_rate_heuristic(self,x,t):
		t_relax = (self.t_4-self.t_0) / 13 #~2500 a
		h0 = self.local_height(x,self.t_0)
		uy_max = -b_sub * (self.local_height(x,self.t_1) - h0)
		uy_med = uy_max + b_reb * (self.local_height(x,self.t_2) - self.local_height(x,self.t_3))
		# time point when the ice has locally completely retreated
		t_startPG = self.t_3 - (self.t_3-self.t_2) * (x-self.x_0) / self.L_max

		vy = 0.0
		if (self.t_0 < t <= self.t_1): #immediate deflection
			vy = -b_sub * self.local_height_rate(x,t)
		if (self.t_1 < t <= self.t_2): #constant subsidence
			vy = 0.0
		if (self.t_2 < t <= self.t_3): #restrained rebound
			vy = -b_reb * self.local_height_rate(x,t)
		if (t > t_startPG): #postglacial rebound (retarded)
			dt = t - t_startPG
			# process starts immediately after local post-glaciation
			vy = - uy_med * exp(-dt/t_relax) / t_relax
		return vy

	def local_deflection_heuristic(self,x,t):
		t_relax = (self.t_4-self.t_0) / 13 #~2500 a
		#t_relax = 2500 * 31557600 #s
		h0 = self.local_height(x,self.t_0)
		uy_max = -b_sub * (self.local_height(x,self.t_1) - h0)
		uy_med = uy_max + b_reb * (self.local_height(x,self.t_2) - self.local_height(x,self.t_3))
		# time point when the ice has locally completely retreated
		t_startPG = self.t_3 - (self.t_3-self.t_2) * (x-self.x_0) / self.L_max
		
		uy = 0.0
		if (self.t_0 < t <= self.t_1): #immediate deflection
			uy = -b_sub * (self.local_height(x,t) - h0)
		if (self.t_1 < t <= self.t_2): #constant subsidence
			uy = uy_max
		if (self.t_2 < t <= self.t_3): #restrained rebound
			uy = uy_max + b_reb * (self.local_height(x,self.t_2) - self.local_height(x,t))
		if (t > t_startPG): #postglacial rebound (retarded)
			dt = t - t_startPG
			# process starts immediately after local post-glaciation
			uy = uy_med * exp(-dt/t_relax)
		return uy

	def local_displacement_heuristic(self,x,y,t):
		# TODO: move constants on top
		eps_yy = -0.0005
		H_dom = 150000 #TODO scaling with 20?
		uy_compaction = (H_dom-(-y)) * eps_yy
		uy_deflection = self.local_deflection_heuristic(x,t)
		
		h0 = self.local_height(x,self.t_0)
		uy_max = -b_sub * (self.local_height(x,self.t_1) - h0)
		
		if (abs(uy_max)>0.0):
			uy = (1 + uy_compaction/uy_max) * uy_deflection
		else:
			uy = uy_deflection
		return uy

	# analytical function for the lithosphere deflection due to glacier load
	def local_deflection_elastic(self,x,t):
		B = 1 #m
		T = 7500 #m
		EI = 10e9 * B*T**3 / 12 * 1e7
		xG = self.length(t)
		L  = self.L_dom
		qG = self.height(t) * self.rho_ice * gravity
		qM = 0.9*qG * xG/L
		uy = 0.0
		if (x-self.x_0 < xG):
			x1 = x - self.x_0
			part1 = qM * ( (L-xG)**3 * (x1/6-L/8-xG/24) )
			part2 = (qG-qM)/2 * (x1**4/12 - (L-xG)*L*xG*x1 - xG**3*x1/3 + 2/3*L**3*xG - 1/2*L**2*xG**2 + xG**4/12)
			uy = -(part1 + part2) / EI
		else:
			x2 = x - self.x_0 - xG
			part1 = qM * ( (L-xG)**3 * (x2/6-(L-xG)/8) - x2**4/24 )
			part2 = (qG-qM)*xG/2 * ( x2**2*(xG/2+x2/3) + (L-xG)*( (L-xG)*(2/3*L - xG/6) - L*x2) )
			uy = -(part1 + part2) / EI
		
		return uy

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
		ax.set_title('Glacier evolution') #'Gletschervorstoß'
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
	
	def check_deflection_continuity(self):
		print("Deflections: ")
		#print(self.local_deflection(self.L_dom,self.t_1), "m")
		xG = self.length(self.t_1)
		print('xG = ', xG)
		print('uy(xG-1)= ',self.local_deflection(xG-1,self.t_1), "m")
		print('uy(xG+1)= ',self.local_deflection(xG+1,self.t_1), "m")
		
	def plot_deflection(self):
		tRange = np.linspace(self.t_4,10*self.t_4,11)
		fig,ax = plt.subplots()
		ax.set_title('Crustal deflection')
		for t in tRange:
			xRange = np.linspace(self.x_0, self.x_0 + self.L_dom,110)
			yRange = np.empty(shape=[0])
			for x in xRange:
				if (deflection_mode == 0):
					y = self.local_deflection_heuristic(x,t)
				if (deflection_mode == 1):
					y = self.local_deflection_elastic(x,t)
				yRange = np.append(yRange,y)
			ax.plot(xRange,yRange,label='t=$%.2f $ ' %(t/s_a/1000))
		ax.set_xlabel('$x$ / m')
		ax.set_ylabel('deflection / m')
		ax.grid()
		fig.legend()
		plt.show()
	
	def plot_deflection_rate(self):
		tRange = np.linspace(self.t_3,self.t_4,11)
		fig,ax = plt.subplots()
		ax.set_title('Crustal deflection rate')
		for t in tRange:
			xRange = np.linspace(self.x_0, self.x_0 + self.L_max,110)
			yRange = np.empty(shape=[0])
			for x in xRange:
				if (deflection_mode == 0):
					y = self.local_deflection_rate_heuristic(x,t)
				yRange = np.append(yRange,y)
			ax.plot(xRange,yRange,label='t=$%.2f $ ' %(t))
		ax.set_xlabel('$x$ / m')
		ax.set_ylabel('deflection rate / m/a')
		ax.grid()
		fig.legend()
		plt.show()
