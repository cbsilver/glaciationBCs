# Model of the evolving glacier extensions (length and height)
# parameterized analytical function for glacier geometry
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt
from glaciationBCs import time_control_AREHS as tcr

from glaciationBCs.constants_AREHS import gravity
from glaciationBCs.constants_AREHS import rho_wat
from glaciationBCs.constants_AREHS import rho_ice
from glaciationBCs.constants_AREHS import s_a
from glaciationBCs.constants_AREHS import head_correction

class glacier():
	# class variables: owned by the class itself, static, shared by all class instances
	T_under = 273.15 + 0.5 #K
	fricnum = 0.2
	qf_melt = 6e-3 * 10 / s_a # = 6mm/a
	lT_tran = 200 #m
	lT_thaw = 200 #m

	# constructor
	def __init__(self, L_max, H_max, u_0, t_):
		# instance variables
		self.L_max = L_max
		self.H_max = H_max
		self.u_0 = u_0
		self.t_ = t_
		self.t_prev = 0.

		H_ = [0.0, 0.0, 0.0, 0.0, H_max, H_max, 0.0, 0.0]
		L_ = [0.0, 0.0, 0.0, 0.0, L_max, L_max, 0.0, 0.0]
		self.tcr_h = tcr.time_control(t_, H_)
		self.tcr_l = tcr.time_control(t_, L_)

	def normalstress(self, u, t):
		return -rho_ice * gravity * self.local_height(u, t)

	def tangentialstress(self, u, t):
		return self.fricnum * self.normalstress(u, t)

	def pressure(self, u, t):
		# reduced pressure/head to a reasonable amount
		return -self.normalstress(u,t) * head_correction

	def temperature(self, u, t):
		return self.T_under

	# analytical function for the glacier's shape
	def shape(self, u, t):
		l = self.length(t)
		if l==0:
			return 0
		# local coordinate starting from the North: s = u - u_0 (u_0=u_min)
		# local coordinate starting from the South: s = u_0 - u (u_0=u_max)
		s = self.u_0 - u
		xi = max(0., s/l)
		if xi<=1:
			return max(0., 1. - (xi**2.5))**1.5
		return 0

	def local_height(self, u, t):
		return self.height(t) * self.shape(u, t)

	# piecewise linear laws for the evolution of the glacier's dimensions
	def height(self, t):
		return self.tcr_h.function_value(t)

	def length(self, t):
		return self.tcr_l.function_value(t)

	# analytical function for the glacier meltwater production
	def local_meltwater(self, u ,t):
		# constant flux at a temperate glacier base
		q = self.qf_melt * self.shape(u, t)
		# constant flux at a frozen glacier base
		#q = 0.0
		return q

	# auxiliary functions
	def smoothstep (self, edge0, edge1, u):
		if u < edge0: return 0
		if u >= edge1: return 1
		#Scale/bias into [0..1] range
		xi = (u - edge0) / (edge1 - edge0)

		return xi*xi * (3 - 2*xi)

	def print_max_load(self):
		print("Maximal normal stress due to glacier load: ")
		print(self.normalstress(self.u_0, self.t_[5])/1e6, "MPa")

	def plot_evolving_shape(self):
		tRange = np.linspace(self.t_[6], self.t_[0],11)
		fig,ax = plt.subplots()
		ax.set_title('Glacier evolution') #'Gletschervorschub'
		for t in tRange:
			#uRange = np.linspace(self.u_0, self.u_0 + self.length(t), 110)
			uRange = np.linspace(self.u_0 - self.length(t), self.u_0 , 110)
			fRange = np.empty(shape=[0])
			for u in uRange:
				h = self.local_height(u, t)
				fRange = np.append(fRange, h)
			ax.plot(uRange, fRange, label='t=$%.2f $ ' %(t/s_a))
			#ax.fill_between(uRange, 0, fRange)
		ax.set_xlabel('$u$ / m')
		ax.set_ylabel('height / m')
		ax.grid()
		# fig.savefig("glacier_test.png")
		plt.legend(loc='upper left', bbox_to_anchor = (1.05, 1.0))
		plt.show()

	def plot_evolution(self):
		self.tcr_h.plot_evolution()
		self.tcr_l.plot_evolution()

	def plot_temperature(self, u_min):
		fig,ax = plt.subplots()
		tRange = np.linspace(self.t_[0], self.t_[6],7)
		t = tRange[4]
		t = self.t_[5]
		lg = self.length(t)
		u_max = self.u_0
		ug_tip = u_max - lg
		u_tran0 = ug_tip - self.lT_thaw - self.lT_tran
		u_tran1 = ug_tip - self.lT_thaw
		uRange = np.linspace(u_min, ug_tip, 500)
		fRange = np.empty(shape=[0])
		Tg = +0.5
		Ta = -1.5
		for u in uRange:
			f = self.smoothstep(u_tran0, u_tran1, u)
			fRange = np.append(fRange, f)
		TRange = Ta + fRange * (Tg - Ta)
		ax.plot(uRange, TRange, label='t=$%.2f $ ' %(t/s_a))
		ax.set_title('Top temperature')
		ax.set_xlabel('$u$ / m')
		ax.set_ylabel('$T$ / °C')
		ax.grid()
		# fig.savefig("glacier_test.png")
		plt.legend(loc='upper left', bbox_to_anchor = (1.05, 1.0))
		plt.show()