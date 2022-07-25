# Auxiliary class for piecewise linear time functions
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt

from glaciationBCs.constants_AREHS import *

class time_control():
	# class variables: owned by the class itself, static, shared by all class instances
	stages = {
	# t0->t1
	0 : "interglacial period",	  # interglacial (warm) period, steady state
	# t1->t2
	1 : "permafrost development", # part glacial (cold) period
	# t2->t3
	2 : "permafrost-only period", # part glacial (cold) period
	# t3->t4
	3 : "glacier advance",		  # part glacial (cold) period
	# t4->t5
	4 : "glacier dormancy",		  # part glacial (cold) period
	# t5->t6
	5 : "glacier retreat",		  # interglacial (warm) period
	# t6->... repeat cycle
	}

	# constructor
	def __init__(self, t_, f_):
		# instance variables
		self.t_ = t_
		self.f_ = f_

	def time_modulation(self, t):
		return t%self.t_[6]

	def stage_control(self, t):
		print("t = %.1f years (%d s)" % (t/s_a, t))
		for i in range(6)[:-1]:
			if (self.time_modulation(t) >= self.t_[i]):
				return self.stages[i]
		return "undefined stage"

	def linear_function(self, t, t_S, t_E, f_S, f_E):
		if (f_E==f_S):
			return f_S
		else:
			Df = f_E - f_S
			Dt = t_E - t_S
			return Df/Dt * (t-t_S) + f_S

	def function_value(self, t):
		t_mod = self.time_modulation(t)
		for i in range(6)[::-1]:
			if (t_mod >= self.t_[i]):
				return self.linear_function(t_mod, self.t_[i], self.t_[i+1],
											self.f_[i], self.f_[i+1])

	def plot_evolution(self):
		tRange = np.ravel([np.array(self.t_[:-1])+i*self.t_[6] for i in range(5)])
		fRange = [self.function_value(t) for t in tRange]
		fig,ax = plt.subplots()
		ax.set_title('Temporal evolution')
		ax.plot(tRange / s_a, fRange, "-o")
		ax.set_xlabel('$t$ / years')
		ax.grid()
		plt.show()
