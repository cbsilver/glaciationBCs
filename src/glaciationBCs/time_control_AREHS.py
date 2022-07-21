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
		# TODO
		return t

	def stage_control(self, t):
		print("t = %.1f years (%d s)" % (t/s_a, t))
		for i in range(6):
			if (self.t_[i] < t <= self.t_[i+1]):
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
		for i in range(6):
			if (self.t_[i] <= t <= self.t_[i+1]):
				return self.linear_function(t, self.t_[i], self.t_[i+1],
											   self.f_[i], self.f_[i+1])

	def plot_evolution(self):
		tRange = np.linspace(self.t_[0],self.t_[6],20)
		fRange = np.empty(shape=[0])
		for t in tRange:
			f = self.function_value(t)
			fRange = np.append(fRange,f)
		fig,ax = plt.subplots()
		ax.set_title('Temporal evolution')
		ax.plot(tRange / s_a, fRange)
		ax.scatter(np.array(self.t_) / s_a, np.array(self.f_))
		ax.set_xlabel('$t$ / years')
		ax.grid()
		plt.show()
