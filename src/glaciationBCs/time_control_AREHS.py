# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt

from glaciationBCs.constants_AREHS import *

class time_control():
	# class variables: owned by the class itself, static, shared by all class instances
	# TODO remove initialization phase
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
		if (self.t_[0] < t <= self.t_[1]):
			return self.stages[0]
		if (self.t_[1] < t <= self.t_[2]):
			return self.stages[1]
		if (self.t_[2] < t <= self.t_[3]):
			return self.stages[2]
		if (self.t_[3] < t <= self.t_[4]):
			return self.stages[3]
		if (self.t_[4] < t <= self.t_[5]):
			return self.stages[4]
		if (self.t_[5] < t <= self.t_[6]):
			return self.stages[5]
		return "undefined stage"

	def linear_function(self, t, t_S, t_E, f_S, f_E):
		if (f_E==f_S):
			return f_S
		else:
			Df = f_E - f_S
			Dt = t_E - t_S
			return Df/Dt * (t-t_S) + f_S

	def function_value(self, t):
		# TODO for loop over all stages
		if (self.t_[0] <= t <= self.t_[1]):
			return self.linear_function(t, self.t_[0], self.t_[1], self.f_[0], self.f_[1])
		if (self.t_[1] <  t <= self.t_[2]):
			return self.linear_function(t, self.t_[1], self.t_[2], self.f_[1], self.f_[2])
		if (self.t_[2] <  t <= self.t_[3]):
			return self.linear_function(t, self.t_[2], self.t_[3], self.f_[2], self.f_[3])
		if (self.t_[3] <  t <= self.t_[4]):
			return self.linear_function(t, self.t_[3], self.t_[4], self.f_[3], self.f_[4])
		if (self.t_[4] <  t <= self.t_[5]):
			return self.linear_function(t, self.t_[4], self.t_[5], self.f_[4], self.f_[5])
		if (self.t_[5] <  t <= self.t_[6]):
			return self.linear_function(t, self.t_[5], self.t_[6], self.f_[5], self.f_[6])

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
