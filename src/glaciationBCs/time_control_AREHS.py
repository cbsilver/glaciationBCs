# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class time_control():
	# class variables: owned by the class itself, static, shared by all class instances
	stages = {
	# t0->t1
	0 : "initialization phase", # for steady state / equilibrium geothermal heat flux
	# t1->t2
	1 : "permafrost development", # in the glacial (cold) period
	# t2->t3
	2 : "permafrost-only period", # in the glacial (cold) period
	# t3->t4
    3 : "glacier advance",
    # t4->t5
    5 : "glacier dormancy",
    # t4->t5
    5 : "glacier retreat",
    # t5->t6
    6 : "interglacial period", # warm period
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
		if (self.t_4 < t <= self.t_5):
			return stages[5]
		if (self.t_5 < t <= self.t_6):
			return stages[6]
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
		

