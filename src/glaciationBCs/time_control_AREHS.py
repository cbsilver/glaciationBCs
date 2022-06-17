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
	def __init__(self, t_0, t_1, t_2, t_3, t_4, t_5, t_6):
		# instance variables
		self.t_0 = t_0; #print("t0 = ", t_0)
		self.t_1 = t_1; #print("t1 = ", t_1)
		self.t_2 = t_2; #print("t2 = ", t_2)
		self.t_3 = t_3; #print("t3 = ", t_3)
		self.t_4 = t_4; #print("t4 = ", t_4)
		self.t_5 = t_5; #print("t5 = ", t_5)
		self.t_6 = t_6; #print("t6 = ", t_6)
		
	# def time modulation

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
		if (self.t_4 < t <= self.t_5):
			return stages[5]
		if (self.t_5 < t <= self.t_6):
			return stages[6]
		return "undefined stage"
	
	# TODO virtual method time_factor, to be overload by glacierclass and airclass
	
	# to multiply with process variable
	def time_factor_glacier(self, t):
		if (     0.0 < t <= self.t_3):
			return 0.0
		if (self.t_3 < t <= self.t_4):
			return (t-self.t_3) / (self.t_4-self.t_3)
		if (self.t_4 < t <= self.t_5):
			return 1.0
		if (self.t_5 < t <= self.t_6):
			return 1 - (t-self.t_5) / (self.t_6-self.t_5)
		return 0.0

	# to multiply with process variable TODO
	def time_factor_air(self, t):
		if (     0.0 < t <= self.t_1):
			return 1.0
		if (self.t_1 < t <= self.t_2):
			return 1 - (t-self.t_1) / (self.t_2-self.t_1)
		if (self.t_2 < t <= self.t_5):
			return 1.0
		if (self.t_5 < t <= self.t_6):
			return (t-self.t_5) / (self.t_6-self.t_5)
		return 0.0
