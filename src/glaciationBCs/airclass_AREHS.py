# Climate model of the evolving atmosphere
# parameterized analytical function for (cyclic) temperature evolution
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt
import time_control_AREHS as tcr

from constants_AREHS import s_a

class air():
	# class variables: owned by the class itself, shared by all instances of the class
	pressure  = 0.e3 #Pa
	
	def __init__(self, T_ini, T_min, t_):
		# instance variables
		self.T_ini = T_ini
		self.T_min = T_min
		
		self.t_ = t_		
		T_ = [T_ini, T_ini, T_min, T_min, T_min, T_min, T_ini]
		self.tcr = tcr.time_control(t_, T_)
	
	# linear temperature profile from north to south
	def temperature(self, t):
		return self.tcr.function_value(t)
		
	def plot_evolution(self):
		tRange = np.linspace(self.t_[0],self.t_[6],20)
		fRange = np.empty(shape=[0])
		#fRange = self.temperature(tRange)
		for t in tRange:
			f = self.temperature(t)	
			fRange = np.append(fRange,f)
		fig,ax = plt.subplots()
		ax.set_title('Temporal evolution')
		ax.plot(tRange / s_a, fRange)
		ax.scatter(tRange / s_a,fRange)
		ax.set_xlabel('$t$ / years')
		ax.set_ylabel('temperatur / m')
		ax.grid()
		plt.show()
