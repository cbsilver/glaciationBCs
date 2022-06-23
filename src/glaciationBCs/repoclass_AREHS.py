# Data model of the evolving deep geological repository (thermal field)
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt

from math import exp
from constants_AREHS import s_a

class repo():
	# class variables:
	# -
		
	# constructor
	def __init__(self, BE_Q, BE_z, HA_Q, HA_z, BE_vol, HA_vol, dgr_area, t_inter_BE, t_inter_HA, t_filled):
		# instance variables: owned by instances of the class, can be different for each instance
		# parameters RK-BE
		self.BE_Q = BE_Q
		self.BE_z = BE_z
		
		# parameters RK-HA:
		self.HA_Q = HA_Q
		self.HA_z = HA_z
		
		# total volume:
		self.BE_vol = BE_vol
		self.HA_vol = HA_vol
		
		self.dgr_area = dgr_area

		# times:
		self.t_inter_BE = t_inter_BE
		self.t_inter_HA = t_inter_HA
		self.t_filled = t_filled
		
	def radioactive_heatflow(self,t):

		Q_BE = lambda t: (2/3) *( (self.BE_Q[0]*np.exp(-t*self.BE_z[0])) 
				                + (self.BE_Q[1]*np.exp(-t*self.BE_z[1])) 
				                + (self.BE_Q[2]*np.exp(-t*self.BE_z[2])) 
				                + (self.BE_Q[3]*np.exp(-t*self.BE_z[3])) 
				                + (self.BE_Q[4]*np.exp(-t*self.BE_z[4])) ) * self.BE_vol

		Q_HA = lambda t: (2* 0.182/0.72) * ( 
				                  (self.HA_Q[0]*np.exp(-t*self.HA_z[0])) 
				                + (self.HA_Q[1]*np.exp(-t*self.HA_z[1])) 
				                + (self.HA_Q[2]*np.exp(-t*self.HA_z[2])) 
				                + (self.HA_Q[3]*np.exp(-t*self.HA_z[3])) 
				                + (self.HA_Q[4]*np.exp(-t*self.HA_z[4])) ) * self.HA_vol

		# stepwise filling of RK-BE
		Qsum_BE= lambda t:(0.05*Q_BE(t + self.t_inter_BE)   
				        + 0.05*Q_BE(t-4*s_a  + self.t_inter_BE)
				        + 0.05*Q_BE(t-8*s_a  + self.t_inter_BE) 
				        + 0.05*Q_BE(t-12*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-16*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-20*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-24*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-28*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-32*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-36*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-40*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-44*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-48*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-52*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-56*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-60*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-64*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-68*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-72*s_a + self.t_inter_BE) 
				        + 0.05*Q_BE(t-76*s_a + self.t_inter_BE) )

		Qsum_HA= lambda t:(0.05*Q_HA(t + self.t_inter_HA)   
				        + 0.05*Q_HA(t-4*s_a  + self.t_inter_HA)
				        + 0.05*Q_HA(t-8*s_a  + self.t_inter_HA) 
				        + 0.05*Q_HA(t-12*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-16*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-20*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-24*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-28*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-32*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-36*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-40*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-44*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-48*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-52*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-56*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-60*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-64*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-68*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-72*s_a + self.t_inter_HA) 
				        + 0.05*Q_HA(t-76*s_a + self.t_inter_HA) )

		Qsum = lambda t: Qsum_BE(t) + Qsum_HA(t)
		
		return Qsum(t) 					# heat flow = Wärmestrom

	def radioactive_heatflux(self,t):	# heat flux = Wärmestromdichte!
		# shift time to when the dgr is filled completely
		t = self.t_filled + t
		# TODO: remove sqrt for 3D
		return np.sqrt(self.radioactive_heatflow(t)) / self.dgr_area


	# auxiliary functions
	def print_max_load(self):
		print("Maximal heat flow from repository: ")
		print(self.radioactive_heatflow(self.t_filled), "W")
	
	def plot_evolution(self):
		time_ = np.linspace(self.t_filled , 10000*s_a, 1000)

		plt.figure(figsize=(12,6))
		plt.title('Gesamtwärmeleistung der Abfälle (RK-BE + RK-HA)\n nach vollständiger Einlagerung')
		plt.plot(time_/s_a, self.radioactive_heatflow(time_)/1000, label='RK-BE + RK-HA stufenweise befüllt', color = 'red', lw = 2.5)
		plt.semilogx()
		plt.xlim(10,10000)
		plt.ylim(0,13000)
		plt.xlabel('Zeit [a]')
		plt.ylabel('Leistung [kW]')
		plt.legend(loc='upper right')
		plt.grid(True)
		plt.show()

