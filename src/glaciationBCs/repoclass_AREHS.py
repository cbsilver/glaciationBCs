# Data model of the evolving deep geological repository (thermal field)
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt

from glaciationBCs.constants_AREHS import s_a

class repo():
	# class variables:
	# -

	# constructor
	def __init__(self, BE_Q, BE_z, BE_f, 
					   HA_Q, HA_z, HA_f, 
					   BE_vol, HA_vol, repo_size,
					   t_inter_BE, t_inter_HA, t_filled):
		# instance variables: owned by instances of the class, can be different for each instance
		# parameters RK-BE:
		self.BE_Q = BE_Q
		self.BE_z = BE_z
		self.BE_f = BE_f

		# parameters RK-HA:
		self.HA_Q = HA_Q
		self.HA_z = HA_z
		self.HA_f = HA_f

		# total volume:
		self.BE_vol = BE_vol
		self.HA_vol = HA_vol

		self.repo_size = repo_size

		# times:
		self.t_inter_BE = t_inter_BE
		self.t_inter_HA = t_inter_HA
		self.t_filled = t_filled
		self.t_prev = 0.
		
	def radioactive_heatflow(self,t_sim): # heat flow = Wärmestrom

		Q_BE = lambda t: np.sum(self.BE_Q[:]*np.exp(-t@self.BE_z[np.newaxis,:]))
		Q_HA = lambda t: np.sum(self.HA_Q[:]*np.exp(-t@self.HA_z[np.newaxis,:]))
		# stepwise filling of repository in steps of 4 years
		# time array which masks out containers which are not yet stored
		t_masked = t_sim-np.arange(0, self.t_filled, 4*s_a)
		t_masked = t_masked[t_masked >= 0]
		Qsum_BE = np.sum(Q_BE(t_masked[:, np.newaxis] + self.t_inter_BE))
		Qsum_HA = np.sum(Q_HA(t_masked[:, np.newaxis] + self.t_inter_HA))

		return Qsum_BE*self.BE_vol*self.BE_f + Qsum_HA*self.HA_vol*self.HA_f

	def radioactive_heatflux(self,t): # heat flux = Wärmestromdichte!
		# use t + self.t_filled to disregard initial heat production from filling
		return self.radioactive_heatflow(t) / self.repo_size


	# auxiliary functions
	def print_max_load(self):
		print("Maximal heat flow from repository: ")
		print(self.radioactive_heatflow(self.t_filled), "W")

	def plot_evolution(self):
		time_ = np.linspace(0, 800, 201)*s_a

		plt.figure(figsize=(12,6))
		plt.title('Gesamtwärmeleistung der Abfälle (RK-BE + RK-HA)\n nach vollständiger Einlagerung')
		plt.plot(time_/s_a, [self.radioactive_heatflow(ti)/1000 for ti in time_],
			label='RK-BE + RK-HA stufenweise befüllt', color = 'red', lw = 2.5)
		#plt.xlim(0,10000)
		#plt.ylim(0,13000)
		plt.xlabel('Zeit [a]')
		plt.ylabel('Leistung [kW]')
		plt.legend(loc='upper right')
		plt.grid(True)
		#plt.savefig("repo_test.png")
		# plt.show()
