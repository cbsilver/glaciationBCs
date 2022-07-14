# Data model of the evolving deep geological repository (thermal field)
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt

from math import exp
from glaciationBCs.constants_AREHS import s_a

class repo():
	# class variables:
	# -

	# constructor
	def __init__(self, BE_Q, BE_z, BE_f, HA_Q, HA_z, HA_f, BE_vol, HA_vol, dgr_area, t_inter_BE, t_inter_HA, t_filled):
		# instance variables: owned by instances of the class, can be different for each instance
		# parameters RK-BE
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

		self.dgr_area = dgr_area

		# times:
		self.t_inter_BE = t_inter_BE
		self.t_inter_HA = t_inter_HA
		self.t_filled = t_filled

	def radioactive_heatflow(self,t): # heat flow = Wärmestrom

		Q_BE = lambda t: np.sum(self.BE_Q[0:4]*np.exp(-t@self.BE_z[np.newaxis,0:4]))
		Q_HA = lambda t: np.sum(self.HA_Q[0:4]*np.exp(-t@self.HA_z[np.newaxis,0:4]))
		# stepwise filling of RK-BE
		Qsum_BE = np.sum(Q_BE(t-np.arange(0,80,4)[:,np.newaxis]*s_a + self.t_inter_BE))
		Qsum_HA = np.sum(Q_HA(t-np.arange(0,80,4)[:,np.newaxis]*s_a + self.t_inter_HA))

		return Qsum_BE*self.BE_vol*self.BE_f + Qsum_HA*self.HA_vol*self.HA_f

	def radioactive_heatflux(self,t): # heat flux = Wärmestromdichte!
		# shift time to when the dgr is filled completely
		t = self.t_filled + t
		# TODO: remove sqrt for 3D
		return np.sqrt(self.radioactive_heatflow(t)) / self.dgr_area #TODO line start and endpoint


	# auxiliary functions
	def print_max_load(self):
		print("Maximal heat flow from repository: ")
		print(self.radioactive_heatflow(self.t_filled), "W")

	def plot_evolution(self):
		time_ = np.linspace(self.t_filled , 10000*s_a, 1000)

		plt.figure(figsize=(12,6))
		plt.title('Gesamtwärmeleistung der Abfälle (RK-BE + RK-HA)\n nach vollständiger Einlagerung')
		plt.plot(time_/s_a, [self.radioactive_heatflow(ti)/1000 for ti in time_],
			label='RK-BE + RK-HA stufenweise befüllt', color = 'red', lw = 2.5)
		plt.semilogx()
		plt.xlim(10,10000)
		plt.ylim(0,13000)
		plt.xlabel('Zeit [a]')
		plt.ylabel('Leistung [kW]')
		plt.legend(loc='upper right')
		plt.grid(True)
		plt.savefig("repo_test.png")
		# plt.show()
