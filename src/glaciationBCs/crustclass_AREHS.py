# Data model of the evolving crust (thermal, hydraulic and mechanical farfield)
# Physical units: kg, m, s, K

import numpy as np
import matplotlib.pyplot as plt

from glaciationBCs.constants_AREHS import gravity
from glaciationBCs.constants_AREHS import rho_wat
from glaciationBCs.constants_AREHS import c_p_wat

class crust():
	# class variables:
	V_fluid_max = 1e-11	#m/s

	# constructor
	def __init__(self, q_geo, v_min, v_max, T_ini, T_bot):
		# instance variables: owned by instances of the class, can be different for each instance
		self.q_geo = q_geo
		self.v_min = v_min
		self.v_max = v_max
		self.T_bot = T_bot
		self.T_ini = T_ini

	def geothermal_heatflux(self):
		return [0.0, self.q_geo, 0.0]

	def displacement_below(self):
		return [0.0, 0.0, 0.0]

	def displacement_aside(self):
		return [0.0, 0.0, 0.0]

	def geothermal_temperature(self, v, T_atm):
		# linear profile according to geothermal heatflux
		DT = self.T_bot - T_atm
		Dv = (self.v_min - self.v_max)
		return DT/Dv * (v - self.v_max) + T_atm

	def lateral_heatflux(self, v, T_atm, props):
		# only for salt models
		if props.model_id == 3:
			for i, lv in enumerate(props.south_layer_bounds[:-1]):
				if lv >= v > props.south_layer_bounds[i+1]:
					if props.rock_type_array[i] == "salt":
						return 0.
					break
		# linear profile according to decreasing fluid velocity		
		Dv = (self.v_min - self.v_max)
		V_fluid = self.V_fluid_max * (v - self.v_max) / Dv
		q_heat = V_fluid * (self.T_ini - T_atm) * rho_wat * c_p_wat
		return q_heat

	def hydrostatic_pressure(self, v):
		# linear profile according to gravity
		p_pore = rho_wat * gravity * (self.v_max - v)
		return p_pore

	def lithostatic_stresses(self, v, props):
		heights = np.abs(np.diff(props.south_layer_bounds))
		rho_eff = (1. - props.poro_array) * props.rho_array + \
					(props.poro_array - props.biot_array) * 1000.
		layer_stress = rho_eff * gravity * heights
		total_stress = np.append(0, np.add.accumulate(layer_stress[:-1]))
		stress = [0., 0., 0.]
		for i, (ls, lh, lv, ts, nu, ab) in enumerate(zip(layer_stress, heights, props.south_layer_bounds[:-1], total_stress, props.nu_array, props.biot_array)):
			if lv >= v >= props.south_layer_bounds[i+1]:
				stress = np.array([nu / (1. - nu), 1., nu / (1. - nu)]) * ls/lh * (lv - v) + ts
				break
		return -stress

	def plot_profile(self, T_atm):
		vRange = np.linspace(self.v_min,self.v_max,20)
		#fRange = self.hydrostatic_pressure(vRange)
		fRange = self.lateral_heatflux(vRange, T_atm)
		fig,ax = plt.subplots()
		ax.set_title('Vertical profile')
		ax.plot(fRange, vRange)
		ax.set_xlabel('$p$ / Pa')
		ax.set_ylabel('$y$ / m')
		ax.grid()
		plt.show()

	def plot_lithostatic_stress(self):
		vRange = np.linspace(0,-1000, 100)
		fRange = [1e-6*self.lithostatic_stresses(v, props)[0] for v in vRange]
		fig,ax = plt.subplots()
		ax.set_title('Vertical profile')
		ax.plot(vRange, fRange)
		ax.set_ylabel('$sigma$ / MPa')
		ax.set_xlabel('$y$ / m')
		ax.grid()
		plt.show()

	def plot_profile_evolution(self):
		vRange = np.linspace(self.v_min,self.v_max,20)
		TRange = np.linspace(self.T_ini,self.T_ini-10,10)
		fig,ax = plt.subplots()
		for T_atm in TRange:
			fRange = self.lateral_heatflux(vRange, T_atm)
			ax.plot(fRange, vRange, label='T_atm=$%.2f $ ' %(T_atm))
		ax.set_title('Vertical profile')
		ax.set_xlabel('$q_x$ / W/m²')
		ax.set_ylabel('$v$ / m')
		ax.grid()
		fig.legend()
		plt.show()
