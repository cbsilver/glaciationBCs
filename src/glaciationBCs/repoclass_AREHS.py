# Data model of the evolving deep geological repository (thermal field)
# Physical units: kg, m, s, K

import copy as cp
import numpy as np
import matplotlib.pyplot as plt
import functools

from math import exp

class repo():
	# class variables:
	# TODO inner class farfield/deepfield -> self.farfield.ux_data
	# TODO members: vertical/horizontal_boundary -> boundary.ux
		
	# constructor
	def __init__(self):
		# instance variables: owned by instances of the class, can be different for each instance
				#Aus Jobmann et al., 2017 - Projekt Ansicht:
		#Parameter RK-BE
		self.BE_Q1 = 842.65
		self.BE_z1 = 3.58e-11

		self.BE_Q2 = 1269.66
		self.BE_z2 = 2.95e-10

		self.BE_Q3 = 3895.17
		self.BE_z3 = 9.12e-10

		self.BE_Q4 = 8308.18
		self.BE_z4 = 1.04e-8

		self.BE_Q5 = 42363.74
		self.BE_z5 = 2.62e-8

		#Parameter RK-HA:
		self.HA_Q1 = 7583.16
		self.HA_z1 = 8.28e-9

		self.HA_Q2 = 2412.91
		self.HA_z2 = 6.12e-10

		self.HA_Q3 = 2458.56
		self.HA_z3 = 7.15e-10

		self.HA_Q4 = 2546.25
		self.HA_z4 = 7.15e-10

		self.HA_Q5 = 7231.62
		self.HA_z5 = 8.28e-9

		#Gesamtvolumen:
		self.BE_vol = 7632 #m^3
		self.HA_vol = 1342.8 #m^3

		#Zeiten:
		self.t_zwischen_BE = 23*s_a
		self.t_zwischen_HA = 30*s_a
		self.t_voll = 80*s_a #Ende der stufenweise Einlagerung --> Start der Waermequelle im Modell
		
	def internal_heatsource(self,x,y,t):
		t = t + 80*s_a

		Q_BE = lambda t: (2/3) *( (BE_Q1*np.exp(-t*BE_z1)) 
				                + (BE_Q2*np.exp(-t*BE_z2)) 
				                + (BE_Q3*np.exp(-t*BE_z3)) 
				                + (BE_Q4*np.exp(-t*BE_z4)) 
				                + (BE_Q5*np.exp(-t*BE_z5)) )*BE_vol #(W/m^3 * m^3)


		Q_HA = lambda t: (2* 0.182/0.72) * ( 
				                  (HA_Q1*np.exp(-t*HA_z1)) 
				                + (HA_Q2*np.exp(-t*HA_z2)) 
				                + (HA_Q3*np.exp(-t*HA_z3)) 
				                + (HA_Q4*np.exp(-t*HA_z4)) 
				                + (HA_Q5*np.exp(-t*HA_z5)) )*HA_vol #(W/m^3 * m^3)


		#"Stufenweise" Einlagerung RK-BE
		Qsum_BE= lambda t:(0.05*Q_BE(t + t_zwischen_BE)   
				        + 0.05*Q_BE(t-4*s_a  + t_zwischen_BE)
				        + 0.05*Q_BE(t-8*s_a  + t_zwischen_BE) 
				        + 0.05*Q_BE(t-12*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-16*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-20*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-24*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-28*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-32*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-36*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-40*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-44*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-48*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-52*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-56*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-60*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-64*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-68*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-72*s_a + t_zwischen_BE) 
				        + 0.05*Q_BE(t-76*s_a + t_zwischen_BE) )


		Qsum_HA= lambda t:(0.05*Q_HA(t + t_zwischen_HA)   
				        + 0.05*Q_HA(t-4*s_a  + t_zwischen_HA)
				        + 0.05*Q_HA(t-8*s_a  + t_zwischen_HA) 
				        + 0.05*Q_HA(t-12*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-16*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-20*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-24*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-28*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-32*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-36*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-40*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-44*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-48*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-52*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-56*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-60*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-64*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-68*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-72*s_a + t_zwischen_HA) 
				        + 0.05*Q_HA(t-76*s_a + t_zwischen_HA) )


		Qsum = lambda t: Qsum_BE(t) + Qsum_HA(t)
		
		return (Qsum_BE(t) + Qsum_HA(t)) / 3000.0

