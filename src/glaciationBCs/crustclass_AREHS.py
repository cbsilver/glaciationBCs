# Data model of the evolving crust (displacement farfield)
# Data must be provided (e.g. from some GIA simulation)
# Physical units: kg, m, s, K

import copy as cp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import functools

from math import pi, sin, cos, sinh, cosh, sqrt, exp

# Numerical constants
s_a = 365.25*24*3600 #=31557600 seconds per year
eps = 1.e-8

# helper functions
def closest(lst, K):
	idx = (np.abs(lst - K)).argmin()
	return lst[idx], idx #tuple!
	
def index_range_yfixed(idx_y, Nx):
    iA = Nx * idx_y
    iE = Nx *(idx_y+1) - 1
    return iA, iE
    
def index_range_Yfixed(y, yvals):
    index_range = np.where(yvals==y)[0][:]
    return index_range.min(), index_range.max()

def index_list_Xfixed(x, xvals):
    index_list = np.where(xvals==x)[0][:]
    return index_list

def setup_lineplot_uxuy(coord, value):
	fig, ax = plt.subplots(ncols=2,figsize=(24,6))
	ax[0].set_title('Horizont displacement for ' + coord + "=%.0f"%(value))
	ax[1].set_title('Vertical displacement for ' + coord + "=%.0f"%(value))
	ax[0].set_ylabel('$u_x$ / m')
	ax[1].set_ylabel('$u_y$ / m')
	for i in [0,1]:
	    ax[i].grid()
	return ax

# read external field data from GIA (Glacier Isostatic Adjustment)
def read_data_GIA(datapath='data/'):
	ux_data = pd.read_csv(datapath + 'ux.dat', sep=r'\s{1,}', engine='python')
	uy_data = pd.read_csv(datapath + 'uy.dat', sep=r'\s{1,}', engine='python')
	tSeries = pd.read_csv(datapath + 'ux.dat', sep=r'\s{1,}', engine='python', 
			  header=None, nrows=1, usecols=range(2,68))
	ux_data.info()
	uy_data.info()
	tSeries.info()
	return ux_data, uy_data, tSeries


class crust():
	# class variables:
	# TODO inner class farfield/deepfield -> self.farfield.ux_data
	# TODO members: vertical/horizontal_boundary -> boundary.ux
		
	# constructor
	def __init__(self, datapath='data/'):
		# instance variables: owned by instances of the class, can be different for each instance
		self.datapath = datapath
		self.ux_data = pd.read_csv(self.datapath + 'ux.dat', sep=r'\s{1,}', engine='python')
		self.uy_data = pd.read_csv(self.datapath + 'uy.dat', sep=r'\s{1,}', engine='python')
		self.tSeries = pd.read_csv(self.datapath + 'uy.dat', sep=r'\s{1,}', engine='python', 
				  header=None, nrows=1, usecols=range(2,68))
		self.xvalues = self.ux_data.x[0:].to_numpy(dtype=float)
		self.yvalues = self.uy_data.y[0:].to_numpy(dtype=float)
		self.tvalues = self.tSeries.to_numpy(dtype=float)[0,:]
		
		self.xmin = self.xvalues.min()
		self.xmax = self.xvalues.max()
		self.length = self.xmax - self.xmin
		self.ymin = self.yvalues.min()
		self.ymax = self.yvalues.max()
		self.depth = self.ymax - self.ymin
		self.idLst_vertical = index_list_Xfixed(x=1150000, xvals= self.xvalues)
		self.idRng_horizont = index_range_Yfixed(y=-7500,  yvals= self.yvalues)
		#self.idRng_horizont = index_range_yfixed(3, self.Nx)
		# TODO get numbers automatically
		self.Nx = 231 # from 0 to 230
		self.Ny = 19  # from 0 to 18
		self.dt = 0.5 #ka
		self.dx = 5000 #m
	
	def geothermal_heatflux():
		#TODO
		return 0

	def dy(self,y):
		if (y > -30000.0): return 2500.0
		else: 			   return 10000.0
	
	# plot data along a vertical line for a defined time period
	def xlineplot_evolution_uxuy(self, x, tRange):
		idxLst = index_list_Xfixed(x,self.xvalues)
		xi = x
		yi = self.yvalues[idxLst]
		ax = setup_lineplot_uxuy('x', xi)
		for t in tRange:
		    ux_t = self.ux_data[str(t)]
		    uy_t = self.uy_data[str(t)]
		    ux_tyi = ux_t[idxLst]
		    uy_tyi = uy_t[idxLst]
		    ax[0].plot(ux_tyi, yi, linewidth=2, label=t)
		    ax[1].plot(uy_tyi, yi, linewidth=2, label=t)
		for i in [0,1]:
		    ax[i].set_ylabel('$y$ / m')
		    ax[i].legend()
		plt.show()
	
	# plot data along a horizontal line for a defined time period
	def ylineplot_evolution_uxuy(self, i, tRange):
		iR = index_range_yfixed(i,self.Nx)
		xi = self.xvalues[iR[0]:iR[1]+1]
		yi = self.yvalues[iR[0]]
		ax = setup_lineplot_uxuy('y', yi)
		for t in tRange:
		    ux_t = self.ux_data[str(t)]
		    uy_t = self.uy_data[str(t)]
		    ux_tyi = ux_t[iR[0]:iR[1]+1]
		    uy_tyi = uy_t[iR[0]:iR[1]+1]
		    ax[0].plot(xi, ux_tyi, linewidth=2, label=t)
		    ax[1].plot(xi, uy_tyi, linewidth=2, label=t)
		for i in [0,1]:
		    ax[i].set_xlabel('$x$ / m')
		    ax[i].legend()
		plt.show()
	
	# Interpolation for a constant y value (along fixed horizontal line)
	# TODO rewrite using u vector
	@functools.lru_cache
	def interpolateX_data_uxuy(self,x,y,t):
		#i  = 3 # index for y=-7500m
		#iR = index_range_yfixed(i, self.Nx)
		iR = self.idRng_horizont
		xi = self.xvalues[iR[0]:iR[1]+1]
		yi = self.yvalues[iR[0]]
		
		# first interpolation between two time points closest to t
		t1 = closest(self.tvalues, t)[0]
		t2 = t1 + np.sign(t-t1) * self.dt

		ux_t1 = self.ux_data[str(t1)].to_numpy(dtype=float)
		uy_t1 = self.uy_data[str(t1)].to_numpy(dtype=float)
		ux_t1yi = ux_t1[iR[0]:iR[1]+1]
		uy_t1yi = uy_t1[iR[0]:iR[1]+1]
		
		ux_ttyi = ux_t1yi.copy()
		uy_ttyi = uy_t1yi.copy()
		# secant correction if necessary
		if (np.abs(t1-t2) > eps):
			ux_t2 = self.ux_data[str(t2)].to_numpy(dtype=float)
			uy_t2 = self.uy_data[str(t2)].to_numpy(dtype=float)
			ux_t2yi = ux_t2[iR[0]:iR[1]+1]
			uy_t2yi = uy_t2[iR[0]:iR[1]+1]
			ux_ttyi += (ux_t2yi - ux_t1yi) * (t1-t)/(t1-t2)
			uy_ttyi += (uy_t2yi - uy_t1yi) * (t1-t)/(t1-t2)
		
		# second interpolation between two positions closest to x
		tu = closest(xi, x)
		x1 = tu[0]
		x2 = x1 + np.sign(x-x1) * self.dx
		id1 = tu[1]
		id2 = int(id1 + np.sign(x-x1))  # + because x>0!
		
		ux_ttyix1 = ux_ttyi[id1]
		uy_ttyix1 = uy_ttyi[id1]
		
		ux_ttyixx = ux_ttyix1.copy()
		uy_ttyixx = uy_ttyix1.copy()
		# secant correction if necessary
		if (np.abs(x1-x2) > eps):
			ux_ttyix2 = ux_ttyi[id2]
			uy_ttyix2 = uy_ttyi[id2]
			ux_ttyixx += (ux_ttyix2 - ux_ttyix1) * (x1-x)/(x1-x2)
			uy_ttyixx += (uy_ttyix2 - uy_ttyix1) * (x1-x)/(x1-x2)
		
		return ux_ttyixx, uy_ttyixx
	
	# Interpolation for a constant x value (along fixed vertical line)
	# TODO rewrite using u vector
	@functools.lru_cache
	def interpolateY_data_uxuy(self,x,y,t):
		# index list for x=1150000
		idxLst = self.idLst_vertical
		xi = x
		yi = self.yvalues[idxLst]

		# first interpolation between two time points closest to t
		t1 = closest(self.tvalues, t)[0]
		t2 = t1 + np.sign(t-t1) * self.dt

		ux_t1 = self.ux_data[str(t1)].to_numpy(dtype=float)
		uy_t1 = self.uy_data[str(t1)].to_numpy(dtype=float)
		ux_t1xi = ux_t1[idxLst]
		uy_t1xi = uy_t1[idxLst]
		
		ux_ttxi = ux_t1xi.copy()
		uy_ttxi = uy_t1xi.copy()
		# secant correction if necessary
		if (np.abs(t1-t2) > eps):
			ux_t2 = self.ux_data[str(t2)].to_numpy(dtype=float)
			uy_t2 = self.uy_data[str(t2)].to_numpy(dtype=float)
			ux_t2xi = ux_t2[idxLst]
			uy_t2xi = uy_t2[idxLst]
			ux_ttxi += (ux_t2xi - ux_t1xi) * (t1-t)/(t1-t2)
			uy_ttxi += (uy_t2xi - uy_t1xi) * (t1-t)/(t1-t2)

		# second interpolation between two positions closest to y
		tu = closest(yi, y)
		y1 = tu[0]
		y2 = y1 + np.sign(y-y1) * self.dy(y)
		id1 = tu[1]
		id2 = int(id1 - np.sign(y-y1)) # - because y<0!
		
		ux_tty1xi = ux_ttxi[id1]
		uy_tty1xi = uy_ttxi[id1]
		
		ux_ttyyxi = ux_tty1xi.copy()
		uy_ttyyxi = uy_tty1xi.copy()
		# secant correction if necessary
		if (np.abs(y1-y2) > eps):
			ux_tty2xi = ux_ttxi[id2]
			uy_tty2xi = uy_ttxi[id2]
			ux_ttyyxi += (ux_tty2xi - ux_tty1xi) * (y1-y)/(y1-y2)
			uy_ttyyxi += (uy_tty2xi - uy_tty1xi) * (y1-y)/(y1-y2)
		
		return ux_ttyyxi, uy_ttyyxi

	# check interpolation algorithms
	def check_interpolationX_uxuy(self, t, y = -7500):
		ax = setup_lineplot_uxuy('y', y)
		for x in np.linspace(self.xmin, self.xmax, 80):
			u = self.interpolateX_data_uxuy(x,y,t)
			ux = u[0]
			uy = u[1]
			ax[0].scatter(x, ux, linewidth=2)
			ax[1].scatter(x, uy, linewidth=2)
		for i in [0,1]:
		    ax[i].set_xlabel('$x$ / m')
		    
		self.ylineplot_evolution_uxuy(3,[t])

	def check_interpolationY_uxuy(self,t, x = 1150000):
		ax = setup_lineplot_uxuy('x', x)
		for y in np.linspace(self.ymin, self.ymax, 80):
			u = self.interpolateY_data_uxuy(x,y,t)
			ux = u[0]
			uy = u[1]
			ax[0].scatter(ux, y, linewidth=2)
			ax[1].scatter(uy, y, linewidth=2)
		for i in [0,1]:
		    ax[i].set_ylabel('$y$ / m')
		    
		self.xlineplot_evolution_uxuy(x,[t])
		
	def internal_heatsource(self,x,y,t):
		_jahre = 365.25*24*60*60
		t = t + 80*_jahre
		#Aus Jobmann et al., 2017 - Projekt Ansicht:
		#Parameter RK-BE
		BE_Q1 = 842.65
		BE_z1 = 3.58e-11

		BE_Q2 = 1269.66
		BE_z2 = 2.95e-10

		BE_Q3 = 3895.17
		BE_z3 = 9.12e-10

		BE_Q4 = 8308.18
		BE_z4 = 1.04e-8

		BE_Q5 = 42363.74
		BE_z5 = 2.62e-8

		#Parameter RK-HA:
		HA_Q1 = 7583.16
		HA_z1 = 8.28e-9

		HA_Q2 = 2412.91
		HA_z2 = 6.12e-10

		HA_Q3 = 2458.56
		HA_z3 = 7.15e-10

		HA_Q4 = 2546.25
		HA_z4 = 7.15e-10

		HA_Q5 = 7231.62
		HA_z5 = 8.28e-9

		#Gesamtvolumen:
		BE_vol = 7632 #m^3
		HA_vol = 1342.8 #m^3

		#Zeiten:
		t_zwischen_BE = 23*_jahre
		t_zwischen_HA = 30*_jahre
		t_voll = 80*_jahre #Ende der stufenweise Einlagerung --> Start der Waermequelle im Modell

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
				        + 0.05*Q_BE(t-4*_jahre  + t_zwischen_BE)
				        + 0.05*Q_BE(t-8*_jahre  + t_zwischen_BE) 
				        + 0.05*Q_BE(t-12*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-16*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-20*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-24*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-28*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-32*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-36*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-40*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-44*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-48*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-52*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-56*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-60*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-64*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-68*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-72*_jahre + t_zwischen_BE) 
				        + 0.05*Q_BE(t-76*_jahre + t_zwischen_BE) )


		Qsum_HA= lambda t:(0.05*Q_HA(t + t_zwischen_HA)   
				        + 0.05*Q_HA(t-4*_jahre  + t_zwischen_HA)
				        + 0.05*Q_HA(t-8*_jahre  + t_zwischen_HA) 
				        + 0.05*Q_HA(t-12*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-16*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-20*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-24*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-28*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-32*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-36*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-40*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-44*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-48*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-52*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-56*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-60*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-64*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-68*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-72*_jahre + t_zwischen_HA) 
				        + 0.05*Q_HA(t-76*_jahre + t_zwischen_HA) )


		Qsum = lambda t: Qsum_BE(t) + Qsum_HA(t)
		
		return (Qsum_BE(t) + Qsum_HA(t)) / 3000.0

