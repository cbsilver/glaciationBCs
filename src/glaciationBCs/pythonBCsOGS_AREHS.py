# Collection of python boundary condition (BC) classes for OpenGeoSys
# BCs reflect the external geosphere: cryo-, litho- and atmosphere
# Physical units: depending on parameter set, see below!

import sys
print(sys.version)

import matplotlib
from glaciationBCs import glacierclass_AREHS as glc	# glacier
from glaciationBCs import crustclass_AREHS as crc 	# earth crust
from glaciationBCs import repoclass_AREHS as dgr	# repository
from glaciationBCs import airclass_AREHS as air		# atmosphere

import numpy as np


import OpenGeoSys
from glaciationBCs.constants_AREHS import *

# Nomenclature: BC Process_LocationQuantity_Component
# 					(THM)			(XYZ)

# Process	Dirichlet BC	Neumann BC (normal to boundary)
# T			temperature		heat flux
# H			pressure		hydraulic flux
# M			displacement	momentum flux (stress vector)

#TODO linear distribution for T gradient and p gradient at left boundary

# ---------------------------------------------------------
# Thermal BCs
# ---------------------------------------------------------
class BCT_SurfaceTemperature(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCT_SurfaceTemperature, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(T_ini, T_min, t_)
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_)
		if plotinput:
			self.air.plot_evolution()

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		print(self.air.tcr.stage_control(t))
		"""
		if x-self.glacier.x_0 > self.glacier.length(t) or self.glacier.length(t)==0.0:
			value = self.air.temperature(t)
		else:
			value = self.glacier.temperature(x,t)
		"""
		under_glacier = x-self.glacier.x_0 <= self.glacier.length(t) > 0
		if under_glacier:
			# prescribe fixed temperature underneath the glacier body
			value = self.glacier.temperature(x,t)
		else:
			#linear profile from north to south
			value = self.air.temperature(t)
		
		return (True, value)

class BCT_InitialTemperature(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCT_InitialTemperature, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(T_ini, T_min, t_)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		value = self.air.T_ini
		
		return (True, value)

class BCT_SourceFromRepository(OpenGeoSys.SourceTerm):

	def __init__(self):
		super(BCT_SourceFromRepository, self).__init__()
		# instantiate member objects of the external geosphere
		self.repo = dgr.repo(BE_Q, BE_z, BE_f, HA_Q, HA_z, HA_f, BE_vol, HA_vol, 
							 lrepo, t_inter_BE, t_inter_HA, t_filled)
		if plotinput:
			self.repo.print_max_load()
			self.repo.plot_evolution()

	def getFlux(self, t, coords, primary_vars):
		x, y, z = coords
		
		value = 0
		if ((xrmin <= x <= xrmax) and (yrmin <= y <= yrmax)): #TODO
			print("y = ",y)
			value = self.repo.radioactive_heatflux(t)
		
		derivative = [0.0]
		return (value, derivative)	

class BCT_BottomHeatFlux(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCT_BottomHeatFlux, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(q_geo)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of heat
		x, y, z = coords
		
		# get heat flux component
		value = self.crust.geothermal_heatflux()[1]
		derivative = [0.0, 0.0]
		return (True, value, derivative)	

class BCT_VerticalGradient(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCT_VerticalGradient, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(T_ini, T_min, t_)
		self.crust = crc.crust(q_geo)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords

		# evolving surface temperature
		T_top = self.air.temperature(t)
		# vertical geothermal gradient
		value = self.crust.geothermal_temperature(x,y,z,t,T_top)

		return (True, value)

# ------------------------------------------------------
# Hydraulic BCs
# ------------------------------------------------------
class BCH_SurfacePressure(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCH_SurfacePressure, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(T_ini, T_min, t_)
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords

		print(self.glacier.tcr_h.stage_control(t))
		
		under_glacier = x-self.glacier.x_0 <= self.glacier.length(t) > 0
		if under_glacier:
			# height dependent pressure from glacier
			value = self.glacier.pressure(x,t)
		else:
			# fixed pressure from ambient air
			value = self.air.pressure
		
		return (True, value)

class BCH_SurfaceInflux(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCH_SurfaceInflux, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_)
	
	def getFlux(self, t, coords, primary_vars): #here Neumann BC: hydraulic flux
		x, y, z = coords
		
		under_glacier = x-self.glacier.x_0 <= self.glacier.length(t) > 0
		if under_glacier:
			# get hydraulic flux under glacier
			value = self.glacier.local_meltwater(x,t)
			derivative = [ 0.0, 0.0 ]
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, [ 0.0, 0.0 ])

class BCH_VerticalGradient(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCH_VerticalGradient, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(T_ini, T_min, t_)
		self.crust = crc.crust(q_geo)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords

		# height dependent pressure from the crust
		p_pore = self.crust.hydrostatic_pressure(x,y,z,t)
		# fixed pressure from ambient air
		p_atmo = self.air.pressure

		value = p_pore + p_atmo

		return (True, value)


# --------------------------------------------------------
# Mechanics BCs
# --------------------------------------------------------
class BCM_SurfaceTraction_X(OpenGeoSys.BoundaryCondition):
	
	def __init__(self):
		super(BCM_SurfaceTraction_X, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_)
		if plotinput:
			self.glacier.print_max_load()
			self.glacier.plot_evolution()
			self.glacier.plot_evolving_shape()
		
	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		x, y, z = coords
		
		under_glacier = x-self.glacier.x_0 <= self.glacier.length(t) > 0
		if under_glacier:
			value = self.glacier.tangentialstress(x,t)
			derivative = [ 0.0, 0.0 ]
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, [ 0.0, 0.0 ])

class BCM_SurfaceTraction_Y(OpenGeoSys.BoundaryCondition):
	
	def __init__(self):
		super(BCM_SurfaceTraction_Y, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		x, y, z = coords
		
		under_glacier = x-self.glacier.x_0 <= self.glacier.length(t) > 0
		if under_glacier:
			value = self.glacier.normalstress(x,t)
			derivative = [ 0.0, 0.0,   ]
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, [ 0.0, 0.0,   ])

class BCM_BottomDisplacement_X(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCM_BottomDisplacement_X, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(q_geo)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_x
		value = self.crust.displacement_below(x,y,z,t)[0]
		
		return (True, value)

class BCM_BottomDisplacement_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCM_BottomDisplacement_Y, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(q_geo)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_y
		value = self.crust.displacement_below(x,y,z,t)[1]
		
		return (True, value)

class BCM_LateralDisplacement_X(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCM_LateralDisplacement_X, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(q_geo)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_x
		value = self.crust.displacement_aside(x,y,z,t)[0]
		
		return (True, value)

class BCM_LateralDisplacement_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCM_LateralDisplacement_Y, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(q_geo)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_y
		value = self.crust.displacement_aside(x,y,z,t)[0]
		
		return (True, value)


# ---------------------------------------------
# instantiate the BC objects used by OpenGeoSys
# ---------------------------------------------
# Naming convention:
# bc_Process_(external)origin_boundary_type(_coefficient)

# Atmosphere BCs
bc_T_atmosph_above_Dirichlet = BCT_InitialTemperature()

# Cryosphere BCs
bc_T_glacier_above_Dirichlet = BCT_SurfaceTemperature()
bc_H_glacier_above_Dirichlet = BCH_SurfacePressure()
bc_H_glacier_above_Neumann   = BCH_SurfaceInflux()
bc_M_glacier_above_Neumann_x = BCM_SurfaceTraction_X()
bc_M_glacier_above_Neumann_y = BCM_SurfaceTraction_Y()

# Internal source
bc_T_dgrepo_inside_VolSource = BCT_SourceFromRepository()

# Lithosphere BCs
bc_T_crustal_below_Neumann_y = BCT_BottomHeatFlux()
bc_T_crustal_aside_Dirichlet = BCT_VerticalGradient()
bc_H_crustal_aside_Dirichlet = BCH_VerticalGradient()
bc_M_crustal_south_Dirichlet_x = BCM_LateralDisplacement_X()
bc_M_crustal_south_Dirichlet_y = BCM_LateralDisplacement_Y()
bc_M_crustal_below_Dirichlet_x = BCM_BottomDisplacement_X()
bc_M_crustal_below_Dirichlet_y = BCM_BottomDisplacement_Y()
#bc_M_crustal_north
#bc_M_crustal_aside

