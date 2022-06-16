# Collection of python boundary condition (BC) classes for OpenGeoSys
# BCs reflect the external geosphere: cryo-, litho- and atmosphere
# Physical units: depending on parameter set, see below!

import OpenGeoSys
import glaciationBCs
from glaciationBCs import glacierclass as glc	#glacial objects
from glaciationBCs import crustclass as crc 	#crustal objects
from glaciationBCs import airclass as air		# aerial objects

import numpy as np

s_a = 365.25*24*3600 #=31557600 seconds per year

# Choose parametrization
T_N = 266.15 #K
T_S = 276.15 #K
T_C = 8 #K
Set = "TH"

if (Set=="M"): # units: kg, m, s, K
	L_dom = 120000 #m
	L_max = 0.7*L_dom
	H_max = 200 #m
	x_0 = -0.5*L_dom
	t_0 = 0.00 #s
	t_1 = 1.0000 #s
	t_2 = 2
	t_3 = 3
	t_4 = 4

if (Set=="T"): # units: kg, m, a, K
	L_dom = 1150000 #m
	L_max = 575000 #m
	H_max = 3200 #m
	x_0 = 0.0 #m
	t_0 = 17500 #a
	t_1 = t_0 + 12500 #a
	t_2 = t_1 +  5000 #a
	t_3 = t_2 +  5000 #a
	t_4 = t_3 + 10000 #a

if (Set=="HM"): # units: kg, m, s, K
	L_dom = 1150000 #m
	L_max = 575000 #m
	H_max = 3200 #m
	x_0 = 0.0 #m
	t_0 = 0.0 * 32500 * s_a #a
	t_1 = t_0 + 12500 * s_a #a
	t_2 = t_1 +  5000 * s_a #a
	t_3 = t_2 +  5000 * s_a #a
	t_4 = t_3 + 10000 * s_a #a 

if (Set=="TH"): # units: kg, m, a, K
	L_dom = 1150000 #m
	L_max = 575000 #m
	H_max = 3200 #m
	x_0 = 0.0 #m
	t_0 = 20000000 #a
	t_1 = t_0 + 12500 #a
	t_2 = t_1 +  5000 #a
	t_3 = t_2 +  5000 #a
	t_4 = t_3 + 10000 #a

# Required: Choose vertical scaling factor
y_sfactor = 20 #TODO!

# Optional: Choose path to external data
path2data = '~/Forschung/Simulations/OpenGeoSys/SedimentaryBasinBense/HM/dataGIA/'
plotinput = False

# Nomenclature: BC Process_LocationQuantity_Component
# 					(THM)			(XYZ)

# Process	Dirichlet BC	Neumann BC (normal to boundary)
# T			temperature		heat flux
# H			pressure		hydraulic flux
# M			displacement	momentum flux (stress vector)


# Thermal BCs
# -----------
class BCT_SurfaceTemperature(OpenGeoSys.BoundaryCondition):

	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCT_SurfaceTemperature, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(L_dom, T_N, T_S, T_C, t_0, t_1, t_2, t_3, t_4)
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		print(self.glacier.stagecontrol(t))
		
		if x-self.glacier.x_0 > self.glacier.length(t) or self.glacier.length(t)==0.0:
			#linear profile from north to south
			value = self.air.temperature_profile(x,t)
		else:
			# prescribe fixed temperature underneath the glacier body
			value = self.glacier.temperature(x,t)
		
		return (True, value)

# Hydraulic BCs
# -------------
class BCH_SurfacePressure(OpenGeoSys.BoundaryCondition):

	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCH_SurfacePressure, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(L_dom, T_N, T_S, T_C, t_0, t_1, t_2, t_3, t_4)
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords

		print(self.glacier.stagecontrol(t))
		
		if x-self.glacier.x_0 <= self.glacier.length(t):
			# height dependent pressure from glacier
			value = self.glacier.pressure(x,t)
		else:
			# fixed pressure from ambient air
			value = self.air.pressure
		
		return (True, value)

class BCH_SurfaceHydrohead(OpenGeoSys.BoundaryCondition):

	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCH_SurfaceHydrohead, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(L_dom, T_N, T_S, T_C, t_0, t_1, t_2, t_3, t_4)
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords

		print(self.glacier.stagecontrol(t))
		# get vertical displacement
		u_y = self.glacier.local_deflection_heuristic(x,t)
		
		# head from surface topology
		h_top = y/20 + u_y # scaled!
		
		if x-self.glacier.x_0 <= self.glacier.length(t):
			# height dependent hydraulic head from glacier
			h_ice = self.glacier.hydrohead(x,t)
			value = h_ice + h_top
		else:
			# fixed head from ambient air
			h_air = self.air.hydrohead
			value = h_air + h_top
		
		return (True, value)

class BCH_SurfaceInflux(OpenGeoSys.BoundaryCondition):

	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCH_SurfaceInflux, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
	
	def getFlux(self, t, coords, primary_vars): #here Neumann BC: hydraulic flux
		x, y, z = coords
		
		if x-self.glacier.x_0 <= self.glacier.length(t):
			# get hydraulic flux under glacier
			value = self.glacier.local_meltwater(x,t)
			derivative = [ 0.0, 0.0 ]
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, [ 0.0, 0.0 ])

class BCH_SourceFromDeflection(OpenGeoSys.SourceTerm):

	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCH_SourceFromDeflection, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
		if plotinput: self.glacier.plot_deflection()

	def getFlux(self, t, coords, primary_vars):
		x, y, z = coords
		
		# get subsidence velocity acting as a hydraulic head source term
		value = self.glacier.local_deflection_rate_heuristic(x,t)
		Jac = [0.0, 0.0]
		return (value, Jac)


# Mechanics BCs
# -------------
class BCM_SurfaceTraction_X(OpenGeoSys.BoundaryCondition):
	
	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCM_SurfaceTraction_X, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
		if plotinput: self.glacier.print_max_load()
		if plotinput: self.glacier.plot_evolution()
		
	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		x, y, z = coords
		
		if x-self.glacier.x_0 <= self.glacier.length(t):
			value = self.glacier.tangentialstress(x,t)
			derivative = [ 0.0, 0.0 ]
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, [ 0.0, 0.0 ])

class BCM_SurfaceTraction_Y(OpenGeoSys.BoundaryCondition):
	
	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCM_SurfaceTraction_Y, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		x, y, z = coords
		
		if x-self.glacier.x_0 <= self.glacier.length(t):
			value = self.glacier.normalstress(x,t)
			derivative = [ 0.0, 0.0,   ]
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, [ 0.0, 0.0,   ])

class BCM_BottomDeflection(OpenGeoSys.BoundaryCondition):

	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCM_BottomDeflection, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
		if plotinput: self.glacier.plot_deflection()

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_y
		# scale here with 20 !TODO!
		value = 20 * self.glacier.local_deflection(x,t)
		
		return (True, value)

class BCM_DomainDisplacement(OpenGeoSys.BoundaryCondition):

	def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4):
		super(BCM_DomainDisplacement, self).__init__()
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
		if plotinput: self.glacier.plot_deflection()

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_y
		# scale here with 20 !TODO!
		value = 20 * self.glacier.local_displacement_heuristic(x,y,t)
		
		return (True, value)

class BCM_BottomDisplacement_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, path2data):
		super(BCM_BottomDisplacement_X, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(path2data)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_x
		t_in_ka = t/s_a/1000
		# ?TODO? y_scale = y/20
		value = self.crust.interpolateX_data_uxuy(x,y,t_in_ka)[0]
		
		return (True, value)

class BCM_BottomDisplacement_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self, path2data):
		super(BCM_BottomDisplacement_Y, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(path2data)
		if plotinput:
			idx = 3
			tRange = np.linspace(t_0/s_a/1000, t_1/s_a/1000, 26)
			self.crust.ylineplot_evolution_uxuy(idx, tRange)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_y
		t_in_ka = t/s_a/1000
		# ?TODO? y_scale = y/20
		value = self.crust.interpolateX_data_uxuy(x,y,t_in_ka)[1]
		
		return (True, value)

class BCM_LateralDisplacement_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, path2data):
		super(BCM_LateralDisplacement_X, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(path2data)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_x
		t_in_ka = t/s_a/1000
		y_scale = y/20
		value = self.crust.interpolateY_data_uxuy(x,y_scale,t_in_ka)[0]
		
		return (True, value)

class BCM_LateralDisplacement_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self, path2data):
		super(BCM_LateralDisplacement_Y, self).__init__()
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(path2data)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		x, y, z = coords
		
		# prescribe displacement u_y
		t_in_ka = t/s_a/1000
		y_scale = y/20
		value = self.crust.interpolateY_data_uxuy(x,y_scale,t_in_ka)[1]
		
		return (True, value)


# instantiate the BC objects used by OpenGeoSys
# ---------------------------------------------
# Naming convention:
# bc_Process_(external)origin_boundary_type(_coefficient)

# Cryosphere BCs
bc_T_glacier_above_Dirichlet = BCT_SurfaceTemperature(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_H_glacier_above_Dirichlet = BCH_SurfacePressure(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_H_glacier_above_Dirichlet_head = BCH_SurfaceHydrohead(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_H_glacier_above_VolSource_head = BCH_SourceFromDeflection(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_H_glacier_above_Neumann   = BCH_SurfaceInflux(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_M_glacier_above_Neumann_x = BCM_SurfaceTraction_X(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_M_glacier_above_Neumann_y = BCM_SurfaceTraction_Y(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
#bc_H_glacier_north_Neumann

# Lithosphere BCs
bc_M_crustal_south_Dirichlet_x = BCM_LateralDisplacement_X(path2data)
bc_M_crustal_south_Dirichlet_y = BCM_LateralDisplacement_Y(path2data)
bc_M_crustal_below_Dirichlet_x = BCM_BottomDisplacement_X(path2data)
bc_M_crustal_below_Dirichlet_y = BCM_BottomDisplacement_Y(path2data)
#bc_M_crustal_below_Dirichlet_y = BCM_BottomDeflection(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
#bc_M_glacier_above_Dirichlet_y = BCM_BottomDeflection(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_M_glacier_above_Dirichlet_y = BCM_DomainDisplacement(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
#bc_M_crustal_north
#bc_M_crustal_aside

# just for downward compatibility
bc_thermally_dirichlet = BCT_SurfaceTemperature(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_hydraulic_dirichlet = BCH_SurfacePressure(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_hydraulic_neumann   = BCH_SurfaceInflux(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_mechanics_neumann_x = BCM_SurfaceTraction_X(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_mechanics_neumann_y = BCM_SurfaceTraction_Y(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)
bc_mechanics_dirichlet = BCM_BottomDeflection(L_dom, L_max, H_max, x_0, t_0, t_1, t_2, t_3, t_4)

bc_y = BCM_SurfaceTraction_Y(L_dom, L_max, H_max, x_0, t_0, t_1, 0, 0, 0)
