# Collection of python boundary condition (BC) classes for OpenGeoSys
# BCs reflect the external geosphere: cryo-, litho- and atmosphere
# Physical units: depending on parameter set, see below!

from glaciationBCs.coord_control_AREHS import uvw_bounds, L_max
from glaciationBCs import coord_control_AREHS as uvw# coordinates
from glaciationBCs import glacierclass_AREHS as glc	# glacier
from glaciationBCs import crustclass_AREHS as crc 	# earth crust
from glaciationBCs import repoclass_AREHS as dgr	# repository
from glaciationBCs import airclass_AREHS as air		# atmosphere
import glaciationBCs.constants_AREHS as ac		    # AREHS constants

import pyvista as pv
import numpy as np
import OpenGeoSys

# Nomenclature: BC Process_LocationQuantity_Component
# 					(THM)					(XYZ) TODO uvw

# Process	Dirichlet BC	Neumann BC (normal to boundary)
# T			temperature		heat flux
# H			pressure		hydraulic flux
# M			displacement	momentum flux (stress vector)

# ---------------------------------------------------------
# Thermal BCs
# ---------------------------------------------------------
class BCT_InitialTemperature(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCT_InitialTemperature, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		value = self.air.T_ini

		return (True, value)

class BCT_SurfaceTemperature(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCT_SurfaceTemperature, self).__init__()
		self.uvw = uvw.coord_control(dim)
		u_max = uvw_bounds(dim, bounds, 0)[1]
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(dim, bounds), ac.H_max, u_max, ac.t_)
		if ac.plotinput:
			self.air.plot_evolution()

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		u = self.uvw.assign_coordinates(coords)[0]

		if t != self.air.t_prev:
			print(self.air.tcr.stage_control(t))
			self.air.t_prev = t

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			# prescribe fixed temperature underneath the glacier body
			if ac.T_smoothtransit:
				shape = self.glacier.shape(u,t)
				Tg = self.glacier.temperature(u,t)
				Ta = self.air.temperature(t)
				value = shape * (Tg - Ta) + Ta
			else:
				value = self.glacier.temperature(u,t)
		else:
			# prescribe fixed temperature from the atmoshpere
			value = self.air.temperature(t)

		return (True, value)

class BCT_SourceFromRepository(OpenGeoSys.SourceTerm):

	def __init__(self, dim, repo_path):
		super(BCT_SourceFromRepository, self).__init__()
		self.uvw = uvw.coord_control(dim)
		repo = pv.XMLUnstructuredGridReader(repo_path).read()
		if dim == 2:
			# one division for representative waste amount per line
			# second division for norming per line length
			repo_size = np.sum(repo.compute_cell_sizes().cell_data["Length"])**2
		if dim == 3:
			repo_size = np.sum(repo.compute_cell_sizes().cell_data["Area"])
		# instantiate member objects of the external geosphere
		self.repo = dgr.repo(ac.BE_Q, ac.BE_z, ac.BE_f, ac.HA_Q, ac.HA_z, ac.HA_f, ac.BE_vol, ac.HA_vol,
							 repo_size, ac.t_inter_BE, ac.t_inter_HA, ac.t_filled)
		if ac.plotinput:
			self.repo.print_max_load()
			self.repo.plot_evolution()

	def getFlux(self, t, coords, primary_vars):		
		# prescribe heat flux from radioactive repository
		value = self.repo.radioactive_heatflux(t)
		derivative = [0.0] * len(primary_vars)
		return (value, derivative)

class BCT_BottomHeatFlux(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCT_BottomHeatFlux, self).__init__()
		self.uvw = uvw.coord_control(dim)
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of heat

		# get vertical heat flux component
		value = self.crust.geothermal_heatflux()[1]

		derivative = [0.0] * len(primary_vars)
		return (True, value, derivative)

class BCT_LateralHeatFlux(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCT_LateralHeatFlux, self).__init__()
		self.uvw = uvw.coord_control(dim)
		u_max = uvw_bounds(dim, bounds, 0)[1]
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(dim, bounds), ac.H_max, u_max, ac.t_)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of heat
		u, v, w = self.uvw.assign_coordinates(coords)

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			T_top = self.glacier.temperature(u,t)
		else:
			T_top = self.air.temperature(t)

		# get heat flux component
		q_mag = self.crust.lateral_heatflux(v, T_top)

		u_max = uvw_bounds(self.dim, self.bound, 0)[1]
		at_northern_boundary = (u-ac.eps < u_max < u+ac.eps)
		if at_northern_boundary:
			value = q_mag			# > 0 : OGS heat influx
		else:#southern boundary
			value =-q_mag			# < 0 : OGS heat outflux

		derivative = [0.0] * len(primary_vars)
		return (True, value, derivative)

class BCT_VerticalGradient(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCT_VerticalGradient, self).__init__()
		self.uvw = uvw.coord_control(dim)
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		v = self.uvw.assign_coordinates(coords)[1]

		# evolving surface temperature
		T_top = self.air.temperature(t)
		# vertical geothermal gradient
		value = self.crust.geothermal_temperature(v, T_top)

		return (True, value)

# ------------------------------------------------------
# Hydraulic BCs
# ------------------------------------------------------
class BCH_InitialPressure(OpenGeoSys.BoundaryCondition):

	def __init__(self):
		super(BCH_InitialPressure, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		value = self.air.pressure

		return (True, value)

# TODO: find a better solution like a mixed boundary condition
class BCH_SurfacePressure(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCH_SurfacePressure, self).__init__()
		self.uvw = uvw.coord_control(dim)
		u_max = uvw_bounds(dim, bounds, 0)[1]
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(dim, bounds), ac.H_max, u_max, ac.t_)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		u = self.uvw.assign_coordinates(coords)[0]

		if t != self.glacier.t_prev:
			print(self.glacier.tcr_h.stage_control(t))
			self.glacier.t_prev = t

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		value = self.air.pressure
		if under_glacier:
			# height dependent pressure from glacier
			value += self.glacier.pressure(u,t)
		apply_bc = True if primary_vars[0] > 273.15 else False
		return (apply_bc, value)

class BCH_SurfaceInflux(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCH_SurfaceInflux, self).__init__()
		self.uvw = uvw.coord_control(dim)
		u_max = uvw_bounds(dim, bounds, 0)[1]
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_max(dim, bounds), ac.H_max, u_max, ac.t_)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: hydraulic flux
		u = self.uvw.assign_coordinates(coords)[0]

		derivative = [0.0] * len(primary_vars)

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			# get hydraulic flux under glacier
			value = self.glacier.local_meltwater(u,t)
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, derivative)

class BCH_VerticalGradient(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCH_VerticalGradient, self).__init__()
		self.uvw = uvw.coord_control(dim)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		v = self.uvw.assign_coordinates(coords)[1]

		# height dependent pressure from the crust
		p_pore = self.crust.hydrostatic_pressure(v)
		# fixed pressure from ambient air
		p_atmo = self.air.pressure

		value = p_pore + p_atmo

		return (True, value)

class BCH_VerticalGradientFromInit(OpenGeoSys.BoundaryCondition):

	def __init__(self, boundary_result_path):
		super(BCH_VerticalGradientFromInit, self).__init__()
		
		pvr = pv.get_reader(boundary_result_path)
		mesh = pvr.read()
		if "pressure_interpolated" in mesh.point_data.keys():
			self.init_data = mesh.point_data["pressure_interpolated"]
		elif "pressure" in mesh.point_data.keys():
			self.init_data = mesh.point_data["pressure"]
		else:
			raise Exception("No pressure data in boundary result!")

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		# p_init = p_pore + p_air
		#index = np.argmin(np.linalg.norm(self.init_coords - coords, axis=1))
		#print(index, node_id)
		p_init = self.init_data[node_id]
		return (True, p_init)

class BCH_VerticalGradientExtended(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCH_VerticalGradientExtended, self).__init__()
		self.uvw = uvw.coord_control(dim)
		u_max = uvw_bounds(dim, bounds, 0)[1]
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(dim, bounds), ac.H_max, u_max, ac.t_)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		u, v, w = self.uvw.assign_coordinates(coords)

		# height dependent pressure from the crust
		p_pore = self.crust.hydrostatic_pressure(v)
		# fixed pressure from ambient air/glacier
		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		p_extr = self.air.pressure
		if under_glacier:
			# height dependent pressure from glacier
			p_extr += self.glacier.pressure(u,t)
			
		value = p_pore + p_extr

		return (True, value)


class BCH_VerticalGradientExtendedFromInit(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds, boundary_result_path):
		super(BCH_VerticalGradientExtendedFromInit, self).__init__()
		self.uvw = uvw.coord_control(dim)
		u_max = uvw_bounds(dim, bounds, 0)[1]
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_max(dim, bounds), ac.H_max, u_max, ac.t_)

		pvr = pv.get_reader(boundary_result_path)
		mesh = pvr.read()
		if "pressure_interpolated" in mesh.point_data.keys():
			self.init_data = mesh.point_data["pressure_interpolated"]
		elif "pressure" in mesh.point_data.keys():
			self.init_data = mesh.point_data["pressure"]
		else:
			raise Exception("No pressure data in boundary result!")


	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		u = self.uvw.assign_coordinates(coords)[0]

		# p_init = p_pore + p_air
		#index = np.argmin(np.linalg.norm(self.init_coords - coords, axis=1))
		p_init = self.init_data[node_id]
		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		p_total = p_init
		if under_glacier:
			p_total += self.glacier.pressure(u,t)

		return (True, p_total)

# --------------------------------------------------------
# Mechanics BCs
# --------------------------------------------------------
class BCM_SurfaceTraction_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCM_SurfaceTraction_X, self).__init__()
		self.uvw = uvw.coord_control(dim)
		u_max = uvw_bounds(dim, bounds, 0)[1]
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_max(dim, bounds), ac.H_max, u_max, ac.t_)
		if ac.plotinput:
			self.glacier.print_max_load()
			self.glacier.plot_evolution()
			self.glacier.plot_evolving_shape()

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		u = self.uvw.assign_coordinates(coords)[0]

		derivative = [0.0] * len(primary_vars)

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			value = self.glacier.tangentialstress(u,t)
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, derivative)

class BCM_SurfaceTraction_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCM_SurfaceTraction_Y, self).__init__()
		self.uvw = uvw.coord_control(dim)
		u_max = uvw_bounds(dim, bounds, 0)[1]
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(dim, bounds), ac.H_max, u_max, ac.t_)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		u = self.uvw.assign_coordinates(coords)[0]

		derivative = [0.0] * len(primary_vars)

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		value = -self.air.pressure
		if under_glacier:
			value += self.glacier.normalstress(u,t)

		return (True, value, derivative)

class BCM_BottomDisplacement_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCM_BottomDisplacement_X, self).__init__()
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		# prescribe displacement u_x
		value = self.crust.displacement_below()[0]

		return (True, value)

class BCM_BottomDisplacement_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCM_BottomDisplacement_Y, self).__init__()
		# instantiate member objects of the external geosphere
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		# prescribe displacement u_y
		value = self.crust.displacement_below()[1]

		return (True, value)

class BCM_LateralDisplacement_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCM_LateralDisplacement_X, self).__init__()
		# instantiate member objects of the external geosphere
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		# prescribe displacement u_x
		value = self.crust.displacement_aside()[0]

		return (True, value)

class BCM_LateralDisplacement_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCM_LateralDisplacement_Y, self).__init__()
		# instantiate member objects of the external geosphere
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		# prescribe displacement u_y
		value = self.crust.displacement_aside()[1]

		return (True, value)

class BCM_LateralTraction_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, dim, bounds):
		super(BCM_LateralTraction_X, self).__init__()
		self.uvw = uvw.coord_control(dim)
		# instantiate member objects of the external geosphere
		v_min, v_max = uvw_bounds(dim, bounds, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		# v = self.uvw.assign_coordinates(coords)[1]
		derivative = [0.0] * len(primary_vars)
		# value = self.crust.lithostatic_stresses(v)
		value = 0
		return (True, value, derivative)

