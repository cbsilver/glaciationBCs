# Collection of python boundary condition (BC) classes for OpenGeoSys
# BCs reflect the external geosphere: cryo-, litho- and atmosphere
# Physical units: depending on parameter set, see below!

from glaciationBCs import coord_control_AREHS as uvw# coordinates
from glaciationBCs import glacierclass_AREHS as glc	# glacier
from glaciationBCs import crustclass_AREHS as crc 	# earth crust
from glaciationBCs import repoclass_AREHS as dgr	# repository
from glaciationBCs import airclass_AREHS as air		# atmosphere
import glaciationBCs.constants_AREHS as ac		    # AREHS constants

import OpenGeoSys

# Nomenclature: BC Process_LocationQuantity_Component
# 					(THM)					(XYZ) TODO uvw

# Process	Dirichlet BC	Neumann BC (normal to boundary)
# T			temperature		heat flux
# H			pressure		hydraulic flux
# M			displacement	momentum flux (stress vector)


# Geomodel-specific parameters
def model_uvw(props, axis):
	coord_ctrl = uvw.coord_control(props.dimension)
	u_min, v_min, w_min = coord_ctrl.assign_coordinates(props.coords_min)
	u_max, v_max, w_max = coord_ctrl.assign_coordinates(props.coords_max)
	if axis==0:
		return u_min, u_max 
	if axis==1:		
		return v_min, v_max
	if axis==2:
		return w_min, w_max


def L_max(props):
	u_min, u_max = model_uvw(props, 0)
	return ac.glacial_advance * (u_max - u_min)

# ---------------------------------------------------------
# Thermal BCs
# ---------------------------------------------------------
class BCT_InitialTemperature(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCT_InitialTemperature, self).__init__()
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.dimension = props.dimension

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		value = self.air.T_ini

		return (True, value)

class BCT_SurfaceTemperature(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCT_SurfaceTemperature, self).__init__()
		self.dimension = props.dimension
		self.uvw = uvw.coord_control(props.dimension)
		u_min, u_max = model_uvw(props, 0)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(props), ac.H_max, u_max, ac.t_)
		if ac.plotinput:
			self.air.plot_evolution()

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		u, v, w = self.uvw.assign_coordinates(coords)

		if t != self.air.t_prev:
			print(self.air.tcr.stage_control(t))
			self.air.t_prev = t

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			# prescribe fixed temperature underneath the glacier body
			value = self.glacier.temperature(u,t)
		else:
			value = self.air.temperature(t)

		return (True, value)

class BCT_SourceFromRepository(OpenGeoSys.SourceTerm):

	def __init__(self, props, repo_size):
		super(BCT_SourceFromRepository, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		# instantiate member objects of the external geosphere
		self.repo = dgr.repo(ac.BE_Q, ac.BE_z, ac.BE_f, ac.HA_Q, ac.HA_z, ac.HA_f, ac.BE_vol, ac.HA_vol,
							 repo_size, ac.t_inter_BE, ac.t_inter_HA, ac.t_filled, props.dimension)
		if ac.plotinput:
			self.repo.print_max_load()
			self.repo.plot_evolution()

	def getFlux(self, t, coords, primary_vars):
		u, v, w = self.uvw.assign_coordinates(coords)
		
		# prescribe heat flux from radioactive repository
		value = self.repo.radioactive_heatflux(t)			

		derivative = [0.0] * len(primary_vars)
		return (value, derivative)

class BCT_BottomHeatFlux(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCT_BottomHeatFlux, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		v_min, v_max = model_uvw(props, 1)
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of heat

		# get vertical heat flux component
		value = self.crust.geothermal_heatflux()[1]

		derivative = [0.0] * len(primary_vars)
		return (True, value, derivative)

class BCT_LateralHeatFlux(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCT_LateralHeatFlux, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		self.props = props
		u_min, u_max = model_uvw(props, 0)
		v_min, v_max = model_uvw(props, 1)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(props), ac.H_max, u_max, ac.t_)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of heat
		u, v, w = self.uvw.assign_coordinates(coords)

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			T_top = self.glacier.temperature(u,t)
		else:
			T_top = self.air.temperature(t)

		# get heat flux component
		q_mag = self.crust.lateral_heatflux(v, T_top, self.props)

		u_min, u_max = model_uvw(self.props, 0)
		at_northern_boundary = (u-ac.eps < u_max < u+ac.eps)
		if at_northern_boundary:
			value = q_mag			# > 0 : OGS heat influx
		else:#southern boundary
			value =-q_mag			# < 0 : OGS heat outflux

		derivative = [0.0] * len(primary_vars)
		return (True, value, derivative)

class BCT_VerticalGradient(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCT_VerticalGradient, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		v_min, v_max = model_uvw(props, 1)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		u, v, w = self.uvw.assign_coordinates(coords)

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


class BCH_SurfacePressure(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCH_SurfacePressure, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		u_min, u_max = model_uvw(props, 0)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(props), ac.H_max, u_max, ac.t_)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		u, v, w = self.uvw.assign_coordinates(coords)

		if t != self.glacier.t_prev:
			print(self.glacier.tcr_h.stage_control(t))
			self.glacier.t_prev = t

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			# height dependent pressure from glacier
			value = self.glacier.pressure(u,t)
		else:
			# fixed pressure from ambient air
			value = self.air.pressure

		return (True, value)

class BCH_SurfaceInflux(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCH_SurfaceInflux, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		u_min, u_max = model_uvw(props, 0)
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_max(props), ac.H_max, u_max, ac.t_)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: hydraulic flux
		u, v, w = self.uvw.assign_coordinates(coords)

		derivative = [0.0] * len(primary_vars)

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			# get hydraulic flux under glacier
			value = self.glacier.local_meltwater(u,t)
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, derivative)

class BCH_VerticalGradient(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCH_VerticalGradient, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		u_min, u_max = model_uvw(props, 0)
		v_min, v_max = model_uvw(props, 1)
		# instantiate member objects of the external geosphere
		self.air = air.air(ac.T_ini, ac.T_min, ac.t_)
		self.glacier = glc.glacier(L_max(props), ac.H_max, u_max, ac.t_)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):
		u, v, w = self.uvw.assign_coordinates(coords)

		# height dependent pressure from the crust
		p_pore = self.crust.hydrostatic_pressure(v)
		# fixed pressure from ambient air/glacier
		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			# height dependent pressure from glacier
			p_extr = self.glacier.pressure(u,t)
		else:
			# fixed pressure from ambient air
			p_extr = self.air.pressure

		value = p_pore + p_extr

		return (True, value)


# --------------------------------------------------------
# Mechanics BCs
# --------------------------------------------------------
class BCM_SurfaceTraction_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCM_SurfaceTraction_X, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		u_min, u_max = model_uvw(props, 0)
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_max(props), ac.H_max, u_max, ac.t_)
		if ac.plotinput:
			self.glacier.print_max_load()
			self.glacier.plot_evolution()
			self.glacier.plot_evolving_shape()

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		u, v, w = self.uvw.assign_coordinates(coords)

		derivative = [0.0] * len(primary_vars)

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			value = self.glacier.tangentialstress(u,t)
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, derivative)

class BCM_SurfaceTraction_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCM_SurfaceTraction_Y, self).__init__()
		self.uvw = uvw.coord_control(props.dimension)
		u_min, u_max = model_uvw(props, 0)
		# instantiate member objects of the external geosphere
		self.glacier = glc.glacier(L_max(props), ac.H_max, u_max, ac.t_)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		u, v, w = self.uvw.assign_coordinates(coords)

		derivative = [0.0] * len(primary_vars)

		under_glacier = self.glacier.u_0-u <= self.glacier.length(t) > 0
		if under_glacier:
			value = self.glacier.normalstress(u,t)
			return (True, value, derivative)
		# no BC => free boundary then (no flux)
		return (False, 0.0, derivative)

class BCM_BottomDisplacement_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCM_BottomDisplacement_X, self).__init__()
		v_min, v_max = model_uvw(props, 1)
		# instantiate member objects of the external geosphere
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		# prescribe displacement u_x
		value = self.crust.displacement_below()[0]

		return (True, value)

class BCM_BottomDisplacement_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCM_BottomDisplacement_Y, self).__init__()
		# instantiate member objects of the external geosphere
		v_min, v_max = model_uvw(props, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		# prescribe displacement u_y
		value = self.crust.displacement_below()[1]

		return (True, value)

class BCM_LateralDisplacement_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCM_LateralDisplacement_X, self).__init__()
		# instantiate member objects of the external geosphere
		v_min, v_max = model_uvw(props, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		# prescribe displacement u_x
		value = self.crust.displacement_aside()[0]

		return (True, value)

class BCM_LateralDisplacement_Y(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCM_LateralDisplacement_Y, self).__init__()
		# instantiate member objects of the external geosphere
		v_min, v_max = model_uvw(props, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getDirichletBCValue(self, t, coords, node_id, primary_vars):

		# prescribe displacement u_y
		value = self.crust.displacement_aside()[1]

		return (True, value)

class BCM_LateralTraction_X(OpenGeoSys.BoundaryCondition):

	def __init__(self, props):
		super(BCM_LateralTraction_X, self).__init__()
		self.props = props
		self.uvw = uvw.coord_control(props.dimension)
		# instantiate member objects of the external geosphere
		v_min, v_max = model_uvw(props, 1)
		self.crust = crc.crust(ac.q_geo, v_min, v_max, ac.T_ini, ac.T_bot)

	def getFlux(self, t, coords, primary_vars): #here Neumann BC: flux of linear momentum
		u, v, w = self.uvw.assign_coordinates(coords)
		derivative = [0.0] * len(primary_vars)
		value = self.crust.lithostatic_stresses(v, self.props)
		return (True, value, derivative)

