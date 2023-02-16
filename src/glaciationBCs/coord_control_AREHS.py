# Auxiliary class for coordinate handling and shifting 2D/3D
# Physical units: kg, m, s, K

from glaciationBCs.constants_AREHS import glacial_advance

# CONVENTION: directions instead of coordinates:
# * u: glacier retreat direction	| u points to the North
# * v: vertical direction (gravity) | v points to the Sky
# * w: lateral direction (only 3D)	| w points to the East
# OLD VERSION:
# * u: glacier advance direction	| u points to the South
# * v: vertical direction (gravity) | v points to the Sky
# * w: lateral direction (only 3D)	| w points to the West

class coord_control():
	# class variables:
	# -

	# constructor
	def __init__(self, dimension):
		# instance variables
		if dimension == 2:
			self.is2D = True
		else:
			self.is2D = False

	def assign_coordinates(self, coords):
		x, y, z = coords
		if self.is2D:
			#return coords
			u, v, w = (x, y, z)
		else: #3D
			u, v, w = (y, z, x)
		return [u, v, w]

# Geomodel-specific parameters
def uvw_bounds(dim, bounds, axis):
	coord_ctrl = coord_control(dim)
	x_min, x_max, y_min, y_max, z_min, z_max = bounds
	u_min, v_min, w_min = coord_ctrl.assign_coordinates((x_min, y_min, z_min))
	u_max, v_max, w_max = coord_ctrl.assign_coordinates((x_max, y_max, z_max))
	if axis==0:
		return u_min, u_max 
	if axis==1:		
		return v_min, v_max
	if axis==2:
		return w_min, w_max


def L_max(dim, bounds):
	u_min, u_max = uvw_bounds(dim, bounds, 0)
	return glacial_advance * (u_max - u_min)