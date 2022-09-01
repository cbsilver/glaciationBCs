# Auxiliary class for coordinate handling and shifting 2D/3D
# Physical units: kg, m, s, K

import numpy as np

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

