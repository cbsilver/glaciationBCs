# Auxiliary class for coordinate handling and shifting 2D/3D
# Physical units: kg, m, s, K

import numpy as np

# CONVENTION: directions instead of coordinates:
# * u: glacier advance direction	| u points to the South
# * v: vertical direction (gravity) | v points to the Sky
# * w: lateral direction (only 3D)	| w points to the West

#TODO 2D -> 3D: coords = swap(coords) # y->x, z->y

class coord_control():
	# class variables:
	# -

	def __init__(self, dim):
		self.dim = dim
		# or "2D" / "3D"

	def assign_coordinates(self, coords):
		x, y, z = coords
		if self.dim == 2: #2D
			u, v, w = (x, y, z)
			#return coords
		if self.dim == 3: #3D
			u, v, w = (y, z, x)
		return [u, v, w]

# u, v, w = assign_coordinates(coords, 2)
