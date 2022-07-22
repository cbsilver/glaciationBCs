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
			u=x
			v=y
			w=z
			#return coords
		if self.dim == 3: #3D
			u=y
			v=z
			w=x		
		return [u, v, w]

# u, v, w = assign_coordinates(coords, 2)

