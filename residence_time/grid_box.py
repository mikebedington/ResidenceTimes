import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import PyFVCOM as pf


class regular_grid():
	def __init__(self, res, grid_lims_file='tamar_est_lims.dat'):
		"""
		An object to represent a regular grid with squares of size res and overall size definied by grid_lims
		
		Parameters
		----------
		res : float
			resolution of grid (in grid coordinates, uses Euclidean geometry so e.g. metres)
		grid_lims_file : str
			path to a file with a saved array of shape 2x2 and defines the overall corners of the regular grid. 
			In the array [0,0] is x_min, [1,0] is x_max, [0,1] is y_min,and [1,1] is y_max

		"""

		self.grid_resolution = res
		grid_lims = np.loadtxt(grid_lims_file)
		assert grid_lims[0,0] < grid_lims[1,0] and grid_lims[0,1] < grid_lims[1,1]
		self.grid_lims = grid_lims		

		self.no_x_boxes = int(np.ceil((grid_lims[1,0] - grid_lims[0,0])/res))
		self.no_y_boxes = int(np.ceil((grid_lims[1,1] - grid_lims[0,1])/res))

		self.total_boxes = self.no_x_boxes * self.no_y_boxes

		x_coords_raw = np.arange(1,self.no_x_boxes+1)*self.grid_resolution
		self.x_coords = np.tile(x_coords_raw, [self.no_y_boxes,1]).T
		self.x_centre_coords = np.tile(x_coords_raw + self.grid_resolution/2, [self.no_y_boxes,1]).T

		y_coords_raw = np.arange(1,self.no_y_boxes+1)*self.grid_resolution
		self.y_coords = np.tile(y_coords_raw, [self.no_x_boxes,1])
		self.y_centre_coords = np.tile(y_coords_raw + self.grid_resolution/2, [self.no_x_boxes,1])
		
		box_no = np.asarray(np.arange(1, self.total_boxes + 1), dtype=int)
		self.box_no = np.reshape(box_no, [self.no_y_boxes, self.no_x_boxes]).T

	def which_box(self, points_array_x, points_array_y):		
		"""
		Finds the index of which box of the regular grid each point of a given array is in

		Parameters
		----------
		points_array_x : array like
			1d array of x coordinates of points to find their 
		points_array_y : array like
			As above but y coordinates

		Returns
		-------
		which_box : array (len(points_array_x))
			Array of a single integer index for the regular grid box which each point falls in
			The indexing is from the bottom left corner going along the x axis first
			-999 is used for points falling outside the regular grid

		"""

		x_boxes = np.floor(points_array_x/self.grid_resolution)
		y_boxes = np.floor(points_array_y/self.grid_resolution)

		which_box = x_boxes + (y_boxes * self.no_x_boxes)

		which_box[x_boxes < 0] = -999
		which_box[x_boxes >= self.no_x_boxes] = -999

		which_box[y_boxes < 0] = -999
		which_box[y_boxes >= self.no_y_boxes] = -999

		return which_box

	def get_centre_coords(self):
		"""
		A gettr for the centres of each box in the regular grid. Not very pythonic but useful to reduce code for plotting.
		"""
		return self.x_centre_coords, self.y_centre_coords, self.box_no

	def get_box_coord(self, box_nos):
		"""
		Gets the coordinates of a given box index (index as in that returned by self.which_box)
		
		Parameters
		----------
		box_nos : array
			1D array of box indices. Indexes are from the bottom left of the grid along x-axis then up y

		Returns
		-------
		box_coords = array (len(box_nos),2)
			Array of coordinats of the centre of each of the box indexs given in the input. -999 is used for duff indices.

		"""
	
		box_y = np.floor(box_nos/self.no_x_boxes)
		box_x = box_nos - box_y*self.no_x_boxes

		box_coords = np.asarray([box_x * self.grid_resolution, box_y * self.grid_resolution]).T + 0.5*self.grid_resolution

		box_coords[box_nos == -999, :] = -999
	
		return np.squeeze(np.asarray(box_coords))

	def plot_boxes(self, box_nos, ax_handle=None, fig=None, color_str='r'):
		"""
		Plots the specified boxes as pyplot patches

		Parameters
		----------
		box_nos : array
			1D array of box indices as e.g. returned by self.which_box.
		ax_handle, fig : optional, pyplot axes and figure objects
			The figure on which to plot the boxes. If none is given a new figure is instantiated
		color_str : optional, str
			Pyplot colour string for the colours of the boxes

		"""
		
		if ax_handle == None:
			fig, ax_handle = plt.subplots()
			ax_handle.set_xlim([0, np.max(self.x_coords)])
			ax_handle.set_ylim([0, np.max(self.y_coords)])

		box_nos = box_nos[box_nos != -999]	
		box_coords = self.get_box_coord(box_nos)

		for this_coord in box_coords:
			this_rect = patches.Rectangle(this_coord + self.grid_resolution/2,self.grid_resolution,self.grid_resolution,linewidth=1,edgecolor=color_str,facecolor=color_str)
			ax_handle.add_patch(this_rect)

		return ax_handle, fig


class fvcom_grid():
	def __init__(self, grid_lims_file, grid_file):
		"""
		An object representing a reduced fvcom grid, cut down from an orignal grid (grid_file) to just that within the limits
		specified in grid_lims_file
		
		Parameters
		----------
		grid_lims_file : str
			path to a file with a saved array of shape 2x2 and defines the overall corners of the reduced grid.
			In the array [0,0] is x_min, [1,0] is x_max, [0,1] is y_min,and [1,1] is y_max

		grid_file : str
			Path to the original fvcom _grd.dat file

		"""

		self.grid_lims_raw = np.loadtxt(grid_lims_file)
		self.estuary_origin = np.asarray([self.grid_lims_raw[0,0], self.grid_lims_raw[0,1]])

		self.raw_triangles, nodes, self.raw_X, self.raw_Y, depth = pf.grid.read_fvcom_mesh(grid_file)
		
		self.X = self.raw_X - self.estuary_origin[0]
		self.Y = self.raw_Y - self.estuary_origin[1]

		self.grid_lims = np.asarray([[0, self.grid_lims_raw[1,0] - self.estuary_origin[0]], [0, self.grid_lims_raw[1,1] - self.estuary_origin[1]]]).T


