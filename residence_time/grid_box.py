import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import PyFVCOM as pf


class regular_grid():
	def __init__(self, res, grid_lims_file='tamar_est_lims.dat'):
		self.grid_resolution = res
		grid_lims = np.loadtxt(grid_lims_file)
		assert grid_lims[0,0] < grid_lims[1,0] and grid_lims[0,1] < grid_lims[1,1]
		self.grid_lims = grid_lims		

		self.no_x_boxes = np.ceil((grid_lims[1,0] - grid_lims[0,0])/res)
		self.no_y_boxes = np.ceil((grid_lims[1,1] - grid_lims[0,1])/res)

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
		
		x_boxes = np.floor(points_array_x/self.grid_resolution)
		y_boxes = np.floor(points_array_y/self.grid_resolution)

		which_box = x_boxes + (y_boxes * self.no_x_boxes)

		which_box[x_boxes < 0] = -999
		which_box[x_boxes >= self.no_x_boxes] = -999

		which_box[y_boxes < 0] = -999
		which_box[y_boxes >= self.no_y_boxes] = -999

		return which_box

	def get_centre_coords(self):

		return self.x_centre_coords, self.y_centre_coords, self.box_no

	def get_box_coord(self, box_nos):
		
		box_y = np.floor(box_nos/self.no_x_boxes)
		box_x = box_nos - box_y*self.no_x_boxes

		box_coords = np.asarray([box_x * self.grid_resolution, box_y * self.grid_resolution]).T + 0.5*self.grid_resolution

		box_coords[box_nos == -999, :] = -999
	
		return np.squeeze(np.asarray(box_coords))

	def plot_boxes(self, box_nos, ax_handle=None, fig=None, color_str='r'):
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
		self.grid_lims_raw = np.loadtxt(grid_lims_file)
		self.estuary_origin = np.asarray([self.grid_lims_raw[0,0], self.grid_lims_raw[0,1]])

		self.raw_triangles, nodes, self.raw_X, self.raw_Y, depth = pf.grid_tools.read_fvcom_mesh(grid_file)
		
		self.X = self.raw_X - self.estuary_origin[0]
		self.Y = self.raw_Y - self.estuary_origin[1]

		self.grid_lims = np.asarray([[0, self.grid_lims_raw[1,0] - self.estuary_origin[0]], [0, self.grid_lims_raw[1,1] - self.estuary_origin[1]]]).T





