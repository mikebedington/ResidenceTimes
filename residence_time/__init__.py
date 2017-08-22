import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import time

import PyFVCOM as pf

import .pylag_reader as pr
import .grid_box as gb

class particle_set():
	def __init__(self, pylag_data_dir, res, grid_lims_file='tamar_est_lims.dat', grid_file='tamar_v2_grd.dat'):
		self.data_reader = pr.pylag_reader(pylag_data_dir, estuary_lims_file=grid_lims_file)
		self.fvcom_grid = gb.fvcom_grid(grid_lims_file, grid_file)
		self.subsample_grid = gb.regular_grid(res, grid_lims_file=grid_lims_file)

	def get_data_file(self, get_files): 
		self.particle_data = self.data_reader.return_particle_coords(files=get_files)
		self.timesteps = list(self.particle_data)

	def plot_tracks(self, ax_handle=None, fig=None):
		if ax_handle == None:
			fig, ax_handle = plt.subplots()
			ax_handle.set_xlim(self.fvcom_grid.grid_lims[:,0])
			ax_handle.set_ylim(self.fvcom_grid.grid_lims[:,1])

		ax_handle.scatter(self.fvcom_grid.X, self.fvcom_grid.Y, color='grey')

		for this_timestep in self.timesteps:
			ax_handle.plot(self.particle_data[this_timestep]['x'], self.particle_data[this_timestep]['y'])	

		return ax_handle, fig

	def get_particle_box_series(self, get_files=None):
		if get_files is not None:
			self.get_data_file(get_files)
		
		for this_timestep in self.timesteps:
			x = self.particle_data[this_timestep]['x']
			y = self.particle_data[this_timestep]['y']
			self.particle_data[this_timestep]['grid_box'] = np.reshape(self.subsample_grid.which_box(x.flatten(), y.flatten()), x.shape)
	
	def get_first_timestep(self):
		return self.particle_data[self.timesteps[0]]



class grid_area():
	def __init__(self, grid_file='tamar_v2_grd.dat', grid_lims_file='tamar_est_lims.dat'):
		self.no_poly_points = int(input('No of polygon points: '))
	
		triangle, nodes, X, Y, Z = pf.grid_tools.read_fvcom_mesh(grid_file)
		plt.figure()
		plt.scatter(X,Y)
		grid_lims = np.loadtxt(grid_lims_file)
		assert grid_lims[0,0] < grid_lims[1,0] and grid_lims[0,1] < grid_lims[1,1]
		plt.xlim(grid_lims[:,0])
		plt.ylim(grid_lims[:,1])

		self.estuary_origin = np.asarray([grid_lims[0,0], grid_lims[0,1]])

		time.sleep(7)
		self.poly_points_raw = np.asarray(plt.ginput(self.no_poly_points))
		self.poly_points = np.asarray([self.poly_points_raw[:,0] - self.estuary_origin[0], self.poly_points_raw[:,1] - self.estuary_origin[1]]).T
		self.path = mplPath.Path(self.poly_points)

		plt.plot(np.append(self.poly_points_raw[:,0], self.poly_points_raw[0,0]), np.append(self.poly_points_raw[:,1], self.poly_points_raw[0,1]), color='red', linewidth=2)
		
		time.sleep(10)
		plt.close()
	
	def add_reg_grid(self, reg_grid):
		x_coords, y_coords, box_no = reg_grid.get_centre_coords()
		self.grid_boxes = box_no.flatten()[self.path.contains_points(np.asarray([x_coords.flatten(), y_coords.flatten()]).T)]
		
	def first_exit_RT(self, particles):
		pass

	def last_exit_RT(self, particles):
		pass



