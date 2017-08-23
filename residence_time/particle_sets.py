import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import time

import PyFVCOM as pf

import residence_time.pylag_reader as pr
import residence_time.grid_box as gb

class particle_set():
	def __init__(self, pylag_data_dir, res, grid_lims_file='tamar_est_lims.dat', grid_file='tamar_v2_grd.dat'):
		self.data_reader = pr.pylag_reader(pylag_data_dir, estuary_lims_file=grid_lims_file)
		self.fvcom_grid = gb.fvcom_grid(grid_lims_file, grid_file)
		self.subsample_grid = gb.regular_grid(res, grid_lims_file=grid_lims_file)

	def get_data_file(self, get_files='all'): 
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

	def add_grid_area(self):
		if not hasattr(self, 'grid_areas'):
			self.grid_areas = {}
			this_area = 1
		else:
			this_area = len(self.grid_areas) + 1

		self.grid_areas[this_area] = grid_area()
		self.grid_areas[this_area].add_reg_grid(self.subsample_grid)	
		
	def residence_times(self, this_area=1):
		for this_time_key, this_time_data in self.particle_data.items():
			in_poly = np.in1d(this_time_data['grid_box'], self.grid_areas[this_area].grid_boxes).reshape(this_time_data['grid_box'].shape)
			 	
			this_time_residences = []
			for this_part_series in in_poly.T:
				this_time_residences.append(parse_occupation_series(this_part_series))
			
			self.particle_data[this_time_key]['residence_area_' + str(this_area)] = np.asarray(this_time_residences)

	def get_all_residence_times(self, residence_area=1):
		return_data = []
		for this_time_key, this_time_data in self.particle_data.items():
			return_data.append(this_time_data['residence_area_' + str(residence_area)])
		return np.asarray(return_data)

	
def parse_occupation_series(bool_series):
	if not np.any(bool_series):
		return 0,0

	else:
		last_present = np.max(np.where(bool_series == True))
		
		first_false = False
		first_leave = np.min(np.where(bool_series == True))
		while first_false == False and first_leave < len(bool_series):
			if bool_series[first_leave] == False:
				first_false = True
			else:
				first_leave += 1

		return last_present, first_leave	


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
		

