import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import time
import pickle as pk

import PyFVCOM as pf

import residence_time.pylag_reader as pr
import residence_time.grid_box as gb

class particle_set():
	def __init__(self, pylag_data_dir, res, grid_lims_file='tamar_est_lims.dat', grid_file='tamar_v2_grd.dat', set_name='particle_set_1'):
		self.data_reader = pr.pylag_reader(pylag_data_dir, estuary_lims_file=grid_lims_file)
		self.grid_lims_file = grid_lims_file
		self.fvcom_grid = gb.fvcom_grid(grid_lims_file, grid_file)
		self.subsample_grid = gb.regular_grid(res, grid_lims_file=grid_lims_file)
		self.set_name = set_name

	def get_data_file(self, get_files='all'): 
		self.particle_data = self.data_reader.return_particle_coords_whole_file(files=get_files)
		self.timesteps = list(self.particle_data)

	def get_data_subset(self, files, subset):
		self.particle_data = self.data_reader.return_particle_coords(files, subset)
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

	def add_grid_area(self, this_area=None):
		if not hasattr(self, 'grid_areas'):
			self.grid_areas = {}
			if this_area is None:
				this_area = 1
		elif this_area is None:
			this_area = len(self.grid_areas) + 1

		self.grid_areas[this_area] = grid_area(grid_lims_file=self.grid_lims_file)
		self.grid_areas[this_area].add_reg_grid(self.subsample_grid)	

	def add_grid_area_by_coord(self, coords, this_area=None):
		if not hasattr(self, 'grid_areas'):
			self.grid_areas = {}
			if this_area is None:
				this_area = 1
		elif this_area is None:
			this_area = len(self.grid_areas) + 1

		self.grid_areas[this_area] = grid_area_coords(coords, grid_lims_file=self.grid_lims_file)
		self.grid_areas[this_area].add_reg_grid(self.subsample_grid)

	def calc_residence_times(self, indices=None, this_area=1):
		if not hasattr(self, 'residence_times'):
			self.residence_times = {}

		for this_time_key, this_time_data in self.particle_data.items():
			in_poly = np.in1d(this_time_data['grid_box'], self.grid_areas[this_area].grid_boxes).reshape(this_time_data['grid_box'].shape)
	
			this_time_residences = []
			for this_part_series in in_poly.T:
				this_time_residences.append(parse_occupation_series(this_part_series))
			
			if indices == None:
				res_array = np.append(np.expand_dims(np.arange(0,len(this_time_residences)), axis=1), np.asarray(this_time_residences), axis=1)
			else:
				res_array = np.append(np.expand_dims(indices, axis=1), np.asarray(this_time_residences), axis=1)
			
			if this_time_key in self.residence_times:
				self.residence_times[this_time_key] = np.append(self.residence_times[this_time_key], res_array, axis=0)
			else:
				self.residence_times[this_time_key] = res_array

	def get_all_residence_times(self, parts_per_loop = 1000):
		if not hasattr(self, 'grid_areas') or 'estuary_area' not in self.grid_areas:
			print('No estuary area defined')
			return
		
		# set up data holder for mean and std first entry and last leave residence times		
		total_parts = self.data_reader.get_particle_nos()
		if total_parts % parts_per_loop:
			print('Particle loop ({}) needs to divide total particles ({})'.format(parts_per_loop, total_parts))
			return
		else:
			part_loops = total_parts/parts_per_loop

		# loop over file and particle subset to update mean and std
		for this_file_ind, this_file in enumerate(self.data_reader.pylag_files_list):
			print('Reading file {}'.format(this_file))
			if total_parts == parts_per_loop:
				print('Retrieving all')
				self.get_data_file(get_files=this_file_ind)
				self.get_particle_box_series()
				self.calc_residence_times(this_area='estuary_area')
				del self.particle_data

			else:
				for this_loop_no in np.arange(0, part_loops):
					print('Loop {} of {}'.format(this_loop_no + 1, part_loops))
					this_indices = np.arange(this_loop_no*parts_per_loop, ((this_loop_no + 1)*parts_per_loop), dtype=int)
					self.get_data_subset(this_file,this_indices)
					self.get_particle_box_series()
					self.calc_residence_times(indices=this_indices, this_area='estuary_area')
					del self.particle_data 



	def save_residence_times(self):
		file_name = 'residence_times_' + self.set_name + '.pk1'
		with open(file_name, 'wb') as f:
			pk.dump(self.residence_times, f, pk.HIGHEST_PROTOCOL)

	def make_grid_index(self, initial_pos_file):
		initial_pos = np.loadtxt(initial_pos_file, skiprows=1)
		self.initial_pos = {'x':initial_pos[:,1] - self.data_reader.estuary_origin[0], 'y':initial_pos[:,2] - self.data_reader.estuary_origin[1], 'z':initial_pos[:,2]}
		self.initial_pos['grid_box'] = self.subsample_grid.which_box(self.initial_pos['x'], self.initial_pos['y'])

	def exclude_parts(self, boolean_exclude):
		for this_key, this_data in self.initial_pos.items():
			self.initial_pos[this_key] = this_data[~boolean_exclude]

	def mean_per_box(self, data_in):
		unique_boxes = np.unique(self.initial_pos['grid_box'])
		box_means = []
		for this_u in unique_boxes:
			this_mean = np.mean(data_in[np.where(self.initial_pos['grid_box'] == int(this_u))])
			box_means.append(this_mean)
		return np.asarray(unique_boxes, dtype=int), np.asarray(box_means)

	def mean_per_grid_area(self, box_data_in):
		grid_area_means = []
		grid_area_box_lists = []
		for this_g_name, this_ga in self.grid_areas.items():
			this_ga_box_list = this_ga.grid_boxes
			records_total = 0
			this_mean = 0
			useful_boxes = this_ga_box_list[np.isin(this_ga_box_list, self.initial_pos['grid_box'])]
			for this_box in useful_boxes:
				this_data = box_data_in[np.where(self.initial_pos['grid_box'] == int(this_box))[0]]
				this_mean = ((this_mean * records_total) + np.sum(this_data)) / (records_total + len(this_data))
				records_total = records_total + len(this_data)

			this_mean = np.tile(this_mean, len(useful_boxes))
			grid_area_means.append(this_mean)
			grid_area_box_lists.append(useful_boxes)			

		grid_area_means = [item for sublist in grid_area_means for item in sublist]
		grid_area_box_lists = [item for sublist in grid_area_box_lists for item in sublist]

		return np.asarray(grid_area_box_lists, dtype=int), np.asarray(grid_area_means)

	def box_plotting_data(self, boxes, data):
		pc_x = self.subsample_grid.x_coords - self.subsample_grid.grid_resolution
		pc_y = self.subsample_grid.y_coords - self.subsample_grid.grid_resolution
		pc_data = np.empty(pc_x.size)
		pc_data[:] = np.NAN
		pc_data[boxes] = data
		pc_data = np.reshape(pc_data, (pc_x.shape[1], pc_x.shape[0])).T
		pc_data = np.ma.masked_where(np.isnan(pc_data), pc_data)
		return pc_x, pc_y, pc_data
		


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
	
		triangle, nodes, X, Y, Z = pf.grid.read_fvcom_mesh(grid_file)
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

	def points_in_area(self, points_array, reorigin=False):
		if reorigin:
			points_array = np.asarray([points_origin[:,0] - self.estuary_origin[0], points_origin[:,1] - self.estuary_origin[1]])

		return np.asarray(self.path.contains_points(points_list))
	

class grid_area_coords(grid_area):
	def __init__(self, coords, grid_file='tamar_v2_grd.dat', grid_lims_file='tamar_est_lims.dat'):
		grid_lims = np.loadtxt(grid_lims_file)
		assert grid_lims[0,0] < grid_lims[1,0] and grid_lims[0,1] < grid_lims[1,1]

		self.estuary_origin = np.asarray([grid_lims[0,0], grid_lims[0,1]])
		self.poly_points = np.asarray([coords[:,0] - self.estuary_origin[0], coords[:,1] - self.estuary_origin[1]]).T
		self.path = mplPath.Path(self.poly_points)
	
