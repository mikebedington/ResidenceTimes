import numpy as np
import os
import netCDF4 as nc

from natsort import natsorted

class pylag_reader:
	def __init__ (self, pylag_data_dir, estuary_lims_file='tamar_est_lims.dat', time_origin=0):
		"""
		Object to read a directory of pylag output and adjust it to a new frame of reference

		Parameters
		----------
		pylag_data_dir : str
			Path with directory 

		estuary_lims_file : str
			path to a file with a saved array of shape 2x2 and defines the overall corners of the estuary subset. 
			In the array [0,0] is x_min, [1,0] is x_max, [0,1] is y_min,and [1,1] is y_max.
			The point [x_min, y_min] will be used as the origin for the estuary reference frame

		time_origin : optional,int
			Can be used to re origin the time stamps 


		"""
		if pylag_data_dir[-1] != '/':
			pylag_data_dir+='/'
		self.data_dir = pylag_data_dir
		self.time_origin = time_origin
		estuary_lims = np.loadtxt(estuary_lims_file)
		assert estuary_lims[0,0] < estuary_lims[1,0] and estuary_lims[0,1] < estuary_lims[1,1] 
		self.estuary_origin = np.asarray([estuary_lims[0,0], estuary_lims[0,1]])
		pylag_files_list_raw = os.listdir(pylag_data_dir)
		self.pylag_files_list = natsorted([this_str for this_str in pylag_files_list_raw if '.nc' in this_str])

	def return_particle_coords_whole_file(self, files='all'):
		"""
		Gets the particle tracks from files in the objects data directory

		Parameters
		----------
		files : optional, str or int
			Which files to read, by default reads all in the data directory ('all'). Otherwise it can return
			the data from a single file identified either by name or by index within the objects file list.

		Returns
		-------
		particle_dict : dict of dicts
			Outer dict has keys corresponding to  the start times of the first particle from each read file, the inner 
			dictionary contains 'x' - x coordinates, 'y' - y coordinates, and 't', times. All returned values are in the
			estuary reference frame

		"""

		if files == 'all':
			files_to_retrieve = self.pylag_files_list
		elif isinstance(files, int):
			files_to_retrieve = [self.pylag_files_list[files]]
		elif isinstance(files, str):
			files_to_retrieve = [files]
		else:
			files_to_retrieve = []
			for this_file_ind in files:
				files_to_retrieve.append(self.pylag_files_list[this_file_ind])
		
		particle_dict = {}		
		tot_files = len(files_to_retrieve)

		for this_ind, this_file in enumerate(files_to_retrieve):
			print('Retrieving file %d of %d' % (this_ind +1, tot_files))
			this_part_dict = self.return_particle_coords(this_file, 'all')
			particle_dict.update(this_part_dict)

		return particle_dict

	def return_particle_coords(self, file_str, part_indices='all'):
		"""
		Get the particle x,y,t coordinates in the estuary reference frame for a single file in the data directory

		Note - this should probably either be private or use a file index rather than a string
		
		Parameters
		----------
		file_str : str
			The name of the pylag output file in the objects data directory to be read
		part_indices : optional, boolean array or integer array or 'all'
			Return only a subset of the particles as defined by an indexing array.	

		Returns
		-------
		part_dict : dict of dicts
			Outer dict has a key of the start time of the first particle in the given file, the inner dictionary
			contains 'x' - x coordinates, 'y' - y coordinates, and 't', times. All returned values are in the estuary
			reference frame

		"""
		this_nc = nc.Dataset(self.data_dir + file_str, 'r')
		if isinstance(part_indices, str) and part_indices == 'all':
			#part_indices = np.arange(0, this_nc.variables['xpos'].shape[1])
			x_adj = this_nc.variables['xpos'][:] - self.estuary_origin[0]
			y_adj = this_nc.variables['ypos'][:] - self.estuary_origin[1]
			z_adj = this_nc.variables['zpos'][:]

		else:
			x_adj = this_nc.variables['xpos'][:, part_indices] - self.estuary_origin[0]
			y_adj = this_nc.variables['ypos'][:, part_indices] - self.estuary_origin[1]
			z_adj = this_nc.variables['zpos'][:, part_indices]
		
		this_t = this_nc.variables['time'][:] - self.time_origin

		part_dict_temp = {'x':x_adj, 'y':y_adj, 'z':z_adj, 't':this_t}
		this_nc.close()
		part_dict = {np.min(part_dict_temp['t']):part_dict_temp}
		return part_dict

	def get_particle_nos(self):
		"""
		Get total number of particles in first pylag file
	
		Returns
		-------
		particle_nos : int
			Number of particles from first pylag file
		"""

		this_nc = nc.Dataset(self.data_dir + self.pylag_files_list[0])
		particle_nos = this_nc.variables['xpos'].shape[1]
		this_nc.close()
		return particle_nos
