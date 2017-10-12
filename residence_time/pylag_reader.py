import numpy as np
import os
import netCDF4 as nc

from natsort import natsorted


class pylag_reader:
	def __init__ (self, pylag_data_dir, estuary_lims_file='tamar_est_lims.dat', time_origin=0):
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
		if files == 'all':
			files_to_retrieve = self.pylag_files_list
		elif isinstance(files, int):
			files_to_retrieve = [self.pylag_files_list[files]]
		else:
			files_to_retrieve = []
			for this_file_ind in files:
				files_to_retrieve.append(self.pylag_files_list[this_file_ind])
		
		particle_dict = {}		
		tot_files = len(files_to_retrieve)

		for this_ind, this_file in enumerate(files_to_retrieve):
			print('Retrieving file %d of %d' % (this_ind +1, tot_files))
			this_part_dict = self.return_particle_coords(this_file, 'all')
			particle_dict[np.min(this_part_dict['t'])] = this_part_dict	
			this_nc.close()		

		return particle_dict

	def return_particle_coords(self, file_str, part_indices):
		this_nc = nc.Dataset(self.data_dir + file_str, 'r')
		if isinstance(part_indices, str) and part_indices == 'all':
			part_indices = np.arange(0, this_nc.variables['xpos'].shape[1])

		x_adj = this_nc.variables['xpos'][:, part_indices] - self.estuary_origin[0]
		y_adj = this_nc.variables['ypos'][:, part_indices] - self.estuary_origin[1]
		z_adj = this_nc.variables['zpos'][:, part_indices]
		this_t = this_nc.variables['time'][:] - self.time_origin

		part_dict_temp = {'x':x_adj, 'y':y_adj, 'z':z_adj, 't':this_t}
		this_nc.close()
		part_dict = {np.min(part_dict_temp['t']):part_dict_temp}
		return part_dict

	def get_particle_nos(self):
		this_nc = nc.Dataset(self.data_dir + self.pylag_files_list[0])
		particle_nos = this_nc.variables['xpos'].shape[1]
		this_nc.close()
		return particle_nos
