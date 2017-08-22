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

	def return_particle_coords(self, files='all'):
		if files == 'all':
			files_to_retrieve = self.pylag_files_list
		elif isinstance(files, int):
			files_to_retrieve = [self.pylag_files_list[files]]
		else:
			files_to_retrieve = []
			for this_file_ind in files:
				files_to_retrieve.append(self.pylag_files_list[this_file_ind])
		
		particle_dict = {}		

		for this_file in files_to_retrieve:
			this_nc = nc.Dataset(self.data_dir + this_file, 'r')
			
			x_adj = this_nc.variables['xpos'][:] - self.estuary_origin[0]
			y_adj = this_nc.variables['ypos'][:] - self.estuary_origin[1]			
			this_t = this_nc.variables['time'][:] - self.time_origin

			this_part_dict = {'x':x_adj, 'y':y_adj, 'z':this_nc.variables['zpos'][:], 't':this_t}
			particle_dict[np.min(this_t)] = this_part_dict	
			this_nc.close()		

		return particle_dict		
