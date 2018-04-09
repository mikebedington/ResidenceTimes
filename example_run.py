import multiprocessing as mp
import numpy as np

import residence_time as rt

data_dir = '/data/euryale2/scratch/mbe/Models_2/PyLAG/output/test_dir'

filestr = 'part_set_1'
parts_set = rt.particle_set(data_dir, 200, grid_lims_file='tamar_est_lims_2.dat', set_name=filestr)
coords = np.loadtxt('estuary_area.txt')
parts_set.add_grid_area_by_coord(coords, this_area='estuary_area')
parts_set.get_all_residence_times()
parts_set.save_residence_times()


with open('part_set_1.pk1', 'rb') as f:
    res_times2 = pk.load(f)

parts_set.make_grid_index('test_out.dat')

res_time_ex = res_times2[1451692800]

bool_exclude = res_time_ex[:,1] == 0
parts_set.exclude_parts(bool_exclude)
res_time_ex = res_time_ex[~bool_exclude,1]
boxes, data = parts_set.mean_per_box(res_time_ex)
pc_x, pc_y, pc_data = parts_set.box_plotting_data(boxes, data)

plt.pcolormesh(pc_x, pc_y, pc_data)


