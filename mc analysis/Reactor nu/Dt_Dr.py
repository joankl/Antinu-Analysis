'''

Function that receives the full_data and returns the Dt and Dr of simulated events usign the evIndex

'''
import numpy as np
import uproot
import pandas as pd

read_file_route = 'data/full_data.csv'
save_file_route = 'data/'

#Read the data file and extract info
data_df = pd.read_csv(read_file_route)

evID = np.array(data_df['evID'])
evIndex = np.array(data_df['evIndex'])
energy = np.array(data_df['energy'])
posx = np.array(data_df['posx'])
posy = np.array(data_df['posy'])
posz = np.array(data_df['posz'])
posr = np.array(data_df['posr'])
time_clock = np.array(data_df['time_clock (mcs)'])


def Dt_Dr(en_cut, r_cut):

	ev_index_condition = (evIndex == 0) | (evIndex == 2)

	#Index of elements in data which verify the ev index condition
	filter_index = np.where(ev_index_condition)[0]

	test_evindex = evIndex[filter_index]  #To test if the selection if done correctly

	#list of indicies of elements to be removed in filter_index. This is used to correct the filter_index list.
	index_to_del = []

	for i in range(len(test_evindex)-1):
		if test_evindex[i] == 0 and test_evindex[i+1] == 0:
			index_to_del.append(i)
		if test_evindex[i] == 2 and test_evindex[i+1] == 2:
			index_to_del.append(i+1)

	#cleaning of the filter_index list
	filter_index = np.delete(filter_index, index_to_del) 

	# Extract observables and Arange by pairs of prompt and delayed. If there is an error due to non-equal division, then set a condition for the filter_index:
	evindex_filter = evIndex[filter_index].reshape((-1,2))
	energy_filter = energy[filter_index].reshape((-1,2))
	posr_filter = posr[filter_index].reshape((-1,2))
	posx_filter = posx[filter_index].reshape((-1,2))
	posy_filter = posy[filter_index].reshape((-1,2))
	posz_filter = posz[filter_index].reshape((-1,2))
	time_filter = time_clock[filter_index].reshape((-1,2))

	# Now we are in conditions to apply cuts on energy and position of the events
	cut_condition = (energy_filter >= en_cut) & (posr_filter <= r_cut)

	mask = np.all(cut_condition, axis = 1)  #evaluates where all the elements of axis 1 (rows) are True

	#cut values of energy and postions
	energy_cut = energy_filter[mask]
	posr_cut = posr_filter[mask]
	posx_cut = posx_filter[mask]
	posy_cut = posy_filter[mask]
	posz_cut = posz_filter[mask]
	time_cut = time_filter[mask]

	#extract the information of prompt and delayed
	#Delta time
	t0 = time_cut[:,0]
	t1 = time_cut[:,1]
	Dt = t1 - t0

	#Delta r
	posx_0 = posx_cut[:,0]
	posy_0 = posy_cut[:,0]
	posz_0 = posz_cut[:,0]

	posx_1 = posx_cut[:,1]
	posy_1 = posy_cut[:,1]
	posz_1 = posz_cut[:,1]

	dx = posx_1 - posx_0
	dy = posy_1 - posy_0
	dz = posz_1 - posz_0

	Dr = np.sqrt((dx**2) + (dy**2) + (dz**2))

	#Now write and save the data
	data = {'Dt (mcs)': Dt,
			'Dr (mm)': Dr
			}

	Dr_Dt_pd_data = pd.DataFrame(data)
	Dr_Dt_pd_data.to_csv(save_file_route + 'Dr_Dt_{en_cut_}MeV_{r_cut_}mm.csv'.format(en_cut_ = en_cut, r_cut_ = r_cut))

	return Dr_Dt_pd_data


