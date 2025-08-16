'''

Main function to extract the variables of interest from the ntuples.root

'''
import numpy as np
import uproot
import pandas as pd
import glob

def pos_r(x, y, z):
	#function to recompute the radius of events due to the posz correction given the x, y and z position
	r = np.sqrt(x**2 + y**2 + z**2)
	return r

route = '/snoplus simulations/real_data_antinu/'
pattern_file = 'Analysis20R_r*'

save_route = 'data/'

#DataFrame Main Structure to be fill
df = pd.DataFrame(columns = ['evIndex', 'evID', 'energy', 'posx', 'posy', 'posz', 'posr', 'time_clock (mcs)'])

for route_i in glob.glob(route+pattern_file):

	print('reading file ' + route_i)
	
	#read ith file
	file_i = uproot.open(route_i)
	output = file_i['output;1']

	#extract variables of interest------------------------------
	#validation info
	scint_fit = np.array(output['scintFit'])

	#event Index
	evIndex = np.array(output['evIndex'])  #evIndex = 0 is prompt, and evIndex > 0 is delayed and tails of prompt and delayed

	#recons. energy
	energy = np.array(output['energy'])

	#recons. position
	posx = np.array(output['posx'])
	posy = np.array(output['posy'])
	posz = np.array(output['posz']) - 184.4
	posr = pos_r(posx, posy, posz)

	#MC info
	mcID = np.array(output['mcIndex'])
	mc_posx = np.array(output['mcPosx'])
	mc_posy = np.array(output['mcPosy'])
	mc_posz = np.array(output['mcPosz'])
	mc_posr = np.array(output['mcPosr'])

	# Time info
	clock_count50 = np.array(output['clockCount50'], dtype = np.int64)

	#ev info
	evID = np.array(output['eventID'])
	nhits = np.array(output['nhits'])

	# Extract valid info -> Valid scint_fit -------------

	evIndex = np.extract(scint_fit, evIndex)

	mc_posx = np.extract(scint_fit, mc_posx)
	mc_posy = np.extract(scint_fit, mc_posy)
	mc_posz = np.extract(scint_fit, mc_posz)

	posx = np.extract(scint_fit, posx)
	posy = np.extract(scint_fit, posy)
	posz = np.extract(scint_fit, posz)
	posr = np.extract(scint_fit, posr)

	time = np.extract(scint_fit, (clock_count50*20)/1000)  #convert ns to μs

	evID = np.extract(scint_fit, evID)

	energy = np.extract(scint_fit, energy)
	#---------------------------------------------------------------
	# Arange data into a ith pandas dataframe and then add to the main file df

	data_i = {'evIndex': evIndex,
				'evID': evID,
				'energy': energy,
				'posx': posx,
				'posy': posy,
				'posz': posz,
				'posr': posr,
				'time_clock (mcs)': time
				}

	df_i = pd.DataFrame(data_i)

	df = pd.concat([df, df_i], ignore_index=True)

df.to_csv(save_route + 'full_data.csv')