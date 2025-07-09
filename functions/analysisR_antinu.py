# -*- coding: utf-8 -*

'''

Main function to extract the variables of interest from the ntuples.root
This Code will save the summary files of the input data files.It will return a pandas data frame in csv format
and save a file.csv per each input file

The R in this function means that it will read the Analysis20R* files
'''

import numpy as np
import uproot
import pandas as pd
import glob
import re

#Funcion para garantizar una lectura ordenada de runs!
def orden_natural(archivo):
	return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(\d+)', archivo)]


route = '/share/neutrino/snoplus/Data/FullFill_2p2/Analysis20R/'
pattern_file = 'Analysis20R*.root'

file_list = glob.glob(route+pattern_file)
file_list.sort(key = orden_natural)

save_route = '/lstore/sno/joankl/anti_nu/real_data/gold_data/all_files/'
out_file_name = 'Analysis20R_nu_real_data_N_'

#DataFrame Main Structure to be fill
df = pd.DataFrame(columns = ['fitValid', 'scintFit', 'runID', 'dcFlag', 'evIndex',
			     'evID', 'triggerWord', 'Nhits', 'NhitsClean', 'correctedNhits', 
			     'neckNhits','energy', 'posx', 'posy', 'posz', 'posr', 'posFOM', 
			      'posFOM2', 'time_clock50 (mcs)', 'ITR', 'beta14', 'alphabeta212', 'alphabeta214'])

for i_dx, route_i in enumerate(file_list):

	print('reading file ' + route_i)
	
	#read ith file
	file_i = uproot.open(route_i)
	#print(file_i.keys())
	output = file_i['output;1']
	#print(output.keys())


	#extract variables of interest------------------------------
        #runID
	run_id = np.array(output['runID'])

	#validation of scintfit
	scint_fit = np.array(output['scintFit'])
	#fitValid
	fit_valid = np.array(output['fitValid'])

	#event Index
	evIndex = np.array(output['evIndex'])  #evIndex = 0 is prompt, and evIndex > 0 is delayed and tails of prompt and delayed
	
	#TriggerWords
	trigger_word = np.array(output['triggerWord'])

	#Hit Info:
	nhits = np.array(output['nhits'])
	nhits_clean = np.array(output['nhitsCleaned'])
	corrected_nhits = np.array(output['correctedNhits'])
	neck_nhits = np.array(output['necknhits'])

	#FOM
	posFOM = np.array(output['posFOM'])
	posFOM2 = np.array(output['posFOM2'])

	#ITR
	itr = np.array(output['itr'])

	#classifiers (?)
	beta14 = np.array(output['beta14'])
	alphabeta212 = np.array(output['alphaBeta212'])
	alphabeta214 = np.array(output['alphaBeta214'])

	#recons. energy
	energy = np.array(output['energy'])

	#recons. position
	posx = np.array(output['posx'])
	posy = np.array(output['posy'])
	posz = np.array(output['posz'])
	posr = np.array(output['posr'])

	#MC info
	#mcID = np.array(output['mcIndex'])
	#mc_posx = np.array(output['mcPosx'])
	#mc_posy = np.array(output['mcPosy'])
	#mc_posz = np.array(output['mcPosz'])
	#mc_posr = np.array(output['mcPosr'])

	# Time info
	clock_count50 = np.array(output['clockCount50'], dtype = np.int64)

	#evID
	evID = np.array(output['eventID'])

	#dcFlag
	dc_flag = np.array(output['dcFlagged'])

	# Extract valid info -> Valid scint_fit and fit_valid -------------

	valid_condition = (scint_fit & fit_valid)

	scint_fit =  np.extract(valid_condition, scint_fit)
	fit_valid = np.extract(valid_condition, fit_valid)

	evIndex = np.extract(valid_condition, evIndex)

	#mc_posx = np.extract(valid_condition, mc_posx)
	#mc_posy = np.extract(valid_condition, mc_posy)
	#mc_posz = np.extract(valid_condition, mc_posz)

	posx = np.extract(valid_condition, posx)
	posy = np.extract(valid_condition, posy)
	posz = np.extract(valid_condition, posz)
	posr = np.extract(valid_condition, posr)

	time = np.extract(valid_condition, (clock_count50*20)/1000)  #convert ns to μs
        
	run_id = np.extract(valid_condition, run_id)
	evID = np.extract(valid_condition, evID)

	energy = np.extract(valid_condition, energy)

	dc_flag = np.extract(valid_condition, dc_flag)
	
	trigger_word = np.extract(valid_condition, trigger_word)
	
	nhits = np.extract(valid_condition, nhits)
	nhits_clean = np.extract(valid_condition, nhits_clean)
	corrected_nhits = np.extract(valid_condition, corrected_nhits)
	neck_nhits = np.extract(valid_condition, neck_nhits)
	
	posFOM = np.extract(valid_condition, posFOM)
	posFOM2 = np.extract(valid_condition, posFOM2)

	itr = np.extract(valid_condition, itr)

	beta14 = np.extract(valid_condition, beta14)
	alphabeta212 = np.extract(valid_condition, alphabeta212)
	alphabeta214 = np.extract(valid_condition, alphabeta214)
	
	#---------------------------------------------------------------
	# Arange data into a ith pandas dataframe and then add to the main file df

	data_i = {'fitValid':fit_valid,
		  'scintFit':scint_fit,
		  'runID': run_id,
		  'dcFlag': dc_flag,
		  'evIndex': evIndex,
		  'evID': evID,
                  'triggerWord': trigger_word,
                  'Nhits': nhits,
                  'NhitsClean': nhits_clean,
                  'correctedNhits': corrected_nhits,
		  'neckNhits': neck_nhits,
		  'energy': energy,
		  'posx': posx,
		  'posy': posy,
		  'posz': posz,
		  'posr': posr,
		  'posFOM': posFOM,
		  'posFOM2': posFOM2,
		  'time_clock50 (mcs)': time,
		  'ITR': itr,
		  'beta14': beta14,
		  'alphabeta212': alphabeta212,
                  'alphabeta214': alphabeta214}
	print(data_i)
	df_i = pd.DataFrame(data_i)

	#df = pd.concat([df, df_i], ignore_index=True)

	df_i.to_parquet(save_route + out_file_name + str(i_dx) +'.parquet')
