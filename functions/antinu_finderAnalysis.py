import numpy as np

import glob
import re
import itertools
import pickle 

from antinu_finder import antinu_finder

def antinu_finderAnalysis(input_file_dir, output_file_dir, file_out_name):

	'''
	This function will receive the numpy arrays of observables and use them to perfomr antinu analysis through the algorithm antinu_finder.
	The function will return two dictionaries: One with the computed energy of prompt and delayed, and the values of the Dt and Dr observables
	'''

	#Anaysis Settings--------------------------------------------------------------------------------------------
	#Primarial observable cuts:
	en_cut_inf = 1
	en_cut_sup = 8.0
	posr_cut_sup = 5500

	#antinu_finder params:
	tau = 251 
	alpha = 9 
	dr_sup_lim = 1000 
	dr_inf_lim = 0 
	energy_delay_inf_cut = 1.8 
	energy_delay_sup_cut = 3.0

	#Useful functions  ------------------------------------------------------------------------------------------
	def posr_cal(x, y, z):
	    dr = np.sqrt(x**2 + y**2 + z**2)
	    return dr

	#ordenar lectura de carpetas:
	def natural_order(file):
	    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', file)]

	#Load data from directories --------------------------------------------------------------------------------

	var_read_name_list = ['energy', 'posx', 'posy', 'posz', 'clockCount50'] #Select observables for the analysis

	flist = [] #File list to be filled with the directories of the observables arrays for the analysis

	for var_i in var_read_name_list:
		flist.append(input_file_dir + var_i + '.npy')

	print(f'Observables List: {flist}')


	# Extract Data Part ------------------------------------------------------------------------------------
	print('construction dictionary of observables')
	data_dict = {var_i: np.load(file_i) for var_i, file_i in zip(var_read_name_list, flist)} #Create the main dictionary where each entry is an observable

	posr = posr_cal(data_dict['posx'], data_dict['posy'], data_dict['posz'])

	#Primary energy and position cuts---
	print('performing primary cuts')
	en_cut_inf = en_cut_inf
	en_cut_sup = en_cut_sup
	posr_cut_sup = posr_cut_sup

	general_condition = (data_dict['energy'] >= en_cut_inf) & (data_dict['energy'] <= en_cut_sup) & (posr <= posr_cut_sup)

	energy_cut = data_dict['energy'][general_condition]
	posr_cut = data_dict['energy'][general_condition]

	posx_cut = data_dict['posx'][general_condition]
	posy_cut = data_dict['posy'][general_condition]
	posz_cut = data_dict['posz'][general_condition]

	clockCount50_cut = data_dict['clockCount50'][general_condition]

	print('antinu finder working ...')
	antinu_dict = antinu_finder(energy_cut, clockCount50_cut, posx_cut, posy_cut, posz_cut, tau, alpha, dr_sup_lim, dr_inf_lim, energy_delay_inf_cut, energy_delay_sup_cut)

    # Save dictionary
	print('Saving dictionary')
	with open(output_file_dir + file_out_name + '.pkl', 'wb') as f:
		pickle.dump(antinu_dict, f)

	print('antinu analysis finished <3')

if __name__ == '__main__':
	input_file_dir = '/lstore/sno/joankl/anti_nu/real_data/gold_data/multiple_jobs_out/analysis/output_0/'
	output_file_dir = '/lstore/sno/joankl/anti_nu/real_data/gold_data/antinu_dict/'
	file_out_name = 'antinu_dict0'
	antinu_finderAnalysis(input_file_dir, output_file_dir, file_out_name)