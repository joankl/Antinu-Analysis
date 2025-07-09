import numpy as np
import uproot
import glob
import re
import pickle 

from antinu_finder import antinu_finder

'''
This code employs all the analysis step until now in a way that are analized events within a range where the clockCount50 doesnt resets. This code read an initial file,
and then evaluates if the clockCount50 (cc50) resets. If not resets, then new files are readen an the information is added  to the previous readen file until the cc50 breaks.
When the  cc50 breaks, the antinu_finder algorithm is applied to this information, and once is done, the code will repeate this process over all the Analysis.root files.
This code employs valid condition selection, datacleanning, observable corrections and transformations and depends on the antinu_finder algorithm to save the
observables of interest (energy of prompt and delayed, Dt and Dr of delayed coincidences)
'''

#File Directories:--------------------------------------------------------------------------------------------

main_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/Analysis20/FF_phase/'
file_name = 'Analysis20_r*.root'

save_dir = '/lstore/sno/joankl/anti_nu/real_data/gold_data/antinu_dict/'

var_read_name_list = ['scintFit', 'fitValid', 'dcFlagged', 'nhits', 'energy', 'posx', 'posy', 'posz', 'clockCount50']

flist = glob.glob(main_dir + file_name)

#Anaysis Settings--------------------------------------------------------------------------------------------
#Base cuts: Datacleanning and min nhits:
nhits_min = 20
mask_cut = 0x2100000042C2

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

#Usefull functions ------------------------------------------------------------------------------------------
#ordenar lectura de carpetas:
def natural_order(file):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', file)]

#compute the norm of a vector
def posr_cal(x, y, z):
    dr = np.sqrt(x**2 + y**2 + z**2)
    return dr

#Main Analysis Function----------------------------------------------------------------------------------------   

# Sort the filelist name by alphanumeric order
flist.sort(key = natural_order)
total_files = len(flist)

print(f'Loaded {total_files} files')


#all_cc50 = np.array([])

# Reading the files

file_counter = 0
loop_counter = 0

while file_counter < total_files:

    if loop_counter == 0:
        decrease_index = []

    file_i = uproot.open(flist[file_counter]) 
    print(f'reading {flist[file_counter]}')     
    data_i = file_i['output;1']

    #Initially, create a new temporal dictionary with the observables of interest for each file_i in directory readen
    if len(decrease_index) == 0:
        #Condition to avoid the situation of the decrease index in the last element of the cc_50
        if loop_counter > 0:
            file_counter += 1
        else:
            dic_var_temp = {var_i: np.array(data_i[var_i]) for var_i in var_read_name_list}

    if len(decrease_index) > 0:
        dic_var_temp = {var_i: np.array(data_i[var_i])[decrease_index[0] + 1:] for var_i in var_read_name_list}
    
    print('Performing Valid and Cleaning Cuts')
    #Valid cuts, Data cleanning and cuts-------------------------------------------------
    valid_condition = (dic_var_temp['scintFit'] & dic_var_temp['fitValid'])
    nhits_min = nhits_min
    nhits_condition = (dic_var_temp['nhits'] >= nhits_min)
    mask_cut = mask_cut
    dcflag_condition = ((int(mask_cut) & dic_var_temp['dcFlagged']) == int(mask_cut))
    
    general_condition = valid_condition & nhits_condition & dcflag_condition

    for var_i in var_read_name_list:
        dic_var_temp[var_i] = np.extract(general_condition, dic_var_temp[var_i])
    #-------------------------------------------------------------------------------------
  
    #Now, for this file evaluates if the clockCount50 decrease within this file. We will use the np.diff which compute the difference between consecutive values within and array. 
    #When a diff<0 means that the cc50 breaks here

    cc50_i = dic_var_temp['clockCount50'].astype(np.float64)
    #print(cc50_i)
    cc50_dif = np.diff(cc50_i) #Compute difference between elements
    decrease_index = np.where(cc50_dif<0)[0] #evaluate if exist the index within cc50_dif which verifies the decrease of the clockCount50
    #print('decrease_index list:', decrease_index)

    #Case where the initial file readen has a cc50 break
    if len(decrease_index) > 0:
        print('decrease index within the same file found!')

        #limit the observables by this index and go to analysis
        for var_i in var_read_name_list:
            dic_var_temp[var_i] = dic_var_temp[var_i][:decrease_index[0]]
        
    
    #Case where within the file there is not a break in the clockCount50. Within this bucle, observables records will grow from reading new files until a clockCount50 break is found. 
    while len(decrease_index) == 0:
        #Entering here means that within the file the cc50 is always increasing
        #Read the next file and its observables, append the new values to the pre-existing observables in dictionary and keep evaluating if cc50 decrease.
        
        #Read the next file
        file_step_counter += 1
        new_file_i = uproot.open(flist[file_step_counter])
        new_data_i = new_file_i['output;1']
        print(f'Concatenating file {flist[file_step_counter]} in while')

        #create new dictionary to evaluate its valid cuts and data cleanning
        new_dic_var_temp = {var_i: np.array(new_data_i[var_i]) for var_i in var_read_name_list}
        
        #Valid cuts, Data cleanning and cuts-------------------------------------------------
        valid_condition = (new_dic_var_temp['scintFit'] & new_dic_var_temp['fitValid'])
        nhits_min = nhits_min
        nhits_condition = (new_dic_var_temp['nhits'] >= nhits_min)
        mask_cut = mask_cut
        dcflag_condition = ((int(mask_cut) & new_dic_var_temp['dcFlagged']) == int(mask_cut))
        
        general_condition = valid_condition & nhits_condition & dcflag_condition
    
        for var_i in var_read_name_list:
            new_dic_var_temp[var_i] = np.extract(general_condition, new_dic_var_temp[var_i])
        #-------------------------------------------------------------------------------------
        #print(f'new cc50: {new_dic_var_temp['clockCount50']}')
        #print(f'new data of file {flist[file_step_counter]} with filtered data of size: {len(new_dic_var_temp['energy'])}')

        # verify if there is a cc50 break in the new file readen:
        new_cc50 = new_dic_var_temp['clockCount50']
        #print(f'new_cc50 {new_cc50}')
        new_cc50_dif = np.diff(new_cc50)
        decrease_index = np.where(new_cc50_dif<0)[0]
        #print(f'new decrease index with elements index {new_decrease_index}') 

        # Case where the new file has not a break. just append all values to the previuous readen files in dic_var_temp
        if len(decrease_index) == 0:
            #print(f'No break of index in new file was found')
            for var_i in var_read_name_list:
                dic_var_temp[var_i] = np.append(dic_var_temp[var_i], new_dic_var_temp[var_i])
           
        # Case where there is a cc50 break in the new readen file, then just append the observable info of this file util the break index and break the while since a break was found
        if len(decrease_index) > 0:
            #print(f'break in index {new_decrease_index[0]} of new file has found')
            for var_i in var_read_name_list:
                dic_var_temp[var_i] = np.append(dic_var_temp[var_i], new_dic_var_temp[var_i][:decrease_index[0]])

    #all_cc50 = np.append(all_cc50, cc50_i)

    #get the index where cc_50 break to start reading the present file observables from exactly this index. Index is relative to the new_file_i.
    #break_index = np.where(np.diff(new_dic_var_temp['clockCount50'])<0)[0]
    #print(break_index[0])

    #Case where a cc50 break is found and we are out of the while. The analysis is performed using the last updated dic_var_temp
    if len(decrease_index) > 0:

        print('Going to analysis')
        # Once a Break is found, perform the analysis:
        #Data correction and transf:
        dic_var_temp['posz'] -= 184.4
        dic_var_temp['clockCount50'] = (dic_var_temp['clockCount50'] * 20)/ 1000 # ns -> μs
        
        posr = posr_cal(dic_var_temp['posx'], dic_var_temp['posy'], dic_var_temp['posz'])

        print(f'Performing Primal Cuts with energy: [{en_cut_inf}, {en_cut_sup}] (MeV) and posr <= {posr_cut_sup} (mm)')
        #Primary Energy and position cuts
        en_cut_inf = en_cut_inf
        en_cut_sup = en_cut_sup
        posr_cut_sup = posr_cut_sup
        
        general_condition = (dic_var_temp['energy'] >= en_cut_inf) & (dic_var_temp['energy'] <= en_cut_sup) & (posr <= posr_cut_sup)

        energy = dic_var_temp['energy'][general_condition]
        posx = dic_var_temp['posx'][general_condition]
        posy = dic_var_temp['posy'][general_condition]
        posz = dic_var_temp['posz'][general_condition]
        clockCount50 = dic_var_temp['clockCount50'][general_condition]

        print('antinu finder working ...')
        #Apply the antinu finder
        antinu_dic = antinu_finder(energy, clockCount50, posx, posy, posz, tau, alpha, dr_sup_lim, dr_inf_lim, energy_delay_inf_cut, energy_delay_sup_cut)
        print('antinu finished!')

        file_name = f'antinu_dic_out_{i_dx}'
        print(f'saving {file_name} file')

        # Save dictionary
        with open(save_dir + file_name + '.pkl', 'wb') as f:
            pickle.dump(antinu_dic, f)


    loop_counter += 1

