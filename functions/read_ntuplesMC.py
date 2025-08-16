'''
Function to read MC data of format .ntuplet.root and save observables of interest.
This function also implement the antinu finder algorithm for MC data based on the evIndex (evIndex == 0 are prompts) and (evIndex == 2 are delays).
This Function will Return the Number of initial events, the number of events after the applied cuts (to compute efficiency), energy of delay, energy of prompt, dt and dr.
'''

import numpy as np
import uproot
import pandas as pd
import glob
import re
import os

def extract_data(read_dir, file_txt_dir, save_dir):

    print('Intializating Functions')

    # Analysis Cuts settings ----------------------------------
    nhits_min = 20
    mask_cut = 0x2100000042C2
    energy_inf_cut = 1.0
    posr_cut = 5500
    #----------------------------------------------------------
    
    #Ensure that save directory exiss by creating it
    os.makedirs(save_dir, exist_ok=True)

    # ---------------- Subfunciones auxiliares ----------------
    def read_files_txt(file_txt_dir):
        """Lee un archivo .txt y devuelve una lista de nombres de archivo."""
        with open(file_txt_dir, "r", encoding="utf-8") as f:
            file_names = [line.strip() for line in f if line.strip()]
        return file_names

    def orden_natural(archivo):
        """Función para ordenar archivos naturalmente por número de run/subrun."""
        return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(\d+)', archivo)]

    def extract_subrunID(filename_list):
        """Extrae el subrunID de cada archivo en la lista."""
        subrunID_list = np.array([], dtype=np.int16)
        for filename_i in filename_list:
            base_filename = os.path.basename(filename_i)
            match = re.search(r'_s(\d+)_', base_filename)
            if match:
                subrun_number = int(match.group(1))
                subrunID_list = np.append(subrunID_list, subrun_number)
        return subrunID_list

    def posr_cal(x, y, z):
        """Recalcula posr usando la corrección en z."""
        return np.sqrt(x**2 + y**2 + z**2)

    # ---------------- Parámetros y variables ----------------
    var_read_name_list = [
        'evIndex', 'scintFit', 'fitValid', 'dcFlagged',
        'nhits', 'energy', 'posx', 'posy', 'posz',
        'clockCount50'
        ]

    var_save_name_list = [var + '_sv' for var in var_read_name_list]

    # Diccionario para almacenar los datos acumulados
    data_dict = {var: np.array([]) for var in var_save_name_list}
    #data_dict['n_init_events'] = np.array([])
    #data_dict['n_final_events'] = np.array([])

    #Datos de Interes para prompt-delay analysis
    energy_prompt = np.array([]) # in MeV
    energy_delay = np.array([])  # in MeV

    dt = np.array([]) # in mcs
    dr = np.array([]) # in mm


    # ---------------- Lectura de archivos ----------------
    file_name_list = read_files_txt(file_txt_dir)
    file_name_list.sort(key=orden_natural)

    print(f'list name of file to read {file_name_list}')

    for i_dx, file_name_i in enumerate(file_name_list):
        full_dir_file_i = os.path.join(read_dir, file_name_i)
        print(f'Leyendo directoria y archivo: {full_dir_file_i}')

        try:
            file_i = uproot.open(full_dir_file_i)
            output = file_i['output;1']
            print(f'Branches of file {file_i.keys()}')
        except Exception as e:
            print(f"Error al abrir archivo {file_name_i}: {e}")
            continue

        # Variables temporales para este archivo
        temp_vars = {} #temporal dictionary

        for var in var_read_name_list:
            if var == 'posz':
                temp_vars[var] = np.array(output[var]) - 184.4
            elif var == 'dcFlagged':
                temp_vars[var] = np.array(output[var]).astype(np.uint64)
            elif var == 'clockCount50':
                temp_vars[var] = (np.array(output[var]) * 20 / 1000).astype(np.uint64)  # ns -> μs
            else:
                temp_vars[var] = np.array(output[var])

        posr = posr_cal(temp_vars['posx'], temp_vars['posy'], temp_vars['posz'])
        #print('posr', type(posr))
        print(f"Events in file: {len(temp_vars['evIndex'])}")

        #Save the Initial number of prompt_events + delay_events: 
        #1) look for evIndex = 0 (prompt) and evIndex = 2 (delay)
        #2) Save the total len as the initil number of events

        #ev_index_condition = (temp_vars['evIndex'] == 0) | (temp_vars['evIndex'] == 2)
        #n_initial_evs = len(temp_vars['evIndex'][ev_index_condition])
        #data_dict['n_init_events'] = np.append(data_dict['n_init_events'], n_initial_evs)

        # Condiciones de validez
        valid_condition = (temp_vars['scintFit'] & temp_vars['fitValid'])
        nhits_condition = (temp_vars['nhits'] >= nhits_min)
        #dcflag_condition = ((int(mask_cut) & temp_vars['dcFlagged']) == int(mask_cut))
        posr_condition = (posr <= posr_cut)
        energy_condition = (temp_vars['energy'] >= energy_inf_cut)
        general_condition = valid_condition & nhits_condition & posr_condition & energy_condition

        #print(np.where(general_condition == True))        

        print('Performing validity and general cuts ...')
        # Select data that satisfy the general_condition
        for var in var_read_name_list:
            data_filtered = np.extract(general_condition, temp_vars[var])
            data_dict[var + '_sv'] = np.append(data_dict[var + '_sv'], data_filtered)

        print('Cuts performed')
        
        #Select the data_filtered to perform prompt and delay analysis
        evIndex_filtered = data_dict['evIndex_sv']
        energy_filtered = data_dict['energy_sv']
        posx_filtered = data_dict['posx_sv']
        posy_filtered = data_dict['posy_sv']
        posz_filtered = data_dict['posz_sv']
        posr_filtered = np.extract(general_condition, posr)
        time_filtered = data_dict['clockCount50_sv']

        print('posr_filtered shape', type(posr_filtered))

        #Filter the Prompt and the Delay events by evIndex = 0 and = 2, respectively
        #prompt_evIndex_condition = (evIndex_filtered == 0 & energy_filtered > 0)
        #delay_evIndex_condition = (evIndex_filtered == 2 & energy_filtered > 0)
        print('Constructing prompt and delayed matchs through evIndex')
        ev_index_condition = (evIndex_filtered == 0) | (evIndex_filtered == 2) 

        filter_index = np.where(ev_index_condition)[0]  # Indices which verifies the evIndex condition.

        #This array shows how we have consecutive eventIndex with 0's or 2's values, which indicates that we could perform a wrong evaluation of dt and dr. We should remove this isue
        test_evindex = evIndex_filtered[filter_index]

        #Remove repeated consecutive values by saving its index
        index_to_del = []  #list with the index of the elements to be removed from filter_index to correctly select our observables
        for i in range(len(test_evindex)-1):
            if test_evindex[i] == 0 and test_evindex[i+1] == 0:
                index_to_del.append(i)
            if test_evindex[i] == 2 and test_evindex[i+1] == 2:
                index_to_del.append(i+1)

        filter_index = np.delete(filter_index, index_to_del)  # cleaning of the filter_index list

        #Evaluar if the filter_index has always a pair of values which correspond to prompt and delay pairs. If not then remover the last index
        if (len(filter_index)%2  != 0):
            print('Removing the last element of filter_index due to unpair of elements.')
            filter_index = filter_index[:-1]

        # Extract observables and Arange by pairs of prompt[:,0] and delayed[:,1].
        evindex_filter_index = evIndex_filtered[filter_index].reshape((-1,2))
        energy_filter_index = energy_filtered[filter_index].reshape((-1,2))
        #posr_filter_index = posr_filtered[filter_index].reshape((-1,2))
        posx_filter_index = posx_filtered[filter_index].reshape((-1,2))
        posy_filter_index = posy_filtered[filter_index].reshape((-1,2))
        posz_filter_index = posz_filtered[filter_index].reshape((-1,2))
        time_filter_index = time_filtered[filter_index].reshape((-1,2))

        #Extract information of the prompt and delay events
        print('Extracting the Prompt and Delay quantities')
        #Energies
        energy_0 = energy_filter_index[:,0]
        energy_1 = energy_filter_index[:,1]

        #Dt
        t0 = time_filter_index[:,0]
        t1 = time_filter_index[:,1]
        Dt = t1 - t0

        #Dr
        posx_0 = posx_filter_index[:,0]
        posy_0 = posy_filter_index[:,0]
        posz_0 = posz_filter_index[:,0]

        posx_1 = posx_filter_index[:,1]
        posy_1 = posy_filter_index[:,1]
        posz_1 = posz_filter_index[:,1]

        dx = posx_1 - posx_0
        dy = posy_1 - posy_0
        dz = posz_1 - posz_0

        Dr = np.sqrt((dx**2) + (dy**2) + (dz**2))


        # Append the observable values to the external file loop empty arrays and save at the same time:
        energy_prompt = np.append(energy_prompt, energy_0)
        energy_delay = np.append(energy_delay, energy_1)
        dt = np.append(dt, Dt)
        dr = np.append(dr, Dr)

        #Save the final number of total events after cuts:
        #data_dict['n_final_events'] = np.append(data_dict['n_final_events'], len(energy_prompt) + len(energy_delay))

        # ---------- Save the information ----------------
        print('Saving information')
        #np.save(save_dir + 'n_init_evs.npy', data_dict['n_init_events'])
        #np.save(save_dir + 'n_final_evs.npy', data_dict['n_final_events'])
        np.save(save_dir + 'energy_prompt.npy', energy_prompt)
        np.save(save_dir + 'energy_delay.npy', energy_delay)
        np.save(save_dir + 'dt.npy', dt)
        np.save(save_dir + 'dr.npy', dr)

    return print('Extraction concluded!')

'''
if __name__ == "__main__":
    read_dir = '/share/neutrino/snoplus/MonteCarlo/FullFill_2p2_709/ScintFit_2p2ReactoribdRun/'
    file_txt_dir = '/lstore/sno/joankl/anti_nu/MonteCarlo/reactor_nu_ibd_data/file_name_list/sublist_0.txt'
    save_dir = '/lstore/sno/joankl/anti_nu/MonteCarlo/reactor_nu_ibd_data/out_results/'
    extract_data(read_dir, file_txt_dir, save_dir)
'''
