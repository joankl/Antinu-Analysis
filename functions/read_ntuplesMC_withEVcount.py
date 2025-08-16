'''
Function to read MC data of format .ntuplet.root and save observables of interest.
This function also implement the antinu finder algorithm for MC data based on the evIndex (evIndex <= 0 are all the generated events)
This Function will Return the Number of initial events, the number of events after the applied cuts (to compute efficiency), energy of delay, energy of prompt, dt and dr.
'''

import numpy as np
import uproot
import pandas as pd
import glob
import re
import os

import numpy as np
import uproot
import re
import os

def extract_data(read_dir, file_txt_dir, save_dir):
    """
    Reads MC data (.ntuplet.root) and applies coincidence-finding and cuts.
    Now counts the number of valid coincidence pairs per file, before and after cuts.
    Saves energy_prompt, energy_delay, dt, dr, plus counts per file.
    """

    print('Initializing Functions')

    # ------------------- Cuts settings -------------------
    nhits_min = 20
    mask_cut = 0x2100000042C2
    energy_inf_cut = 1.0
    posr_cut = 5500
    # ------------------------------------------------------

    os.makedirs(save_dir, exist_ok=True)

    # ---------------- Helper functions -------------------
    def read_files_txt(file_txt_dir):
        """Reads a .txt file containing the list of ROOT files."""
        with open(file_txt_dir, "r", encoding="utf-8") as f:
            return [line.strip() for line in f if line.strip()]

    def orden_natural(archivo):
        """Natural sort key."""
        return [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', archivo)]

    def posr_cal(x, y, z):
        """Recalculates posr with z-correction."""
        return np.sqrt(x**2 + y**2 + z**2)

    # --------------- Variables to Read and Save -----------------------------

    var_read_name_list = [
        'evIndex', 'scintFit', 'fitValid', 'dcFlagged',
        'nhits', 'energy', 'posx', 'posy', 'posz',
        'clockCount50'
    ]

    # Arrays to store per-file coincidence counts
    initial_events_per_file = []
    final_pairs_per_file = []

    # Arrays to store observables across all files
    energy_prompt_all = np.array([])
    energy_delay_all = np.array([])
    dt_all = np.array([])
    dr_all = np.array([])

    # ---------------- File list ----------------
    file_name_list = read_files_txt(file_txt_dir)
    file_name_list.sort(key=orden_natural)
    print(f'Files to read: {file_name_list}')

    for file_name_i in file_name_list:
        full_dir_file_i = os.path.join(read_dir, file_name_i)
        print(f'Reading: {full_dir_file_i}')

        try:
            file_i = uproot.open(full_dir_file_i)
            output = file_i['output;1']
        except Exception as e:
            print(f"Error opening file {file_name_i}: {e}")
            continue

        # Temporary vars for this file
        temp_vars = {}
        for var in var_read_name_list:
            if var == 'posz':
                temp_vars[var] = np.array(output[var]) - 184.4
            elif var == 'dcFlagged':
                temp_vars[var] = np.array(output[var]).astype(np.uint64)
            elif var == 'clockCount50':
                temp_vars[var] = (np.array(output[var]) * 20 / 1000).astype(np.uint64)
            else:
                temp_vars[var] = np.array(output[var])

        posr = posr_cal(temp_vars['posx'], temp_vars['posy'], temp_vars['posz'])

        # -------- Count all the MC generated events --------
        evIndex = temp_vars['evIndex']
        filter_index = np.where((evIndex <= 0))[0]

        n_init_events = len(filter_index)
        initial_events_per_file.append(n_init_events)
        print(f'Initial MC generated events = {n_init_events}')

        # -------- Apply cuts --------
        valid_condition = (temp_vars['scintFit'] & temp_vars['fitValid'])
        nhits_condition = (temp_vars['nhits'] >= nhits_min)
        posr_condition = (posr <= posr_cut)
        energy_condition = (temp_vars['energy'] >= energy_inf_cut)
        general_condition = valid_condition & nhits_condition & posr_condition & energy_condition

        # Filtered variables after cuts
        evIndex_f = np.extract(general_condition, evIndex)
        energy_f = np.extract(general_condition, temp_vars['energy'])
        posx_f = np.extract(general_condition, temp_vars['posx'])
        posy_f = np.extract(general_condition, temp_vars['posy'])
        posz_f = np.extract(general_condition, temp_vars['posz'])
        time_f = np.extract(general_condition, temp_vars['clockCount50'])

        # -------- Count final valid coincidence pairs (after cuts) --------
        filter_index_final = np.where((evIndex_f == 0) | (evIndex_f == 2))[0]

        # Remove retriggers
        idx_to_remove_final = []
        test_evindex_final = evIndex_f[filter_index_final]
        for i in range(len(test_evindex_final)-1):
            if test_evindex_final[i] == 0 and test_evindex_final[i+1] == 0:
                idx_to_remove_final.append(i)
            if test_evindex_final[i] == 2 and test_evindex_final[i+1] == 2:
                idx_to_remove_final.append(i+1)
        filter_index_final = np.delete(filter_index_final, idx_to_remove_final)

        if len(filter_index_final) % 2 != 0:
            filter_index_final = filter_index_final[:-1]

        n_final_pairs = len(filter_index_final) // 2
        final_pairs_per_file.append(n_final_pairs)

        # -------- Extract observables for final pairs --------
        energy_pairs = energy_f[filter_index_final].reshape((-1, 2))
        posx_pairs = posx_f[filter_index_final].reshape((-1, 2))
        posy_pairs = posy_f[filter_index_final].reshape((-1, 2))
        posz_pairs = posz_f[filter_index_final].reshape((-1, 2))
        time_pairs = time_f[filter_index_final].reshape((-1, 2))

        energy_prompt_all = np.append(energy_prompt_all, energy_pairs[:, 0])
        energy_delay_all = np.append(energy_delay_all, energy_pairs[:, 1])

        Dt = time_pairs[:, 1] - time_pairs[:, 0]
        dt_all = np.append(dt_all, Dt)

        dx = posx_pairs[:, 1] - posx_pairs[:, 0]
        dy = posy_pairs[:, 1] - posy_pairs[:, 0]
        dz = posz_pairs[:, 1] - posz_pairs[:, 0]
        Dr = np.sqrt(dx**2 + dy**2 + dz**2)
        dr_all = np.append(dr_all, Dr)

    # -------- Save results --------
    np.save(os.path.join(save_dir, 'initial_pairs.npy'), np.array(initial_events_per_file))
    np.save(os.path.join(save_dir, 'final_pairs.npy'), np.array(final_pairs_per_file))
    np.save(os.path.join(save_dir, 'energy_prompt.npy'), energy_prompt_all)
    np.save(os.path.join(save_dir, 'energy_delay.npy'), energy_delay_all)
    np.save(os.path.join(save_dir, 'dt.npy'), dt_all)
    np.save(os.path.join(save_dir, 'dr.npy'), dr_all)

    print('Extraction concluded!')
    print(f"Initial valid coincidence pairs per file: {initial_events_per_file}")
    print(f"Final valid coincidence pairs per file: {final_pairs_per_file}")

if __name__ == "__main__":
    read_dir = '/share/neutrino/snoplus/MonteCarlo/FullFill_2p2_709/ScintFit_2p2ReactoribdRun/'
    file_txt_dir = '/lstore/sno/joankl/anti_nu/MonteCarlo/reactor_nu_ibd_data/file_name_list/sublist_0.txt'
    save_dir = '/lstore/sno/joankl/anti_nu/MonteCarlo/reactor_nu_ibd_data/proof/'
    extract_data(read_dir, file_txt_dir, save_dir)
