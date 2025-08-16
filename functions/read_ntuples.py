'''                                                                                                                                                               
Function to read the ntuples.root of real data. It will extract the variables of interest and save it in separeted numpy arrays (faste to save and lighter than files.csv!)
'''

import numpy as np
import uproot
import re
import os

def extract_data(read_dir, file_txt_dir, save_dir):
    # Ensure save directory exists
    os.makedirs(save_dir, exist_ok=True)

    # ---------------- Helper functions ----------------
    def read_files_txt(file_txt_dir):
        """Read a .txt file and return a list of file names."""
        with open(file_txt_dir, "r", encoding="utf-8") as f:
            return [line.strip() for line in f if line.strip()]

    def orden_natural(archivo):
        """Sort files naturally by run/subrun number."""
        return [int(texto) if texto.isdigit() else texto.lower()
                for texto in re.split(r'(\d+)', archivo)]

    def extract_subrunID(filename_list):
        """Extract subrunID from file names."""
        subrunID_list = np.array([], dtype=np.int16)
        for filename_i in filename_list:
            base_filename = os.path.basename(filename_i)
            match = re.search(r'_s(\d+)_', base_filename)
            if match:
                subrunID_list = np.append(subrunID_list, int(match.group(1)))
        return subrunID_list

    def posr_cal(x, y, z):
        """Recalculate posr using z correction."""
        return np.sqrt(x**2 + y**2 + z**2)

    # ---------------- Variables ----------------
    var_read_name_list = [
        'runID', 'eventID', 'scintFit', 'fitValid', 'dcFlagged', 'triggerWord',
        'nhits', 'nhitsCleaned', 'correctedNhits', 'necknhits', 'posFOM',
        'posFOM2', 'itr', 'energy', 'posx', 'posy', 'posz',
        'clockCount50', 'beta14', 'alphaBeta212', 'alphaBeta214'
    ]

    var_save_name_list = [var + '_sv' for var in var_read_name_list]
    data_dict = {var: np.array([]) for var in var_save_name_list}
    subrunID_sv = np.array([], dtype=np.int16)

    # ---------------- Read file list ----------------
    file_name_list = read_files_txt(file_txt_dir)
    file_name_list.sort(key=orden_natural)

    # ---------------- Process each file ----------------
    for file_name_i in file_name_list:
        full_dir_file_i = os.path.join(read_dir, file_name_i)
        print(f"Reading file: {full_dir_file_i}")

        try:
            file_i = uproot.open(full_dir_file_i)
            output = file_i['output;1']
        except Exception as e:
            print(f"Error opening file {file_name_i}: {e}")
            continue

        # Temporary dictionary for this file
        temp_vars = {}
        for var in var_read_name_list:
            if var == 'posz':
                temp_vars[var] = np.array(output[var]) - 184.4
            elif var == 'dcFlagged':
                temp_vars[var] = np.array(output[var]).astype(np.uint64)
            elif var == 'clockCount50':
                temp_vars[var] = (np.array(output[var]) * 20 / 1000).astype(np.uint64)  # ns -> μs
            else:
                temp_vars[var] = np.array(output[var])

        # Skip empty files
        if 'runID' not in temp_vars or len(temp_vars['runID']) == 0:
            print(f"Empty file or missing runID: {file_name_i}")
            continue

        # Extract subrun ID
        subrunID_number = extract_subrunID([full_dir_file_i])
        if len(subrunID_number) == 0:
            print(f"Could not extract subrunID from {file_name_i}")
            continue
        subrunID_array = np.ones(len(temp_vars['runID'])) * subrunID_number[0]

        # Apply cuts
        valid_condition = (temp_vars['scintFit'] & temp_vars['fitValid'])
        nhits_condition = (temp_vars['nhits'] >= 20)
        mask_cut = 0x2100000042C2
        dcflag_condition = ((int(mask_cut) & temp_vars['dcFlagged']) == int(mask_cut))
        general_condition = valid_condition & nhits_condition & dcflag_condition

        # Accumulate filtered data
        subrunID_sv = np.append(subrunID_sv, subrunID_array[general_condition].astype(np.int16))
        for var in var_read_name_list:
            data_filtered = np.extract(general_condition, temp_vars[var])
            data_dict[var + '_sv'] = np.append(data_dict[var + '_sv'], data_filtered)

    # ---------------- Save all accumulated data once ----------------
    np.save(os.path.join(save_dir, 'subrunID'), subrunID_sv)
    for var_sv in var_save_name_list:
        np.save(os.path.join(save_dir, var_sv.replace('_sv', '')), data_dict[var_sv])

    print("Extraction concluded!")


if __name__ == "__main__":
    read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/Analysis20R/Bronze/'
    file_txt_dir = '/lstore/sno/joankl/anti_nu/real_data/bronze_data/file_name_list/name_list_analysis/sublist_0.txt'
    save_dir = '/lstore/sno/joankl/anti_nu/real_data/bronze_data/proof'
    extract_data(read_dir, file_txt_dir, save_dir)


        

