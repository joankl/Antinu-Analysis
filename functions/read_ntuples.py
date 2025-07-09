'''                                                                                                                                                               
Function to read the ntuples.root of real data. It will extract the variables of interest and save it in separeted numpy arrays (faste to save and lighter than files.csv!)
'''

import numpy as np
import uproot
import pandas as pd
import glob
import re
import os

def extract_data(read_dir, file_txt_dir, save_dir):
    
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
        'runID', 'eventID', 'scintFit', 'fitValid', 'dcFlagged', 'triggerWord',
        'nhits', 'nhitsCleaned', 'correctedNhits', 'necknhits', 'posFOM',
        'posFOM2', 'itr', 'energy', 'posx', 'posy', 'posz',
        'clockCount50', 'beta14', 'alphaBeta212', 'alphaBeta214'
    ]

    var_save_name_list = [var + '_sv' for var in var_read_name_list]

    # Diccionario para almacenar los datos acumulados
    data_dict = {var: np.array([]) for var in var_save_name_list}
    subrunID_sv = np.array([], dtype=np.int16)

    # ---------------- Lectura de archivos ----------------
    file_name_list = read_files_txt(file_txt_dir)
    file_name_list.sort(key=orden_natural)

    for i_dx, file_name_i in enumerate(file_name_list):
        full_dir_file_i = os.path.join(read_dir, file_name_i)
        print(f'Leyendo archivo: {full_dir_file_i}')

        try:
            file_i = uproot.open(full_dir_file_i)
            output = file_i['output;1']
        except Exception as e:
            print(f"Error al abrir archivo {file_name_i}: {e}")
            continue

        # Variables temporales para este archivo
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

        # Verifica que runID exista antes de usarlo
        if 'runID' not in temp_vars or len(temp_vars['runID']) == 0:
            print(f"Archivo vacío o sin runID: {file_name_i}")
            continue

        print(f"Eventos en archivo: {len(temp_vars['runID'])}")

        # Subrun ID por archivo
        subrunID_number = extract_subrunID([full_dir_file_i])
        if len(subrunID_number) == 0:
            print(f"No se pudo extraer subrunID de {file_name_i}")
            continue

        subrunID_array = np.ones(len(temp_vars['runID'])) * subrunID_number[0]

        # Condiciones de validez
        valid_condition = (temp_vars['scintFit'] & temp_vars['fitValid'])
        nhits_min = 20
        nhits_condition = (temp_vars['nhits'] >= nhits_min)
        mask_cut = 0x2100000042C2
        dcflag_condition = ((int(mask_cut) & temp_vars['dcFlagged']) == int(mask_cut))
        general_condition = valid_condition & nhits_condition & dcflag_condition

        # Aplicar condiciones y acumular datos
        subrunID_sv = np.append(subrunID_sv, subrunID_array[general_condition].astype(np.int16))

        np.save(os.path.join(save_dir, 'subrunID'), subrunID_sv)
       
        for var in var_read_name_list:
            data_filtered = np.extract(general_condition, temp_vars[var])
            data_dict[var + '_sv'] = np.append(data_dict[var + '_sv'], data_filtered)

    # ---------------- Guardado final ----------------
        for var_sv in var_save_name_list:
            np.save(os.path.join(save_dir, var_sv.replace('_sv', '')), data_dict[var_sv])


    return print('extraction concluded!')

        

