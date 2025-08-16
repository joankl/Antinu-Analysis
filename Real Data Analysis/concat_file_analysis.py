'''
Code to concatenate all the output files.csv in the folder all_files
'''

import pandas as pd
#import polars as pl
import glob
import re

# Function to read in order the files
def orden_natural(archivo):
        return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(/d+)', archivo)]

C:\snoplus simulations\real_data_antinu\bronze_data\all_files

#input info
read_file_dir = 'C:/snoplus simulations/real_data_antinu/bronze_data/all_files' # directory of files to be concatenated
in_file_name = 'Analysis20_nu_real_data_N_1.csv' #input file name pattern

#Output info
out_file_dir = 'C:/snoplus simulations/real_data_antinu/bronze_data/'  # directory to save the file
out_file_name = 'Bronze_data_Analysis20.csv'  # output file name

# list of files to be readen:
files = glob.glob(read_file_dir + in_file_name)
files.sort(key = orden_natural)

print('Number of input files:', len(files))

#DataFrame Main Structure to be fill
df = pd.DataFrame(columns = ['fitValid', 'scintFit', 'runID', 'dcFlag', 'evIndex',
                                 'evID', 'triggerWord', 'NHits', 'NhitsClean', 'correctedNhits',
                                 'neckNhits','energy', 'posx', 'posy', 'posz', 'posr', 'posFOM',
                                 'posFOM2', 'time_clock50 (mcs)', 'ITR', 'beta14', 'alphabeta212', 'alphabeta214'])

#read and concatenate using pd.dataframe
for file_i in files:
	print('reading file ' + file_i)
	df = pd.concat([pd.read_csv(file_i), df], ignore_index = True)
	df.to_csv(out_file_dir + out_file_name, index = False)

print('files concatenated sucessfully!')

#read and concatenate using polars (more efficient)
#df = pl.concat([pl.read_csv(file_i) for file_i in files])
#df.write_csv(out_file_dir + out_file_name)

print('file ' + out_file_name + ' created sucessfully!')
