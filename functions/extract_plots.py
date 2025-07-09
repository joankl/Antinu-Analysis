# Code to extract plots from csv filtered data for antinu analysis

import numpy as np 
import pandas as pd 

import matplotlib.pyplot as plt
import seaborn as sn 

#read file
data_route = ''
#save figs
save_route = ''

file = pd.read_csv(data_route)

#cuts:
en_cut_inf = 0
en_cut_sup = 8
pos_r_cut = 5500
condition_cut = (file['energy'] >= en_cut_inf) & (file['energy'] <= en_cut_sup) & (file['pos_r'] <= pos_r_cut)
file_cut = file[condition_cut].reset_index()

#Extract Observables
energy = file_cut['energy']
posx_ev = file_cut['posx']
posy_ev = file_cut['posy']
posz_ev = file_cut['posz']
time_clock = file_cut['time_clock (mcs)']

#Plots
sn.histplot()
