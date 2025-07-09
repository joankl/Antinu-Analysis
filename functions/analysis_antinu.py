'''
This Analysis.py is used to read the files.csv obtained by the summary_antinu_info.py and saved in all_files directory to extract useful information and concatenated information. 
This will return the energy, reconstructed position and time_clock50 of all event in these files.csv.
The code also contains the logic to obtain the Dt and Dr results of suspected antinu events and save this information and the energy of the prompt and delayed as well.
However, since the files.csv are readen in a loop over runs (because each file correspond to a subrun), a prompt event may have a delayed in a separated file, which is 
not considered here.
'''

#Read Files
read_route = 'C:/snoplus simulations/real_data_antinu/bronze_data/summary_files/'
pattern_file = 'Analysis20_nu_real_data_N_*.csv'

#Save file route
save_route = ''

#Function to read the files by runID order which is important when computing the timeclock!
def orden_natural(archivo):
	return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(/d+)', archivo)]

file_list = glob.glob(read_route+pattern_file)
file_list.sort(key = orden_natural)

#Adjust cuts -----------------------------------------------

#Main cuts on data:
en_cut_inf = 0.5  #To avoid accidentals and retriggers
en_cut_sup = 8    # Max energy to be observed
pos_r_cut = 5500  #Select FV
dc_flag_cut = 0x2100000042C2

#antinu event finder algorithm settings:
# Values selected by looking at the delta_t and delta_r distribution of MC reactor nu data and the fit of the curves
alpha = 8 

#Expected values dr
tau = 251
dr_exp = 357
dr_unc = 200 #choose a large interval on uncertainty otherwise we are rejecting a lot of events!

dr_sup_lim = 1000 #choose a large interval on interval otherwise we are rejecting a lot of events!
dr_inf_lim = 0

#Energy cut on delayed events
energy_delay_inf_cut = 1.8
energy_delay_sup_cut = 3
# ----------------------------------------------------------

#Data to save within the cuts:
all_energy_spec = [] #Energy spectrum of all data
all_timeclock50 = [] #Time clock 50 of all data
all_ev_recons_xyz = np.empty((0, 3)) #all events reconstructed position coordinates (x,y,z)
energy_prompt = [] #energy of suspected prompt
energy_delay = []  #energy of suspected delay
delta_t = [] #time interval between suspected prompt and delay
delta_r = [] #space interval between suspected prompt and delay.

#Read files and do analysis
for i_dx, route_i in enumerate(file_list):
    print('reading file:', route_i)

    #Extract the data from the ith file:
    file_i = pd.read_csv(route_i)  #All data in pandas Dataframe format

    #Cut to be applied in general data
    condition_cut = (file_i['energy'] > en_cut_inf) & (file_i['energy'] <= en_cut_sup) & (file_i['posr'] <= pos_r_cut) & (file_i['dcFlag'] == int(dc_flag_cut))

    # Extract the file with cuts
    file_cut = file_i[condition_cut].reset_index()

    #Extract the observables of interest
    energy = file_cut['energy'].to_numpy()
    posx = file_cut['posx'].to_numpy()
    posy = file_cut['posy'].to_numpy()
    posz = file_cut['posz'].to_numpy()
    pos_vec_xyz = np.array([posx, posy, posz])
    time_clock = file_cut['time_clock50 (mcs)'].to_numpy()

    #Save the overall energy spectrum
    all_energy_spec.append(energy)
    all_ev_recons_xyz = np.concatenate((all_ev_recons_xyz, pos_vec_xyz.T))

    #Algorithm to look for antinu events: First take a suspected prompt and look for a compatible delayed
    '''
    The algorithm will work as follow:
    (1) take an event with time_clock[i] and pos[i];
    (2) evaluate the difference between events dr = dr[i+1] - dr[i]
    (3) if dr lies in dr_exp +/- dr_unc or dr_exp lies in interval [dr_sup_lim, dr_inf_lim] -> save event observables 
    
    This logic must be applied along all the events and within a while loop until a size dt = (alpha)(tau) is reached.
    The code will use a candidate prompt and then find a delay candidate.
    '''
    
    N_ev = energy.shape[0]  # Number of events in file_i over which iterate
    
    for ev_i in range(N_ev):
    
        delay_index = ev_i + 1
        #print('In event', ev_i)
    
        #Suspected promt event observables
        t1 = time_clock[ev_i]
        posx1 = posx[ev_i]
        posy1 = posy[ev_i]
        posz1 = posz[ev_i]
        energy_1 = energy[ev_i]
        
        try:
            #Suspected delayed
            t2 = time_clock[delay_index]
            posx2 = posx[delay_index]
            posy2 = posy[delay_index]
            posz2 = posz[delay_index]
            energy_2 = energy[delay_index]
            
            #Evaluate the dr and dt between events
            dx = posx2 - posx1
            dy = posy2 - posy1
            dz = posz2 - posz1
            dr = np.sqrt(dx**2 + dy**2 + dz**2)
            dt = t2 - t1
            
        except IndexError:
            continue
           
        #Start looking for a delay event within a windows alpha times tau
        while (dt > 0) and (dt <= alpha*tau):
    
            # If the dr coincide with the expected values, then save the observables of interest and break the while
            #if (dr >= dr_exp - dr_unc) and (dr <= dr_exp + dr_unc) and (dt > 0) and (dt <= alpha*tau) and (energy_2 >= energy_delay_inf_cut) and (energy_2 <= energy_delay_sup_cut):
            if (dr >= dr_inf_lim) and (dr <= dr_sup_lim) and (dt > 0) and (dt <= alpha*tau) and (energy_2 >= energy_delay_inf_cut) and (energy_2 <= energy_delay_sup_cut):
                print('Pair Found')
                print(dr)
                delta_t.append(dt)
                delta_r.append(dr)
                energy_prompt.append(energy_1)
                energy_delay.append(energy_2)
    
                break
            # Evolve the delay index if the condition is not fulfilled
            delay_index += 1
    
            # Look for the next Suspected delayed
            try:
                t2 = time_clock[delay_index]
                posx2 = posx[delay_index]
                posy2 = posy[delay_index]
                posz2 = posz[delay_index]
                energy_2 = energy[delay_index]
        
                #Evaluate the dr and dt between events
                dx = posx2 - posx1
                dy = posy2 - posy1
                dz = posz2 - posz1
                dr = np.sqrt(dx**2 + dy**2 + dz**2)
                dt = t2 - t1
            except IndexError:
                break

        
    #Transform to numpy array and save relevant data
    print(f'saving file information with index file {i_dx}')
    np.save(save_route + 'all_energy.npy', np.array(all_energy_spec).reshape((-1,))) #All Energy
    np.save(save_route + 'all_ev_recons_xyz.npy', all_ev_recons_xyz)
    np.save(save_route + 'energy_prompt.npy', np.array(energy_prompt)) #prompt energy
    np.save(save_route + 'energy_delay.npy', np.array(energy_delay)) #delay energy
    np.save(save_route + 'dt.npy', np.array(delta_t)) #delta_t
    np.save(save_route + 'dr.npy', np.array(delta_r)) #delta_r