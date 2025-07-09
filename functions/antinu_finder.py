import numpy as np

def antinu_finder(energy, clockCount50, posx, posy, posz, tau, alpha, dr_sup_lim, dr_inf_lim, energy_delay_inf_cut, energy_delay_sup_cut):

	'''
	*Parameters:
	- energy (MeV), clockCount50 (mcs) and pos{xyz} (mm) of the events;
	- tau: half-time of neutron capture in the experiment. Value extracted from MC ~ 251 (mcs);
	- alpha: multiplicative factor of tau. It defines the temporal windows size to find a supected delayed for a suspected promtp. Value from MC ~ 9;
	- dr_{sup, inf}_lim: superior and inferior limit to restrict  the distant between supsected prompt and delayed. Values in interval [0, 1000] (mm)
	'''


	'''
	The algorithm will work as follow:
	(1) take an event with time_clock[i] and pos[i];
	(2) evaluate the difference between events dr = dr[i+1] - dr[i]
	(3) if dr lies in dr_exp +/- dr_unc or dr_exp lies in interval [dr_sup_lim, dr_inf_lim] -> save event observables 

	This logic must be applied along all the events and within a while loop until a size dt = (alpha)(tau) is reached.
	The code will use a candidate prompt and then find a delay candidate.
	'''

	#Data to save
	delta_t = []
	delta_r = []

	energy_prompt = []
	energy_delay = []

	N_ev = energy.shape[0] #number of events

	for ev_i in range(N_ev):

	    delay_index = ev_i + 1
	    #print('In event', ev_i)

	    #Suspected promt event observables
	    t1 = clockCount50[ev_i]
	    posx1 = posx[ev_i]
	    posy1 = posy[ev_i]
	    posz1 = posz[ev_i]
	    energy_1 = energy[ev_i]
	    
	    try:
	        #Suspected delayed
	        t2 = clockCount50[delay_index]
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

	    #print('intial prompt_i: ', ev_i)
	    #print('intial delay_i: ', delay_index)
	       
	    #Start looking for a delay event within a windows alpha times tau
	    while (dt > 0) and (dt <= alpha*tau):
	        
	        #print(f'prompt_i within dt interval. index = {ev_i}, t1 = {t1}')
	        #print(f'delay_i within dt interval. index = {delay_index}, t2 = {t2}')
	        #print(f'dt within interval = {dt}')

	        # If the dr coincide with the expected values, then save the observables of interest and break the while
	        #if (dr >= dr_exp - dr_unc) and (dr <= dr_exp + dr_unc) and (dt > 0) and (dt <= alpha*tau) and (energy_2 >= energy_delay_inf_cut) and (energy_2 <= energy_delay_sup_cut):
	        if (dr >= dr_inf_lim) and (dr <= dr_sup_lim) and (dt > 0) and (dt <= alpha*tau) and (energy_2 >= energy_delay_inf_cut) and (energy_2 <= energy_delay_sup_cut):
	            print('Pair Found')
	            print(f'dr: {dr} and dt: {dt}')
	            delta_t.append(dt)
	            delta_r.append(dr)
	            energy_prompt.append(energy_1)
	            energy_delay.append(energy_2)

	            break
	        # Evolve the delay index if the condition is not fulfilled
	        delay_index += 1

	        # Look for the next Suspected delayed
	        try:
	            t2 = clockCount50[delay_index]
	            #print(f'actual prompt index = {ev_i}')
	            #print(f'new delay candidate: index = {delay_index}, t2 ={t2}')
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
	            #print(f'new dt = {dt}')
	        except IndexError:
	            break

	energy_prompt = np.array(energy_prompt)
	energy_delay = np.array(energy_delay)
	delta_t = np.array(delta_t)
	delta_r = np.array(delta_r)

	dict_out = {'all_energy': energy, 
				'energy_prompt': energy_prompt,
				'energy_delay': energy_delay,
				'delta_t': delta_t,
				'delta_r': delta_r
				}