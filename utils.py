'''
utils.py

Utility class for holding extra methods from mesa.py, nuh5p.py


'''

class Utils():
	
	#List of element names and stable elements for mppnp.py
	elements_names = ['Neutron','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
	    'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
	    'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',
	    'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te', 'I','Xe','Cs','Ba','La','Ce',
	    'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',
	    'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu']

    	stable_el = [['Neutron','999'],['H',1, 2],['He', 3, 4],['Li', 6, 7],['Be', 9],['B', 10, 11],['C', 12, 13],['N', 14, 15],['O', 16, 17, 18],['F', 19],['Ne', 20, 21, 22],['Na', 23],['Mg', 24, 25, 26],['Al', 27],['Si', 28, 29, 30],['P', 31],['S', 32, 33, 34, 36],['Cl', 35, 37],['Ar', 36, 38, 40],['K', 39, 40, 41],['Ca', 40, 42, 43, 44, 46, 48],['Sc', 45],['Ti', 46, 47, 48, 49, 50],['V', 50, 51],['Cr', 50, 52, 53, 54],['Mn', 55],['Fe', 54, 56, 57, 58],['Co', 59],['Ni', 58, 60, 61, 62, 64],['Cu', 63, 65],['Zn', 64, 66, 67, 68, 70],['Ga', 69, 71],['Ge', 70, 72, 73, 74, 76],['As', 75],['Se', 74, 76, 77, 78, 80, 82],['Br', 79, 81],['Kr', 78, 80, 82, 83, 84, 86],['Rb', 85, 87],['Sr', 84, 86, 87, 88],['Y', 89],['Zr', 90, 91, 92, 94, 96],['Nb', 93],['Mo', 92, 94, 95, 96, 97, 98, 100],['Tc',999],['Ru', 96, 98, 99, 100, 101, 102, 104],['Rh', 103],['Pd', 102, 104, 105, 106, 108, 110],['Ag', 107, 109],['Cd', 106, 108, 110, 111, 112, 113, 114, 116],['In', 113, 115],['Sn', 112, 114, 115, 116, 117, 118, 119, 120, 122, 124],['Sb', 121, 123],['Te', 120, 122, 123, 124, 125, 126, 128, 130],['I', 127],['Xe', 124, 126, 128, 129, 130, 131, 132, 134, 136],['Cs', 133],['Ba', 130, 132, 134, 135, 136, 137, 138],['La', 138, 139],['Ce', 136, 138, 140, 142],['Pr', 141],['Nd', 142, 143, 144, 145, 146, 148, 150],['Pm',999],['Sm', 144, 147, 148, 149, 150, 152, 154],['Eu', 151, 153],['Gd', 152, 154, 155, 156, 157, 158, 160],['Tb', 159],['Dy', 156, 158, 160, 161, 162, 163, 164],['Ho', 165],['Er', 162, 164, 166, 167, 168, 170],['Tm', 169],['Yb', 168, 170, 171, 172, 173, 174, 176],['Lu', 175, 176],['Hf', 174, 176, 177, 178, 179, 180],['Ta', 180, 181],['W', 180, 182, 183, 184, 186],['Re', 185, 187],['Os', 184, 186, 187, 188, 189, 190, 192],['Ir', 191, 193],['Pt', 190, 192, 194, 195, 196, 198],['Au', 197],['Hg', 196, 198, 199, 200, 201, 202, 204],['Tl', 203, 205],['Pb', 204, 206, 207, 208],['Bi', 209],['Th', 232],['U',235,238]] 

def close_wins(win_min,win_max):
	''' close all windows in a certain window number range
		
	win_min/max  minumum and maximum window number to close
	'''

	for i in range(win_min,win_max+1):
		close(i)

def species_list(what_list):
	''' provide default lists of elements to plot
	what_list   string name of species lists provided
	CNONe       C,N, O and some other light elements
	s-process   s-process indicators
	'''
	if what_list is "CNONe":
		list_to_print = ['H-1','He-4','C-12','O-16','Ne-20']
	elif what_list is "sprocess":
		list_to_print = ['Fe-56','Ge-70','Zn-70','Se-76','Kr-80','Kr-82','Kr-86','Sr-88','Ba-138','Pb-208']
	elif what_list is "burn_stages":
		list_to_print = ['H-1','He-4','C-12','O-16','Ne-20','Si-28']	
	elif what_list is "list_marco_1":
		list_to_print = ['C-12','O-16','Ne-20','Ne-22','Na-23','Fe-54','Fe-56','Zn-70','Ge-70','Se-76','Kr-80','Kr-82','Sr-88','Y-89','Zr-96','Te-124','Xe-130','Xe-134','Ba-138']	
			
	return list_to_print
	
def symbol_list(what_list):
	''' provide default symbol lists
	what_list   string name of symbol lists provided
	list1, list2
	'''
	if what_list is "list1":
		symbol=['ro','bo','ko','go','mo'\
			,'r-','b-','k-','g-','m-','r--','b--','k--'\
			,'g--','r1']
		#symbol=['r+','ro','r-']
	elif what_list is "list2":
		symbol=['r-','b--','g-.','k:','md','.','o','v','^','<','>','1','2',\
			'3','4','s','p','*','h','H','+']
	elif what_list is "lines1":
		symbol=['r-','b--','g-.','k:','c-','m-','b-','g--','k-.','c:','m-','b-']
			
	return symbol
		
def make_list(default_symbol_list,len_list_to_print):
	''' provide the list of symbols to use according for the list of species/arrays to plot.
	default_symbol_list = list of symbols that the user choose to use.  	    	
	len_list_to_print   = len of list of species/arrays to print.
	'''

	import scipy as sc

	symbol_used = []
	for i in range(len_list_to_print):
		symbol_used.append(default_symbol_list[sc.mod(i,len(default_symbol_list))])
			
	return symbol_used
		
def solar(filename_solar,solar_factor):
    ''' read solar abundances from filename_solar.
    solar_factor is the correction factor to apply, in case filename_solar is not solar,
    but some file used to get initial abundances at metallicity lower than solar. However, notice 
    that this is really  rude, since alpha-enahncements and things like that are not properly considered.
    Only H and He4 are not multiplied. So, for publications PLEASE use proper filename_solar at...solar, and 
    use solar_factor = 1. Marco 	 '''

    import numpy as np	

    f0=open(filename_solar)
    sol=f0.readlines()
    f0.close 
    sol[0].split("         ")

    # Now read in the whole file and create a hashed array:
    global names_sol
    names_sol=[]
    global z_sol
    z_sol=[]    
    yps=np.zeros(len(sol))
    mass_number=np.zeros(len(sol))
    for i in range(len(sol)):
        z_sol.append(int(sol[i][1:3]))
        names_sol.extend([sol[i].split("         ")[0][4:]])
        yps[i]=float(sol[i].split("         ")[1]) * solar_factor
        mass_number[i]=int(names_sol[i][2:5])
	if mass_number[i] == 1 or mass_number[i] == 4:
		yps[i] = yps[i]/solar_factor
    #  convert 'h   1' in prot, not needed any more??
    #names_sol[0] = 'prot '
    
    
    # now zip them together:
    global solar_abundance
    solar_abundance={}
    for a,b in zip(names_sol,yps):
        solar_abundance[a] = b



    global z_bismuth
    z_bismuth = 83
    global z_for_elem
    z_for_elem = []
    global index_stable
    index_stable = []
    global solar_elem_abund
    solar_elem_abund = np.zeros(z_bismuth)

    i_for_stable = 1
    i_for_unstable = 0
    for i in range(z_bismuth):
        z_for_elem.append(int(i+1))
    	# the only elements below bismuth with no stable isotopes are Tc and Pm
    	if i+1 == 43 or i+1 == 61:
        	index_stable.append(i_for_unstable) 
    	else:
        	index_stable.append(i_for_stable)

    for i in range(z_bismuth):
        dummy = 0.
        for j in range(len(solar_abundance)):
            if z_sol[j] == i+1:
                dummy = dummy + float(solar_abundance[names_sol[j]])
    	solar_elem_abund[i] = dummy

