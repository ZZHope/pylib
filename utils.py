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


def   convert_specie_naming_from_h5_to_ppn(isotope_names):
    	''' read isotopes names from h5 files, and convert them according to standard scheme used inside ppn and mppnp.
    	Also Z and A are recalculated, for these species. Isomers are excluded for now, since there were recent changes 
    	in isomers name. As soon as the isomers names are settled, than Z and A provided here will be obsolete, and can be
    	changed by usual Z and A. 	 '''

	import numpy as np

        spe_rude1 = []
        spe_rude2 = []
        spe_rude3 = []
        for i in range(len(isotope_names)):
            spe_rude1.append(isotope_names[i].split('-')[0])
            spe_rude2.append(isotope_names[i].split('-')[1])
        # spe_rude1 is elem name and spe_rude2 is mass number.
        #print spe_rude1,spe_rude2
        k = 0
        for i in range(len(spe_rude1)):
            try: 
            	if int(spe_rude2[i]) < 10:
            	    spe_rude3.append(str(spe_rude1[i][0:2])+str('  ')+str(spe_rude2[i][0:3]))    
            	elif int(spe_rude2[i]) >= 10 and int(spe_rude2[i]) < 100 :    
            	    spe_rude3.append(str(spe_rude1[i][0:2])+str(' ')+str(spe_rude2[i][0:3]))    
            	elif int(spe_rude2[i]) >= 100 :    
            	    spe_rude3.append(str(spe_rude1[i][0:2])+str(spe_rude2[i][0:3]))    
            except ValueError:
                k = k+1
                None

        global spe
        spe = []
	global n_array
        n_array = []
        for i in range(len(spe_rude3)):
            if len(str(spe_rude1[i])) == 1:
                spe.append(str(spe_rude3[i][0:1])+str(' ')+str(spe_rude3[i][1:4]))        
            else:
                spe.append(spe_rude3[i]) 
	    n_array.append(i)
        if spe[0]=='Ne  1':
		spe[0] = 'N   1' 
                                          
        # spe_rude is the isotope name, in agreement with what we use in ppn, etc.
        # need to do this to can use other functions without changing them drastically.

                 
        
        # here I skip isomers...
        global amass_int
        amass_int=np.zeros(len(spe_rude2)) 
        for i in range(len(spe_rude2)-k):   
            amass_int[i]=int(spe_rude2[i])
            #print amass_int


        # here I have to create an array for the atomic numbers.
        # I need to this when I calculate and plot element abundances 
        
        global znum_int
        znum_int=np.zeros(len(spe)) 
        
        for i in range(len(spe)):
            if str(spe[i][0:2]) == 'H ':
                znum_int[i] = 1        
            elif str(spe[i][0:2]) == 'He':
                znum_int[i] = 2
            elif str(spe[i][0:2]) == 'Li':
                znum_int[i] = 3
            elif str(spe[i][0:2]) == 'Be':
                znum_int[i] = 4
            elif str(spe[i][0:2]) == 'B ':
                znum_int[i] = 5
            elif str(spe[i][0:2]) == 'C ':
                znum_int[i] = 6
            elif str(spe[i][0:2]) == 'N ':
                znum_int[i] = 7
            elif str(spe[i][0:2]) == 'O ':
                znum_int[i] = 8
            elif str(spe[i][0:2]) == 'F ':
                znum_int[i] = 9
            elif str(spe[i][0:2]) == 'Ne':
                znum_int[i] = 10
            elif str(spe[i][0:2]) == 'Na':
                znum_int[i] = 11
            elif str(spe[i][0:2]) == 'Mg':
                znum_int[i] = 12
            elif str(spe[i][0:2]) == 'Al':
                znum_int[i] = 13
            elif str(spe[i][0:2]) == 'Si':
                znum_int[i] = 14
            elif str(spe[i][0:2]) == 'P ':
                znum_int[i] = 15
            elif str(spe[i][0:2]) == 'S ':
                znum_int[i] = 16
            elif str(spe[i][0:2]) == 'Cl':
                znum_int[i] = 17
            elif str(spe[i][0:2]) == 'Ar':
                znum_int[i] = 18
            elif str(spe[i][0:2]) == 'K ':
                znum_int[i] = 19
            elif str(spe[i][0:2]) == 'Ca':
                znum_int[i] = 20
            elif str(spe[i][0:2]) == 'Sc':
                znum_int[i] = 21
            elif str(spe[i][0:2]) == 'Ti':
                znum_int[i] = 22
            elif str(spe[i][0:2]) == 'V ':
                znum_int[i] = 23
            elif str(spe[i][0:2]) == 'Cr':
                znum_int[i] = 24
            elif str(spe[i][0:2]) == 'Mn':
                znum_int[i] = 25
            elif str(spe[i][0:2]) == 'Fe':
                znum_int[i] = 26
            elif str(spe[i][0:2]) == 'Co':
                znum_int[i] = 27
            elif str(spe[i][0:2]) == 'Ni':
                znum_int[i] = 28
            elif str(spe[i][0:2]) == 'Cu':
                znum_int[i] = 29
            elif str(spe[i][0:2]) == 'Zn':
                znum_int[i] = 30
            elif str(spe[i][0:2]) == 'Ga':
                znum_int[i] = 31
            elif str(spe[i][0:2]) == 'Ge':
                znum_int[i] = 32
            elif str(spe[i][0:2]) == 'As':
                znum_int[i] = 33
            elif str(spe[i][0:2]) == 'Se':
                znum_int[i] = 34
            elif str(spe[i][0:2]) == 'Br':
                znum_int[i] = 35
            elif str(spe[i][0:2]) == 'Kr':
                znum_int[i] = 36
            elif str(spe[i][0:2]) == 'Rb':
                znum_int[i] = 37
            elif str(spe[i][0:2]) == 'Sr':
                znum_int[i] = 38
            elif str(spe[i][0:2]) == 'Y ':
                znum_int[i] = 39
            elif str(spe[i][0:2]) == 'Zr':
                znum_int[i] = 40
            elif str(spe[i][0:2]) == 'Nb':
                znum_int[i] = 41
            elif str(spe[i][0:2]) == 'Mo':
                znum_int[i] = 42
            elif str(spe[i][0:2]) == 'Tc':
                znum_int[i] = 43
            elif str(spe[i][0:2]) == 'Ru':
                znum_int[i] = 44
            elif str(spe[i][0:2]) == 'Rh':
                znum_int[i] = 45
            elif str(spe[i][0:2]) == 'Pd':
                znum_int[i] = 46
            elif str(spe[i][0:2]) == 'Ag':
                znum_int[i] = 47
            elif str(spe[i][0:2]) == 'Cd':
                znum_int[i] = 48
            elif str(spe[i][0:2]) == 'In':
                znum_int[i] = 49
            elif str(spe[i][0:2]) == 'Sn':
                znum_int[i] = 50
            elif str(spe[i][0:2]) == 'Sb':
                znum_int[i] = 51
            elif str(spe[i][0:2]) == 'Te':
                znum_int[i] = 52
            elif str(spe[i][0:2]) == 'I ':
                znum_int[i] = 53
            elif str(spe[i][0:2]) == 'Xe':
                znum_int[i] = 54
            elif str(spe[i][0:2]) == 'Cs':
                znum_int[i] = 55
            elif str(spe[i][0:2]) == 'Ba':
                znum_int[i] = 56
            elif str(spe[i][0:2]) == 'La':
                znum_int[i] = 57
            elif str(spe[i][0:2]) == 'Ce':
                znum_int[i] = 58
            elif str(spe[i][0:2]) == 'Pr':
                znum_int[i] = 59
            elif str(spe[i][0:2]) == 'Nd':
                znum_int[i] = 60
            elif str(spe[i][0:2]) == 'Pm':
                znum_int[i] = 61
            elif str(spe[i][0:2]) == 'Sm':
                znum_int[i] = 62
            elif str(spe[i][0:2]) == 'Eu':
                znum_int[i] = 63
            elif str(spe[i][0:2]) == 'Gd':
                znum_int[i] = 64
            elif str(spe[i][0:2]) == 'Tb':
                znum_int[i] = 65
            elif str(spe[i][0:2]) == 'Dy':
                znum_int[i] = 66
            elif str(spe[i][0:2]) == 'Ho':
                znum_int[i] = 67
            elif str(spe[i][0:2]) == 'Er':
                znum_int[i] = 68
            elif str(spe[i][0:2]) == 'Tm':
                znum_int[i] = 69
            elif str(spe[i][0:2]) == 'Yb':
                znum_int[i] = 70
            elif str(spe[i][0:2]) == 'Lu':
                znum_int[i] = 71
            elif str(spe[i][0:2]) == 'Hf':
                znum_int[i] = 72
            elif str(spe[i][0:2]) == 'Ta':
                znum_int[i] = 73
            elif str(spe[i][0:2]) == 'W ':
                znum_int[i] = 74
            elif str(spe[i][0:2]) == 'Re':
                znum_int[i] = 75
            elif str(spe[i][0:2]) == 'Os':
                znum_int[i] = 76
            elif str(spe[i][0:2]) == 'Ir':
                znum_int[i] = 77
            elif str(spe[i][0:2]) == 'Pt':
                znum_int[i] = 78
            elif str(spe[i][0:2]) == 'Au':
                znum_int[i] = 79
            elif str(spe[i][0:2]) == 'Hg':
                znum_int[i] = 80
            elif str(spe[i][0:2]) == 'Tl':
                znum_int[i] = 81
            elif str(spe[i][0:2]) == 'Pb':
                znum_int[i] = 82
            elif str(spe[i][0:2]) == 'Bi':
                znum_int[i] = 83
            elif str(spe[i][0:2]) == 'Po':
                znum_int[i] = 84
            elif str(spe[i][0:2]) == 'At':
                znum_int[i] = 85

	if spe[0] == 'N   1':
		znum_int[0] = 0

        # here the index to connect name and atomic numbers.
        global index_atomic_number
        index_atomic_number = {}    
        for a,b in zip(spe,znum_int):
            index_atomic_number[a]=b


def   define_zip_index_for_species(names_ppn_world,number_names_ppn_world):
    	''' This just give back cl, that is the original index as it is read from files from a data file. 	 '''

        #connect the specie number in the list, with the specie name
        global cl
        cl={}
        for a,b in zip(names_ppn_world,number_names_ppn_world):
            cl[a] = b  



def element_abund_marco(i_decay,stable_isotope_list,stable_isotope_identifier,mass_fractions_array_not_decayed,mass_fractions_array_decayed):
    ''' Given an array of isotopic abundances not decayed and a similar array of isotopic abundances not decayed, 
    here elements abundances, and production factors for elements are calculated'''


    # this way is done in a really simple way. May be done better for sure, in a couple of loops.
    # I keep this, since I have only to copy over old script. Falk will probably redo it.
    
    import numpy as np
    #import utils as u

    global elem_abund
    elem_abund = np.zeros(z_bismuth)
    global elem_abund_decayed
    elem_abund_decayed = np.zeros(z_bismuth)
    global elem_prod_fac
    elem_prod_fac = np.zeros(z_bismuth)
    global elem_prod_fac_decayed
    elem_prod_fac_decayed = np.zeros(z_bismuth)
    
        
    # notice that elem_abund include all contribution, both from stables and unstables in
    # that moment.
    for i in range(z_bismuth):
        dummy = 0.
        for j in range(len(spe)):
            if znum_int[j] == i+1 and stable_isotope_identifier[j] > 0.5:
                dummy = dummy + float(mass_fractions_array_not_decayed[j])
    	elem_abund[i] = dummy


    for i in range(z_bismuth):
        if index_stable[i] == 1:
            elem_prod_fac[i] = float(elem_abund[i]/solar_elem_abund[i])
        elif index_stable[i] == 0:
            elem_prod_fac[i] = 0.    


    if i_decay == 2:
        for i in range(z_bismuth):
            dummy = 0.
            for j in range(len(mass_fractions_array_decayed)):
                if znum_int[cl[stable_isotope_list[j].capitalize()]] == i+1:
                    #print znum_int[cl[stable[j].capitalize()]],cl[stable[j].capitalize()],stable[j]
                    dummy = dummy + float(mass_fractions_array_decayed[j])
       	    elem_abund_decayed[i] = dummy


        for i in range(z_bismuth):
            if index_stable[i] == 1:
                elem_prod_fac_decayed[i] = float(elem_abund_decayed[i]/solar_elem_abund[i])
            elif index_stable[i] == 0:
                elem_prod_fac_decayed[i] = 0.    

