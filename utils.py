'''
utils.py

Utility class for holding extra methods from mesa.py, nuh5p.py


'''

import numpy as np
import scipy as sc
import ascii_table as att
from scipy import optimize
import matplotlib.pyplot as pl

class data_fitting():
	'''
	Wrapper for the scipy method optimize.leastsq

	Typically you have some data y(x) and you want to
	fit some function to this data.

	in order to test this create some sample data y(x)
	x=arange(0,100)
	y=(((random_sample(100)+0.5)*50.)-50)+x # this is the data to fit

	This class provides two simple fit functions:
	linear and power law
	Feel free to add more, or provide your custom function
	as an argument. For example, in order to provide an exponential
	fit function, first define the function:
	def ff(coef,x):
	    return coef[0]*sc.exp(coef[1]*x)+coef[2]
	and then call the data_fitting instance:
	f=utils.data_fitting(ff,coef=(1,0,0))

        Once you have initialized this class, the instance provides two methods:
        fit and plot to check the fit:
        f.fit(x,y)
	f.plot()
	The fit coefficients are stored in self.fcoef.

	(FH)
	'''

	def __init__(self,func='linear',coef=(1,1)):
		'''
		func       'linear', 'powerlaw' or a custom function
		coef       a guess for the list of coeffiecients, for
		           'powerlaw' coef must have three enries, if
			   you provide your own function, provide as many
			   coef entries as your function needs
		'''		
		if func is 'linear':
			print "Information: 'linear' fit needs coef list with 2 entries"
			print " -> will use default: coef = "+str(coef)
			if len(coef) is not 2:
				print "Warning: you want a linear fit but you have not"
				print "         provided a guess for coef with the"
				print "         right length (2)."
				print " -> I will continue and assume coef=(1,1)"
				coef = (1,1)
			def ff(coef,x):
				return coef[0]*x + coef[1]
		        self.func_name = func
		elif func is 'powerlaw':
			print "Information: 'powerlaw' fit needs coef list with 3 entries"
			print " -> will use default: coef = "+str(coef)
			if len(coef) is not 3:
				print "Warning: you want a power law fit but you have"
				print "         not provided a guess for coef with the"
				print "         right length (3)."
				print " -> I will continue and assume coef=(1,1,1)"
				coef = (1,1,1)
			def ff(coef,x):
				return coef[0]*x**coef[1] + coef[2]
		        self.func_name = func
		else:
			print "Information: You provide a fit function yourself. I trust"
			print "             you have provided a matching guess for the "
			print "             coefficient list!"
			ff = func
			self.func_name = func.__name__


		# we want to determine the coefficients that
		# fit the power law to the data
		# this is done by finding the minimum to a
		# residual function:
		# func(params) = ydata - f(xdata, params)
		# therefore we define a residual function
		def fres(coef,y,ff,x):
			return y-ff(coef,x)

		self.residual = fres
		self.coef     = coef
		self.func     = ff

	def fit(self,x,y,dcoef='none'):
		'''
		performs the fit

		x,y     matching data arrays that define a numerical function
		        y(x), this is the data to be fitted
		dcoef   optionally you can provide a different guess for the
		        coefficients
			
		output:
		self.fcoef   contains the fitted coefficients

		returns:
		ierr         values between 1 and 4 signal success
		'''
                self.x = x
                self.y = y

		if dcoef is not 'none':
			coef = dcoef
		else:
			coef = self.coef
			
		fcoef=optimize.leastsq(self.residual,coef,args=(y,self.func,x))
		self.fcoef = fcoef[0].tolist()
		return fcoef[1]

	def plot(self,ifig=1):
		'''
		plot the data and the fitted function

		ifig  figure window number
		'''

		if len(self.coef) is not len(self.fcoef):
			print "Warning: the fitted coefficient list is not same"
			print "         length as guessed list - still I will try ..."

		pl.figure(ifig)
		pl.plot(self.x,self.y, label='data')
		pl.plot(self.x,self.func(self.fcoef,self.x),label=self.func_name)
		pl.legend()
		
class constants():
	mass_sun=1.9891e+33
	mass_sun_unit='g'
	one_year=31558149.984
	avogadro=6.02214179e23
	avogadro_unit='mol^-1'

class Utils():
	'''This private class contains utilities that are used by methods,
	mostly in the ppn and mppnp classes. Users what normally not use these
	methods directly. Things go here when it can be imagined that they may
	be used not in immediate conjunction with plotting. Otherwise they
	would go into the superclass data_plot.'''

	#elements_names is marked for deletion (FH, Oct2011) and
	#should be replaced with self.stable_names from self._stable_names()
	elements_names = ['Neutron','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
	    'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
	    'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',
	    'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te', 'I','Xe','Cs','Ba','La','Ce',
	    'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',
	    'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu']

    	stable_el = [['Neutron','999'],['H',1, 2],['He', 3, 4],['Li', 6, 7],['Be', 9],
		     ['B', 10, 11],['C', 12, 13],['N', 14, 15],['O', 16, 17, 18],['F', 19],
		     ['Ne', 20, 21, 22],['Na', 23],['Mg', 24, 25, 26],['Al', 27],['Si', 28, 29, 30],
		     ['P', 31],['S', 32, 33, 34, 36],['Cl', 35, 37],['Ar', 36, 38, 40],['K', 39, 40, 41],
		     ['Ca', 40, 42, 43, 44, 46, 48],['Sc', 45],['Ti', 46, 47, 48, 49, 50],['V', 50, 51],
		     ['Cr', 50, 52, 53, 54],['Mn', 55],['Fe', 54, 56, 57, 58],['Co', 59],
		     ['Ni', 58, 60, 61, 62, 64],['Cu', 63, 65],['Zn', 64, 66, 67, 68, 70],['Ga', 69, 71],
		     ['Ge', 70, 72, 73, 74, 76],['As', 75],['Se', 74, 76, 77, 78, 80, 82],['Br', 79, 81],
		     ['Kr', 78, 80, 82, 83, 84, 86],['Rb', 85, 87],['Sr', 84, 86, 87, 88],['Y', 89],
		     ['Zr', 90, 91, 92, 94, 96],['Nb', 93],['Mo', 92, 94, 95, 96, 97, 98, 100],
		     ['Tc',999],['Ru', 96, 98, 99, 100, 101, 102, 104],['Rh', 103],
		     ['Pd', 102, 104, 105, 106, 108, 110],['Ag', 107, 109],
		     ['Cd', 106, 108, 110, 111, 112, 113, 114, 116],['In', 113, 115],
		     ['Sn', 112, 114, 115, 116, 117, 118, 119, 120, 122, 124],['Sb', 121, 123],
		     ['Te', 120, 122, 123, 124, 125, 126, 128, 130],['I', 127],
		     ['Xe', 124, 126, 128, 129, 130, 131, 132, 134, 136],['Cs', 133],
		     ['Ba', 130, 132, 134, 135, 136, 137, 138],['La', 138, 139],['Ce', 136, 138, 140, 142],
		     ['Pr', 141],['Nd', 142, 143, 144, 145, 146, 148, 150],['Pm',999],
		     ['Sm', 144, 147, 148, 149, 150, 152, 154],['Eu', 151, 153],
		     ['Gd', 152, 154, 155, 156, 157, 158, 160],['Tb', 159],
		     ['Dy', 156, 158, 160, 161, 162, 163, 164],['Ho', 165],
		     ['Er', 162, 164, 166, 167, 168, 170],['Tm', 169],['Yb', 168, 170, 171, 172, 173, 174, 176],
		     ['Lu', 175, 176],['Hf', 174, 176, 177, 178, 179, 180],['Ta', 180, 181],
		     ['W', 180, 182, 183, 184, 186],['Re', 185, 187],['Os', 184, 186, 187, 188, 189, 190, 192],
		     ['Ir', 191, 193],['Pt', 190, 192, 194, 195, 196, 198],['Au', 197],
		     ['Hg', 196, 198, 199, 200, 201, 202, 204],['Tl', 203, 205],['Pb', 204, 206, 207, 208],
		     ['Bi', 209],['Th', 232],['U',235,238]] 

	def _stable_names(self):
		'''This private method extracts the element names from stable_el. Note
		that stable_names is a misnomer as stable_el also contains unstable
		element names with a number 999 for the *stable* mass numbers. (FH
		wonders who would come up up with this name ...)'''
		stable_names=[]
		for i in range(len(self.stable_el)):
			stable_names.append(self.stable_el[i][0])
		self.stable_names=stable_names

	def _process_abundance_vector(self,a,z,isomers,yps):
		'''This private method takes as input one vector definition and
		processes it, including sorting by charge number and
		mass number. It returns the processed input variables
		plus an element and isotope vector and a list of
		isomers.'''
		
		tmp=[]
		isom=[]
		for i in xrange(len(a)):
			if z[i]!=0 and isomers[i]==1: #if its not 'NEUt and not an isomer'
				tmp.append([self.stable_names[int(z[i])]+'-'+str(int(a[i])),yps[i],z[i],a[i]])
			elif isomers[i]!=1: #if it is an isomer
				if yps[i]==0:
					isom.append([self.stable_names[int(z[i])]+'-'+str(int(a[i]))+'-'+str(int(isomers[i]-1)),1e-99])
				else:
					isom.append([self.stable_names[int(z[i])]+'-'+str(int(a[i]))+'-'+str(int(isomers[i]-1)),yps[i]])	
		tmp.sort(self.compar)
		tmp.sort(self.comparator)
		abunds=[]
		isotope_to_plot=[]
		z_iso_to_plot=[]
		a_iso_to_plot=[]
		el_iso_to_plot=[]
		for i in xrange(len(tmp)):
			isotope_to_plot.append(tmp[i][0])
			abunds.append(tmp[i][1])
			z_iso_to_plot.append(int(tmp[i][2]))
			a_iso_to_plot.append(int(tmp[i][3]))
			el_iso_to_plot.append(self.stable_names[int(tmp[i][2])])

		return a_iso_to_plot,z_iso_to_plot,abunds,isotope_to_plot,el_iso_to_plot,isom

	def compar(self,x, y):
		'''
		simple comparator method
		'''
		
		indX=0
		indY=0
		
		a= int(x[0].split('-')[1])
		
		b= int(y[0].split('-')[1])
		

		if a>b:
			return 1
		if a==b:
			return 0
		if a<b:
			return -1
	
	def comparator(self,x, y):
		'''
		simple comparator method
		'''
		
		indX=0
		indY=0
		for i in xrange(len(self.stable_names)):
			if self.stable_names[i] == x[0].split('-')[0]:
				indX=i
			if self.stable_names[i] == y[0].split('-')[0]:
				indY=i

		if indX>indY:
			return 1
		if indX==indY:
			return 0
		if indX<indY:
			return -1

	def _read_isotopedatabase(self,ffname='isotopedatabase.txt'):
		'''This private method reads the isotopedatabase.txt file in sldir run
		dictory and returns z, a, elements, the cutoff mass
		for each species that delineate beta+ and beta- decay
		and the logical in the last column. Also provides
		charge_from_element dictionary according to
		isotopedatabase.txt.''' 
		name=self.sldir+ffname
		z_db, a_db, el_db, stable_a_db,logic_db=\
		    np.loadtxt(name,unpack=True,dtype='str') 
		z_db=np.array(z_db,dtype='int')
		a_db=np.array(a_db,dtype='int')
		stable_a_db=np.array(stable_a_db,dtype='int')

		# charge number for element name from dictionary in isotopedatabase.txt
		charge_from_element_name={}
		for name in self.stable_names:
			if name=='Neutron' or name=='Neut' or name=='NEUT' or name=='N-1':
				name='nn'
			try:
				zz=z_db[np.where(el_db==name)][0]
				charge_from_element_name[name]=zz
			except IndexError:
				print name+" does not exist in this run"
		return z_db, a_db, el_db, stable_a_db,logic_db,charge_from_element_name
	
	def decay_indexpointer(self):
		'''This private method provides decay indexpointers which allow to
		instantaneously decay an abundance vector an abundance
		vector instantaneously. These are:

		self.decay_idp                        points in the iso_to_plot (i.e. the 
		                                      undecayed abundance vector 
		                                      index space) to the decay target
		self.idp_to_stables_in_isostoplot     points to the stable isotopes in the 
		                                      undecayed abundance vector index space
		For an application example see ppn.py-abu_vector-_getcycle.
		'''
		a_iso_to_plot   =self.a_iso_to_plot
		isotope_to_plot =self.isotope_to_plot
		z_iso_to_plot   =self.z_iso_to_plot
		el_iso_to_plot  =self.el_iso_to_plot
		abunds          =self.abunds
		isom            =self.isom

		z_db, a_db, el_db, stable_a_db,logic_db,charge_from_element_name=\
		    self._read_isotopedatabase()
		# find out which species  beta+ and which beta- decay: 
		beta=np.sign(stable_a_db-a_db) # if a species is unstable and if beta < 0 => beta- decay
		                            # else beta > 0 => beta+ decay

		# now we need an index array on the scale of the abundance
		# distribution to be plotted that points to itself for stable species,
		# and to the stable element to which it decays in case of an unstable
		# species
		decay_index_pointer=np.zeros(len(isotope_to_plot), dtype='int')-1
		idp_to_stables_in_isostoplot=[]
		for i in range(len(isotope_to_plot)):
		    element_name=isotope_to_plot[i].split('-')[0] 
		    try:
			stable_a=stable_a_db[np.where(el_db==element_name)][0] # 4th column for that element in isotopedatabase.txt
		    except IndexError:
			print "Can't find element "+element_name+" in isotopedatabase.txt"
		    if a_iso_to_plot[i] <= 209 and stable_a <=209:  # Bi209 is last stable element        
			stable_mass_numbers=self.stable_el[self.stable_names.index(element_name)][1:]
			iso_db_index_range_el=np.where(el_db==element_name)
			beta_for_this_species=\
			    beta[iso_db_index_range_el][np.where(a_db[iso_db_index_range_el]==a_iso_to_plot[i])]
			if beta_for_this_species == 0:  # if there are no stable species for an element (Tc,Pm) the cutoff specifies 
			    beta_for_this_species = -1  # the lowest mass beta- isotope
			if a_iso_to_plot[i] in stable_mass_numbers:
			    # print isotope_to_plot[i]+" is stable"
			    decay_index_pointer[i]=i
			    idp_to_stables_in_isostoplot.append(i)
			elif  a_iso_to_plot[i]==8: # Be8 -> He4
			   decay_index_pointer[i]=isotope_to_plot.index('He-4')
			else: # beta decay
			    found_decay_target=False
			    i_search=-1*beta_for_this_species
			    while not found_decay_target:
				try:
				    try_target_el=self.stable_names[charge_from_element_name[element_name]+i_search]
				except TypeError:
				    print "Maybe information about species "+isotope_to_plot[i]+" is not available in isotopedatabase.txt"
				    decay_index_pointer[i]=-1
				    break
				# print try_target_el
				try:
				    stable_mass_numbers=self.stable_el[self.stable_names.index(try_target_el)][1:]
				except ValueError:
				    print "Can not find decay target for "+isotope_to_plot[i]
				if a_iso_to_plot[i] in stable_mass_numbers:
				    ind_range=np.where(np.array(el_iso_to_plot)==try_target_el)[0]
				    if a_iso_to_plot[i] in np.array(a_iso_to_plot)[ind_range]:
					this_ind=\
					    ind_range[np.where(np.array(a_iso_to_plot)[ind_range]==a_iso_to_plot[i])[0]]
					# print isotope_to_plot[i]+" is unstable and decays to "+isotope_to_plot[this_ind]
					decay_index_pointer[i]=this_ind
				    else:
					print "It seems unstable species "+isotope_to_plot[i]+" wants to decay to " \
					    +try_target_el+"-"+str(a_iso_to_plot[i])+", however this species is not in this run." \
					    +" This points to an inconsistency in the network build. Here we will ignore the abundance of " \
					    +isotope_to_plot[i]+'.'
					decay_index_pointer[i]=-1                            
				    found_decay_target=True
				else:
				    i_search += -1*beta_for_this_species
		if self.debug:
			print "Decay rules:"
			for i in range(len(isotope_to_plot)):
				if decay_index_pointer[i]>= 0:
					print isotope_to_plot[i]+" -> "+isotope_to_plot[decay_index_pointer[i]]					
		ind_tmp=idp_to_stables_in_isostoplot
                #ind_tmp=utils.strictly_monotonic(decay_index_pointer)  # this would do the same, but the method above is more straight forward

		self.decay_idp=decay_index_pointer
		self.idp_to_stables_in_isostoplot=ind_tmp

class iniabu(Utils):
	'''
	This class in the utils package reads an abundance
	distribution file of the type iniab.dat. It then provides you
	with methods to change some abundances, modify, normalise and
	eventually write out the final distribution in a format that
	can be used as an initial abundance file for ppn. This class
	also contains a method to write initial abundance files for a
	MESA run, for a given MESA netowrk.
	'''
	# clean variables that we will use in this class
	
	filename = ''
    
	def __init__(self,filename):
		'''
		Init method will read file of type iniab.dat, as they are for
		example found in the frames/mppnp/USEPP directory.

		An instance of this class will have the following data arrays
		z      charge number
		a      mass number
		abu    abundance
		names  name of species
		habu   a hash array of abundances, referenced by species name
		hindex hash index returning index of species from name

		E.g. if x is an instance then x.names[4] gives you the
		name of species 4, and x.habu['c 12'] gives you the
		abundance of C12, and x.hindex['c 12'] returns
		4. Note, that you have to use the species names as
		they are provided in the iniabu.dat file.

		Example - generate modified input file ppn calculations: 

		import utils
		p=utils.iniabu('iniab1.0E-02.ppn_asplund05')
		sp={}
		sp['h   1']=0.2
		sp['c  12']=0.5
		sp['o  16']=0.2
		p.set_and_normalize(sp)
		p.write('p_ini.dat','header for this example')

		p.write_mesa allows you to write this NuGrid initial abundance
		file into a MESA readable initial abundance file.
		'''
		f0=open(filename)
		sol=f0.readlines()
		f0.close 

		# Now read in the whole file and create a hashed array:
		names=[]
		z=[]    
                yps=np.zeros(len(sol))
                mass_number=np.zeros(len(sol))
		for i in range(len(sol)):
		    z.append(int(sol[i][1:3]))
		    names.extend([sol[i].split("         ")[0][4:]])
		    yps[i]=float(sol[i].split("         ")[1])
		    try:
			mass_number[i]=int(names[i][2:5])
		    except ValueError:
			print "WARNING:"
			print "This initial abundance file uses an element name that does"
			print "not contain the mass number in the 3rd to 5th position."
			print "It is assumed that this is the proton and we will change"
			print "the name to 'h   1' to be consistent with the notation used"
			print "in iniab.dat files"
			names[i]='h   1'
		    mass_number[i]=int(names[i][2:5])
		# now zip them together:
                hash_abu={}
                hash_index={}
		for a,b in zip(names,yps):
		   hash_abu[a] = b

		for i in range(len(names)):
		   hash_index[names[i]] = i

		self.z=z
		self.abu=yps
		self.a=mass_number
		self.names=names
		self.habu=hash_abu
		self.hindex=hash_index

	def write(self,outfile='initial_abundance.dat',header_string='initial abundances for a PPN run'):
		'''
		Write initial abundance file (intended for use with ppn)

                outfile          name of output file
		header_string    srting with header line
		'''
		dcols=['Z', 'species','mass fraction']
		data=[self.z,self.names,self.abu]
		hd=[header_string]
		att.write(outfile,hd,dcols,data)

	def write_mesa(self,mesa_isos_file='isos.txt',add_excess_iso='fe56',outfile='xa_iniabu.dat',header_string='initial abundances for a MESA run',header_char='!'):
		'''
		Write initial abundance file

		mesa_isos_file   list with isos copied from mesa network 
                                 definition file in mesa/data/net_data/nets
	        add_excess_iso   add 1.-sum(isos in mesa net) to this isotope
                outfile          name of output file
		header_string    srting with header line

		Useage example:
		import utils
		!ls ~/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/   # find ppn initial abundance file
		!cat ~/mesa/data/net_data/nets/agb.net                   # find isos needed in mesa net
		!cat > isos.txt                                          # paste needed isos into file             
		help(utils.iniabu)                                       # check documentation of method
		x=utils.iniabu('path_to_here/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn')
		x.write_mesa?
		x.write_mesa(add_excess_iso='ne22',header_string='mppnp/USEEPP/iniab2.0E-02GN93.ppn for mesa/agb.net',outfile='xa_2.0E-02GN93.mesa')
		'''


		f=open('isos.txt')
		a=f.readlines()
		isos=[]
		for i in range(len(a)):
			isos.append(a[i].strip().rstrip(','))

		mesa_names=[]
		abus=[]
		for i in range(len(self.z)):
			b=self.names[i].split()
			a=''
			a=a.join(b)
			if a in isos:
				mesa_names.append(a)
				abus.append(self.abu[i])
			# mesa_names.append(elements_names[int(x.z[i])].lower()+str(int(x.a[i])))
				
		for i in range(len(isos)):
			if isos[i] not in mesa_names:
				mesa_names.append(isos[i])
				abus.append(0.0)

		excess=1.-np.sum(np.array(abus))
		abus=np.array(abus)
		abus[mesa_names.index(add_excess_iso)] += excess

		dcols=['','']
		data=[mesa_names,abus]
		hd=[header_string]
		att.write(outfile,hd,dcols,data,header_char=header_char)

	def set_and_normalize(self,species_hash):
		'''
		species_hash is a hash array in which you provide
		abundances referenced by species names that you want
		to set to some particular value; all other species are
		then normalised so that that the total sum is 1

		You can set up the argument array for this method for
		example in the following way:		
		In [117]: sp={}
		In [118]: sp['he  4']=0.2
		In [119]: sp['h   1']=0.5
		'''
		sum_before = 1.
		for i in range(len(species_hash)):
			sum_before -=  self.abu[self.hindex[species_hash.keys()[i]]]
		print "sum_before = "+str(sum_before)
		normalization_factor=(1.-sum(species_hash.values()))/sum_before
		print "normalizing the rest witih factor "+str(normalization_factor)
		self.abu *= normalization_factor
		for i in range(len(species_hash)):
			self.abu[self.hindex[species_hash.keys()[i]]]=species_hash.values()[i]
		for name in self.habu:
			self.habu[name]=self.abu[self.hindex[name]]

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
		list_to_print = ['H-1','He-4','C-12','N-14','O-16','Ne-20']
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
	list1, list2, lines1, lines2
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
		symbol=['b--','k--','r--','c--','m--','g--','b-','k-','r-','c-','m-','g-','b.','b-.','k-.','r-.','c-.','m-.','g-.','b:','k:','r:','c:','m:','g:']
	elif what_list is "lines2":
		symbol=['g:','r-.','k-','b--','k--','b.','r--','b-','c--','m--','g--','r-','c-','m-','g-','k-.','c-.','m-.','g-.','k:','r:','c:','m:','b-.','b:']
	return symbol
		
def make_list(default_symbol_list,len_list_to_print):
	''' provide the list of symbols to use according for the list of species/arrays to plot.
	default_symbol_list = list of symbols that the user choose to use.  	    	
	len_list_to_print   = len of list of species/arrays to print.
	'''
	
	symbol_used = []
	for i in range(len_list_to_print):
		symbol_used.append(default_symbol_list[sc.mod(i,len(default_symbol_list))])
			
	return symbol_used
		
def strictly_monotonic(bb):
    '''bb is an index array which may have numerous double or triple
    occurrences of indices, such as for example the
    decay_index_pointer. This method removes all entries <= -, then
    all dubliates and finally returns a sorted list of indices.
    '''
    cc=bb[np.where(bb>=0)]
    cc.sort()
    dc=cc[1:]-cc[:-1] # subsequent equal entries have 0 in db
    dc=np.insert(dc,0,1) # the first element is always unique (the second occurence is the dublicate)
    dc_mask=np.ma.masked_equal(dc,0)
    return np.ma.array(cc,mask=dc_mask.mask).compressed()

def solar(filename_solar,solar_factor):
    ''' read solar abundances from filename_solar. solar_factor is
    the correction factor to apply, in case filename_solar is not
    solar, but some file used to get initial abundances at metallicity
    lower than solar. However, notice that this is really rude, since
    alpha-enahncements and things like that are not properly
    considered.  Only H and He4 are not multiplied. So, for
    publications PLEASE use proper filename_solar at...solar, and use
    solar_factor = 1. Marco'''

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
        try:
		mass_number[i]=int(names_sol[i][2:5])
	except ValueError:
		print "WARNING:"
		print "This initial abundance file uses an element name that does"
		print "not contain the mass number in the 3rd to 5th position."
		print "It is assumed that this is the proton and we will change"
		print "the name to 'h   1' to be consistent with the notation used in"
		print "iniab.dat files"
		names_sol[i]='h   1'
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



    z_bismuth = 83
    global solar_elem_abund
    solar_elem_abund = np.zeros(z_bismuth)


    for i in range(z_bismuth):
        dummy = 0.
        for j in range(len(solar_abundance)):
            if z_sol[j] == i+1:
                dummy = dummy + float(solar_abundance[names_sol[j]])
    	solar_elem_abund[i] = dummy


def convert_specie_naming_from_h5_to_ppn(isotope_names):
	''' read isotopes names from h5 files, and convert them
    	according to standard scheme used inside ppn and mppnp.  Also
    	Z and A are recalculated, for these species. Isomers are
    	excluded for now, since there were recent changes in isomers
    	name. As soon as the isomers names are settled, than Z and A
    	provided here will be obsolete, and can be changed by usual Z
    	and A.  '''

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


def stable_specie():
    ''' provide the list of stable species, and decay path feeding stables '''


    import numpy as np


    stable_raw=[]
    stable_raw = ['H   1', 'H   2',\
    'HE  3', 'HE  4',\
    'LI  6', 'LI  7',\
    'BE  9',\
    'B  10', 'B  11',\
    'C  12', 'C  13',\
    'N  14', 'N  15',\
    'O  16', 'O  17', 'O  18',\
    'F  19',\
    'NE 20', 'NE 21', 'NE 22',\
    'NA 23',\
    'MG 24', 'MG 25', 'MG 26',\
    'AL 27',\
    'SI 28', 'SI 29', 'SI 30',\
    'P  31',\
    'S  32', 'S  33', 'S  34', 'S  36',\
    'CL 35', 'CL 37',\
    'AR 36', 'AR 38', 'AR 40',\
    'K  39', 'K  40', 'K  41',\
    'CA 40', 'CA 42', 'CA 43', 'CA 44', 'CA 46', 'CA 48',\
    'SC 45',\
    'TI 46', 'TI 47', 'TI 48', 'TI 49', 'TI 50',\
    'V  50', 'V  51',\
    'CR 50', 'CR 52', 'CR 53', 'CR 54',\
    'MN 55',\
    'FE 54', 'FE 56', 'FE 57', 'FE 58',\
    'CO 59',\
    'NI 58', 'NI 60', 'NI 61', 'NI 62', 'NI 64',\
    'CU 63', 'CU 65',\
    'ZN 64', 'ZN 66', 'ZN 67', 'ZN 68', 'ZN 70',\
    'GA 69', 'GA 71',\
    'GE 70', 'GE 72', 'GE 73', 'GE 74', 'GE 76',\
    'AS 75',\
    'SE 74', 'SE 76', 'SE 77', 'SE 78', 'SE 80', 'SE 82',\
    'BR 79', 'BR 81',\
    'KR 78', 'KR 80', 'KR 82', 'KR 83', 'KR 84', 'KR 86',\
    'RB 85', 'RB 87',\
    'SR 84', 'SR 86', 'SR 87', 'SR 88',\
    'Y  89',\
    'ZR 90', 'ZR 91', 'ZR 92', 'ZR 94', 'ZR 96',\
    'NB 93',\
    'MO 92', 'MO 94', 'MO 95', 'MO 96', 'MO 97', 'MO 98', 'MO100',\
    'RU 96', 'RU 98', 'RU 99', 'RU100', 'RU101', 'RU102', 'RU104',\
    'RH103',\
    'PD102', 'PD104', 'PD105', 'PD106', 'PD108', 'PD110',\
    'AG107', 'AG109',\
    'CD106', 'CD108', 'CD110', 'CD111', 'CD112', 'CD113', 'CD114', 'CD116',\
    'IN113',  'IN115',\
    'SN112', 'SN114', 'SN115', 'SN116', 'SN117', 'SN118', 'SN119', 'SN120', 'SN122', 'SN124',\
    'SB121', 'SB123',\
    'TE120', 'TE122', 'TE123', 'TE124', 'TE125', 'TE126', 'TE128', 'TE130',\
    'I 127',\
    'XE124', 'XE126', 'XE128', 'XE129', 'XE130', 'XE131', 'XE132', 'XE134', 'XE136',\
    'CS133',\
    'BA130', 'BA132', 'BA134', 'BA135', 'BA136', 'BA137', 'BA138',\
    'LA138', 'LA139',\
    'CE136', 'CE138', 'CE140', 'CE142',\
    'PR141',\
    'ND142', 'ND143', 'ND144', 'ND145', 'ND146', 'ND148', 'ND150',\
    'SM144', 'SM147', 'SM148', 'SM149', 'SM150', 'SM152', 'SM154',\
    'EU151', 'EU153',\
    'GD152', 'GD154', 'GD155', 'GD156', 'GD157', 'GD158', 'GD160',\
    'TB159',\
    'DY156', 'DY158', 'DY160', 'DY161', 'DY162', 'DY163', 'DY164',\
    'HO165',\
    'ER162', 'ER164', 'ER166', 'ER167', 'ER168', 'ER170',\
    'TM169',\
    'YB168', 'YB170', 'YB171', 'YB172', 'YB173', 'YB174', 'YB176',\
    'LU175', 'LU176',\
    'HF174', 'HF176', 'HF177', 'HF178', 'HF179', 'HF180',\
    'TA180', 'TA181',\
    'W 180', 'W 182', 'W 183', 'W 184', 'W 186',\
    'RE185', 'RE187',\
    'OS184', 'OS186', 'OS187', 'OS188', 'OS189', 'OS190', 'OS192',\
    'IR191', 'IR193',\
    'PT190', 'PT192', 'PT194', 'PT195', 'PT196', 'PT198',\
    'AU197',\
    'HG196', 'HG198', 'HG199', 'HG200', 'HG201', 'HG202', 'HG204',\
    'TL203', 'TL205',\
    'PB204', 'PB206', 'PB207', 'PB208',\
    'BI209']

    jj=-1
    global count_size_stable
    count_size_stable=[]
    global stable
    stable=[]
    global jdum
    jdum=np.zeros(len(stable_raw))
    global jjdum
    jjdum=np.zeros(len(spe))
    for i in range(len(stable_raw)):
       dum_str = stable_raw[i]
       for j in range(len(spe)):
               if stable_raw[i].capitalize() == spe[j]:
                   stable.append(stable_raw[i]) 
                   jdum[i]=1
		   jjdum[j]=1
                   jj=jj+1
                   count_size_stable.append(int(jj))
    #print stable
    # back_ind is an index to go back, to use the order of stable
    # useful for example for decayed yields.
    global back_ind
    back_ind={}     
    for a,b in zip(stable,count_size_stable):
         back_ind[a]=b       
    #print 'in stable:',back_ind['SE 74']   
    # definition of decay paths
    global decay_raw
    decay_raw=[]
    decay_raw=[['H   1'],\
    ['H   2'],\
    ['HE  3'],\
    ['HE  4','B   8'],\
    ['LI  6'],\
    ['LI  7','BE  7'],\
    ['BE  9'],\
    ['B  10','BE 10'],\
    ['B  11','C  11','BE 11'],\
    ['C  12'],\
    ['C  13','N  13','O  13'],\
    ['N  14','C  14','O  14'],\
    ['N  15','C  15','O  15','F  15'],\
    ['O  16'],\
    ['O  17','F  17'],\
    ['O  18','F  18','NE 18'],\
    ['F  19','O  19','NE 19'],\
    ['NE 20','F  20','NA 20'],\
    ['NE 21','F  21','NA 21'],\
    ['NE 22','NA 22','MG 22','F  22'],\
    ['NA 23','MG 23','NE 23'],\
    ['MG 24','AL 24','NA 24','NE 24',],\
    ['MG 25','NA 25','AL 25'],\
    ['MG 26','SI 26','AL 26','NA 26','AL*26'],\
    ['AL 27','SI 27','MG 27'],\
    ['SI 28','AL 28','MG 28'],\
    ['SI 29','P  29','AL 29','MG 29'],\
    ['SI 30','S  30','P  30','AL 30','MG 30'],\
    ['P  31','S  31','SI 31'],\
    ['S  32','P  32','SI 32'],\
    ['S  33','CL 33','P  33','SI 33'],\
    ['S  34','CL 34','P  34'],\
    ['S  36','P  36'],\
    ['CL 35','S  35','AR 35'],\
    ['CL 37','S  37','P  37','AR 37','K  37'],\
    ['AR 36','CL 36'],\
    ['AR 38','CL 38','S  38','P  38','K  38'],\
    ['AR 40','CL 40','S  40'],\
    ['K  39','AR 39','CL 39','S  39','CA 39'],\
    ['K  40'],\
    ['K  41','AR 41','CL 41','S  41','CA 41','SC 41'],\
    ['CA 40','SC 40'],\
    ['CA 42','K  42','AR 42','CL 42','S  42','SC 42','TI 42'],\
    ['CA 43','K  43','AR 43','CL 43','SC 43','TI 43','V  43'],\
    ['CA 44','K  44','AR 44','CL 44','SC 44','TI 44'],\
    ['CA 46','K  46','AR 46'],\
    ['CA 48','K  48','AR 48'],\
    ['SC 45','CA 45','K  45','AR 45','CL 45','TI 45','V  45'],\
    ['TI 46','SC 46','V  46','CR 46'],\
    ['TI 47','SC 47','CA 47','K  47','AR 47','V  47','CR 47'],\
    ['TI 48','SC 48','V  48','CR 48'],\
    ['TI 49','SC 49','CA 49','K  49','V  49','CR 49','MN 49'],\
    ['TI 50','SC 50','CA 50','K  50'],\
    ['V  50'],\
    ['V  51','CR 51','TI 51','SC 51','CA 51','MN 51'],\
    ['CR 50','MN 50'],\
    ['CR 52','MN 52','FE 52','V  52','TI 52','SC 52','CA 52'],\
    ['CR 53','MN 53','FE 53','V  53','TI 53','SC 53'],\
    ['CR 54','MN 54','V  54','TI 54','SC 54'],\
    ['MN 55','FE 55','CR 55','V  55','TI 55','CO 55'],\
    ['FE 54','CO 54'],\
    ['FE 56','NI 56','CO 56','MN 56','CR 56'],\
    ['FE 57','NI 57','CO 57','MN 57','CR 57'],\
    ['FE 58','CO 58','MN 58','CR 58'],\
    ['CO 59','FE 59','MN 59','CR 59','NI 59','CU 59'],\
    ['NI 58','CU 58'],\
    ['NI 60','CO 60','FE 60','MN 60','CR 60','CU 60','ZN 60'],\
    ['NI 61','CO 61','FE 61','MN 61','CU 61','ZN 61'],\
    ['NI 62','CO 62','FE 62','CU 62','ZN 62'],\
    ['NI 64','CO 64','FE 64','CU 64'],\
    ['CU 63','NI 63','CO 63','FE 63','MN 63','ZN 63','GA 63'],\
    ['CU 65','NI 65','CO 65','FE 65','ZN 65','GA 65','GE 65'],\
    ['ZN 64','CU 64','GA 64','GE 64'],\
    ['ZN 66','CU 66','NI 66','CO 66','FE 66','GA 66','GE 66'],\
    ['ZN 67','CU 67','NI 67','CO 67','FE 67','GA 67','GE 67','AS 77'],\
    ['ZN 68','NI 68','CO 68','GA 68','GE 68','CU 68','AS 68','SE 68'],\
    ['ZN 70','CU 70','NI 70','CO 70'],\
    ['GA 69','ZN 69','CU 69','NI 69','GE 69','AS 69','SE 69'],\
    ['GA 71','ZN 71','CU 71','NI 71','GE 71','AS 71','SE 71','BR 71'],\
    ['GE 70','GA 70','AS 70','SE 70','BR 70'],\
    ['GE 72','GA 72','ZN 72','CU 72','NI 72','AS 72','SE 72','BR 72','KR 72'],\
    ['GE 73','GA 73','ZN 73','CU 73','NI 73','AS 73','SE 73','BR 73','KR 73'],\
    ['GE 74','GA 74','ZN 74','CU 74','NI 74','AS 74'],\
    ['GE 76','GA 76','ZN 76','CU 76'],\
    ['AS 75','GE 75','GA 75','ZN 75','CU 75','SE 75','BR 75','KR 75','RB 75'],\
    ['SE 74','AS 74','BR 74','KR 74'],\
    ['SE 76','AS 76','BR 76','KR 76','RB 76','SR 76'],\
    ['SE 77','AS 77','GE 77','BR 77','GA 77','ZN 77','KR 77','RB 77','SR 77'],\
    ['SE 78','AS 78','GE 78','GA 78','ZN 78','BR 78'],\
    ['SE 80','AS 80','GE 80','GA 80','ZN 80'],\
    ['SE 82','AS 82','GE 82','GA 82'],\
    ['BR 79','SE 79','AS 79','GE 79','GA 79','ZN 79','KR 79','RB 79','SR 79','Y  79'],\
    ['BR 81','SE 81','KR 81','AS 81','GE 81','GA 81','RB 81','SR 81','Y  81','ZR 81'],\
    ['KR 78','RB 78','SR 78','Y  78'],\
    ['KR 80','BR 80','RB 80','SR 80','ZR 80'],\
    ['KR 82','BR 82','RB 82','SR 82','Y  82','ZR 82'],\
    ['KR 83','BR 83','SE 83','AS 83','GE 83','RB 83','SR 83','Y  83','ZR 83','NB 83'],\
    ['KR 84','BR 84','SE 84','AS 84','GE 84','RB 84'],\
    ['KR 86','BR 86','SE 86','AS 86'],\
    ['RB 85','KR 85','SR 85','KR*85','BR 85','SE 85','AS 85','Y  85','ZR 85','NB 85','MO 85'],\
    ['RB 87','KR 87','BR 87','SE 87','AS 87'],\
    ['SR 84','Y  84','ZR 84','NB 84','MO 84'],\
    ['SR 86','RB 86','Y  86','ZR 86','NB 86','MO 86'],\
    ['SR 87','Y  87','ZR 87','NB 87','MO 87','TC 87'],\
    ['SR 88','RB 88','KR 88','BR 88','SE 88','Y  88','ZR 88','NB 88','MO 88','TC 88'],\
    ['Y  89','SR 89','RB 89','KR 89','BR 89','ZR 89','NB 89','MO 89','TC 89','RU 89'],\
    ['ZR 90','Y  90','SR 90','NB 90'],\
    ['ZR 91','Y  91','SR 91','NB 91'],\
    ['ZR 92','Y  92','SR 92','NB 92'],\
    ['ZR 94','Y  94'],\
    ['ZR 96'],\
    ['NB 93','ZR 93','MO 93'],\
    ['MO 92','TC 92','RU 92'],\
    ['MO 94','NB 94','TC 94','RU 94'],\
    ['MO 95','NB 95','ZR 95','Y  95'],\
    ['MO 96','NB 96','TC 96'],\
    ['MO 97','NB 97','ZR 97','TC 97'],\
    ['MO 98','NB 98'],\
    ['MO100'],\
    ['RU 96','RH 96','PD 96'],\
    ['RU 98','TC 98','RH 98','PD 98'],\
    ['RU 99','TC 99','MO 99'],\
    ['RU100','TC100','RH100'],\
    ['RU101','TC101','RH101','MO101'],\
    ['RU102','MO102','TC102'],\
    ['RU104','TC104'],\
    ['RH103','RU103','TC103','PD103'],\
    ['PD102','AG102','CD102'],\
    ['PD104','RH104'],\
    ['PD105','RH105','RU105','AG105','AG105'],\
    ['PD106','RH106','RU106','AG106'],\
    ['PD108'],\
    ['PD110'],\
    ['AG107','PD107','CD107'],\
    ['AG109','PD109','CD109'],\
    ['CD106','IN106','SN106'],\
    ['CD108','AG108','IN108','SN108'],\
    ['CD110','AG110'],\
    ['CD111','AG111','PD111','IN111'],\
    ['CD112','AG112','PD112'],\
    ['CD113','AG113'],\
    ['CD114'],\
    ['CD116'],\
    ['IN113','SN113','SB113','TE113'],\
    ['IN115','CD115'],\
    ['SN112','IN112','SB112','TE112'],\
    ['SN114','IN114','SB114','TE114'],\
    ['SN115','SB115','TE115','I 115'],\
    ['SN116','IN116'],\
    ['SN117','IN117','SB117','CD117'],\
    ['SN118','SB118'],\
    ['SN119','SB119'],\
    ['SN120','SB120'],\
    ['SN122'],\
    ['SN124'],\
    ['SB121','SN121','TE121'],\
    ['SB123','SN123'],\
    ['TE120','I 120','XE120'],\
    ['TE122','SB122'],\
    ['TE123','I 123'],\
    ['TE124','SB124','I 124'],\
    ['TE125','SB125','I 125','SN125'],\
    ['TE126','SB126','SN126'],\
    ['TE128','SB128','SN128','I 128'],\
    ['TE130'],\
    ['I 127','TE127','XE127'],\
    ['XE124','CS124','BA124'],\
    ['XE126','CS126','BA126'],\
    ['XE128','I 128'],\
    ['XE129','I 129','TE129','SB129'],\
    ['XE130','I 130'],\
    ['XE131','I 131','TE131','CS131'],\
    ['XE132','I 132','TE132','CS132'],\
    ['XE134','I 134'],\
    ['XE136'],\
    ['CS133','XE133','I 133','BA133'],\
    ['BA130','LA130','CE130','PR130'],\
    ['BA132','LA132','CE132','PR132'],\
    ['BA134','CS134'],\
    ['BA135','CS135','XE135','I 135'],\
    ['BA136','CS136'],\
    ['BA137','CS137','XE137','LA137'],\
    ['BA138','CS138','XE138'],\
    ['LA138'],\
    ['LA139','BA139','CS139','CE139'],\
    ['CE136','PR136','ND136'],\
    ['CE138','PR138','ND138','PM138','SM138'],\
    ['CE140','LA140','BA140'],\
    ['CE142','LA142','BA142'],\
    ['PR141','CE141','LA141','BA141'],\
    ['ND142','PR142'],\
    ['ND143','PR143','CE143','LA143','PM143'],\
    ['ND144','PR144','CE144','LA144','PM144'],\
    ['ND145','PR145','CE145','PM145'],\
    ['ND146','PR146','CE146'],\
    ['ND148'],\
    ['ND150'],\
    ['SM144','EU144','GD144'],\
    ['SM147','PM147','ND147'],\
    ['SM148','PM148'],\
    ['SM149','PM149','ND149'],\
    ['SM150','PM150'],\
    ['SM152','PM152','ND152','EU152'],\
    ['SM154','PM154'],\
    ['EU151','SM151','PM151','ND151'],\
    ['EU153','SM153','PM153','GD153'],\
    ['GD152','EU152','TB152','DY152','HO152'],\
    ['GD154','EU154'],\
    ['GD155','EU155','SM155','TB155'],\
    ['GD156','EU156','SM156'],\
    ['GD157','EU157','SM157','TB157'],\
    ['GD158','EU158','SM158'],\
    ['GD160'],\
    ['TB159','GD159','EU159'],\
    ['DY156','HO156','ER156','TM156','YB156'],\
    ['DY158','HO158','ER158','TM158','YB158'],\
    ['DY160','TB160'],\
    ['DY161','TB161','GD161'],\
    ['DY162','TB162','GD162'],\
    ['DY163','TB163','HO163'],\
    ['DY164','TB164'],\
    ['HO165','DY165','ER165'],\
    ['ER162','TM162','YB162','LU162'],\
    ['ER164','TM164','YB164','LU164'],\
    ['ER166','HO166','DY166'],\
    ['ER167','HO167','DY167'],\
    ['ER168','HO168','DY168','TB161','TM168'],\
    ['ER170'],\
    ['TM169','ER169','HO169','YB169'],\
    ['YB168','LU168','HF168','TA168','W 168'],\
    ['YB170','TM170'],\
    ['YB171','TM171','ER171','LU171'],\
    ['YB172','TM172','ER172','LU172'],\
    ['YB173','TM173','ER173','LU173'],\
    ['YB174','TM174','ER174','LU174'],\
    ['YB176'],\
    ['LU175','YB175','TM175','HF175'],\
    ['LU176'],\
    ['HF174','TA174','W 174','RE174','OS174'],\
    ['HF176'],\
    ['HF177','LU177','YB177'],\
    ['HF178','LU178','YB178'],\
    ['HF179','LU179','TA179'],\
    ['HF180','LU180'],\
    ['TA180'],\
    ['TA181','HF181','W 181'],\
    ['W 180','RE180','OS180','IR180','PT180'],\
    ['W 182','TA182','HF182'],\
    ['W 183','TA183','HF183','RE183'],\
    ['W 184','TA184','HF184','RE184'],\
    ['W 186','TA186'],\
    ['RE185','W 185','TA185'],\
    ['RE187','W 187'],\
    ['OS184','IR184','PT184','AU184'],\
    ['OS186','RE186'],\
    ['OS187'],\
    ['OS188','RE188','W 188','IR188'],\
    ['OS189','RE189','W 189','IR189'],\
    ['OS190','RE190','W 190','IR190'],\
    ['OS192'],\
    ['IR191','OS191','RE191','PT191'],\
    ['IR193','OS193','PT193'],\
    ['PT190','AU190','HG190','TL190','PB190'],\
    ['PT192','IR192'],\
    ['PT194','IR194','OS194'],\
    ['PT195','IR195','OS195','AU195'],\
    ['PT196','IR196','AU196'],\
    ['PT198'],\
    ['AU197','PT197','HG197'],\
    ['HG196','TL196','PB196'],\
    ['HG198','AU198'],\
    ['HG199','AU199','PT199'],\
    ['HG200','AU200','PT200'],\
    ['HG201','AU201','PT201'],\
    ['HG202','AU202','PT202'],\
    ['HG204'],\
    ['TL203','HG203','PB203'],\
    ['TL205','HG205','PB205'],\
    ['PB204','TL204'],\
    ['PB206','TL206','HG206','PO210'],\
    ['PB207','TL207','HG207','BI207'],\
    ['PB208','TL208','HG208','BI208'],\
    ['BI209','PB209','TL209']]
    #print decay_raw    


def give_zip_element_z_and_names(element_name):
    ''' create 2 indexes that, given the name of the element/specie, give the atomic number.	 '''

    import numpy as np	

    global z_bismuth
    z_bismuth = 83
    global z_for_elem
    z_for_elem = []
    global index_stable
    index_stable = []

    i_for_stable = 1
    i_for_unstable = 0
    for i in range(z_bismuth):
        z_for_elem.append(int(i+1))
    	# the only elements below bismuth with no stable isotopes are Tc and Pm
    	if i+1 == 43 or i+1 == 61:
        	index_stable.append(i_for_unstable) 
    	else:
        	index_stable.append(i_for_stable)

        dummy_index = np.zeros(len(element_name))
        for i in range(len(element_name)):
            if element_name[i] == 'Neutron':
                dummy_index[i] = 0        
            elif element_name[i] == 'H':
                dummy_index[i] = 1        
            elif element_name[i] == 'He':
                dummy_index[i] = 2
            elif element_name[i] == 'Li':
                dummy_index[i] = 3
            elif element_name[i] == 'Be':
                dummy_index[i] = 4
            elif element_name[i] == 'B':
                dummy_index[i] = 5
            elif element_name[i] == 'C':
                dummy_index[i] = 6
            elif element_name[i] == 'N':
                dummy_index[i] = 7
            elif element_name[i] == 'O':
                dummy_index[i] = 8
            elif element_name[i] == 'F':
                dummy_index[i] = 9
            elif element_name[i] == 'Ne':
                dummy_index[i] = 10
            elif element_name[i] == 'Na':
                dummy_index[i] = 11
            elif element_name[i] == 'Mg':
                dummy_index[i] = 12
            elif element_name[i] == 'Al':
                dummy_index[i] = 13
            elif element_name[i] == 'Si':
                dummy_index[i] = 14
            elif element_name[i] == 'P':
                dummy_index[i] = 15
            elif element_name[i] == 'S':
                dummy_index[i] = 16
            elif element_name[i] == 'Cl':
                dummy_index[i] = 17
            elif element_name[i] == 'Ar':
                dummy_index[i] = 18
            elif element_name[i] == 'K':
                dummy_index[i] = 19
            elif element_name[i] == 'Ca':
                dummy_index[i] = 20
            elif element_name[i] == 'Sc':
                dummy_index[i] = 21
            elif element_name[i] == 'Ti':
                dummy_index[i] = 22
            elif element_name[i] == 'V':
                dummy_index[i] = 23
            elif element_name[i] == 'Cr':
                dummy_index[i] = 24
            elif element_name[i] == 'Mn':
                dummy_index[i] = 25
            elif element_name[i] == 'Fe':
                dummy_index[i] = 26
            elif element_name[i] == 'Co':
                dummy_index[i] = 27
            elif element_name[i] == 'Ni':
                dummy_index[i] = 28
            elif element_name[i] == 'Cu':
                dummy_index[i] = 29
            elif element_name[i] == 'Zn':
                dummy_index[i] = 30
            elif element_name[i] == 'Ga':
                dummy_index[i] = 31
            elif element_name[i] == 'Ge':
                dummy_index[i] = 32
            elif element_name[i] == 'As':
                dummy_index[i] = 33
            elif element_name[i] == 'Se':
                dummy_index[i] = 34
            elif element_name[i] == 'Br':
                dummy_index[i] = 35
            elif element_name[i] == 'Kr':
                dummy_index[i] = 36
            elif element_name[i] == 'Rb':
                dummy_index[i] = 37
            elif element_name[i] == 'Sr':
                dummy_index[i] = 38
            elif element_name[i] == 'Y':
                dummy_index[i] = 39
            elif element_name[i] == 'Zr':
                dummy_index[i] = 40
            elif element_name[i] == 'Nb':
                dummy_index[i] = 41
            elif element_name[i] == 'Mo':
                dummy_index[i] = 42
            elif element_name[i] == 'Tc':
                dummy_index[i] = 43
            elif element_name[i] == 'Ru':
                dummy_index[i] = 44
            elif element_name[i] == 'Rh':
                dummy_index[i] = 45
            elif element_name[i] == 'Pd':
                dummy_index[i] = 46
            elif element_name[i] == 'Ag':
                dummy_index[i] = 47
            elif element_name[i] == 'Cd':
                dummy_index[i] = 48
            elif element_name[i] == 'In':
                dummy_index[i] = 49
            elif element_name[i] == 'Sn':
                dummy_index[i] = 50
            elif element_name[i] == 'Sb':
                dummy_index[i] = 51
            elif element_name[i] == 'Te':
                dummy_index[i] = 52
            elif element_name[i] == 'I':
                dummy_index[i] = 53
            elif element_name[i] == 'Xe':
                dummy_index[i] = 54
            elif element_name[i] == 'Cs':
                dummy_index[i] = 55
            elif element_name[i] == 'Ba':
                dummy_index[i] = 56
            elif element_name[i] == 'La':
                dummy_index[i] = 57
            elif element_name[i] == 'Ce':
                dummy_index[i] = 58
            elif element_name[i] == 'Pr':
                dummy_index[i] = 59
            elif element_name[i] == 'Nd':
                dummy_index[i] = 60
            elif element_name[i] == 'Pm':
                dummy_index[i] = 61
            elif element_name[i] == 'Sm':
                dummy_index[i] = 62
            elif element_name[i] == 'Eu':
                dummy_index[i] = 63
            elif element_name[i] == 'Gd':
                dummy_index[i] = 64
            elif element_name[i] == 'Tb':
                dummy_index[i] = 65
            elif element_name[i] == 'Dy':
                dummy_index[i] = 66
            elif element_name[i] == 'Ho':
                dummy_index[i] = 67
            elif element_name[i] == 'Er':
                dummy_index[i] = 68
            elif element_name[i] == 'Tm':
                dummy_index[i] = 69
            elif element_name[i] == 'Yb':
                dummy_index[i] = 70
            elif element_name[i] == 'Lu':
                dummy_index[i] = 71
            elif element_name[i] == 'Hf':
                dummy_index[i] = 72
            elif element_name[i] == 'Ta':
                dummy_index[i] = 73
            elif element_name[i] == 'W':
                dummy_index[i] = 74
            elif element_name[i] == 'Re':
                dummy_index[i] = 75
            elif element_name[i] == 'Os':
                dummy_index[i] = 76
            elif element_name[i] == 'Ir':
                dummy_index[i] = 77
            elif element_name[i] == 'Pt':
                dummy_index[i] = 78
            elif element_name[i] == 'Au':
                dummy_index[i] = 79
            elif element_name[i] == 'Hg':
                dummy_index[i] = 80
            elif element_name[i] == 'Tl':
                dummy_index[i] = 81
            elif element_name[i] == 'Pb':
                dummy_index[i] = 82
            elif element_name[i] == 'Bi':
                dummy_index[i] = 83
            elif element_name[i] == 'Po':
                dummy_index[i] = 84
            elif element_name[i] == 'At':
                dummy_index[i] = 85

	#if spe[0] == 'N   1':
	#	znum_int[0] = 0

        # here the index to connect name and atomic numbers.
        global index_z_for_elements
        index_z_for_elements = {}    
        for a,b in zip(element_name,dummy_index):
            index_z_for_elements[a]=b


