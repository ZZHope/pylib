"""
Analyze and visualize the output of single-zone ppn simulations

Contains:
class xtime - for analyzing x-time.data time evolutions
class abu_vector - for analyzing iso_massfxxxx.DAT time evolutions

PPN Assumptions:
The first non white space character in a header line is a '#'.
Header attributes are separated from their value by white space or by white space 
	surrounding an equals sign.
An Header attribute is separated by the previous Header attribute by white space
	or a line break.
There are only 6 data columns. The first being the number, second being Z, third being A, Fourth isomere state
	A fifth abundance_yps and finally the element name.
The first Five columns consist purely of numbers, no strings are allowed. 
Of the values in the final column, the name column, the first two are letters 
specifying the element name, and the rest are spaces or numbers (in that strict
order), except for the element names: Neut and Prot
All the profile files in the directory have the same cycle attributes.	
The cycle numbers of the 'filename'+xxxxx start at 0.
PPN files allways end in .DAT and are not allowed any '.'
The can not be any blank lines in the data files.
No cycle numbers are skipped, ie if cycle 0 and 3 are present, 1 and 2 must be 
	here aswell.
"""

from numpy import *
from data_plot import *
import matplotlib
from matplotlib.pylab import *
from data_plot import *
from utils import *
import os


class xtime(DataPlot):
    ''' read and plot x-time.dat output files
    Example:
    In [1]: import ppn
    In [3]: f=ppn.xtime()
    There are 1099 species found.
    There are 19 time steps found.

    In [4]: f.col
    f.col_num  f.col_tot  f.cols     

    In [4]: f.cols[0:10]
    Out[4]: 
    ['age',
    't9',
    'rho',

    'sum_yps',
    'NEUT',
    'PROT',
    'H   2',
    'HE  3',
    'HE  4',
    'BE  7']

    In [5]: f.plot('HE  4')
    age HE  4
    
    In [6]: f.plot('C  12')
    age C  12
    
    In [7]: f.plot('PROT')
    '''
    sldir  = '' #Standard Directory
    data = []
    cols=[]	#list of column attribute names
    col_num = []#Dict of column attribute names and their associated values
    xdat = []
    ydat = []
    
    def __init__(self,sldir='./',fname='x-time.dat'):
        ''' read x-time.dat file
        input:
		sldir - the directory of the pecified filename
		fname - specify alternative filename of file of type x-time.dat
	output:
		A xtime instance
        '''
        self.sldir= sldir 
		
	if not os.path.exists(sldir):  # If the path does not exist
		print 'error: Directory, '+sldir+ ' not found'
		print 'Now returning None'
		return None
	else:
		self.data, self.col_num, self.cols, self.col_tot, self.ilines =\
		self._readFile(fname,sldir)
    
   
    def _readFile(self,fname,sldir):
    	        '''
    	        Private method that reads in the data file and organizes it within 
    	        this object.
    	        '''
		if sldir.endswith('/'):
			fname = str(sldir)+str(fname)
		else:
			fname = str(sldir)+'/'+str(fname)
	   	f=open(fname,'r')
	
	# read header line
		line=f.readline()
		cols  = []
		ispec = 0
		for i in range(1,len(line.split('|'))):
		    col = line.split('|')[i].strip()
		    if '-' in col:
			ispec += 1
			col   = col.split('-')[1]
		    cols.append(col)
		    col_num={}
		col_tot = len(cols)

		print 'number of species: ', str(ispec)
		print 'number of cols: ', str(col_tot)


		col_num={}
		for a,b in zip(cols,range(col_tot)):
			col_num[a]=b
	
	# read remainder of the file
		lines=f.readlines()
		data=[]
		for i in range(len(lines)):
		    v=lines[i].split()
		    vv=array(v,dtype='float')
		    data.append(vv)
		ilines=i
		print "There are "+str(ilines)+" time steps found."
		return data,col_num,cols,col_tot,ilines
    
    def get(self,col_str):
        '''get one data column with the data

        col_str - one of the column strings in self.cols
        '''
        data_column=zeros(self.ilines)
        for i in range(self.ilines):
            data_column[i]=self.data[i][self.col_num[col_str]]
        return data_column
    
    def plot_xtime(self,y,x='t_y',label='default',labelx=None, labely=None , 
    	title=None, shape='.',logx=False, logy=True, base=10):
        '''make a simple plot of two columns against each other.
        An example would be instance.plot_xtime('PB206', label='PB206 vs t_y'
        Recomend using the plot function DataPlot.plot() it has more functionality 
        Y         - column on Y-axis
        X         - column on X-axis, defaults to "t_y"
        label	  - Legend label
        labelX: The label on the X axis
	labelY: The label on the Y axis
	title: The Title of the Graph
	logX: A boolean of weather the user wants the x axis logarithmically
	logY: A boolean of weather the user wants the Y axis logarithmically
	base: The base of the logarithm. Default = 10
	shape: What shape and colour the user would like their plot in.
	Please see 
	http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
	for all possable choices
        '''
        if label is 'default':
            lab_str=y
        else:
            lab_str=label
        
        try:
        	self.get(x)
        except KeyError:
        	x='age'
	
        DataPlot.plot(self,x,y,legend=lab_str,labelx=labelx, labely=labely, 
        	      title=title, shape=shape,logx=logx, logy=logy, base=base)
        '''
        print X,Y
        xdat=self.get(X)
        ydat=self.get(Y)
        self.xdat = xdat
        self.ydat = ydat
        

        plot(xdat,log10(ydat),label=lab_str)
        legend()
        '''
        
class abu_vector(DataPlot,Utils):
	'''
	Class for reading selem00xxxx.DAT files
	Example run through for cycle 0:
	>>> import ppn
	>>> p=ppn.abu_vector('./run/')
	39 cycle numbers found in ./run/
	Ranging from 0 to 38

	To find the cycle attributes:
	>>> p.cattrs
	['mod', 'dzeit', 'agej', 't9', 'rho', 'densn', 'densp', 'densa']
	
	To find the data column attributes
	>>> p.dcols
	['NUM', 'Z', 'A', 'ISOM', 'ABUNDANCE_MF', 'ISOTP']
	>>> p.get('Z',0)
	array([1, 2, 2, 4, 5, 3, 6, 6, 7, 7, 6, 7, 8, 8, 8, 9, 9, 9])
	>>> p.get('ABUNDANCE_MF',0)
	array([  1.43722000e-10,   1.00000000e-99,   9.81499000e-01,
         4.08738000e-20,   1.00000000e-99,   2.06944000e-21,
         3.42800000e-04,   9.62307000e-05,   1.05081000e-12,
         1.25518000e-02,   1.90131000e-08,   1.42230000e-07,
         4.98449000e-05,   4.80246000e-07,   4.24345000e-12,
         9.85201000e-17,   6.30866000e-16,   9.12726000e-11])
         
        or if the user wants the data from the first 3 cycles:
        >>> p.get('ABUNDANCE_MF',[0,1,2])
        [array([  1.43722000e-10,   ...,   9.81499000e-01,],
        array([  1.43722000e-10,   ...,   9.81499000e-01,],
        array([  1.43722000e-10,   ...,   9.81499000e-01,]
        >>> p.getElement('C 14',0)
        array([  1.50000000e+01,   6.00000000e+00,   1.40000000e+01,
         1.90131000e-08])
        >>> p.plot('abundance_yps', 'Z',0)
	plots data
	>>> p.iso_abund(0)
	Plots an isotope abundance distribution
	>>> p.abu_chart(0)
	Plots an isotope abundance chart
	
	One note about the plot functions, if instead of a single cycle the user 
	inputs a list of cycles, the method will then, instead of plotting them, 
	will then save a .png for each cycle. Also if you just want a singular 
	plot saved, the user can input their cycle, in a list like [0]. And that 
	will save their plot.
	'''
	sldir = ''  #Standard Directory
	inputdir = '' # A copy of Standard Directory which never changes 
	cattrs={} # cycle attributes
	dcols=[]  # list of the column attributes
	index=0   # index of were column data begins in the file
	files=[]  # list of files
	isotopes=[]# list of isotopes 
	def __init__(self,sldir='./', filenames='iso_massf'):
		''' 
		initial method of this class
		Input:
			fname - The .DAT file the user is looking at.
			sldir - where fname exists
			file-names - the default file-names of the abundance vectors
				    Defaults to iso_massf
		Output: 
			A PPN instance
		
		'''
                self.debug=False
                self._stable_names() # provides in addition to stable_el from 
                                     # utils also just the stable element names
		self.sldir = sldir
		self.inputdir = ''
		self.startdir = os.getcwd()
		self.cattrs=[]
		self.dcols=[]
		self.files=[]
		self.isotopes=[]
		if not os.path.exists(sldir):  # If the path does not exist
			print 'error: Directory, '+sldir+ ' not found'
			print 'Now returning None'
			return None
		f=os.listdir(sldir) # reads the directory
		for file in f:  
			# Removes any files that are not ppn files
			filelength=len(filenames)+4
		    	if filenames in file and 'DAT' in file and '~' not in file \
                                and '#' not in file and len(file)>filelength \
                                and 'restart' not in file:
		    		self.files.append(file)
		self.files.sort()
		
		if len(self.files)==0: 
			# If there are no Files in the directory
		    	print 'Error: no '+filenames+ ' named files exist in Directory'
		    	print 'Now returning None'
		    	return None	
		fname=self.files[len(self.files)-1]
		self.cattrs,self.dcols,self.index=self._readFile(fname,sldir)
		
                indexp_cyc2filels={}  # created index pointer from mod (cycle
                i = 0                 # name) to index in files array
                for file in self.files:
                    mod=self.get('mod',fname=file,numtype='file')
                    indexp_cyc2filels[mod] = i
                    i += 1
                self.indexp_cyc2filels = indexp_cyc2filels

		for i in xrange(len(self.files)):
			self.files[i]=self.sldir+self.files[i]
		print str(len(self.files))+' cycle numbers found in '+sldir
		print 'Ranging from 0 to '+str(len(self.files)-1)
		self.isotopes=self.get('ISOTP',self.files[0],numtype='file')
		
	
	def getCycleData(self,attri,fname,numtype='cycNum'): 
		"""
		In this method a column of data for the associated cycle 
		attribute is returned
		Input: 
		attri: The name of the attribute we are looking for.
		Fname: the name of the file we are getting the data from or
			the cycle number found in the filename.
		numtype: Determines whether fname is the name of a file or, the 
			 Cycle number. If it is 'file' it will then  interpret 
			 it as a file
			 if it is 'cycNum' it will then  interpret it as a cycle 
			 number
		"""
		
		fname=self.findFile(fname,numtype) 
		
		if self.inputdir == '':
			self.inputdir = self.sldir      # This chunk of code changes into the directory where fname is,		
		os.chdir(self.inputdir)				  # and appends a '/' to the directory title so it accesses the
		self.sldir=os.getcwd() + '/'		  # file correctly
		
		f=open(fname,'r')
		lines=f.readlines()
		
		if self.inputdir != './': 				#This chunk of code changes back into the directory you started in.
			os.chdir(self.startdir)				
			self.sldir = self.inputdir

			
		for i in range(len(lines)):
			lines[i]=lines[i].strip()
			
		for i in range(len(lines)):
			if lines[i].startswith('#'):
				lines[i]=lines[i].strip('#')
				tmp=lines[i].split()
				tmp1=[]
				for j in range(len(tmp)):
					if tmp[j] != '=' or '':
						tmp1.append(tmp[j])
				tmp=tmp1
				for j in range(len(tmp)):
					if tmp[j]== attri:
						try:
                                                    if '.' in tmp[j+1]:
                                                        return float(tmp[j+1])
                                                    else:
                                                        return int(tmp[j+1])
						except ValueError:
							return str(tmp[j+1])

			elif lines[i].startswith('H'):
				continue
			else:
				print 'This cycle attribute does not exist'
				print 'Returning None'
				return None
		
		
			
		
	def getColData(self,attri,fname,numtype='cycNum'):
		"""
		In this method a column of data for the associated column 
		attribute is returned
		Input: 
		attri: The name of the attribute we are looking for.
		Fname: the name of the file we are getting the data from or
			the cycle number found in the filename.
		numtype: Determines whether fname is the name of a file or, the 
			 Cycle number. If it is 'file' it will then  interpret 
			 it as a file
			 if it is 'cycNum' it will then  interpret it as a cycle 
			 number
		
		"""
		fname=self.findFile(fname,numtype)
		f=open(fname,'r')
		for i in range(self.index+1):
			f.readline()
		lines=f.readlines()
		for i in range(len(lines)):
			lines[i]=lines[i].strip()
			lines[i]=lines[i].split()
		index=0
		data=[]
		
		while index < len (self.dcols):
			if attri== self.dcols[index]:
				break
			index+=1
			
		for i in range(len(lines)):
			
			if index==5 and len(lines[i])==7:
				data.append(str(lines[i][index].capitalize())+'-'\
				+str(lines[i][index+1]))
			elif index==5 and len(lines[i])!=7:
				tmp=str(lines[i][index])
				if tmp[len(tmp)-1].isdigit():
					tmp1=tmp[0]+tmp[1]
					tmp1=tmp1.capitalize()
					tmp2=''
					for j in range(len(tmp)):
						if j == 0 or j == 1:
							continue
						tmp2+=tmp[j]
					data.append(tmp1+'-'+tmp2)
				elif tmp=='PROT':
					data.append('H-1')
				elif tmp==('NEUT'or'NEUTR'or'nn'or'N   1'or'N-1'):
					data.append('N-1')
				else:
					data.append(tmp)
			elif index==0:
				data.append(int(lines[i][index]))
			else:
				data.append(float(lines[i][index]))
		
		return array(data)
		
	def getElement(self,attri,fname,numtype='cycNum'):
		'''
		In this method instead of getting a particular column of data,
		the program gets a paticular row of data for a paticular 
		element name.
		Input: attri is the name of the attribute we are looking for.
		       A complete list of them can be obtaind by calling,
		       'get('element_name')'
		Fname: the name of the file we are getting the data from or
			the cycle number found in the filename.
		numtype: Determines whether fname is the name of a file or, the 
			 Cycle number. If it is 'file' it will then  interpret 
			 it as a file
			 if it is 'cycNum' it will then  interpret it as a 
			 cycle number
		Output:
			A numpy array of the four ellement attributes, number, Z, A
			and abundance, in that order
		Warnig
		'''
		element=[] #Variable for holding the list of element names
		number=[]  #Variable for holding the array of numbers
		z=[]	   #Variable for holding the array of z
		a=[]	   #Variable for holding the array of a
		abd=[]	   #Variable for holding the array of Abundance
		data=[]	   #variable for the final list of data
		
		fname=self.findFile(fname,numtype)
		f=open(fname,'r')
		for i in range(self.index+1):
			f.readline()
		lines=f.readlines()
		for i in range(len(lines)):
			lines[i]=lines[i].strip()
			lines[i]=lines[i].split()
		index=0
		data=[]
		
		while index < len (self.dcols):
			if attri== self.dcols[index]:
				break
			index+=1
		
		element=self.get(self.dcols[5],fname,numtype)
		number=[]
		z=[]
		a=[]
		isom=[]
		abd=[]
		for i in range(len(lines)):
			number.append(int(lines[i][0]))
			z.append(float(lines[i][1]))
			isom.append(float(lines[i][2]))
			abd.append(float(lines[i][1]))
		index=0 #Variable for determing the index in the data columns
		
		
		while index < len(element):
			if attri == element[index]:				
				break
			
			index+=1

		data.append(number[index])
		data.append(z[index])
		data.append(a[index])
		data.append(isom[index])
		data.append(abd[index])
		
		return array(data)



        def get(self,attri,fname=None,numtype='cycNum',decayed=False): 
		'''In this method all data for an entire cycle (basically
		the content of an iso_massfnnnn.DAT file) or a column
		of data for the associated attribute is returned. 

		1: Input attri is string attribute: 
		attri: The name of the attribute we are looking for or the cycle. 
		fname: The name of the file we are getting the data from or
			the cycle number found in the filename. Or a List of 
			cycles or filenames.  If this is None, the data from all
			cycles is returned.
		numtype: Determines whether fname is the name of a file or, the 
			 Cycle number. If it is 'file' it will then  interpret 
			 it as a file
			 if it is 'cycNum' it will then  interpret it as a cycle 
			 number
                decayed  not supported in attri mode
		Output: Data in the form of a numpy array
		
		2: Input attri is just one integer cycle number (cycle arrays are not supported):
                decayed   boolean: instantaneously decay abundance distribution
                Output: the following varibales will be added to the instance
                a_iso_to_plot      mass number of plotted range of species"
                isotope_to_plot    corresponding list of isotopes"
                z_iso_to_plot      corresponding charge numbers"
                el_iso_to_plot     corresponding element names"
                abunds             corresponding abundances"
                isom               list of isomers with their abundances
		'''
                if isinstance(attri,int):
                    print "Calling get method in cycle mode, adding a_iso_to_plot, z.. el.. isotope.. isotope... to instance"
                    self._getcycle(attri,decayed)
                elif isinstance(attri,str):
                    data=self._getattr(attri,fname,numtype)
                    return data

	def _getcycle(self,cycle,decayed=False):
            ''' Private method for getting a cycle, called from get.
            '''
            yps=self.get('ABUNDANCE_MF', cycle)
            z=self.get('Z', cycle) #charge
            a=self.get('A', cycle) #mass
            isomers=self.get('ISOM', cycle)
            
            a_iso_to_plot,z_iso_to_plot,abunds,isotope_to_plot,el_iso_to_plot,isom=\
                self._process_abundance_vector(a,z,isomers,yps)

            self.a_iso_to_plot=a_iso_to_plot
            self.isotope_to_plot=isotope_to_plot
            self.z_iso_to_plot=z_iso_to_plot
            self.el_iso_to_plot=el_iso_to_plot
            self.abunds=np.array(abunds)
            self.isom=isom

            if decayed:
                try:
                    self.decay_idp
                except AttributeError:
                    print "WARNING: decayed in _getcycle ignores isomers " \
                        "and will decay alpha-unstable p-rich nuclei as if they were beta+ stable."
                    print "Initialising decay index pointers ...."
                    self.decay_indexpointer() # provides self.decay_idp and
                ind_tmp=self.idp_to_stables_in_isostoplot                 

                isotope_decay=array(isotope_to_plot)[ind_tmp]
                z_iso_decay=array(z_iso_to_plot)[ind_tmp]
                a_iso_decay=array(a_iso_to_plot)[ind_tmp]
                el_iso_decay=array(el_iso_to_plot)[ind_tmp]
                abunds_decay=zeros(len(ind_tmp), dtype='float64')
                for i in xrange(len(isotope_to_plot)):
                    idp=where(isotope_decay==isotope_to_plot[self.decay_idp[i]])[0] # points from
                    # i on isotope_to_plot scale to decay target_on_decayed array scale
                    abunds_decay[idp] += abunds[i]

                if self.debug:
                    print "Decayed array:"
                    for i in xrange(len(ind_tmp)):
                        print isotope_decay[i], z_iso_decay[i], a_iso_decay[i], el_iso_decay[i], abunds_decay[i]

                self.a_iso_to_plot=a_iso_decay
                self.isotope_to_plot=isotope_decay
                self.z_iso_to_plot=z_iso_decay
                self.el_iso_to_plot=el_iso_decay
                self.abunds=abunds_decay


	def _getattr(self,attri,fname=None,numtype='cycNum'):
            ''' Private method for getting an attribute, called from get.
            '''
            if str(fname.__class__)=="<type 'list'>":
                isList=True
            else:
                isList=False
		
            data=[]
            if fname==None:
                fname=self.files
                numtype='file'
                isList=True
            if isList:
                for i in xrange(len(fname)):
                    if attri in self.cattrs:
                        data.append(self.getCycleData(attri,fname[i],numtype))
                    elif attri in self.dcols:
                        data.append(self.getColData(attri,fname[i],numtype))
                    elif attri in self.get('ISOTP',fname,numtype):
                        data.append(self.getElement(attri,fname[i],numtype))
                    else:
                        print 'Attribute '+attri+ ' does not exist'
                        print 'Returning none'
                        return None
                    
            else:
                if attri in self.cattrs:
                    return self.getCycleData(attri,fname,numtype)
                elif attri in self.dcols:
                    return self.getColData(attri,fname,numtype)
                elif attri in self.get('ISOTP',fname,numtype):
                    return self.getElement(attri,fname,numtype)
                else:
                    print 'Attribute '+attri+ ' does not exist'
                    print 'Returning none'
                    return None
		
            return data
			
	def _readFile(self,fname,sldir):
		'''
		private method that reads in and organizes the .DAT file
		Loads the data of the .DAT File into the variables cattrs and cols.
		In both these cases they are dictionaries, but in the case of cols,
		it is a dictionary of numpy array exect for the element , 
		element_name where it is just a list
		
		'''
		cattrs=[]
		if sldir.endswith(os.sep): 
			#Making sure fname will be formated correctly
			fname = str(sldir)+str(fname)
		else:
			fname = str(sldir)+os.sep+str(fname)
			self.sldir+=os.sep
		f=open(fname,'r')
		lines=f.readlines()
		for i in range(len(lines)):
			lines[i]=lines[i].strip()
			
		
		cols=lines[0].strip('H')
		cols=cols.strip()
		cols=cols.split()
		for i in range(len(lines)):
			if lines[i].startswith('#'): 
				# if it is a cycle attribute line
				lines[i]=lines[i].strip('#')
				tmp=lines[i].split()
				tmp1=[]
				for j in range(len(tmp)):
					if tmp[j] != '=' or '':
						tmp1.append(tmp[j])
				tmp=tmp1
				
				j=0
				while j <len(tmp):
					cattrs.append(tmp[j])
					j+=2
					
			
			elif not lines[i].startswith('H'):
				index = i-1
				break
		
		
		return cattrs,cols, index
		
	def findFile(self,fname,numtype):
		"""
		Function that finds the associated file for FName when Fname is 
		time or NDump
		input:
		numType designates how this function acts and how it interprets 
			FName
		if numType is 'file', this function will get the desired 
			attribute from that file
		if numType is 'cycNum'this function will get the desired 
			attribute from that file
			with fname's model number
		Fname the name of the file we are looking or
		"""
		numType=numtype.upper() 
		if numType == 'FILE':
                    #do nothing
                    return fname
		elif numType == 'CYCNUM':
                    try:
                        fname = int(fname)
                    except ValueError:
                        print 'Improper choice:'+ str(fname)
                        print 'Reselecting as 0'
                        fname = 0
                        print 'Using '+self.files[fname]
                try:
                    return self.files[self.indexp_cyc2filels[fname]]
                except IndexError:
                    mods = array(self.get('mod'), dtype=int)
                    if fname not in mods:
                        print 'You seem to try to plot a cycle that is not present: '+str(fname)
                        fname = mods[-1]
                        print 'I will assume you want to plot the last cycle in the run: '+str(fname)
                        print '[I am not 100% sure this escape is debugged. You better do this again with'
                        print 'the correct input.]'
                        return self.files[fname]

