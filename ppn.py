"""
Analyze and visualize the output of single-zone ppn simulations

Contains:
class xtime - for analyzing x-time.data time evolutions
class PPN - for analyzing selem0xxxx.DAT time evolutions

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
        self.sldir = sldir
		
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
		cols=[]
		cols.append(line[1:4])
		cols.append(line[15:17])
		cols.append(line[23:26])
		cols.append('sum_yps')
		cols.append('NEUT')
		cols.append('PROT')
		i=81
		ispec=2
		while True:
		    if line[i:i+5] is '':
			break
		    else:
			cols.append(line[i:i+5])
			ispec+=1
			i+=13
		print "There are "+str(ispec)+" species found."
		col_tot=ispec+4
		
		col_num={}        
		for a,b in zip(cols,range(ispec)):
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
    
    def plot_xtime(self,Y,X='t_y',label='default',labelX=None, labelY=None , 
    	title=None, shape='.',logX=False, logY=True, base=10):
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
            lab_str=Y
        else:
            lab_str=label
        
        try:
        	self.get(X)
        except KeyError:
        	X='age'
	
        DataPlot.plot(self,X,Y,legend=lab_str,labelX=labelX, labelY=labelY, 
        	      title=title, shape=shape,logX=logX, logY=logY, base=base)
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
	Rangeing from 0 to 38

	To find the cycle attributes:
	>>> p.cattrs
	['mod', 'dzeit', 'agej', 't9', 'rho', 'densn', 'densp', 'densa']
	
	To find the data column attributes
	>>> p.dcols
	['NUM', 'Z', 'A', 'ISOM', 'ABUNDNACE_MF', 'ISOTP']
	>>> p.get('Z',0)
	array([1, 2, 2, 4, 5, 3, 6, 6, 7, 7, 6, 7, 8, 8, 8, 9, 9, 9])
	>>> p.get('ABUNDNACE_MF',0)
	array([  1.43722000e-10,   1.00000000e-99,   9.81499000e-01,
         4.08738000e-20,   1.00000000e-99,   2.06944000e-21,
         3.42800000e-04,   9.62307000e-05,   1.05081000e-12,
         1.25518000e-02,   1.90131000e-08,   1.42230000e-07,
         4.98449000e-05,   4.80246000e-07,   4.24345000e-12,
         9.85201000e-17,   6.30866000e-16,   9.12726000e-11])
         
        or if the user wants the data from the first 3 cycles:
        >>> p.get('ABUNDNACE_MF',[0,1,2])
        [array([  1.43722000e-10,   ...,   9.81499000e-01,],
        array([  1.43722000e-10,   ...,   9.81499000e-01,],
        array([  1.43722000e-10,   ...,   9.81499000e-01,]
        >>> p.getElement('C 14',0)
        array([  1.50000000e+01,   6.00000000e+00,   1.40000000e+01,
         1.90131000e-08])
        >>> p.plot('abundance_yps', 'Z',0)
	plots data
	>>> p.iso_abund(0)
	Plots an isotope abundence distribution
	>>> p.abu_chart(0)
	Plots an isotope abundence chart
	
	One note about the plot functions, if instead of a single cycle the user 
	inputs a list of cycles, the method will then, instead of plotting them, 
	will then save a .png for each cycle. Also if you just want a singular 
	plot saved, the user can input their cycle, in a list like [0]. And that 
	will save their plot.
	'''
	sldir=''  #Standard Directory
	cattrs={} # cycle attributes
	dcols=[]  # list of the column attributes
	index=0   # index of were column data begins in the file
	files=[]  # list of files
	def __init__(self,sldir='./', filenames='iso_massf'):
		''' 
		initial method of this class
		Input:
			fname - The .DAT file the user is looking at.
			sldir - where fname exists
			filenames - the default finenames of the abundence vectors
				    Defaults to iso_massf
		Output: 
			A PPn instance
		
		'''
		self.sldir = sldir
		self.cattrs=[]
		self.dcols=[]
		self.files=[]	
		if not os.path.exists(sldir):  # If the path does not exist
			print 'error: Directory, '+sldir+ ' not found'
			print 'Now returning None'
			return None
		f=os.listdir(sldir) # reads the directory
		for i in range(len(f)):  
			# Removes any files that are not ppn files
			filelength=len(filenames)+4
		    	if filenames in f[i] and 'DAT' in f[i] and '~' not in f[i] and len(f[i])>filelength and 'restart' not in f[i]:
		    		self.files.append(f[i])
		
		self.files.sort()
		
		if len(self.files)==0: 
			# If there are no selem Files in thes Directory
		    	print 'Error: no '+filenames+ ' named files exist in Directory'
		    	print 'Now returning None'
		    	return None
		fname=self.files[len(self.files)-1]
		self.cattrs,self.dcols,self.index=self._readFile(fname,sldir)
		
		for i in xrange(len(self.files)):
			self.files[i]=self.sldir+self.files[i]
		print str(len(self.files))+' cycle numbers found in '+sldir
		print 'Rangeing from 0 to '+str(len(self.files)-1)
	
	def getCycleData(self,attri,fname,numType='cycNum'):
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
		
		fname=self.findFile(fname,numType)
		f=open(fname,'r')
		lines=f.readlines()
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
							return float(tmp[j+1])
						except ValueError:
							return str(tmp[j+1])

			elif lines[i].startswith('H'):
				abc=1+1 #do nothing
			else:
				print 'This cycle attribute does not exist'
				print 'Returning None'
				return None
		
	def getColData(self,attri,fname,numType='cycNum'):
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
		fname=self.findFile(fname,numType)
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
				elif tmp==('NEUT'or'NEUTR'or'nn'or'N   1'):
					data.append('NEUT')
				else:
					data.append(tmp)
			elif index==0:
				data.append(int(lines[i][index]))
			else:
				data.append(float(lines[i][index]))
		
		return array(data)
		
	def getElement(self,attri,fname,numType='cycNum'):
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
		
		fname=self.findFile(fname,numType)
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
		
		element=self.get(self.dcols[5],fname,numType)
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
			
	def get(self,attri,fname=None,numtype='cycNum'):
		'''
		In this method a column of data for the associated attribute is
		returned. If fname is a list or None a list of each cycles in 
		fname or all cycles is returned
		Input: 
		attri: The name of the attribute we are looking for.
		Fname: The name of the file we are getting the data from or
			the cycle number found in the filename. Or a List of 
			cycles or filenames.  If this is None, the data from all
			cycles is returned.
		numtype: Determines whether fname is the name of a file or, the 
			 Cycle number. If it is 'file' it will then  interpret 
			 it as a file
			 if it is 'cycNum' it will then  interpret it as a cycle 
			 number
		Output: Data in the form of a numpy array
		
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
		
	def findFile(self,fname,numType):
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
		numType=numType.upper()
		if numType=='FILE':
			
			#do nothing
			return self.sldir+fname
		
		elif numType=='CYCNUM':
			try:
				fname=int(fname)
			except ValueError:
				print 'Improper choice:'+ str(fname)
				print 'Reselecting as 0'
				fname=0
			#print 'Using '+self.files[fname]
			return self.files[fname]
