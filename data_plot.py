"""
data_plot.py

SuperClass module for the YProfile and Table Classes
It contains numerous plot function for the YProfile and Table Classes

If one in the future wants their class to inherit this superclasses methods
this is what is required:
A. Place 'from data_table import *' at the top of the module
B. If the class is defined like 'class MyClass:', change that to 
   'class MyClass(DataTable):' 
C. To properly use DataTable's methods properly one will need these methods:
	a get(atri) that returns a numpy array of Data, or a 
	list of numpy arrays of data.  The arguments of this function would need to be
	atri which is the name of the data one is looking for. 

"""
from numpy import *
from math import *
import matplotlib.pylab as pyl
import matplotlib.pyplot as pl
from matplotlib.mpl import colors,cm
from matplotlib.patches import Rectangle, Arrow
from matplotlib.collections import PatchCollection
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
import os
import threading
import time
import sys

class DataPlot:
	
	def classTest(self):
		'''
		Determines what the type of class instance the subclass is, so
		we can dynically determine the behaviour of methods.
		
		This method NEEDS to be modified if any names of files or classes
		are changed
		'''
		
		c=str(self.__class__)
		tmp=''
		if 'ppm.yprofile' == c:
			tmp='YProfile'
		elif 'ascii_table.ascii_table' == c:
			tmp='AsciiTable'
		elif 'mppnp.se' == c:
			tmp='se'
		elif 'mesa.mesa_profile' == c:
			tmp='mesa_profile'
		elif 'mesa.star_log' == c:
			tmp='mesa.star_log'
		elif 'ppn.xtime' == c:
			tmp='xtime'
		elif 'ppn.abu_vector' == c:
			tmp='PPN'
		
		return tmp
		
	def which(self, program):
	    '''
	    Mimics which in the unix shell
	    '''
	    def is_exe(fpath):
		return os.path.exists(fpath) and os.access(fpath, os.X_OK)
	
	    fpath, fname = os.path.split(program)
	    if fpath:
		if is_exe(program):
		    return program
	    else:
		for path in os.environ["PATH"].split(os.pathsep):
		    exe_file = os.path.join(path, program)
		    if is_exe(exe_file):
			return exe_file
	
	    return None

	def logarithm(self,tmpX,tmpY,logX,logY,base):
		logXER=False
		logYER=False
		for i in range(len(tmpX)):
				if tmpX[i]<=0. and logX:
					print 'We can not log a number less than or equal to zero'
					print 'Attempting to remove incompatible values from X'
					logXER=True
				if tmpY[i]<=0. and logY:
					print 'We can not log a number less than or equal to zero'
					print 'Attempting to remove incompatible values from Y'
					logYER=True
		tmX=[]
		tmY=[]
		
		if logXER:
			for i in range(len(tmpX)):
				if tmpX[i]>0.:
					tmX.append( tmpX[i])
					tmY.append(tmpY[i])
			tmpX=tmX
			tmpY=tmY
		elif logYER:
			for i in range(len(tmpY)):
				if tmpY[i]>0.:
					tmX.append( tmpX[i])
					tmY.append(tmpY[i])
			tmpX=tmX
			tmpY=tmY
		
		tmX=tmpX
		tmY=tmpY
			
		if logX:
			tmX=tmpX
			try:
				for i in range(len(tmpX)):
					tmX[i]=log(tmpX[i],base)
			except ValueError:
				print 'We can not log a number less than or equal to zero'
				print 'Attempting to remove incompatible values from X'
				logXER=True
		if logY:
			tmY=tmpY
			try:
				for i in range(len(tmpY)):
					tmY[i]=log(tmpY[i],base)
			except ValueError:
				print 'We can not log a number less than or equal to zero'
				print 'Attempting to remove incompatible values from Y'
				logYER=True
		
		if logX:
			tmpX=tmX
		if logY:
			tmpY=tmY
		
		return tmpX,tmpY
		
	def sparse(self,x,y,sparse):
			"""
			Method that removes every non sparse th element.  For example 
			if this argument was 5, This method would plotthe 0th, 5th, 10th
			... elements.
			Input:
			x: list of x values, of lenthe j
			y: list of y values, of lenthe j
			sparse: Argument that skips every so many data points
			"""
			tmpX=[]
			tmpY=[]
			
			for i in range(len(x)):
				if sparse == 1:
					return x,y
				if (i%sparse)==0:
					tmpX.append(x[i])
					tmpY.append(y[i])
			return tmpX, tmpY
			
	def plotMulti(self,atrix,atriy, cyclist,title,legend=None,labelx=None, labely=None,logx=False, logy=False, base=10,sparse=1,pdf=False,limits=None,):
		'''
		Method for plotting multiple plots and saving it to multiple pngs 
		or PDFs
		Input:
		atrix: The name of the attribute you want on the x axis
		atriy: The name of the attribute you want on the Y axis
		cyclist: List of cycles that you would like plotted
		title: The title of the grapgh and the name of the file.
		Legend: A list of legends for each of your cycles, or one legend 
			for all of the cycles
		pdf: A boolean of if the image should be saved to a pdf file.
			xMin,xMax, yMin, YMax:  plot coopdinates.	
		logx: A boolean of weather the user wants the x axis logarithmically
		logy: A boolean of weather the user wants the Y axis logarithmically
		base: The base of the logarithm. Default = 10
		sparse: Argument that skips every so many data points. For 
			example if this argument was 5, This method would plot
			the 0th, 5th, 10th ... elements.
		limits: The length four list of the x and y limits. The order of
			the list is xmin,xmax,ymin,ymax
		'''
		print 'This method may achieve speedup by calling this method from python or ipython, rather than ipython --pylab --q4thread'
		if str(legend.__class__)!="<type 'list'>":# Determines the legend is a list
			legendList=False
		else:
			legendList=True
			
		if legendList and len(cyclist) !=len(legend): #if it is a list, make sure there is an entry for each cycle
			print 'Please input a proper legend, with correct length, aborting plot'
			return None
		for i in xrange(len(cyclist)):
			if legendList:
				self.plot(atrix,atriy,cyclist[i],'ndump',legend[i],labelx,labely,base=base,sparse=sparse, logx=logx,logy=logy,show=False,limits=limits)
			else:
				self.plot(atrix,atriy,cyclist[i],'ndump',legend,labelx,labely,base=base,sparse=sparse, logx=logx,logy=logy,show=False,limits=limits)
			
			pl.title(title)
			if not pdf:
				pl.savefig(title+str(cyclist[i])+'.png', dpi=400)
			else:
				pl.savefig(title+cyclist[i]+'.pdf', dpi=400)
			pl.clf()
		return None
	
	def plot(self,atrix,atriy, fname=None,numtype='ndump',legend=None,labelx=None, labely=None ,
		indexx=None, indexy=None, title=None, shape='.',logx=False, logy=False, base=10,sparse=1, show=True,limits=None):
		"""
		Simple function that plots atriy as a function of atrix
		This method will automatically find and plot the requested data.
		Input:
		atrix, The name of the attribute you want on the x axis
		atriy, The name of the attribute you want on the Y axis
		Fname: Be the filename, Ndump or time, or cycle,  If fname is a 
		       list, this method will then save a png for each cycle in
		       the list. Warning, this must be a list of cycles and not
		       a list of filenames
		numtype: designates how this function acts and how it interprets 
			 fname. Defaults to file
		if numtype is 'file', this function will get the desird 
		attribute from that file
		if numtype is 'NDump' function will look at the cycle with that 
		nDump
		if numtype is 't' or 'time' function will find the _cycle with 
		the closest time stamp 
		labelx: The label on the X axis
		labely: The label on the Y axis
		indexx: Depreciated: If the get method returns a list of lists, 
			indexx would be the list at the index indexx in the list.
		indexy: Depreciated: If the get method returns a list of lists, 
			indexy would be the list at the index indexx in the list.
		shape: What shape and colour the user would like their plot in.
		       Please see 
		       http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
		       for all possable choices
		title: The Title of the Graph
		logx: A boolean of weather the user wants the x axis logarithmically
		logy: A boolean of weather the user wants the Y axis logarithmically
		base: The base of the logarithm. Default = 10
		sparse: Argument that skips every so many data points. For 
			example if this argument was 5, This method would plot
			the 0th, 5th, 10th ... elements.
		show: A boolean of if the plot should be displayed usefull with
			the multiPlot method
		WARNING: Unstable if get returns a list with only one element (x=[0])
		limits: The length four list of the x and y limits. The order of
			the list is xmin,xmax,ymin,ymax
		"""
		
		#Setting the axis labels

		if labelx== None :
			labelx=atrix
		if labely== None :
			labely=atriy
		
		if title!=None:
			title=title
		else:
			title=labely+' vs '+labelx
			
		if str(fname.__class__)=="<type 'list'>":
			self.plotMulti(atrix,atriy,fname,title,legend ,labelx,labely,logx, logy, 10,1,limits=limits)
			return
		tmpX=[]
		tmpY=[]
		singleX=False
		singleY=False
		#Getting data
		plotType=self.classTest()
		if plotType=='YProfile':
			if fname==None:
				fname=len(self.files)-1
		
			listY=self.get(atriy,fname, numtype,resolution='a')
			listX=self.get(atrix,fname, numtype,resolution='a')
		elif plotType=='se':
			if fname==None:
				listY=self.get( atriy)
				listX=self.get(atrix)
			else:
				listY=self.get(fname, atriy)
				listX=self.get(fname, atrix)
		elif plotType=='PPN' :
			if fname==None and atrix not in self.cattrs and atriy not in self.cattrs:
				fname=len(self.files)-1
			if numtype=='ndump':
				numtype='cycNum'
			listY=self.get(atriy,fname,numtype)
			listX=self.get(atrix,fname,numtype)
		elif plotType=='xtime' or plotType=='mesa_profile' or plotType=='AsciiTable' or plotType=='mesa.star_log':
			listY=self.get(atriy)
			listX=self.get(atrix)
		else:
			listY=self.get(atriy)
			listX=self.get(atrix)
		tmpX=[]
		tmpY=[]
		if isinstance(listX[0], basestring) or isinstance(listY[0], basestring):
			for i in range(len(listX)):
				if '*****' == listX[i] or '*****' == listY[i]:
					print 'There seems to be a string of * in the lists'
					print 'Cutting out elements in both the lists that have an index equal to or greater than the index of the location of the string of *'
					break
				tmpX.append(float(listX[i]))
				tmpY.append(float(listY[i]))
			
			listX=tmpX
			listY=tmpY
			
			
		
		
		#Determining if listX is a list or a list of lists
		try:
			j=listX[0][0]
		except:
			singleX = True
			
		if len(listX) == 1:  # If it is a list of lists with one element.
			tmpX=listX[0]
		elif singleX == True:# If it is a plain list of values.
			tmpX=listX
		elif indexx==None and len(listX)>1: # If it is a list of lists of values.
						    # take the largest
			tmpX=listX[0]
			for i in range(len(listX)):
				if len(tmpX)<len(listX[i]):
					tmpX=listX[i]
		elif indexx<len(listX): # If an index is specified, use that index
			tmpX=listX[indexx]
		else:
			print 'Sorry that indexx does not exist, returning None'
			return None
								
		#Determining if listY is a list or a list of lists
		try:
			j=listY[0][0]
		except:
			singleY = True
			
		if len(listY) == 1: # If it is a list of lists with one element.
			#print 'hello'
			tmpY=listY[0]
		elif singleY == True: # If it is a plain list of values.
			#print 'world'
			tmpY=listY
		elif indexy==None and len(listY)>1:# If it is a list of lists of values.
						    # take the largest
			#print 'fourth'
			tmpY=listY[0]
			for i in range(len(listY)):	
				if len(tmpY)<len(listY[i]):
					tmpY=listY[i]
		elif indexy<len(listY): # If an index is specified, use that index
			#print 'sixth'
			tmpY=listY[indexy]
		else:
			print 'Sorry that indexy does not exist, returning None'
			return None
		'''					
		elif indexy==None and len(listY)==1:
			#print 'fifth'
			tmpY=listY
		'''
		
		
		
		
		#Here, if we end up with different sized lists to plot, it
		#searches for a list that is of an equal length
		if len(tmpY)!=len(tmpX):
			found=False
			print "It seems like the lists are not of equal length"
			print "Now attempting to find a compatible list for ListX"
			for i in range(len(listY)):
				if not singleY and len(tmpX)==len(listY[i]):
					tmpY=listY[i]
					found=True

			if not found:
				print "Now attempting to find a compatible list for ListY"
				for i in range(len(listX)):
					
					if not singleX and len(tmpY)==len(listX[i]):
						tmpX=listX[i]
						found=True
			
			if found:
				print "Suitable list found"
			else:
		
				print "There is no suitalble list, returning None"
				return None
		if len(tmpY)!=len(tmpX) and single == True:
			print 'It seems that the selected lists are of different\nsize, now returning none'
			return None
		# Sparse stuff
		tmpX,tmpY=self.sparse(tmpX,tmpY, sparse)
		
		# Logarithm stuff
		if logy or logx:	
			tmpX,tmpY=self.logarithm(tmpX,tmpY,logx,logy,base)
		
		# Here it ensures that if we are plotting ncycle no values of '*' will be plotted
		tmX=[]
		tmY=[]
		for i in range(len(tmpX)):
			tmX.append(str(tmpX[i]))
			tmY.append(str(tmpY[i]))
		
		tmpX=[]
		tmpY=[]
		for i in range(len(tmX)):
			if '*' in tmX[i] or '*' in tmY[i]:
				print 'There seems to be a string of * in the lists'
				print 'Cutting out elements in both the lists that have an index equal to or greater than the index of the location of the string of *'
				break
			tmpX.append(float(tmX[i]))
			tmpY.append(float(tmY[i]))
		listX=tmpX
		listY=tmpY
		
		#Setting the axis labels
		
		if logx:
			labelx='log '+labelx
		if logy:
			labely='log '+labely
		
		if legend!=None:
			legend=legend
		else:
			legend=labely+' vs '+labelx	
		
		
			
		pl.plot(listX,listY,shape,label=legend)
		pl.legend()
		pl.title(title)
		pl.xlabel(labelx)
		pl.ylabel(labely)
		if show:
			pl.show()

		if limits != None and len(limits)==4:
			
			pl.xlim(limits[0],limits[1])
			pl.ylim(limits[2],limits[3])
		
		
	def clear(self, title=True, xlabel=True, ylabel=True):
		'''
		Method for removing the title and/or xlabel and/or Ylabel
		input:
		Title -  boolean of if title will be cleared 
		xlabel - boolean of if xlabel will be cleared 
		ylabel - boolean of if ylabel will be cleared 
		'''
		if title:
			pyl.title('')
		if xlabel:
			pyl.xlabel('')
		if ylabel:
			pyl.ylabel('')
	# From mesa.py
	def xlimrev(self):
		''' reverse xrange'''
		xmax,xmin=pyl.xlim()
		pyl.xlim(xmin,xmax)
		
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
		for i in xrange(len(self.elements_names)):
			if self.elements_names[i] == x[0].split('-')[0]:
				indX=i
			if self.elements_names[i] == y[0].split('-')[0]:
				indY=i

		if indX>indY:
			return 1
		if indX==indY:
			return 0
		if indX<indY:
			return -1
	def abu_chartMulti(self,cyclist, mass_range=None ,ilabel = True,imlabel = True,imagic =  False,boxstable=True,lbound=20,plotaxis=[0,0,0,0],pdf=False,title=None):
		'''
		Method that plots abundence chart and saves those figures to a .png file 
		(by default). Plots a figure for each cycle in the argument cycle
		input:
		cyclist: The list of cycles we are plotting
		ilabel: elemental labels off/on [0/1]
		imlabel: label for isotopic masses off/on [0/1]
		imagic:  turn lines for magic numbers off/on [0/1]
		plotaxis: Set axis limit: If default [0,0,0,0] the complete 
			  range in (N,Z) will be plotted
			  format: What format will this be saved in ['pdf'/'png']
		title: The title of the plots and the saved images
		'''
		print 'This method may achieve speedup calling this method from python or ipython, rather than ipython --pylab --q4thread'
		
		
		if self.which('dvipng')==None:
			print "This method may need the third party program dvipng to operate"
			print 'It is located at http://sourceforge.net/projects/dvipng/'

		for i in xrange(len(cyclist)):
			self.abu_chart( cyclist[i], mass_range ,ilabel,imlabel,imagic,boxstable,lbound,plotaxis,False)
			if title !=None:
				pl.title(title)
			else:
				name='AbuChart'
			if not pdf:
				pl.savefig(name+str(cyclist[i])+'.png', dpi=400)
			else:
				pl.savefig(name+cyclist[i]+'.pdf', dpi=400)
			pl.clf()
		
		return None
	#from mppnp.se
	def abu_chart(self, cycle, mass_range=None ,ilabel = True,imlabel = True,imagic =  False,
		boxstable=True,lbound=20,plotaxis=[0,0,0,0], show=True):
		'''
		Plots an abundence chart
		input:
		cycle: The cycle we are looking in. It it is a list of cycles, 
			this method will then do a plot for each of theses cycles 
			and save them all to a file
		ilabel: elemental labels off/on [False/True] defaults to True
		imlabel: label for isotopic masses off/on [False/True], defaults to True
		imagic:  turn lines for magic numbers off/on [False/True] defaults to False
		boxstable: plot the black boxes around the stable elements, 
			   defaults to true
		lbound: The lower bound of the colour spectrum ploted. Defaults 
			to 20
		plotaxis: Set axis limit: If default [0,0,0,0] the complete 
			range in (N,Z) will be plotted. It equates to 
			[xMin,xMax,Ymin,Ymax]
		format: What format will this be saved in ['pdf'/'png']
		mass_range - a 1x2 array containing the lower and upper mass range.
		    		 If this is an instance of abu_vector this will 
		    		 only plot isotopes that have an atominc mass 
		    		 within this range. This will throw an error if
		    		 this range does not make sence ie [45,2]
			 	if None, it will plot over the entire range
				Defaults to None
		show:  boolean of if the plot should be displayed useful with
		       saving multiple plots using abu_chartMulti		
		'''
		#######################################################################
		#### plot options
		# Set axis limit: If default [0,0,0,0] the complete range in (N,Z) will 
		# be plotted, i.e. all isotopes, else specify the limits in 
		# plotaxis = [xmin,xmax,ymin,ymax] 
		
		#######################################################################
		
		# read data file
		#inpfile = cycle
		#ff = fdic.ff(inpfile)
		
		if str(cycle.__class__)=="<type 'list'>":
			self.abu_chartMulti(cycle, mass_range,ilabel,imlabel,imagic,boxstable,lbound,plotaxis)
			return
		plotType=self.classTest()
		
		if mass_range!=None and mass_range[0]>mass_range[1]:
			print 'Please input a proper mass range'
			print 'Returning None'
			return None
		
		if plotType=='se':
			cycle=self.se.findCycle(cycle)
			nin=zeros(len(self.se.A))
			zin=zeros(len(self.se.Z))
			for i in xrange(len(nin)):
				nin[i]=self.se.A[i]
				zin[i]=self.se.Z[i]
			for i in xrange(len(nin)):
				nin[i]=nin[i]-zin[i]
			yin=self.get(cycle, 'iso_massf')
			isom=self.se.isomeric_states
			
			masses = self.se.get(cycle,'mass')
			if mass_range != None:
				masses = self.se.get(cycle,'mass')
				masses.sort()
			
			if mass_range != None:
				tmpyps=[] 
				masses = self.se.get(cycle,'mass')
				masses = self.se.get(cycle,'mass')
				masses.sort()
				for i in xrange(len(masses)):
					if (masses[i] >mass_range[0] and masses[i]<mass_range[1]) or (masses[i]==mass_range[0] or masses[i]==mass_range[1]):
						
						tmpyps.append(yin[i])
				yin=tmpyps
			
			
			tmp=zeros(len(yin[0]))
			for i in xrange(len(yin)):
				for j in xrange(len(yin[i])):
					tmp[j]+=yin[i][j]
			
			tmp=tmp/len(yin)
			
			yin=tmp
			
		elif plotType=='PPN':
			
			nin=self.get('A',cycle)
			zin=self.get('Z',cycle)
			for i in xrange(len(nin)):
				nin[i]=nin[i]-zin[i]
			yin=self.get('ABUNDNACE_MF',cycle)
			isom=self.get('ISOM',cycle)
			
			if mass_range != None:
				tmpA=[]
				tmpZ=[]
				tmpIsom=[]
				tmpyps=[]
				for i in xrange(len(nin)):
					if (nin[i] >mass_range[0] and nin[i]<mass_range[1]) or (nin[i]==mass_range[0] or nin[i]==mass_range[1]):
						tmpA.append(nin[i])
						tmpZ.append(zin[i])
						tmpIsom.append(isom[i])
						tmpyps.append(yin[i])
				zin=tmpZ
				nin=tmpA
				yin=tmpyps
				isom=tmpIsom
			
		else:
			print 'This method, abu_chart, is not supported by this class'
			print 'Returning None'
			return None
		# in case we call from ipython -pylab, turn interactive on at end again
		turnoff=False
		if not show:   
			try:
				ioff()
				turnoff=True
			except NameError:
				turnoff=False
				
		nnmax = int(max(nin))+1
		nzmax = int(max(zin))+1
		nzycheck = zeros([nnmax,nzmax,3])
		for i in range(len(nin)):
			if isom[i]==1:
				ni = int(nin[i])
				zi = int(zin[i])
				
				nzycheck[ni,zi,0] = 1
				nzycheck[ni,zi,1] = yin[i]
		
		
		
 
		#######################################################################
		# elemental names: elname(i) is the name of element with Z=i 
		
		elname=self.elements_names
		
		#### create plot
		
		## define axis and plot style (colormap, size, fontsize etc.)
		if plotaxis==[0,0,0,0]:
		  xdim=10
		  ydim=6
		else:
		  dx = plotaxis[1]-plotaxis[0]
		  dy = plotaxis[3]-plotaxis[2]
		  ydim = 6
		  xdim = ydim*dx/dy
		  
		
		params = {'axes.labelsize':  15,
			  'text.fontsize':   15,
			  'legend.fontsize': 15,
			  'xtick.labelsize': 15,
			  'ytick.labelsize': 15,
			  'text.usetex': True}
		#pl.rcParams.update(params) #May cause Error, someting to do with tex
		fig=pl.figure(figsize=(xdim,ydim),dpi=100)
		axx = 0.10
		axy = 0.10
		axw = 0.85
		axh = 0.8
		ax=pl.axes([axx,axy,axw,axh])
		
		# color map choice for abundances
		cmapa = cm.jet
		# color map choice for arrows
		cmapr = cm.autumn
		# if a value is below the lower limit its set to white
		cmapa.set_under(color='w')
		cmapr.set_under(color='w')
		# set value range for abundance colors (log10(Y))
		norma = colors.Normalize(vmin=-lbound,vmax=0)
		# set x- and y-axis scale aspect ratio to 1
		ax.set_aspect('equal')
		#print time,temp and density on top
		temp = ' '#'%8.3e' %ff['temp']
		time = ' '#'%8.3e' %ff['time']
		dens = ' '#'%8.3e' %ff['dens']
		
		#May cause Error, someting to do with tex
		'''
		#box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
		#	      + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
		anchored_box = AnchoredOffsetbox(loc=3,
				child=box1, pad=0.,
				frameon=False,
				bbox_to_anchor=(0., 1.02),
				bbox_transform=ax.transAxes,
				borderpad=0.,
				)
		ax.add_artist(anchored_box)
		'''
		## Colour bar plotted
		
		patches = []
		color = []
		
		for i in range(nzmax):
		    for j in range(nnmax):
		      if nzycheck[j,i,0]==1:
			xy = j-0.5,i-0.5
			
			rect = Rectangle(xy,1,1,)
			
			# abundance 
			yab = nzycheck[j,i,1]
			if yab == 0:
				
				yab=1e-99
				
			
			col =log10(yab)
		
			patches.append(rect)
			color.append(col)
		
		
		p = PatchCollection(patches, cmap=cmapa, norm=norma)
		p.set_array(array(color))
		p.set_zorder(1)
		ax.add_collection(p)
		cb = pl.colorbar(p)
		  
		# colorbar label
		cb.set_label('log$_{10}$(Y)')
		  
		# plot file name
		graphname = 'abundance-chart'+str(cycle)
		  
		# Add black frames for stable isotopes
		'''
		f = open('stable.dat')
		
		head = f.readline()
		stable = []
		'''
		if boxstable:
			for i in xrange(len(self.stable_el)):
				if i == 0:
					continue
				
				
				tmp = self.stable_el[i]
				try:
					zz= self.elements_names.index(tmp[0]) #charge
				except:
					continue
						
				for j in xrange(len(tmp)):
					if j == 0:
						continue
					
					nn = int(tmp[j]) #atomic mass
					nn=nn-zz
					
					xy = nn-0.5,zz-0.5
					rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=3.)
					rect.set_zorder(2)
					ax.add_patch(rect)
		
		
		
		
		# decide which array to take for label positions
		iarr = 0
		
		# plot element labels
		if ilabel:
		  for z in range(nzmax):
		    try:
		      nmin = min(argwhere(nzycheck[:,z,iarr]))[0]-1
		      ax.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
		    except ValueError:
		      continue
		      
		# plot mass numbers
		if imlabel:
		  for z in range(nzmax):
		     for n in range(nnmax):
			a = z+n
			if nzycheck[n,z,iarr]==1:
			  ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='small',clip_on=True)
		
		# plot lines at magic numbers
		if imagic:
		  ixymagic=[2, 8, 20, 28, 50, 82, 126]
		  nmagic = len(ixymagic)
		  for magic in ixymagic:
		    if magic<=nzmax:
		      try:
			xnmin = min(argwhere(nzycheck[:,magic,iarr]))[0]
			xnmax = max(argwhere(nzycheck[:,magic,iarr]))[0]
			line = ax.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
		      except ValueError:
			dummy=0
		    if magic<=nnmax:
		      try:
			yzmin = min(argwhere(nzycheck[magic,:,iarr]))[0]
			yzmax = max(argwhere(nzycheck[magic,:,iarr]))[0]
			line = ax.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
		      except ValueError:
			dummy=0
		
		# set axis limits
		if plotaxis==[0,0,0,0]:
		  
		  xmax=max(nin)
		  ymax=max(zin)
		  ax.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
		else:
		  ax.axis(plotaxis)
		
		# set x- and y-axis label
		ax.set_xlabel('neutron number')
		ax.set_ylabel('proton number')
		pl.title('Abu Chart for cycle '+str(int(cycle)))
		#fig.savefig(graphname)
		print graphname,'is done'
		if show:
			pl.show()
		if turnoff:
			ion()
		return
		
	def iso_abundMulti(self,cyclist, stable=False,amass_range=None,mass_range=None,
		ylim=[1e-13,10],shape='o',ref=-1,title=None,pdf=False):
		'''
		Method that plots figures and saves those figures to a .png file 
		(by default). Plots a figure for each cycle in the argument cycle
		input:
		 cycle       - a string/integer of the cycle of interest.
		 stable     - a boolean of whether to filter out the unstables.
		    		Defaults to False
		 amass_range -a 1x2 array containing the lower and upper Atomic
		    		 mass range. If not None this method will only 
		    		 plot the isomers with an atomic mass that
		    		 appears with in this range.
		 mass_range - a 1x2 array containing the lower and upper mass range.
		    		 If this is an instance of abu_vector this will 
		    		 only plot isotopes that have an atominc mass 
		    		 within this range. This will throw an error if
		    		 this range does not make sence ie [45,2]
			 	if None, it will plot over the entire range
				Defaults to None	     
		 ylim - A 1x2 array containing the lower and upper Y limits.
		    	   Defaults to 1e-13 and 10
		 ref  - reference cycle, If it is not -1, this method will 
		    plot the abundences of cycle devided by the .
		    	   The default is -1, it will do nothing
		 shape -The Shape of the dataplots, will default to circles
		 pdf - Boolean of if the output images will be pdf files
		 Title - Title of the plot and the image file
		 
		 NOTE, This method will achieve a 33% speed up and a 25% sppedup
		 when calling this method from python or ipython, rather than
		 ipython --pylab --q4thread
		'''
		print 'This method will achieve a 33% speedup and a 25% speedupwhen calling this method from python or ipython, rather than ipython --pylab --q4thread'
			
		
		for i in xrange(len(cyclist)):
			self.iso_abund(cyclist[i],stable,amass_range,mass_range,ylim,shape,ref,False)
			if title !=None:
				pl.title(title)
			else:
				name='IsoAbund'
			if not pdf:
				pl.savefig(name+str(cyclist[i])+'.png', dpi=400)
			else:
				pl.savefig(name+cyclist[i]+'.pdf', dpi=400)
			pl.clf()
		
		return None
		
	def read(self,fileName):
		'''
		Method for reading very simple data files for refrence cycles.
		The data files are just 3 columns, Charge, followed by 
		isotope name followed by the abundence, each column seperated by
		any number of spaces. Also there can be nospaces in the atomic 
		mass or abundence. An example would look like:
			1.0   H 1     0.728698   
			1.0   H 2     2.77841e-12
		Returns a dicionary with the keys, z, name, abu
		input:
		fileName: Fully specified path of filename
		'''
		z=[]
		name=[]
		abu=[]
		f=open(fileName,'r')
		lines=f.readlines()
		for i in xrange(len(lines)):
			lines[i]=lines[i].strip()
			lines[i]=lines[i].split()
		for i in xrange(len(lines)):
			z.append(lines[i][0])
			if len(lines[i])==4:
				name.append(lines[i][1]+' '+lines[i][2])
				abu.append(lines[i][3])
			elif lines[i][1]=='PROT':
				name.append('h 1')
				abu.append(lines[i][2])
			elif len(lines[i])==3:
				name.append(lines[i][1][0]+lines[i][1][1]+' '+lines[i][1][2]+lines[i][1][3]+lines[i][1][4])
				abu.append(lines[i][2])
			else:
				print 'Incompatable file, returning None'
				return None
		
		dat={'z':z, 'name':name,'abu':abu}
		return dat
		
	def iso_abund(self,cycle, stable=False,amass_range=None,mass_range=None,ylim=[1e-13,10],shape='o',ref=-1,show=True,given_abund_array=None):
		''' plot the abundance of all the chemical species
		inputs:
		    
		    cycle       - a string/integer of the cycle of interest.
		    		  It it is a list of cycles, this method will 
		    		  then do a plot for each of theses cycles and 
				  save them all to a file
		    stable     - a boolean of whether to filter out the unstables.
		    		Defaults to False
		    amass_range -a 1x2 array containing the lower and upper Atomic
		    		 mass range. If not None this method will only 
		    		 plot the isomers with an atomic mass that
		    		 appears with in this range.
		    mass_range - a 1x2 array containing the lower and upper mass range.
		    		 If this is an instance of abu_vector this will 
		    		 only plot isotopes that have an atominc mass 
		    		 within this range. This will throw an error if
		    		 this range does not make sense ie [45,2]
			 	if None, it will plot over the entire range
				Defaults to None	     
		    ylim - A 1x2 array containing the lower and upper Y limits.
		    	   Defaults to 1e-13 and 10
		    ref  - reference cycle, If it is not -1, this method will 
		    	   plot the abundances of cycle devided by the 
		    	   abundances of the refrence
		    	   Cycle. If any abundence in the reference cycle is zero, 
		    	   it will then interpret that abundence of 1e-99.
		    	   If ref is a string then this is interpreted as 
		    	   a filename.  It will then use the abundance column as
		    	   a reference cycle. The standards for reading this file 
		    	   can be found in self.read()'s docstring, ie in the 
		    	   method read() found above.
		    	   The default is -1, it will do nothing.
                           If ref=-2, than the array given_abund_array is plotted.
		    Shape -The Shape of the dataplots, will default to circles 
		    given_abund_array - array of abundances that can be given from script comparison_ppn, and if ref=-2 are plotted.
		    	  
		'''
		elem_list = []
		elem_index = []
		masses = []
		plotType=self.classTest()
		if str(cycle.__class__)=="<type 'list'>":
			self.iso_abundMulti(cycle, stable,amass_range,mass_range,ylim,shape,ref)
			return
			
		if str(ref.__class__)=="<type 'str'>":
			stringRef=True
			fileName=ref
			ref=-1
		else:
			stringRef=False
		if mass_range!=None and mass_range[0]>mass_range[1]:
			print 'Please input a proper mass range'
			print 'Returning None'
			return None
		if amass_range!=None and amass_range[0]>amass_range[1]:
			print 'Please input a proper Atomic mass range'
			print 'Returning None'
			return None
		if plotType=='se':
			isotope_to_plot = self.se.isotopes
			cycle=self.se.findCycle(cycle)
			z=self.se.Z #charge
			a=self.se.A #mass
			isomers=self.se.isomeric_states
			abunds = self.se.get(cycle,'iso_massf')
			if ref >-1:
				ref=self.se.findCycle(ref)
				abundsRef=self.se.get(ref,'iso_massf')
			
			if stringRef:
				abundsRef=self.read(fileName)
				tmpypsRef=zeros(len(z))
				for i in xrange(len(z)):
					for j in xrange(len(abundsRef['z'])):
						if isomers[i]==1 and z[i]==int(abundsRef['z'][j]) and a[i]==int(abundsRef['name'][j].split()[1]):
							tmpypsRef[i]=abundsRef['abu'][j] 
							break
						
				abundsRef=tmpypsRef
			masses = self.se.get(cycle,'mass')
			if mass_range == None:
			    print 'Using default mass range'
			    mass_range = [min(masses),max(masses)]    
			masses.sort()
			mass_range.sort()
			if amass_range != None:
				tmpA=[]
				tmpZ=[]
				tmpIso=[]
				tmpIsom=[]
				tmpyps=[]
				if ref >-1 or stringRef:
					tmpRef=[]
				
				for i in xrange(len(abunds)):
					tmpyps.append([])
					if ref >-1:
						tmpRef.append([])
				for i in xrange(len(a)):
					if (a[i] >=amass_range[0] and a[i]<=amass_range[1]):
						tmpA.append(a[i])
						tmpZ.append(z[i])
						tmpIso.append(isotope_to_plot[i])
						tmpIsom.append(isomers[i])
						for j in xrange(len(abunds)):
							tmpyps[j].append(abunds[j][i])
						if ref >-1 :
							for j in xrange(len(abundsRef)):
							
								tmpRef[j].append(abundsRef[j][i])
						if stringRef:
							tmpRef.append(abundsRef[i])
						
				isotope_to_plot=tmpIso
				z=tmpZ
				a=tmpA
				isomers=tmpIsom
				abunds=tmpyps
				if ref >-1 or stringRef:
					abundsRef=tmpRef
					
			
			isom=[]
			tmp=[]
			for i in range (len(isotope_to_plot)):
				if z[i]!=0 and isomers[i]==1: #if its not 'NEUt and not an isomer'
					tmp.append(self.elements_names[int(z[i])]+'-'+str(int(a[i])))
				elif isomers[i]!=1: #if it is an isomer
					isom.append(self.elements_names[int(z[i])]+'-'+str(int(a[i]))+'-'+str(int(isomers[i]-1)))	
			isotope_to_plot=tmp
			#tmp.sort(self.compar)
			#tmp.sort(self.comparator)			
			#for i in xrange(len(tmp)):
			#	isotope_to_plot.append(tmp[i][0])
			
				
			
			
		elif plotType=='PPN':
			isotope_to_plot = self.get('ISOTP', cycle)
			z=self.get('Z', cycle) #charge
			a=self.get('A', cycle) #mass
			isomers=self.get('ISOM', cycle)
			if ref > -2:
				yps=self.get('ABUNDNACE_MF', cycle)
			if ref == -2:
				yps=given_abund_array
			if ref >-1:
				ypsRef=self.get('ABUNDNACE_MF', ref)
			if stringRef:
				ypsRef=self.read(fileName )
				tmpypsRef=zeros(len(yps))
				for i in xrange(len(z)):
					for j in xrange(len(ypsRef['z'])):
						if isomers[i]==1 and z[i]==int(float((ypsRef['z'][j]))) and a[i]==int(float(ypsRef['name'][j].split()[1])):
							tmpypsRef[i]=ypsRef['abu'][j]
							break
							'''
						elif isomers[i]==1 and j ==len(ypsRef['z'])-1:
							print 'Reference Mismatch, Isotopes in cycle differ from isotepes in '+fileName
							print 'Returning None'
							return None
							'''
							
							
				ypsRef=tmpypsRef
				
			if amass_range != None:
				tmpA=[]
				tmpZ=[]
				tmpIso=[]
				tmpIsom=[]
				tmpyps=[]
				tmpypsRef=[]
				for i in xrange(len(a)):
					if (a[i] >=amass_range[0] and a[i]<=amass_range[1]):
						tmpA.append(a[i])
						tmpZ.append(z[i])
						tmpIso.append(isotope_to_plot[i])
						tmpIsom.append(isomers[i])
						tmpyps.append(yps[i])
						if ref >-1 or stringRef:
							tmpypsRef.append(ypsRef[i])
				isotope_to_plot=tmpIso
				z=tmpZ
				a=tmpA
				isomers=tmpIsom
				yps=tmpyps
				ypsRef=tmpypsRef
				
			if ref >-1 and len(yps)!=len(ypsRef):
				print 'Refrence Cycle mismatch, Aborting plot'
				return None
			
				
			if ref >-1 or stringRef:
				for i in xrange(len(yps)):
					if ypsRef[i]!=0:
						yps[i]=yps[i]/ypsRef[i]
					else:
						
						yps[i]=yps[i]/1e-99
			tmp1=[]
			tmp=[]
			isom=[]
			
			for i in range (len(isotope_to_plot)):
				if z[i]!=0 and isomers[i]==1: #if its not 'NEUt and not an isomer'
					tmp.append([self.elements_names[int(z[i])]+'-'+str(int(a[i])),yps[i]])
				elif isomers[i]!=1: #if it is an isomer
					
					if yps[i]==0:
						
						isom.append([self.elements_names[int(z[i])]+'-'+str(int(a[i]))+'-'+str(int(isomers[i]-1)),1e-99])
					else:
						isom.append([self.elements_names[int(z[i])]+'-'+str(int(a[i]))+'-'+str(int(isomers[i]-1)),yps[i]])	
					
			
			tmp.sort(self.compar)
			tmp.sort(self.comparator)
			
			abunds=[]
			isotope_to_plot=[]
			for i in xrange(len(tmp)):
				isotope_to_plot.append(tmp[i][0])
				abunds.append(tmp[i][1])
			
		else:
			print 'This method, iso_abund, is not supported by this class'
			print 'Returning None'
			return None
		#    Check the inputs
		#if not self.se.cycles.count(str(cycle)):
		#    print 'You entered an cycle that doesn\'t exist in this dataset:', cycle
		#    print 'I will try and correct your format.'
		
		#    print cyc_len, len(str(cycle))
		#
		#    while len(str(cycle)) < cyc_len:
		#        cycle = '0'+str(cycle)
		#        print cycle
		
		#    if not self.se.cycles.count(str(cycle)):
		#        print 'I was unable to correct your cycle.  Please check that it exists in your dataset.'
		
		print 'Using The following conditions:'
		if mass_range != None:
			print '\tmass_range:', mass_range[0], mass_range[1]
		if amass_range != None:
			print '\tAtomic mass_range:', amass_range[0], amass_range[1]
		print '\tcycle:', cycle
		print '\tplot only stable:',stable
		
		
		elem_list = ['' for x in xrange(len(self.elements_names)) ]
		elem_index = [[-1] for x in xrange(len(self.elements_names))]
		#elem_list = self.elements_names
		
		for elem in elem_index:
		    elem.remove(-1)
		for x in xrange(len(isotope_to_plot)):
		    
		    if self.elements_names.count(isotope_to_plot[x].split('-')[0]):
			#print self.se.isotopes[x].split('-')[0], x
			elem_index[self.elements_names.index(isotope_to_plot[x].split('-')[0])].append(x)
			if elem_list.count(isotope_to_plot[x].split('-')[0]) == 0:
				elem_list[self.elements_names.index(isotope_to_plot[x].split('-')[0])] += \
								     isotope_to_plot[x].split('-')[0]
		
		
		for x in xrange(len(elem_index)):
		    numbers = []
		    for y in xrange(len(elem_index[x])):
			numbers.append(isotope_to_plot[elem_index[x][y]][1])
	
		    list3 = zip(elem_index[x], numbers)
		    list3.sort()
		    list3.sort()
		    for y in xrange(len(list3)):
			elem_index[x][y] = list3[y][0]        
		
		
		
	       
		
		index = xrange(len(isotope_to_plot))
		if stable:
		    
		    #    loop through the index list
		    for t in xrange(15):
			removed = False
			#print 'new step'
			for ind in index:    # and the indice of isotope
			    try:
				for iso in elem_index[ind]:
				    for elem in self.stable_el:    
					#    Loop through the list of stable isotopes
					
					if isotope_to_plot[iso].split('-')[0] == elem[0]:    
					    #    Is it the same element?                                    
					    matched = elem.count(int(isotope_to_plot[iso].split('-')[1]))    
					    #    Does the stable el list contain he isotope
					    #print matched, self.se.isotopes[iso].split('-')[0], elem[0]
					    if matched == 0:
						try:
						    #print 'removing', self.elem_index[ind], iso
						    #print 'removed', elem_index[ind], iso
						    elem_index[ind].remove(iso)        
						    removed = True
						except:
						    print 'failed to remove: ', isotope_to_plot[iso]
			    except IndexError:
				removed = True
				#print iso, ind
			if removed == False:
			    #print 'done popping'
			    break                
	
		
		#print elem_list
		#print elem_index
	
		abund_plot = []
		mass_num = []
		
		
		if plotType=='se':
			
			while len(abunds) == 1:
				abunds = abunds[0]
			for j in xrange(len(index)):    #    Loop through the elements
			    temp = []
		
			    try:
				x =  masses[0]
				del x
			    except IndexError:
				masses = [masses]
				
				
			    try:        
				for k in xrange(len(elem_index[index[j]])):    #    Loop through the isotopes
				    abundance = 0
				    
				    #print 'elemIndexKs',  self.elem_index[index[j]][k]                
				    try:
					for l in xrange(len(abunds)):
					    #print mass_range[0], masses[l]  , mass_range[1], masses[l]
					    try:
						if mass_range[0] <=  masses[l]  and mass_range[1] >=  masses[l] :    
						    #    Only collect data in between ranges.        
						    #print abunds[i][l][self.elem_index[index[j]][k]]
						    #if abunds[i][l][self.elem_index[index[j]][k]] < 1e-20:
						    #    abundance += 1e-20
						
						    try:
							abundance += (abunds[l][elem_index[index[j]][k]])*abs(masses[l+1]-masses[l])
						    except IndexError:    #  The last step requires us to interpolate to the next highest step
							abundance += (abunds[l][elem_index[index[j]][k]])*abs(round(masses[l])-masses[l])
							#print abundance
						    
						    	    
					    except IndexError:
						None#print 'end of the line'
					#print abundance/abs(ranges[1]-ranges[0])
					if abundance ==0:
						    abundance = 1e-99
					temp.append(abundance/abs(mass_range[1]-mass_range[0]))
					#print abundance,temp
				    except AttributeError:
					print i
				#mass_num.append(temp)
			    except IndexError:
				None
				#print j, index[j], elem_index
			    #time_step.append(temp)
			    #print time_step
			    abund_plot.append(temp)
			
			if ref > -1:
				abund_plotRef=[]	
				while len(abundsRef) == 1:
					abundsRef = abundsRef[0]
				for j in xrange(len(index)):    #    Loop through the elements
				    temp = []
			
				    try:
					x =  masses[0]
					del x
				    except IndexError:
					masses = [masses]
					
					
				    try:        
					for k in xrange(len(elem_index[index[j]])):    #    Loop through the isotopes
					    abundance = 0
					    
					    #print 'elemIndexKs',  self.elem_index[index[j]][k]                
					    try:
						for l in xrange(len(abundsRef)):
						    #print mass_range[0], masses[l]  , mass_range[1], masses[l]
						    try:
							if mass_range[0] <=  masses[l]  and mass_range[1] >=  masses[l] :    
							    #    Only collect data in between ranges.        
							    #print abunds[i][l][self.elem_index[index[j]][k]]
							    #if abunds[i][l][self.elem_index[index[j]][k]] < 1e-20:
							    #    abundance += 1e-20
							
							    try:
								abundance += (abundsRef[l][elem_index[index[j]][k]])*abs(masses[l+1]-masses[l])
							    except IndexError:    #  The last step requires us to interpolate to the next highest step
								abundance += (abundsRef[l][elem_index[index[j]][k]])*abs(round(masses[l])-masses[l])
								#print abundance
						    except IndexError:
							None#print 'end of the line'
						#print abundance/abs(ranges[1]-ranges[0])
						if abundance ==0:
						    abundance = 1e-99
						temp.append(abundance/abs(mass_range[1]-mass_range[0]))
						#print abundance,temp
					    except AttributeError:
						print i
					#mass_num.append(temp)
				    except IndexError:
					None
					#print j, index[j], elem_index
				    #time_step.append(temp)
				    #print time_step
				    abund_plotRef.append(temp)
				    
				if len(abund_plot)!=len(abund_plotRef):
					print 'Error 1'
				for i in xrange(len(abund_plot)):
					if len(abund_plot[i])!=len(abund_plotRef[i]):
						print 'Error 2'
				    	    	print abund_plot[i],abund_plotRef[i]
				    	    
				    	for j in xrange(len(abund_plot[i])):
						
				    	    	if abund_plotRef[i][j]!=0:
				    	    	    	abund_plot[i][j]=abund_plot[i][j]/abund_plotRef[i][j]
				    	    	else:
				    	    		
				    	    	    	abund_plot[i][j]=abund_plot[i][j]/1e-99
			elif stringRef:
				for i in xrange(len(elem_index)):
					for j in xrange(len(elem_index[i])):
						if abundsRef[elem_index[i][j]]!=0:
							abund_plot[i][j]=abund_plot[i][j]/abundsRef[elem_index[i][j]]
						else:
							abund_plot[i][j]=abund_plot[i][j]/1e-99
			
	        elif plotType=='PPN':
	        	
	        	for i in xrange(len(elem_list)):
	        		tmp=[]
	        		yps=abunds
	        		if elem_list[i]!='':
	        			
	        			for j in xrange(len(elem_index[i])):
	        				
	        				#print self.get(isotope_to_plot[index], cycle)[3]
	        				tmp.append(yps[elem_index[i][j]])
	        				
	        		abund_plot.append(tmp)
		#print cycle
		
		#temp3 = []
		mass_num = []
		
		
		#for j in xrange(len(index)):
		for j in xrange(len(elem_index)):
		    
		    temp = []
		    #temp2 = []    
		    try:
			
			for k in xrange(len(elem_index[j])):
				
				temp.append(float(isotope_to_plot[elem_index[j][k]].split('-')[1]))
	
			mass_num.append(temp)
		    #temp3.append(temp2)
		    except IndexError:
			None
			#print len(elem_index), index[j],j
		#ref_mass = temp3
		#del temp2             
		
		#23
		plot_type = ['-','--','-.']
		pl_index = 0
		#6
		#colors = ['g','r','c','m','k']
		colors = ['k','k','k','k','k']
		cl_index = 0
		
		l1 = []
		l2 = []        
		for j in xrange(len(abund_plot)):        #Loop through the elements of interest
		    #    Process the line
		    #print 'processing line'
		    for l in xrange(len(abund_plot[j])):
                        #print mass_num[j][l] 	
                        print abund_plot[j][l] 	 		
			if abund_plot[j][l] == 0:
			    abund_plot[j][l] = 1e-99
			    
			    
		    
		    try:
			l1.append(pl.semilogy(mass_num[j],abund_plot[j],str(colors[cl_index]+plot_type[pl_index])))
			cl_index+=1
			pl_index+=1
			if pl_index > 2:
			    pl_index = 0
			if cl_index > 4:
			    cl_index = 0
		    except IndexError:
			None
			#print 'div by zero'
		    
			
		    try:
		    	
		    	tmpList=[0,0]
		    	for i in xrange(len(abund_plot[j])):
		    		
		    		coordinates=[mass_num[j][i],abund_plot[j][i]]
		    		
		    		
				if coordinates[1]<=ylim[1] and coordinates[1]>=ylim[0]:
					if coordinates[1]>tmpList[1]:
						if amass_range !=None:
							if coordinates[0]>=amass_range[0]-.5 and coordinates[0]<=amass_range[1]+.5:
								tmpList=coordinates
								
						elif plotType=='PPN':
							tmpList=coordinates
							
						elif plotType!='PPN':
							tmpList=coordinates
							
				#print self.parent.elem_list[self.index[j]],coordinates[0],coordinates[1]
			coordinates=tmpList
			
			if coordinates != [0,0]:
				
				pl.text(coordinates[0],coordinates[1], elem_list[j])
		
		    except ValueError:
			None
			#print 'Empty var:  ', abund_plot[j]
		    except IndexError:
			None
			#print 'out of bounds: ', len(abund_plot), j
			    
		    try:
			pl.semilogy(mass_num[j],abund_plot[j],'k'+shape)
		    except OverflowError:
			None
			#print 'div by zero', len(mass_num[j]), len(abund_plot[j])
		    except IndexError:
			None
			#print 'out of bounds II: ', len(mass_num),len(abund_plot),j
		pl_index = 0
		
		cl_index = 0
		
		for j in xrange(len(isom)):
			#pl.semilogy(isom[j][0].split('-')[1],isom[j][1],'r'+shape)#,str(colors[cl_index]+plot_type[pl_index])))
			cl_index+=1
			pl_index+=1
			if pl_index > 2:
			    pl_index = 0
			if cl_index > 4:
			    cl_index = 0
			coordinates=[int(isom[j][0].split('-')[1]),isom[j][1]]
			name=isom[j][0]

			#pl.text(coordinates[0],coordinates[1], name.split('-')[0]+'m'+name.split('-')[2])
		if plotType=='se':
			title = str('Abundance of Isotopes over range %4.2f' %mass_range[0]) + str('-%4.2f' %mass_range[1]) +\
				str(' for cycle %d' %int(cycle))
			if amass_range !=None:
				pl.xlim([amass_range[0]-.5,amass_range[1]+.5])
		else:
			
			if amass_range ==None:
				title = str('Abundance of Isotopes for Cycle '+str(cycle))
			else:
				#title = str('Abundance of Isotopes with A between '+str(amass_range[0])+' and '+str(amass_range[1])+' for Cycle '+str(cycle))
				title = str('')
				pl.xlim([amass_range[0]-.5,amass_range[1]+.5])
			
		pl.ylim(ylim)
		pl.title(title)
		pl.xlabel('Mass Number')
		if ref>-1:
			#pl.ylabel('Relative Abundance / Refrence Abundance of Cycle '+str(int(ref)))
		        pl.ylabel('Relative Abundance') 
		elif stringRef:
			pl.ylabel('Relative Abundance / Refrence Abundance of '+str(fileName))
		else:
			pl.ylabel('Relative Abundance')
		pl.grid()
		if show:
			pl.show()
		return
	
	def plotprofMulti(self,ini,end,delta,what_specie,xlim1,xlim2,ylim1,ylim2):
	
		''' create a movie with mass fractions vs mass coordinate 
		between xlim1 and xlim2, ylim1 and ylim2. Only works with instances of se 
	
		ini          - initial model
		end          - final model 
		delta        - sparsity factor of the frames
		what_specie  - array with species in the plot
		xlim1, xlim2 - mass coordinate range 
		ylim1, ylim2 - mass fraction coordinate range
	  
		'''
		plotType=self.classTest()
		if plotType=='se':
			for i in range(ini,end+1,delta):
			    step = int(i)
			    print step
			    for j in range(len(what_specie)):
				self.plot_prof_1(step,what_specie[j],xlim1,xlim2,ylim1,ylim2)
			    #          
			    filename = str('%03d' % step)+'_test.png'             
			    pl.savefig(filename, dpi=400) 
			    print 'wrote file ', filename
			    #
			    pl.clf()
			
		else:
			print 'This method is not supported for '+str(self.__class__)
			return
	# From mesa_profile
    	def plot_prof_1(self,species,keystring,xlim1,xlim2,ylim1,ylim2, show=False):
		''' plot one species for cycle between xlim1 and xlim2 
		    Only works with instances of se and mesa _profile
		
		species      - which species to plot
		keystring    - label that appears in the plot or in the cas of se,
			       A cycle or list of cycles
		xlim1, xlim2 - mass coordinate range                                                 
		ylim1, ylim2 - mass fraction coordinate range '''
		plotType=self.classTest()
		if plotType=='se':
			tot_mass=self.se.get(keystring,'total_mass')    
			age=self.se.get(keystring,'age')    
			mass=self.se.get(keystring,'mass')    
			Xspecies=self.se.get(keystring,'iso_massf',species)
			
			mod=keystring
		elif plotType=='mesa_profile':
			tot_mass=self.header_attr['star_mass'] 
			age=self.header_attr['star_age'] 
			mass=self.get('mass')
			mod=self.header_attr['model_number']
			Xspecies=self.get(species)
		else:
			print 'This method is not supported for '+str(self.__class__)
			return
		x,y=self.logarithm(Xspecies,mass,True,False,10)
		print x
		pl.plot(y,x,'-',label=str(species))
		pl.xlim(xlim1,xlim2)
		pl.ylim(ylim1,ylim2)
		pl.legend()
	
		pl.xlabel('$Mass$ $coordinate$', fontsize=20)
		pl.ylabel('$X_{i}$', fontsize=20)
		pl.title('Mass='+str(tot_mass)+', Time='+str(age)+' years, cycle='+str(mod))
		if show:
			pl.show()
	
	# From mesa.star_log
	
def flux_chart(file_name,plotaxis,plot_type):
	'''
	Plots a chart with fluxes
	input:
	file_name: name of the file of fluxes we are looking at.
        plotaxis: [xmin,xmax,ymin,ymax], where on x axis there is neutron number and on y axis there is Z. 
        plot_types: 0 for standard flux plot, 1 if fluxes focused on one specie.
        Note: the script is terribly slow, need to be improved. For now I put here in data_plot:
        [1]: import data_plot
        [2]: data_plot.flux_chart('file_name',[xmin,xmax,ymin,ymax],int)
	The pdf is created, but an error bumped up and the gui is empty. To avoid this, I had to set 'text.usetex': False. See below.
        Also, for the same reason no label in x axys is written using 'text.usetex': True.  
        Note also that the GUI works really slow with this plot. so, we need to optimize from the graphic point of view.
        This need to be included in ppn.py I think, and set in multi option too, in case we want to read more flux files at the same time. 
        Finally, you need to have stable.dat to read in to make it work....
        '''

	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.mpl import colors,cm
	from matplotlib.patches import Rectangle, Arrow
	from matplotlib.collections import PatchCollection
	from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
	import sys

    	f = open(file_name)
    	lines = f.readline()
        lines = f.readlines()
    	f.close()
        
        # starting point of arrow
        coord_x_1 = []
        coord_y_1 = []
        # ending point of arrow (option 1)
        coord_x_2 = []
        coord_y_2 = []
        # ending point of arrow (option 2) 
        coord_x_3 = []
        coord_y_3 = []
        # fluxes
        flux_read = []
        flux_log10 = []
 
        single_line = [] 
        for i in range(len(lines)):
                single_line.append(lines[i].split())
 		coord_y_1.append(float(single_line[i][1]))
                coord_x_1.append(float(single_line[i][2])-coord_y_1[i])
 		coord_y_2.append(float(single_line[i][5]))
                coord_x_2.append(float(single_line[i][6])-coord_y_2[i])
 		coord_y_3.append(float(single_line[i][7]))
                coord_x_3.append(float(single_line[i][8])-coord_y_3[i])
                flux_read.append(float(single_line[i][9]))
                flux_log10.append(np.log10(flux_read[i]+1.0e-99))    
    
        #print coord_x_1[72],coord_y_1[72] 
    	print 'file read!'

        if plot_type == 1:
		I_am_the_target = [23.-11.,11.]
                #I_am_the_target = [13.-7.,7.]
		#I_am_the_target = [22.-10.,10.]
                # = [79.-34.,34.]
        # here below need for plotting 
	# plotaxis = [xmin,xmax,ymin,ymax] 
	#plotaxis=[1,20,1,20]
	#plotaxis=[0,0,0,0]

	# elemental labels off/on [0/1]
	ilabel = 1

	# label for isotopic masses off/on [0/1]
	imlabel = 1

	# turn lines for magic numbers off/on [0/1]
	imagic = 0

	# flow is plotted over "prange" dex. If flow < maxflow-prange it is not plotted
	prange = 8.

        ############################################# 

        max_flux = max(flux_log10)
        
  	nzmax = int(max(max(coord_y_1),max(coord_y_2),max(coord_y_3)))+1
  	nnmax = int(max(max(coord_x_1),max(coord_x_2),max(coord_x_3)))+1

        nzycheck = np.zeros([nnmax,nzmax,3])   
        coord_x_out = np.zeros(len(coord_x_2))
        coord_y_out = np.zeros(len(coord_y_2))
        for i in range(len(flux_log10)):
  		nzycheck[coord_x_1[i],coord_y_1[i],0] = 1
  		nzycheck[coord_x_1[i],coord_y_1[i],1] = flux_log10[i]
  		if coord_x_2[i] >= coord_x_3[i]:
                        coord_x_out[i] = coord_x_2[i] 
                        coord_y_out[i] = coord_y_2[i] 
    			nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
    			nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10[i]
                elif coord_x_2[i] < coord_x_3[i]:               
                        coord_x_out[i] = coord_x_3[i] 
                        coord_y_out[i] = coord_y_3[i] 
    			nzycheck[coord_x_out[i],coord_y_out[i],0] = 1
    			nzycheck[coord_x_out[i],coord_y_out[i],1] = flux_log10[i]
    		if flux_log10[i]>max_flux-prange:
      			nzycheck[coord_x_1[i],coord_y_1[i],2] = 1
      			nzycheck[coord_x_out[i],coord_y_out[i],2] = 1

	#######################################################################
	# elemental names: elname(i) is the name of element with Z=i 
	elname=	('none','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe',
        'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb',
        'Te', 'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
        'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu')


	#### create plot

	## define axis and plot style (colormap, size, fontsize etc.)
	if plotaxis==[0,0,0,0]:
  		xdim=10
  		ydim=6
	else:
  		dx = plotaxis[1]-plotaxis[0]
  		dy = plotaxis[3]-plotaxis[2]
  		ydim = 6
  		xdim = ydim*dx/dy
  
	format = 'pdf'
        # note that I had to set 'text.usetex': False, to avoid Exception in Tkinter callback.
        # and to make the GUI work properly. Why? some missing package?
	params = {'axes.labelsize':  15,
          'text.fontsize':   15,
          'legend.fontsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'text.usetex': False}
	plt.rcParams.update(params)
	fig=plt.figure(figsize=(xdim,ydim),dpi=100)
	axx = 0.10
	axy = 0.10
	axw = 0.85
	axh = 0.8
	ax=plt.axes([axx,axy,axw,axh])

	# color map choice for abundances
	cmapa = cm.jet
	# color map choice for arrows
	cmapr = cm.autumn
	# if a value is below the lower limit its set to white
	cmapa.set_under(color='w')
	cmapr.set_under(color='w')
	# set value range for abundance colors (log10(Y))
	norma = colors.Normalize(vmin=-20,vmax=0)
	# set x- and y-axis scale aspect ratio to 1
	ax.set_aspect('equal')
	#print time,temp and density on top
	#temp = '%8.3e' %ff['temp']
	#time = '%8.3e' %ff['time']
	#dens = '%8.3e' %ff['dens']

	#box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
        #      + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
	#anchored_box = AnchoredOffsetbox(loc=3,
        #        child=box1, pad=0.,
        #        frameon=False,
        #        bbox_to_anchor=(0., 1.02),
        #        bbox_transform=ax.transAxes,
        #        borderpad=0.,
        #        )
	#ax.add_artist(anchored_box)

	# Add black frames for stable isotopes
	f = open('stable.dat')

	head = f.readline()
	stable = []

	for line in f.readlines():
  		tmp = line.split()
  		zz = int(tmp[2])
  		nn = int(tmp[3])
  		xy = nn-0.5,zz-0.5
  		rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=3.)
  		rect.set_zorder(2)
  		ax.add_patch(rect)

 	apatches = []
  	acolor = []
  	m = 0.8/prange
  	vmax=np.ceil(max(flux_log10))
  	vmin=max(flux_log10)-prange
  	b=-vmin*m+0.1
  	normr = colors.Normalize(vmin=vmin,vmax=vmax)
  	ymax=0.
  	xmax=0.

  	for i in range(len(flux_log10)):
    		x = coord_x_1[i]
    		y = coord_y_1[i]
    		dx = coord_x_out[i]-coord_x_1[i]
    		dy = coord_y_out[i]-coord_y_1[i]                
                if plot_type == 0: 
    			if flux_log10[i]>=vmin:
      				arrowwidth = flux_log10[i]*m+b
      				arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				if xmax<x:
        				xmax=x
      				if ymax<y:
        				ymax=y
      				acol = flux_log10[i]
      				apatches.append(arrow)
      				acolor.append(acol)
    		elif plot_type == 1:
			if x==I_am_the_target[0] and y==I_am_the_target[1] and flux_log10[i]>=vmin:
      				arrowwidth = flux_log10[i]*m+b
      				arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				if xmax<x:
        				xmax=x
      				if ymax<y:
        				ymax=y
      				acol = flux_log10[i]
      				apatches.append(arrow)
      				acolor.append(acol)
			if x+dx==I_am_the_target[0] and y+dy==I_am_the_target[1] and flux_log10[i]>=vmin:
      				arrowwidth = flux_log10[i]*m+b
      				arrow = Arrow(x,y,dx,dy, width=arrowwidth)
     				if xmax<x:
        				xmax=x
      				if ymax<y:
        				ymax=y
      				acol = flux_log10[i]
      				apatches.append(arrow)
      				acolor.append(acol)



        	xy = x-0.5,y-0.5
        	rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
        	rect.set_zorder(2)
        	ax.add_patch(rect)
        	xy = x+dx-0.5,y+dy-0.5
        	rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=1.)
        	rect.set_zorder(2)
        	ax.add_patch(rect)

  
  	a = PatchCollection(apatches, cmap=cmapr, norm=normr)
  	a.set_array(np.array(acolor))
  	a.set_zorder(3)
  	ax.add_collection(a)
  	cb = plt.colorbar(a)

  	# colorbar label
  	cb.set_label('log$_{10}$(f)')
  
  	# plot file name
  	graphname = 'flow-chart.'+format

	# decide which array to take for label positions
	iarr = 2

	# plot element labels
  	for z in range(nzmax):
    		try:
      			nmin = min(np.argwhere(nzycheck[:,z,iarr-2]))[0]-1
      			ax.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
    		except ValueError:
      			continue
      
	# plot mass numbers
	if imlabel==1:
  		for z in range(nzmax):
     			for n in range(nnmax):
        			a = z+n
				if nzycheck[n,z,iarr-2]==1:
          				ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='small',clip_on=True)

	# plot lines at magic numbers
	if imagic==1:
  		ixymagic=[2, 8, 20, 28, 50, 82, 126]
  		nmagic = len(ixymagic)
  		for magic in ixymagic:
    			if magic<=nzmax:
      				try:
        				xnmin = min(np.argwhere(nzycheck[:,magic,iarr-2]))[0]
        				xnmax = max(np.argwhere(nzycheck[:,magic,iarr-2]))[0]
        				line = ax.plot([xnmin,xnmax],[magic,magic],lw=3.,color='r',ls='-')
      				except ValueError:
        				dummy=0
    			if magic<=nnmax:
      				try:
        				yzmin = min(np.argwhere(nzycheck[magic,:,iarr-2]))[0]
        				yzmax = max(np.argwhere(nzycheck[magic,:,iarr-2]))[0]
        				line = ax.plot([magic,magic],[yzmin,yzmax],lw=3.,color='r',ls='-')
      				except ValueError:
        				dummy=0

	# set axis limits
	if plotaxis==[0,0,0,0]:
		ax.axis([-0.5,xmax+0.5,-0.5,ymax+0.5])
	else:
  		ax.axis(plotaxis)

	# set x- and y-axis label
	ax.set_xlabel('neutron number')
	ax.set_ylabel('proton number')

	fig.savefig(graphname)
	print graphname,'is done'

	plt.show()
