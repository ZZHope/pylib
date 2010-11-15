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
import os
import threading

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
		if 'ppm.YProfile' == c:
			tmp='YProfile'
		elif 'ascii_table.AsciiTable' == c:
			tmp='AsciiTable'
		elif 'mppnp.h5TPlotTools' == c:
			tmp='h5TPlotTools'
		elif 'mesa.mesa_profile' == c:
			tmp='mesa_profile'
		elif 'mesa.star_log' == c:
			tmp='mesa.star_log'
		elif 'ppn.xtime' == c:
			tmp='xtime'
		elif 'ppn.ppn_profile' == c:
			tmp='PPN'
		
		return tmp
	
	def logarithm(self,tmpX,tmpY,logX,logY,base):
		logXER=False
		logYER=False
		for i in range(len(tmpX)):
				if tmpX[i]<=0 and logX:
					print 'We can not log a number less than or equal to zero'
					print 'Attempting to remove incompatible values from X'
					logXER=True
				if tmpY[i]<=0 and logY:
					print 'We can not log a number less than or equal to zero'
					print 'Attempting to remove incompatible values from Y'
					logYER=True
		tmX=[]
		tmY=[]
		
		if logXER:
			for i in range(len(tmpX)):
				if tmpX[i]>0:
					tmX.append( tmpX[i])
					tmY.append(tmpY[i])
			tmpX=tmX
			tmpY=tmY
		elif logYER:
			for i in range(len(tmpY)):
				if tmpY[i]>0:
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
	
	def plot(self,atriX,atriY, FName=None,numType='ndump',legend=None,labelX=None, labelY=None ,
		indexX=None, indexY=None, title=None, shape='.',logX=False, logY=False, base=10):
		"""
		Simple function that plots atriY as a function of atriX
		This method will automatically find and plot the requested data.
		Input:
		atriX, The name of the attribute you want on the x axis
		atriY, The name of the attribute you want on the Y axis
		Fname: Be the filename, Ndump or time, Defaults to the last NDump
		numType: designates how this function acts and how it interprets FName
			 Defaults to file
		if numType is 'file', this function will get the desird attribute from that file
		if numType is 'NDump' function will look at the cycle with that nDump
		if numType is 't' or 'time' function will find the _cycle with the closest time stamp 
		labelX: The label on the X axis
		labelY: The label on the Y axis
		indexX: Depreciated: If the get method returns a list of lists, indexX
			would be the list at the index indexX in the list.
		indexY: Depreciated: If the get method returns a list of lists, indexY
			would be the list at the index indexX in the list.
		shape: What shape and colour the user would like their plot in.
		       Please see http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
		       for all possable choices
		title: The Title of the Graph
		logX: A boolean of weather the user wants the x axis logarithmically
		logY: A boolean of weather the user wants the Y axis logarithmically
		base: The base of the logarithm. Default = 10
		WARNING: Unstable if get returns a list with only one element (x=[0])
		"""
		
		tmpX=[]
		tmpY=[]
		singleX=False
		singleY=False
		#Getting data
		plotType=self.classTest()
		if plotType=='YProfile':
			if FName==None:
				FName=len(self.files)-1
		
			listY=self.get(atriY,FName, numType,resolution='a')
			listX=self.get(atriX,FName, numType,resolution='a')
		elif plotType=='h5TPlotTools':
			print 'Not supported yet for h5TPlotTools'
			return None
		else :
			listY=self.get(atriY)
			listX=self.get(atriX)
		
		
		
			
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
		elif indexX==None and len(listX)>1: # If it is a list of lists of values.
						    # take the largest
			tmpX=listX[0]
			for i in range(len(listX)):
				if len(tmpX)<len(listX[i]):
					tmpX=listX[i]
		elif indexX<len(listX): # If an index is specified, use that index
			tmpX=listX[indexX]
		else:
			print 'Sorry that indexX does not exist, returning None'
			return None
								
		'''
		elif indexX==None and len(listX)==1:
			tmpX=listX
		'''
		
		
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
		elif indexY==None and len(listY)>1:# If it is a list of lists of values.
						    # take the largest
			#print 'fourth'
			tmpY=listY[0]
			for i in range(len(listY)):	
				if len(tmpY)<len(listY[i]):
					tmpY=listY[i]
		elif indexY<len(listY): # If an index is specified, use that index
			#print 'sixth'
			tmpY=listY[indexY]
		else:
			print 'Sorry that indexY does not exist, returning None'
			return None
		'''					
		elif indexY==None and len(listY)==1:
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
					#print listY[i]
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
		
		# Logarithm stuff
		if logY or logX:	
			tmpX,tmpY=self.logarithm(tmpX,tmpY,logX,logY,base)
		
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

		if labelX== None :
			labelX=atriX
		if labelY== None :
			labelY=atriY
		
		if logX:
			labelX='log '+labelX
		if logY:
			labelY='log '+labelY
		
		if title!=None:
			title=title
		else:
			title=labelY+' vs '+labelX
		
		if legend!=None:
			legend=legend
		else:
			legend=labelY+' vs '+labelX	
			
		pyl.plot(listX,listY,shape,label=legend)
		pyl.legend()
		pyl.title(title)
		pyl.xlabel(labelX)
		pyl.ylabel(labelY)
		pyl.show()
		
		
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
	
	# From mesa_profile
	
		
	
    	def plot_prof_1(self,species,keystring,xlim1,xlim2,ylim1,ylim2):
	
		''' plot one species for cycle between xlim1 and xlim2                               
		
		species      - which species to plot
		keystring    - label that appears in the plot
		xlim1, xlim2 - mass coordinate range                                                 
		ylim1, ylim2 - mass fraction coordinate range '''
		plotType=self.classTest()
		if plotType=='h5TPlotTools':
			tot_mass=self.se.get(keystring,'total_mass')    
			age=self.se.get(keystring,'age')    
			mass=self.se.get(keystring,'mass')    
			Xspecies=self.se.get(keystring,'iso_massf',species)
			mod=keystring
		else:
			tot_mass=self.header_attr['star_mass'] 
			age=self.header_attr['star_age'] 
			mass=self.get('mass')
			mod=self.header_attr['model_number']
			Xspecies=self.get(species)
		pyl.plot(mass,log10(Xspecies),'-',label=keystring)
		pyl.xlim(xlim1,xlim2)
		pyl.ylim(ylim1,ylim2)
		pyl.legend()
	
		pl.xlabel('$Mass$ $coordinate$', fontsize=20)
		pl.ylabel('$X_{i}$', fontsize=20)
		pl.title('Mass='+str(tot_mass)+', Time='+str(age)+' years, cycle='+str(mod))
	
	# From mesa.star_log
	def hrd(self):
		''' plot an HR diagram '''
	
		pyl.plot(self.data[:,self.cols['log_Teff']-1],self.data[:,self.cols['log_L']-1],label = "M="+str(self.header_attr['initial_mass'])+", Z="+str(self.header_attr['initial_z']))
		pyl.legend()
		pyl.xlabel('log Teff')
		pyl.ylabel('log L')
    
    	def hrd_key(self,key_str):
		''' plot an HR diagram 
		
		key_str    a label string'''
	
		pyl.plot(self.data[:,self.cols['log_Teff']-1],self.data[:,self.cols['log_L']-1],label = key_str)
		pyl.legend()
		pyl.xlabel('log Teff')
		pyl.ylabel('log L')
    
    	def kippenhahn(self,num_frame,xax):
		''' Kippenhahn plot as a function of time or model
		
		num_frame    number of frame to plot this plot into
		xax          string that is either model or time to indicate what is to 
			     be used on the x-axis'''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction.'+\
			  ' needs to be "time" or "model"'
		
	    
		h1_boundary_mass  = self.get('h1_boundary_mass')
		he4_boundary_mass = self.get('he4_boundary_mass')
		star_mass         = self.get('star_mass')
		mx1_bot           = self.get('conv_mx1_bot')*star_mass
		mx1_top           = self.get('conv_mx1_top')*star_mass
		mx2_bot           = self.get('conv_mx2_bot')*star_mass
		mx2_top           = self.get('conv_mx2_top')*star_mass
		surface_c12       = self.get('surface_c12')
		surface_o16       = self.get('surface_o16')
	
		COratio=(surface_c12*4.)/(surface_o16*3.)
	
		pyl.plot(xaxisarray,COratio,'-k',label='CO ratio')
		pyl.ylabel('C/O ratio')
		pyl.legend(loc=4)
	
	
		pyl.twinx()
		pyl.plot(xaxisarray,h1_boundary_mass,label='h1_boundary_mass')
		pyl.plot(xaxisarray,he4_boundary_mass,label='he4_boundary_mass')
		pyl.plot(xaxisarray,mx1_bot,',r',label='conv bound')
		pyl.plot(xaxisarray,mx1_top,',r')
		pyl.plot(xaxisarray,mx2_bot,',r')
		pyl.plot(xaxisarray,mx2_top,',r')
		pyl.ylabel('mass coordinate')
		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')
		pyl.plot(xaxisarray,star_mass,label='star_mass')
		pyl.legend(loc=6)

    	def t_lumi(self,num_frame,xax):
		''' Luminosity evolution as a function of time or model
		
		num_frame    number of frame to plot this plot into
		xax          string that is either model or time to indicate what is 
			     to be used on the x-axis'''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction. needs to be "time" or "model"'
		
	    
		logLH   = self.get('log_LH')
		logLHe  = self.get('log_LHe')
	
		pyl.plot(xaxisarray,logLH,label='L_(H)')
		pyl.plot(xaxisarray,logLHe,label='L(He)')
		pyl.ylabel('log L')
		pyl.legend(loc=2)
	
	
		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')

    	def t_surf_parameter(self,num_frame,xax):
		''' Surface parameter evolution as a function of time or model
		
		num_frame    number of frame to plot this plot into
		xax          string that is either model or time to indicate what is 
			     to be used on the x-axis'''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction. needs to be "time" or "model"'
		
	    
		logL    = self.get('log_L')
		logTeff    = self.get('log_Teff')
	
		pyl.plot(xaxisarray,logL,'-k',label='log L')
		pyl.plot(xaxisarray,logTeff,'-k',label='log Teff')
		pyl.ylabel('log L, log Teff')
		pyl.legend(loc=2)
	
	
		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')

