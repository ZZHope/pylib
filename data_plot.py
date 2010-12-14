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
					print x
					print y
					return x,y
				if (i%sparse)==0:
					tmpX.append(x[i])
					tmpY.append(y[i])
			return tmpX, tmpY
			
	def plotMulti(self,atriX,atriY, cycList,title,legend=None ,logX=False, logY=False, base=10,sparse=1,pdf=False,xMin=None,xMax=None,yMin=None,Ymax=None):
		'''
		Method for superimposing multiple plots and saving it to a png or PDF
		Input:
		atriX: The name of the attribute you want on the x axis
		atriY: The name of the attribute you want on the Y axis
		cycList: List of cycles that you would like plotted
		title: The title of the grapgh and the name of the file.
		Legend: A list of legends for each of your cycles, or one legend 
			for all of the cycles
		pdf: A boolean of if the image should be saved to a pdf file.
			xMin,xMax, yMin, YMax:  plot coopdinates.	
		logX: A boolean of weather the user wants the x axis logarithmically
		logY: A boolean of weather the user wants the Y axis logarithmically
		base: The base of the logarithm. Default = 10
		sparse: Argument that skips every so many data points. For 
			example if this argument was 5, This method would plot
			the 0th, 5th, 10th ... elements.
		'''
		if str(legend.__class__)!="<type 'list'>":# Determines the legend is a list
			legendList=False
		else:
			legendList=True
			
		if legendList and len(cycList) !=len(legend): #if it is a list, make sure there is an entry for each cycle
			print 'Please input a proper legend, with correct length, aborting plot'
			return None
		for i in xrange(len(cycList)):
			if legendList:
				self.plot(atriX,atriY,cycList[i], legend=legend[i],base=base,sparse=sparse, logX=logX,logY=logY,show=False)
			else:
				self.plot(atriX,atriY,cycList[i], legend=legend,base=base,sparse=sparse, logX=logX,logY=logY,show=False)
		
			pl.title(title)
			
			
			
			if xMin!=None and xMax!=None and yMin!=None and yMax!=None:
				pl.xlim(xMin,xMax)
				pl.ylim(yMin,yMax)
			
			if not pdf:
				pl.savefig(title+str(cycList[i])+'.png', dpi=100)
			else:
				pl.savefig(title+cycList[i]+'.pdf', dpi=100)
			pl.clf()
		return None
	
	def plot(self,atriX,atriY, FName=None,numType='ndump',legend=None,labelX=None, labelY=None ,
		indexX=None, indexY=None, title=None, shape='.',logX=False, logY=False, base=10,sparse=1, show=True):
		"""
		Simple function that plots atriY as a function of atriX
		This method will automatically find and plot the requested data.
		Input:
		atriX, The name of the attribute you want on the x axis
		atriY, The name of the attribute you want on the Y axis
		Fname: Be the filename, Ndump or time, or cycle
		numType: designates how this function acts and how it interprets 
			 FName. Defaults to file
		if numType is 'file', this function will get the desird 
		attribute from that file
		if numType is 'NDump' function will look at the cycle with that 
		nDump
		if numType is 't' or 'time' function will find the _cycle with 
		the closest time stamp 
		labelX: The label on the X axis
		labelY: The label on the Y axis
		indexX: Depreciated: If the get method returns a list of lists, 
			indexX would be the list at the index indexX in the list.
		indexY: Depreciated: If the get method returns a list of lists, 
			indexY would be the list at the index indexX in the list.
		shape: What shape and colour the user would like their plot in.
		       Please see 
		       http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
		       for all possable choices
		title: The Title of the Graph
		logX: A boolean of weather the user wants the x axis logarithmically
		logY: A boolean of weather the user wants the Y axis logarithmically
		base: The base of the logarithm. Default = 10
		sparse: Argument that skips every so many data points. For 
			example if this argument was 5, This method would plot
			the 0th, 5th, 10th ... elements.
		show: A boolean of if the plot should be displayed> usefull with
			the multiPlot method
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
		elif plotType=='se':
			print 'This method is not supported for '+str(self.__class__)
			return
		elif plotType=='PPN' :
			if FName==None:
				FName=len(self.files)-1
			if numType=='ndump':
				numType='cycNum'
			listY=self.get(atriY,FName,numType)
			listX=self.get(atriX,FName,numType)
		elif plotType=='xtime' or plotType=='mesa_profile' or plotType=='AsciiTable' or plotType=='mesa.star_log':
			listY=self.get(atriY)
			listX=self.get(atriX)
		else:
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
		# Sparse stuff
		
		tmpX,tmpY=self.sparse(tmpX,tmpY, sparse)
		
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
			
		pl.plot(listX,listY,shape,label=legend)
		pl.legend()
		pl.title(title)
		pl.xlabel(labelX)
		pl.ylabel(labelY)
		if show:
			pl.show()
		
		
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
			
	#from mppnp.se
	def abu_chart(self, cycle, mass_range=None ,ilabel = 1,imlabel = 1,imagic = 0,plotaxis=[0,0,0,0],format = 'pdf'):
		'''
		Plots an abundence chart
		input:
		cycle: The cycle we are looking in
		ilabel: elemental labels off/on [0/1]
		imlabel: label for isotopic masses off/on [0/1]
		imagic:  turn lines for magic numbers off/on [0/1]
		plotaxis: Set axis limit: If default [0,0,0,0] the complete 
			  range in (N,Z) will be plotted
			  format: What format will this be saved in ['pdf'/'png']
		mass_range - a 1x2 array containing the lower and upper mass range.
		    		 If this is an instance of abu_vector this will 
		    		 only plot isotopes that have an atominc mass 
		    		 within this range. This will throw an error if
		    		 this range does not make sence ie [45,2]
			 	if None, it will plot over the entire range
				Defaults to None
				
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
		
		plotType=self.classTest()
		
		
		if mass_range!=None and mass_range[0]>mass_range[1]:
			print 'Please input a proper mass range'
			print 'Returning None'
			return None
		
		if plotType=='se':
			nin=self.se.A
			zin=self.se.Z
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
		  
		format = 'pdf'
		params = {'axes.labelsize':  15,
			  'text.fontsize':   15,
			  'legend.fontsize': 15,
			  'xtick.labelsize': 15,
			  'ytick.labelsize': 15,
			  'text.usetex': True}
		pl.rcParams.update(params)
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
		norma = colors.Normalize(vmin=-20,vmax=0)
		# set x- and y-axis scale aspect ratio to 1
		ax.set_aspect('equal')
		#print time,temp and density on top
		temp = ' '#'%8.3e' %ff['temp']
		time = ' '#'%8.3e' %ff['time']
		dens = ' '#'%8.3e' %ff['dens']
		
		box1 = TextArea("t : " + time + " s~~/~~T$_{9}$ : " + temp + "~~/~~$\\rho_{b}$ : " \
			      + dens + ' g/cm$^{3}$', textprops=dict(color="k"))
		anchored_box = AnchoredOffsetbox(loc=3,
				child=box1, pad=0.,
				frameon=False,
				bbox_to_anchor=(0., 1.02),
				bbox_transform=ax.transAxes,
				borderpad=0.,
				)
		ax.add_artist(anchored_box)
		
		## only abundance plotted
		
		patches = []
		color = []
		
		for i in range(nzmax):
		    for j in range(nnmax):
		      if nzycheck[j,i,0]==1:
			xy = j-0.5,i-0.5
			rect = Rectangle(xy,1,1,ec='k')
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
		graphname = 'abundance-chart.'+format
		  
		# Add black frames for stable isotopes
		'''
		f = open('stable.dat')
		
		head = f.readline()
		stable = []
		'''
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
				
				xy = nn-0.5,zz-0.5
				rect = Rectangle(xy,1,1,ec='k',fc='None',fill='False',lw=3.)
				rect.set_zorder(2)
				ax.add_patch(rect)
		
		
		
		
		# decide which array to take for label positions
		iarr = 0
		
		# plot element labels
		if ilabel==1:
		  for z in range(nzmax):
		    try:
		      nmin = min(argwhere(nzycheck[:,z,iarr]))[0]-1
		      ax.text(nmin,z,elname[z],horizontalalignment='center',verticalalignment='center',fontsize='medium',clip_on=True)
		    except ValueError:
		      continue
		      
		# plot mass numbers
		if imlabel==1:
		  for z in range(nzmax):
		     for n in range(nnmax):
			a = z+n
			if nzycheck[n,z,iarr]==1:
			  ax.text(n,z,a,horizontalalignment='center',verticalalignment='center',fontsize='small',clip_on=True)
		
		# plot lines at magic numbers
		if imagic==1:
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
		
		#fig.savefig(graphname)
		print graphname,'is done'
		
		pl.show()
		return
		
	def iso_abund(self,  cycle, stable=False,mass_range=None,ylim=[1e-13,10],shape='o',ref=-1):
		''' plot the abundance of all the chemical species
		inputs:
		    
		    cycle       - a string/integer of the cycle of interest.
		    stable     - a boolean of whether to filter out the unstables.
		    		Defaults to False
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
		    Shape -The Shape of the dataplots, will default to circles 
		    	  
		'''
		elem_list = []
		elem_index = []
		masses = []
		plotType=self.classTest()
		if mass_range!=None and mass_range[0]>mass_range[1]:
			print 'Please input a proper mass range'
			print 'Returning None'
			return None
		if plotType=='se':
			isotope_to_plot = self.se.isotopes
			
			z=self.se.Z #charge
			a=self.se.A #mass
			isomers=self.se.isomeric_states
			abunds = self.se.get(cycle,'iso_massf')
			masses = self.se.get(cycle,'mass')
			if mass_range == None:
			    print 'Using default mass range'
			    mass_range = [min(masses),max(masses)]    
			masses.sort()
			mass_range.sort()
			isom=[]
			tmp=[]
			for i in range (len(isotope_to_plot)):
				if z[i]!=0 and isomers[i]==1: #if its not 'NEUt and not an isomer'
					tmp.append(self.elements_names[int(z[i])]+'-'+str(int(a[i])),)
				elif isomers[i]!=1: #if it is an isomer
					isom.append(self.elements_names[int(z[i])]+'-'+str(int(a[i]))+'-'+str(int(isomers[i]-1)))	
						
			
			
		elif plotType=='PPN':
			isotope_to_plot = self.get('ISOTP', cycle)
			z=self.get('Z', cycle) #charge
			a=self.get('A', cycle) #mass
			isomers=self.get('ISOM', cycle)
			yps=self.get('ABUNDNACE_MF', cycle)
			if ref >-1:
				ypsRef=self.get('ABUNDNACE_MF', ref)
			if mass_range != None:
				tmpA=[]
				tmpZ=[]
				tmpIso=[]
				tmpIsom=[]
				tmpyps=[]
				tmpypsRef=[]
				for i in xrange(len(a)):
					if (a[i] >mass_range[0] and a[i]<mass_range[1]) or (a[i]==mass_range[0] or a[i]==mass_range[1]):
						tmpA.append(a[i])
						tmpZ.append(z[i])
						tmpIso.append(isotope_to_plot[i])
						tmpIsom.append(isomers[i])
						tmpyps.append(yps[i])
						if ref >-1:
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
			
			for i in xrange(len(yps)):
				if ref >-1:
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
		if  plotType=='se' or mass_range != None:
			print '\tmass_range:', mass_range[0], mass_range[1]
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
					if abundance < 1e-20:
					    abundance = 1e-20
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
		colors = ['g','r','c','m','k']
		cl_index = 0
		
		l1 = []
		l2 = []        
		
		#print abund_plot,mass_num
		for j in xrange(len(abund_plot)):        #Loop through the elements of interest
		    #    Process the line
		    #print 'processing line'
		    for l in xrange(len(abund_plot[j])):
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
						if mass_range !=None and plotType=='PPN' :
							if coordinates[0]>=mass_range[0]-.5 and coordinates[0]<=mass_range[1]+.5:
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
			pl.semilogy(mass_num[j],abund_plot[j],'b'+shape)
		    except OverflowError:
			None
			#print 'div by zero', len(mass_num[j]), len(abund_plot[j])
		    except IndexError:
			None
			#print 'out of bounds II: ', len(mass_num),len(abund_plot),j
		pl_index = 0
		
		cl_index = 0
		
		for j in xrange(len(isom)):
			pl.semilogy(isom[j][0].split('-')[1],isom[j][1],'r'+shape)#,str(colors[cl_index]+plot_type[pl_index])))
			cl_index+=1
			pl_index+=1
			if pl_index > 2:
			    pl_index = 0
			if cl_index > 4:
			    cl_index = 0
			coordinates=[int(isom[j][0].split('-')[1]),isom[j][1]]
			name=isom[j][0]

			pl.text(coordinates[0],coordinates[1], name.split('-')[0]+'m'+name.split('-')[2])
		if plotType=='se':
			title = str('Abundance of Isotopes over range %4.2f' %mass_range[0]) + str('-%4.2f' %mass_range[1]) +\
				str(' for cycle %d' %int(cycle))
		else:
			if ref==-1:
				if mass_range ==None:
					title = str('Abundance of Isotopes for Cycle '+str(cycle))
				else:
					title = str('Abundance of Isotopes with A between '+str(mass_range[0])+' and '+str(mass_range[1])+' for Cycle '+str(cycle))
					pl.xlim([mass_range[0]-.5,mass_range[1]+.5])
			else:
				if mass_range ==None:
					title = str('Abundance of Isotopes for Cycle '+str(cycle))
				else:
					title = str('Abundance of Isotopes with A between '+str(mass_range[0])+' and '+str(mass_range[1])+' for Cycle '+str(cycle))
					pl.xlim([mass_range[0]-.5,mass_range[1]+.5])
		pl.ylim(ylim)
		pl.title(title)
		pl.xlabel('Mass Number')
		pl.ylabel('Relative Abundance')
		pl.grid()
		pl.show()
		return
	
	def plotprofMulti(self,ini,end,delta,what_specie,xlim1,xlim2,ylim1,ylim2):
	
		''' create a movie with mass fractions vs mass coordinate 
		between xlim1 and xlim2, ylim1 and ylim2. 
	
		ini          - initial model
		end          - final model 
		delta        - sparsity factor of the frames
		what_specie  - array with species in the plot
		xlim1, xlim2 - mass coordinate range 
		ylim1, ylim2 - mass fraction coordinate range
	  
		'''
		plotType=self.classTest()
		if plotType=='se'or plotType=='mesa_profile':
			for i in range(ini,end+1,delta):
			    step = int(i)
			    print step
			    for j in range(len(what_specie)):
				self.plot_prof_1(what_specie[j],step,xlim1,xlim2,ylim1,ylim2, False)
			    #          
			    filename = str('%03d' % step)+'_test.png'             
			    pl.savefig(filename, dpi=100) 
			    print 'wrote file ', filename
			    #
			    pl.clf()
		else:
			print 'This method is not supported for '+str(self.__class__)
			return
	# From mesa_profile
    	def plot_prof_1(self,species,keystring,xlim1,xlim2,ylim1,ylim2, show=True):
	
		''' plot one species for cycle between xlim1 and xlim2                               
		
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
		pl.plot(y,x,'-',label=str(keystring))
		pl.xlim(xlim1,xlim2)
		pl.ylim(ylim1,ylim2)
		pl.legend()
	
		pl.xlabel('$Mass$ $coordinate$', fontsize=20)
		pl.ylabel('$X_{i}$', fontsize=20)
		pl.title('Mass='+str(tot_mass)+', Time='+str(age)+' years, cycle='+str(mod))
		if show:
			pl.show()
	
	# From mesa.star_log
	

