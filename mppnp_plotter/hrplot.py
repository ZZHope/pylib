#!/usr/bin/python
import os
import sys
from string import *
import numpy as np
import scipy as scipy
from PyQt4 import QtCore as qc
import time as t

import matplotlib.pyplot as graph
from matplotlib.figure import Figure
from matplotlib import mpl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import h5file


#################################################################################################################
#
#	This set of scripts controlls the actual plotting.  It allows for the plot processing, and then the passing
#	of the plot object from this script to the h5Plot.py script
#
################################################################################################################




#	The HR Plot is a common, base plot for this system.
#	It gnerates an interactive plot containing the temp
#	vs Luminosity data.  Poinst are selected generating 
#	the aforementioned plots.

class hrplot(qc.QThread):
	
	# Initialize the class.   Get the data
	def __init__(self, teff, logL, sp_controll, sparsity_factor, parent):
		
		# Initialize the thread, import global variables.
		qc.QThread.__init__(self,None)

		self.teff = teff
		self.logL = logL
		#self.textEditor = textEdit
		self.spc = sp_controll
		self.sparsity = sparsity_factor
		self.parent = parent
			
	def __del__(self):
		print 'Deleting HR plot'
		self.teff = []
		self.logL = []
		#self.textEditor = None
		self.spc = None
		self = None
		#self.fig = None
		self.canvas = None
		self.fig.close()
		print str(self.quit())
		#print 'finished ', qc.QThread.finished()
		
	#  Generate a figure, canvas and plot
	def run(self): 
			
		#	Build the figure
		self.fig = Figure((10.0, 8.0), dpi=100,frameon=True)
		self.canvas = FigureCanvas(self.fig)
		try:
			xbounds = [max(max(self.teff)),min(min(self.teff))]
		except:
			xbounds = [max(self.teff),min(self.teff)]
		
		#	Decide on the scale for labels based on simple logic (a star is not going to be 10 K or less!!!)
		if xbounds[0] < 10:
			xlabel = 'Temperature [log K]'
			ylabel = 'Luminosity [log L]'
		else:
			xlabel = 'Temperature [K]'
			ylabel = 'Luminosity [L]' 

		self.axes = self.fig.add_subplot(111,title='HR Plot',xlabel=xlabel, ylabel=ylabel)#,xlim=xbounds-removed due to changes in datatype
		self.canvas.mpl_connect('pick_event', self.on_pick)
		#self.canvas.mpl_connect('pick_event',h5s.do_hr_subplot)
		self.axes.grid()
		#print 'self.teff', self.teff
		#print 'self.logl', self.logL
		self.axes.plot(self.teff,self.logL) #,'o',picker=3
		self.axes.plot(self.teff,self.logL,'o',picker=3) #
		self.axes.set_xlim(xbounds)
		
		self.canvas.draw()		
		#del self.teff , self.logL, xbounds
		if self.parent.write.checkState() == 2:
			print 'writing to file'
			self.parent.write_to_file(self.teff, self.logL, 'teff', 'logL', 'HR',)
			self.parent.write.setCheckState(0)
			print self.parent.write.checkState()

	
	#	This function controls the interactivity with the plot.  
	def on_pick(self, event):
		#	Extract pick event data
		thisline = event.artist
		xdata = thisline.get_xdata()
		ydata = thisline.get_ydata()
		ind = event.ind*self.sparsity
		try:
			ind = ind[0]
		except IndexError:
			None
		print ind
		#print 'Axes Coords: ', event.artist.get_axes(), event.get_x(), event.get_y()
		
		self.spc.setText(str(ind))
		
		
		
		
#	This class is designed to plot cattrs from a list of cycles.  
#	cattrs categories are selected by the user.  The data is provided 
#	for this function, which then is plotted

class subplot(qc.QThread):
	def __init__(self, parent):
		qc.QThread.__init__(self, None)
		self.parent = parent
		print 'doing plot'
		
	def __del__(self):
		print 'Deleting subplot'
		print str(self.quit())
		self.fig.close()
		self.canvas = None
		#print 'finished ', qc.QThread.finished()
		self = None
		print self

	#	Does the plotting.  
	#	BUilds, figure, canvas, plot and legend
	def run(self,x_data,y_data,xval,yval,cycle, indices, isos, logx, logy):
		self.fig = Figure((10.0, 8.0),frameon=True)
		self.canvas = FigureCanvas(self.fig)
		
		print 'subplot', logx, logy
		
		#print 'using range: ', ranger
		
		if len(cycle) < 5:
			title  = str(x_data) + ' vs ' +str(y_data + ' of cycle ' + str(cycle))
		else: 
			title = str(x_data) + ' vs ' +str(y_data)

		
		self.axes = self.fig.add_subplot(111,title=title,xlabel=x_data, ylabel=y_data)
		
#		if ranger != []:
#			self.axes.ylim(ranger)
		#self.canvas.mpl_connect('pick_event', self.on_pick)
		#self.canvas.mpl_connect('pick_event',h5s.do_hr_subplot)
		self.axes.grid()
		l1=[]
		#for x in range(len(xval)):
		#	l1.append(self.axes.plot(xval[x],yval[x]))
		if indices == []:
			if logx == True and logy == True:
				l1.append(self.axes.loglog(xval,yval))
			elif logx == True and logy == False:
				l1.append(self.axes.semilogx(xval,yval))
			elif logx == False and logy == True:
				l1.append(self.axes.semilogy(xval,yval))
			else:
				l1.append(self.axes.plot(xval,yval))
			leg = self.axes.legend(l1,cycle,loc=0)
		else:
			
			isotopes = []
			#print len(xval), len(yval), len(yval[0])
			for index in indices:
				ytemp = []
				try:
					for x in xrange(len(yval)):
						ytemp.append(yval[x][index])
					if logx == True and logy == True:
						l1.append(self.axes.loglog(xval,ytemp))
					elif logx == True and logy == False:
						l1.append(self.axes.semilogx(xval,ytemp))
					elif logx == False and logy == True:
						l1.append(self.axes.semilogy(xval,ytemp))
					else:
						l1.append(self.axes.plot(xval,ytemp))
				except (ValueError, AssertionError):
					for x in xrange(len(yval[0])):
						ytemp.append(yval[x][index])
					if logx == True and logy == True:
						l1.append(self.axes.loglog(xval[0],ytemp))
					elif logx == True and logy == False:
						l1.append(self.axes.semilogx(xval[0],ytemp))
					elif logx == False and logy == True:
						l1.append(self.axes.semilogy(xval[0],ytemp))
					else:
						l1.append(self.axes.plot(xval[0],ytemp))				
				isotopes.append(isos[index])
			leg = self.axes.legend(l1,isotopes)
		leg.get_frame().set_fill(False)		

	

		self.canvas.draw()		
	#	ranger = []

		if self.parent.write.checkState() == 2:
			print 'writing to file'
			self.parent.write_to_file(xval,yval,x_data,y_data,'Subplot'+str(cycle))
			self.parent.write.setCheckState(0)
			print self.parent.write.checkState()
		
		
		
		del x_data,y_data,xval,yval,cycle
		
		
#	The abundance plot uses the abundance data from an h5 file (passed in from outside) for 
#	a list of cycles.  
class abundance_plot(qc.QThread):
	values = []
	cycs = []
	isos = []
	#  Initialize the thread and set the class data
	def __init__(self,yvals,cycles,isotopes, parent):
		qc.QThread.__init__(self,None)
		self.values = yvals
		self.cycs = cycles
		self.isos = isotopes
		self.parent = parent
		
	def __del__(self):
		print 'Deleting abundance plot'
		self.values = []
		self.cycs = []
		self.isotopes = []
		self.fig.close()
		self.canvas = None
		print str(self.quit())
		#print 'finished ', qc.QThread.finished()
		self = None
		print self

	#  Perform the plot
	def run(self):
		# Find the limits
		maxY = np.log10(max(max(max(self.values))))
		minY = np.log10(min(min(min(self.values))))
		
		#  Build the figure, canvas and subplot. 
		self.fig = Figure((10.0, 8.0),frameon=True)
		self.canvas = FigureCanvas(self.fig)
		#title  = str(x_data) + ' vs ' +str(y_data + ' of cycle ' + str(cycle))
		self.axes = self.fig.add_subplot(111,title='Abundance plot',xlabel='Mass Coord', ylabel='Abundance', ylim=[minY,maxY])

		self.axes.grid()
		l1=[]

		#  Loop through the list of cycles, isotopes, and mass coordinates.
		#  Plots them.
		#  Since the legend is finicky I could only plot the cycles (if muliple are selected) 
		#  Or the Isotopes themselves.
		#  FIXME- Make the the legend better
		if len(self.isos) == 1:
			for x in range(len(self.values)):
				for y in range(len(self.values[x])):
					coos = len(self.values[x][y])
					coords = scipy.linspace(0,coos-1,coos)

					l1.append(self.axes.semilogy(coords,self.values[x][y]))

			leg1 = self.axes.legend(l1,self.cycs,loc=0)
		elif len(self.isos) > 1:	
			for y in range(len(self.values[0])):
				coos = len(self.values[0][y])
				coords = scipy.linspace(0,coos-1,coos)
				
				l1.append(self.axes.semilogy(coords,self.values[0][y]))

			leg1 = self.axes.legend(l1,self.isos,loc=0)
		
		#  Add a legend and draw the canvas
		leg1.get_frame().set_fill(False)		
		
		self.canvas.draw()		

		if self.parent.write.checkState() == 2:
			print 'writing to file'
			self.parent.write_to_file(self.coords, self.values, 'Mass_coord', 'Abundance', 'Abundance_Plot')
			self.parent.write.setCheckState(0)
			print self.parent.write.checkState()
		
		
		
		del self.values , self.cycs, self.isos, coords, l1
		
#	Kippenham plots are nasty little beasts.   This performs a KH plot

class KH_Plot(qc.QThread):
	
	a_age = []
	a_mass = []
	X = []
	Z = []
	
	def __init__(self,all_age,conv, scale, checkstate ,parent):
		qc.QThread.__init__(self,None)
		self.picker = False #	This variable is used to control the picker, allowing selection between points.
		self.a_age = all_age
		self.conv = conv
		self.scale = scale
		self.parent = parent
		self.checkstate = checkstate
		
		
	def __del__(self):
		print 'destroying Kippenhahn Diagram'
		self.a_age = []
		self.conv = []
		self.fig.close()
		self.canvas = None
		print str(self.quit())
		#print 'finished ', qc.QThread.finished()
		self = None
		print self
		
	def run(self):
		
		print 'start'
		self.fig = Figure((10.0, 10.0),frameon=True)#, figsize=(9,9)
		self.canvas = FigureCanvas(self.fig)
		self.axes1 = self.fig.add_subplot(111,title='Kippenhahn Diagram',xlabel='Age/Cycle', ylabel='Mass Coord')
		self.canvas.mpl_connect('pick_event', self.on_pick)
		
		print 'cavas built'
		self.axes1.grid()
		
		l1=[]	
		plot = []
		
		print len(self.conv)
		max_len = len(self.conv[0])
		
		for x in xrange(len(self.conv)):
			if len(self.conv[x]) > max_len:
				max_len = len(self.conv)
		#print max_len
		
		#print self.conv


		Z = []
		
		self.conv = np.transpose(self.conv)
		
		
		
		for x in xrange(len(self.conv[0])):
			#print len(self.a_age), len(self.conv[x]), len(self.conv)
			self.axes1.plot(self.a_age,self.conv[x],'g.',picker=1) #, alpha=0.25
		
		self.axes1.set_ylim([0,self.conv.max()])
		self.axes1.grid()

		if self.checkstate == 2:
			print 'writing to file'
			self.parent.write_to_file(self.a_age, self.conv, 'age', 'convective_indicator', 'Kippenhahn Plot')
			self.parent.write.setCheckState(0)
			print self.parent.write.checkState()


		
		
		#print 'after loopsssssssssssssss'
		del	self.a_age, self.conv
		print self.canvas
		self.canvas.draw()	
		self.mass_coord = []	#	A range of mass coordinates for the picker to populate
		
		
	#	This function controls the interactivity with the plot.  
	def on_pick(self, event):
		#	Extract pick event data
		thisline = event.artist
		xdata = thisline.get_xdata()
		ydata = thisline.get_ydata()		
		ind = event.ind
		if len(ind) > 1:
			ind = ind[0] 
		
		self.picked_cycle = self.parent.h5s.cycles[int(ind*self.scale)]
		print 'plotting cycle: ', self.picked_cycle
		self.parent.sub_plot_controller.setText(str(int(ind*self.scale)))
			
		
class SP_KH_Plot(qc.QThread):
	elements = []
	limits = []
	colors = ['b','g','r','c','m','y']	#	Hopefully all the colors you could need
	linestyle = ['.','-.','--','-','o','v','^','<','>','s','p']
	
	
	def __init__(self, variables, limits, cycle_list, age,Z, parent):
		qc.QThread.__init__(self,None)
		self.variables = variables #	Index of elements of interest
		self.limits = limits #	list of limits.  Indices match self.elements
		self.cycle_list = cycle_list
		self.age = age
		self.Z = Z
		self.parent = parent
		
	def __del__(self):
		print 'special kH delete event'
		self.fig.close()
		#self = None
	
	def run(self):
		#	Define the plot objects
		self.fig = Figure((10.0, 10.0),frameon=True)#, figsize=(9,9)
		self.canvas = FigureCanvas(self.fig)
		self.axes1 = self.fig.add_subplot(111,title='Kippenhahn Diagram',xlabel='Cycle/Age', ylabel='Mass Coord')
		#self.axes2 = self.fig.add_subplot(122,title='ColorBar')
		print 'cavas built'
		self.axes1.grid()	
		
		#print self.elements
		l1 = []
		l2 = []
		step = 0
		#linewidths = []
		for x in xrange(len(self.limits)):
			lineweight = 1
			for y in xrange(len(self.limits[x])):
				#if x == 0:
				#	l2.append(self.axes1.plot(self.age,self.Z[x+y+step],self.colors[x], linewidth=lineweight))
				#	linewidths.append(lineweight)
				if y == 0:
					l1.append(self.axes1.plot(self.cycle_list,self.Z[x+y+step],self.colors[x]+'o', linewidth=lineweight)) #+self.linestyle[y]
				else:
					self.axes1.plot(self.cycle_list,self.Z[x+y+step],self.colors[x]+'o',linewidth=lineweight)
				if step%2 == 1:
					lineweight+=1
					
			step+=1
		
		#self.axes1.set_xticklabels(self.age)
		leg1 = self.axes1.legend(l1, self.variables)
	#	leg2 = self.axes1.legend(l2, self.limits)
	#	for x in xrange(len(self.Z)):
			#print self.Z[x]
	#		self.axes1.plot(self.age,self.Z[x])
		
		
		if self.parent.write.checkState() == 2:
			print 'writing to file'
			self.parent.write_to_file(None, abundance, '', 'Abundance', 'Special_Kippenhahn_Plot')
			self.parent.write.setCheckState(0)
			print self.parent.write.checkState()


		

		print self.canvas
		self.canvas.draw()	


class AB_Plot(qc.QThread): #
	
	a_age = []
	a_mass = []
	X = []
	Z = []	
			
	def __init__(self,all_age,all_mass,x,z, parent):
		#Thread.__init__(self)
		qc.QThread.__init__(self,None)
		self.a_age = all_age
		self.a_mass = all_mass
		self.X = x
		self.Z = z
		self.parent = parent
		
	def __del__(self):
		print 'deleting AB Thread'
		self.a_age = []
		self.a_mass = []
		self.X = []
		self.Z = []
		self.fig.close()
		self.canvas = None
		print str(self.quit())
		#print 'finished ', qc.QThread.finished()
		self = None
		print self
		
	def run(self):
		print 'start'
		self.fig = Figure((10.0, 10.0),frameon=True)#, figsize=(9,9)
		CB = self.fig.colorbar
		self.canvas = FigureCanvas(self.fig)
		self.axes1 = self.fig.add_subplot(111,title='Abundance Plot',xlabel='Age', ylabel='Shell Number')
		#self.axes2 = self.fig.add_subplot(122,title='ColorBar')
		
		self.axes1.grid()
		l1=[]	

		plot = []
		lim = len(self.Z[0])

		if len(self.X) < len(self.Z):
			diff = len(self.Z) - len(self.X)
			while diff > 0:
				dx = self.X[len(self.X)-1] -self.X[len(self.X)-2]
				np.append(self.X, (self.X[len(self.X)-1]+dx))
				diff-=1
		
		if len(self.a_age) < len(self.Z):
			diff = len(self.Z)-len(self.a_age)
			while diff > 0:
				dx = self.a_age[len(self.a_age)-1]-self.a_age[len(self.a_age)-2]
				self.a_age.append(self.a_age[len(self.a_age)-1]+dx)
				diff-=1	
		
		print len(np.transpose(self.Z)),len(np.transpose(self.Z)[0]),len(self.a_age), len(self.X)
		l1 = (self.axes1.contourf(np.transpose(self.Z),alpha=0.5))#, cmap=cmap    pcolormesh looks good  
		self.canvas.draw()		
		if self.parent.write.checkState() == 2:
			print 'writing to file'
			self.parent.write_to_file(None, self.Z, '', 'Abundance', 'Special_Kippenhahn_Plot')
			self.parent.write.setCheckState(0)
			print self.parent.write.checkState()
		
		#del self.a_age, self.a_mass, self.X, self.Z 		
		return self		
	
class generic_plot(qc.QThread):
	"""This preforms a generic plot, currently (04-26-11) not working"""
	def __init__(self,x_parm, y_parm, x_data, y_data, x_log, y_log,  plot_type, parent):
		qc.QThread.__init__(self,None)
		self.xparm = x_parm
		self.yparm = y_parm
		self.xdata = x_data
		self.ydata = y_data
		self.x_log = x_log
		self.y_log = y_log
		self.plottype = plot_type
		self.parent = parent
		self.write = self.parent.write.checkState()
		
	def __del__(self):
		del self.xparm, self.yparm, self.xdata, self.ydata, self.plottype
		self.fig.close()
		print str(self.quit())
		self = None
		
	
	def run(self):
		print self.parent.write
		print "in generic_plot"
		
		print self.xparm, self.yparm
		
		#print 'start ', np.shape(self.xdata), len(self.ydata[0])
		
		l1 = []		
		self.fig = Figure((10.0, 10.0),frameon=True)#, figsize=(9,9)
		self.canvas = FigureCanvas(self.fig)
		self.axes1 = self.fig.add_subplot(111,title=str(self.xparm) + str(' vs ') + str(self.yparm), xlabel=str(self.xparm), ylabel=str(self.yparm))
		self.axes1.grid()
		
		while len(self.xdata) == 1:
			print 'trimming x'
			self.xdata = self.xdata[0]
	#	print self.ydata
		#if len(self.ydata) == 1:
		#	self.ydata = self.ydata[0]
		
		while len(self.ydata) == 1:
			print 'trimming y'
			self.ydata = self.ydata[0]
		
		print "xdata" + str(len(self.xdata))
		print self.xdata
		print "ydata"+ str(len(self.ydata))
		print self.ydata
		
		if self.xparm == 'Cycle':
			for x in xrange(len(self.xdata)):
				#print self.xdata[x]
				if len(self.xdata) != 1:
					self.xdata[x] = int(self.xdata[x])
				else:
					self.xdata[0][x] = int(self.xdata[0][x])
		
		if self.yparm.count('Cycle'):
			ind = self.yparm.index('Cycle')
			for x in xrange(len(self.y_data[ind])):
				self.ydata[ind][x] = int(self.ydata[ind][x])
		
		#print len(self.xdata[0]), len(self.ydata[0])
		if self.plottype == 0:
			for i in xrange(len(self.yparm)):					
				#
				try:
					if self.x_log == 0 and self.y_log == 0:
						l1.append(self.axes1.plot(self.xdata[0],self.ydata[i]))
					elif self.x_log == 0 and self.y_log == 2:
						l1.append(self.axes1.semilogy(self.xdata[0],self.ydata[i]))
					elif self.x_log == 2 and self.y_log == 0:
						l1.append(self.axes1.semilogx(self.xdata[0],self.ydata[i]))
					else:
						l1.append(self.axes1.loglog(self.xdata[0],self.ydata[i]))
						
				except AssertionError:
					if self.x_log == 0 and self.y_log == 0:
						l1.append(self.axes1.plot(self.xdata,self.ydata[i]))
					elif self.x_log == 0 and self.y_log == 2:
						l1.append(self.axes1.semilogy(self.xdata,self.ydata[i]))
					elif self.x_log == 2 and self.y_log == 0:
						l1.append(self.axes1.semilogx(self.xdata,self.ydata[i]))
					else:
						l1.append(self.axes1.loglog(self.xdata,self.ydata[i]))					
						
						
				except IndexError:
					try:
						if self.x_log == 0 and self.y_log == 0:
							l1.append(self.axes1.plot(self.xdata,self.ydata[i]))
						elif self.x_log == 0 and self.y_log == 2:
							l1.append(self.axes1.semilogy(self.xdata,self.ydata[i]))
						elif self.x_log == 2 and self.y_log == 0:
							l1.append(self.axes1.semilogx(self.xdata,self.ydata[i]))
						else:
							l1.append(self.axes1.loglog(self.xdata,self.ydata[i]))
							
					except IndexError:
						if self.x_log == 0 and self.y_log == 0:
							l1.append(self.axes1.plot(self.xdata,self.ydata))
						elif self.x_log == 0 and self.y_log == 2:
							l1.append(self.axes1.semilogy(self.xdata,self.ydata))
						elif self.x_log == 2 and self.y_log == 0:
							l1.append(self.axes1.semilogx(self.xdata,self.ydata))
						else:
							l1.append(self.axes1.loglog(self.xdata,self.ydata))


																						
		else:	
			for i in xrange(len(self.yparm)):					
				#print 'zip', len(self.ydata[0])
				try:
					while len(self.ydata) == 1:
						self.ydata = self.ydata[0]
				except IndexError:
					print 'failed to trim'

				try:
					for j in xrange(len(self.ydata[i])):
						print 'ydata', self.ydata[i]
						if self.ydata[i][j] < 1e-20:
							self.ydata[i][j] = 1e-20

				except TypeError:
					for x in xrange(len(self.ydata)):
						if self.ydata[x] < 1e-20:
							self.ydata[x] = 1e-20
					if self.x_log == 0 and self.y_log == 0:
						l1.append(self.axes1.plot(self.xdata[i],self.ydata))
					elif self.x_log == 0 and self.y_log == 2:
						l1.append(self.axes1.semilogy(self.xdata[i],self.ydata))
					elif self.x_log == 2 and self.y_log == 0:
						l1.append(self.axes1.semilogx(self.xdata[i],self.ydata))
					else:
						l1.append(self.axes1.loglog(self.xdata[i],self.ydata))						
					
				
				try:

					if self.x_log == 0 and self.y_log == 0:
						l1.append(self.axes1.plot(self.xdata[i],self.ydata[i]))
					elif self.x_log == 0 and self.y_log == 2:
						l1.append(self.axes1.semilogy(self.xdata[i],self.ydata[i]))
					elif self.x_log == 2 and self.y_log == 0:
						l1.append(self.axes1.semilogx(self.xdata[i],self.ydata[i]))
					else:
						l1.append(self.axes1.loglog(self.xdata[i],self.ydata[i]))							
						
						
				except IndexError:
					#print 'two'
					try:
						for j in xrange(len(self.ydata[i])):
							if self.ydata[i][j] < 1e-20:
								self.ydata[i][j] = 1e-20
						if self.x_log == 0 and self.y_log == 0:
							l1.append(self.axes1.plot(self.xdata[0], self.ydata[i]))  
						elif self.x_log == 0 and self.y_log == 2:
							l1.append(self.axes1.semilogy(self.xdata[0],self.ydata[i]))
						elif self.x_log == 2 and self.y_log == 0:
							l1.append(self.axes1.semilogx(self.xdata[0],self.ydata[i]))
						else:
							l1.append(self.axes1.loglog(self.xdata[0],self.ydata[i]))							
							
							
					except ValueError:
						#print 'two_b'
						while len(self.ydata[i]) == 1:
						#	print 'trimming before run'
							self.ydata[i] = self.ydata[i][0]
						for j in xrange(len(self.ydata[i])):
							if self.ydata[i][j] < 1e-20:
								self.ydata[i][j] = 1e-20
						
						
						if self.x_log == 0 and self.y_log == 0:
							l1.append(self.axes1.plot(self.xdata,self.ydata[i][0]))
						elif self.x_log == 0 and self.y_log == 2:
							l1.append(self.axes1.semilogy(self.xdata,self.ydata[i][0]))
						elif self.x_log == 2 and self.y_log == 0:
							l1.append(self.axes1.semilogx(self.xdata,self.ydata[i][0]))
						else:
							l1.append(self.axes1.loglog(self.xdata,self.ydata[i][0]))							
							
							
				except ValueError:
					#print len(self.xdata), len(self.ydata[0][i])
					#print 'three'
					try:
						try:
							
							for j in xrange(len(self.ydata[0][i])):
								if self.ydata[0][i][j] < 1e-20:
									self.ydata[0][i][j] = 1e-20
							
							if self.x_log == 0 and self.y_log == 0:
								l1.append(self.axes1.plot(self.xdata,self.ydata[0][i]))
							elif self.x_log == 0 and self.y_log == 2:
								l1.append(self.axes1.semilogy(self.xdata,self.ydata[0][i]))
							elif self.x_log == 2 and self.y_log == 0:
								l1.append(self.axes1.semilogx(self.xdata,self.ydata[0][i]))
							else:
								l1.append(self.axes1.loglog(self.xdata,self.ydata[0][i]))							
						except (IndexError, TypeError):
							print 'index error', self.ydata[i]
							#try:
							for j in xrange(len(self.ydata[i])):
								if self.ydata[i][j] < 1e-20:
									self.ydata[i][j] = 1e-20
							#except:
								#None
							if self.x_log == 0 and self.y_log == 0:
								l1.append(self.axes1.plot(self.xdata,self.ydata[i]))
							elif self.x_log == 0 and self.y_log == 2:
								l1.append(self.axes1.semilogy(self.xdata,self.ydata[i]))
							elif self.x_log == 2 and self.y_log == 0:
								l1.append(self.axes1.semilogx(self.xdata,self.ydata[i]))
							else:
								l1.append(self.axes1.loglog(self.xdata,self.ydata[i]))								
						
							
					except ValueError:
						#print 'three_b'
					#	while len(self.ydata[i]) == 1:
					#		self.ydata[i] = self.ydata[i][0]
					#		print 'trimming before run y', len(self.ydata[i])
					#	try:
					#		while len(self.xdata[i]) == 1:
					#			self.xdata[i] = self.xdata[i][0]
					#			print 'trimming before run x', len(self.xdata[i])						
					#	except (TypeError,IndexError):
					#		print 'skipped trimming x'
						#np.transpose(self.ydata)
						
						while len(self.ydata[i]) ==1 :
							self.ydata[i] = self.ydata[i][0]
						while len(self.xdata) == 1:
							self.xdata = self.xdata[0]
						
						#try:
					#		for j in xrange(len(self.ydata[i])):
					#			if self.ydata[i][j] < 1e-20:
					#				self.ydata[i][j] = 1e-20
					#	except ValueError:
							None
							#for j in xrange(len(self.ydata[i])):
							#	print self.ydata[i][j]
						
						
						
						while len(self.ydata) == 1:
							print len(self.ydata), len(self.ydata[i]), len(self.xdata)
							self.ydata = self.ydata[0]
						#print len(self.ydata), len(self.ydata[i]), len(self.xdata)
						if self.x_log == 0 and self.y_log == 0:
							l1.append(self.axes1.plot(self.xdata,self.ydata[i]))
						elif self.x_log == 0 and self.y_log == 2:
							l1.append(self.axes1.semilogy(self.xdata,self.ydata[i]))
						elif self.x_log == 2 and self.y_log == 0:
							l1.append(self.axes1.semilogx(self.xdata,self.ydata[i]))
						else:
							l1.append(self.axes1.loglog(self.xdata,self.ydata[i]))
		
		l2 = []
		try:
			for x in xrange(len(self.yparm)):
				l2.append(l1[0][x])
		except IndexError:
			print l1[0]
		#print 'l2', l2
		#print 'l1', l1
			for x in xrange(len(self.yparm)):
				#print l1[x], self.yparm[x]
				try:
					l2.append(l1[x])
				except ValueError:
					print self.yparm[x]
				except IndexError:
					print self.yparm
		try:
			leg1 = self.axes1.legend(l1[:len(self.yparm)], self.yparm)
		except (IndexError,ValueError, ZeroDivisionError):
			print 'Legend failed.'
		
		print l2
		print l1
		#print self.ydata	
		
		self.canvas.draw()		
				
				
		
		if self.parent.write.checkState() == 2:
			print 'writing to file'
			self.parent.write_to_file(self.xdata, self.ydata, str(self.xparm), str(self.yparm), 'Generic_plot')
			self.parent.write.setCheckState(0)
			print self.parent.write.checkState()
			
class IA_Plot(qc.QThread):

	def __init__(self,parent, abund_plot, isotope_to_plot, index, rel_abunPlot, ranges, ref_mass):
		print 'at init'
		self.parent = parent
		self.abund_plot = abund_plot
		self.isotope_to_plot = isotope_to_plot
		self.index = index
		self.rel_abunPlot = rel_abunPlot
		self.ranges = ranges
		self.ref_mass = ref_mass
		print len(self.ref_mass), len(ref_mass)
		print 'initiated'
	
	def __del__(self):
		self.fig.close()
		del self.parent,self.abund_plot,self.isotope_to_plot,self.index,self.rel_abunPlot,self.ranges,self.ref_mass

	def run(self):
		
		l1 = []		
		print 'running'
		self.fig = Figure((10.0, 10.0),frameon=True)#, figsize=(9,9)
		self.canvas = FigureCanvas(self.fig)
		self.axes1 = self.fig.add_subplot(111,title='Isotopic Abundance over range: '+ str(self.ranges[0]) + ' To '+str(self.ranges[1]), xlabel='Mass Number', ylabel='Abundance')
		self.axes1.grid()
		
		
		print 'plot started::::::::::right here', self.rel_abunPlot
		
	
		atom_num = []
		if self.rel_abunPlot < 2:
			
			for i in xrange(len(self.parent.cycle_to_plot)):	#Loop through the cycles
				l1 = []
				color = 100000
				colors = []
				
				
				for j in xrange(len(self.index)):		#Loop through the elements of interest
					#	Process the line
					#print 'processing line'
					print self.abund_plot[i][j]
					for l in xrange(len(self.abund_plot[i][j])):
						
						if self.abund_plot[i][j][l] < 1e-20:
							#print 'fixing',  self.abund_plot[i][j][l]	
							self.abund_plot[i][j][l] = 1e-20
					
					
					if i == 0:		
						try:
							
							l1.append(self.axes1.semilogy(self.parent.mass_num[j],self.abund_plot[i][j]))
						except OverflowError:
							print 'div by zero'
						
						
						try:
							coordinates = [self.parent.mass_num[j][self.abund_plot[i][j].index(max(self.abund_plot[i][j]))],max(self.abund_plot[i][j]) ]	#	Get the coordinates of the highest abundance of each element
							self.axes1.text(coordinates[0],coordinates[1], self.parent.elem_list[self.index[j]])
							#print 'coordinates', coordinates
							#print self.parent.elem_list[self.index[j]],coordinates[0],coordinates[1]
						except ValueError:
							None#print 'Empty var p2:  ', self.abund_plot[i][j]
					else:
						self.axes1.semilogy(self.parent.mass_num[j],self.abund_plot[i][j])
						#print '2',self.parent.mass_num[j] ,self.abund_plot[i][j]
					try:
						self.axes1.semilogy(self.parent.mass_num[j],self.abund_plot[i][j],'bo')
						#print '1', self.parent.mass_num[j], self.abund_plot[i][j]
					except OverflowError:
						None#print 'div by zero', len(self.parent.mass_num[j]), len(self.abund_plot[i][j])
			self.axes1.set_ylim([1e-13,1])
			
			#leg1 = self.axes1.legend(l1, self.isotope_to_plot)
		
		else:
			
			print 'doing relative abundance plot'
			linestyles = ['-' , '--' , '-.' , ':' , 'None' , ' ' , '']
			
			for i in xrange(len(self.parent.cycle_to_plot)):
				l2 = []

				
				
				for j in xrange(len(self.index)):
			
					numerator = np.array(self.abund_plot[i][j])
					denominator = np.array(self.ref_mass[j])
						
					#print 'trouble', len(numerator), len(denominator)
					self.axes1.semilogy(self.parent.mass_num[j],(numerator/denominator),linestyle=linestyles[i])
					self.axes1.semilogy(self.parent.mass_num[j],(numerator/denominator),'bo')
					if i == 0:
						try:
							#		Get the coordinates of the highest abundance of each element
							coordinates = [self.parent.mass_num[j][self.abund_plot[i][j].index(max(self.abund_plot[i][j]))],max(numerator/denominator) ]
							self.axes1.text(coordinates[0],coordinates[1], self.parent.elem_list[self.index[j]])
							#print self.parent.elem_list[self.index[j]], coordinates[0]+0.5,coordinates[1]*1.1
						except ValueError:
							None#print 'Empty var:  ', self.abund_plot[i][j]	, self.parent.mass_num[j], self.parent.elem_list[self.index[j]]	
		
			self.axes1.set_ylim([1e-23,10])		
		#	leg1 = self.axes1.legend(l1,self.parent.cycle_to_plot[1:],loc=0)
		
		
		self.canvas.draw()
		


