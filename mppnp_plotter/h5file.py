#!/usr/bin/python
import os
import numpy as np
import time as t
from PyQt4 import QtCore as qc
import functools
import traceback
try:
	import h5py as mrT
except:
	os.system('HDF5_DISABLE_VERSION_CHECK=1')
	import h5py as mrT
try:
	from ascii_table import *
	
except ImportError:
	print 'No module ascii_table'
	print 'Please checkout ascii_table.py svn://forum.astro.keele.ac.uk/utils/pylib and add to python path'
	
debug=True

class h5FileHolder(qc.QThread):
	
	h5sStarted=[]#A list of booleans of if the thread in h5files has had its 
    			#start and join methods called
    	preprocName='h5Preproc.txt' #filename of the preprocessor file
    	
	
	#	This is a list of isotopes to match to.
	isos = ['Neutron','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
    'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
    'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',
    'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te', 'I','Xe','Cs','Ba','La','Ce',
    'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',
    'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At']
    	
    	def findCycle(self,cycNum):
    	    '''
    	    Method that looks through the self.cycles and returns the nearest 
    	    cycle:
    	    
    	    Input:
    	    cycNum: int of the cycle desired cycle 
    	    
    	    '''
    	    
    	    cycNum=int(cycNum)
    	    i=0
    	    
    	    while i < len(self.cycles):
    	    	    if cycNum < int(self.cycles[i]):
    	    	    	    break
    	    	    i+=1
    	    
    	    if i ==0:
    	    	    return self.cycles[i]
    	    elif i == len(self.cycles):
    	    	    return self.cycles[i-1]
	    lower=int(self.cycles[i-1])
	    higher=int(self.cycles[i])
	    
	    if higher- cycNum >= cycNum-lower:
	    	    return self.cycles[i-1]
	    else:
	    	    return self.cycles[i]


   	#	Upon initialization of an h5fileholder, the h5 files are defined (in their own wrapper-see below) and 
	#	and gathers some important data from the files.  
	def __init__(self,filename,filepath, textEdit):
		qc.QThread.__init__(self,None)
		self.textEdit = textEdit

		#	Cleanup the object just in case
		self.h5files = []
		self.h5s = [] #not resetting !!
		self.cycles = []
		self.ages =	 []
		self.hattrs = []
		self.cattrs = []
		self.Tables = []
		self.dcols =	 []
		self.filepaths = []
		self.isotopes = []
		self.isomeric_states = []
		self.A = []
		self.Z = []
		self.done=[]
		self.h5sStarted=[]
		self.filepath=filepath #Variable for the directory
		self.temp = []
		preprocExists=False		#Boolean of if the preprocessor file exists
		
		print filename
		for name in filename:
			self.filepaths.append(name)
		#print 'Opening Files ' + str(self.filepaths)
		t0 = t.time()
		
		#preprocessor stuff
		if self.filepath.endswith(os.sep):
			preprocName = str(self.filepath)+self.preprocName
		else:
			preprocName = str(self.filepath)+os.sep+self.preprocName
		
		
		
		self.preprocExists=os.path.exists(preprocName)
		for i in xrange(len(filename)):
			self.h5sStarted.append(False)
			
		self.h5s.append(h5File(self.filepaths[0],True, True))
		self.connect(self.h5s[-1], qc.SIGNAL('finished()'), self.continue_h5s)
		
		self.h5s[-1].start()
		self.h5sStarted[0]=True
		print "stated:"
		print self.h5sStarted
	def continue_h5s(self):	#tis is equivilent to the run in h5t.py
		if debug:
			print "continue_h5s"
		self.cycles.extend(self.h5s[0].cycle)
		self.ages.extend(self.h5s[0].age)
		
		self.hattrs = self.h5s[0].hattr
		self.cattrs =   self.h5s[0].cattr
		self.Tables    =   self.h5s[0].Table
		self.dcols     =   self.h5s[0].dcol
		self.cycle_header = self.h5s[0].cycle_header	#This string handles the name of the cycle
		new = self.h5s[0].new	#	This boolean handles the changes to the cycle nomenclature format
		
		#preproc stuff
		if self.filepath.endswith(os.sep):
			b = str(self.filepath)+self.preprocName
			
		else:
			b = str(self.filepath)+os.sep+self.preprocName
		
		if self.preprocExists:
			preprocTable=ascii_table(self.preprocName,self.filepath)
			if int(preprocTable.hattrs[0])<len(self.h5files):
				self.preprocExists=False
				print 'A File was added, rewriteing preprocessor file'
        	
        	if self.preprocExists:
			for i in xrange(len(self.h5files)):
				
				if self.h5files[i]+'-cyc' not in preprocTable.dcols and self.preprocExists:
					print 'A File was renamed, rewriteing preprocessor file'
					
					self.preprocExists=False
        
        	if not self.preprocExists and os.path.exists(b):
        		os.system('rm '+b)
        		
		# create list of isotopes stored in this h5 file
		try:
			for x in xrange(len(self.Tables[0])):
				self.isotopes.append([self.isos[int(self.Tables[1][x])],str(int(self.Tables[0][x]))])		
			#print 'new file: ', new, self.isotopes
		except IndexError:
			print self.Tables, self.h5s[0].Table 
		self.t1 = t.time()
		print "continueing"
		if len(self.filepaths) > 1:
			
			for x in range(len(self.filepaths)-1):
				self.h5s.append(h5File(self.filepaths[x+1],False, new))
				if not self.preprocExists:
					self.h5sStarted[x+1]=True
					self.h5s[-1].start()
					self.connect(self.h5s[-1], qc.SIGNAL('finished()'), self.add_data)
		
					
		
		if not self.preprocExists:
			print "all done?"
			self.connect(self, qc.SIGNAL('finished()'), self.all_done)
		else:
			self.all_done()
		
		
	def all_done(self):
		if debug:
			print "all done"
		if not self.preprocExists:
			for x in xrange(len(self.h5s)):
				print len(self.h5s[x].cycle)
				print len(self.h5s[x].age)
			for x in xrange(len(self.h5s)-1):
				self.cycles.extend(self.h5s[x+1].cycle)
				self.ages.extend(self.h5s[x+1].age)
			print len(self.ages)
			print len(self.cycles)
			header=[str(len(self.h5s)),'This is a preprocessor file for the directory: '+str(self.filepath),\
			'At the time of the creation of this file there were '+str(len(self.h5files))+\
			' h5 files.']
			self.cycles = sorted(self.cycles, cmp=self.numeric_compare)
			#self.cycles = sorted(self.cycles, cmp=self.numeric_compare)
			"""
			for cycle in self.cycles:
				print cycle
			try:
				#self.ages = self.get(self.cycles,'age',1)
			except IndexError:
				print 'enountered error fetching age'
			"""
			#self.ages = sorted(self.ages, cmp=self.numeric_compare)		
	
			self.textEdit.append('File search complete.  You may begin plotting')
			
			print 'Writeing preprocessor files'
			data=[]
			dcols=[]
			length=0
			for i in xrange(len(self.h5s)):
				dcols.append(self.h5s[i].filename+'-cyc')
				dcols.append(self.h5s[i].filename+'-age')
				data.append(self.h5s[i].cycle)
				data.append(self.h5s[i].age)
				if len(self.h5s[i].cycle)>length:
					length=len(self.h5s[i].cycle)
				if len(self.h5s[i].age)>length:
					length=len(self.h5s[i].age)
			
			for i in xrange(len(data)):
				for j in xrange(length-len(data[i])):
					data[i].append(3.14159265) #identifier number
			
			write(self.preprocName,header,dcols,data,sldir=self.filepath)
		else:
			print 'Reading preprocessor files'
			preprocTable=ascii_table(self.preprocName,self.filepath)
			
			for i in xrange(len(self.h5s)-1):
				print self.h5s[i+1].filename
				dat=preprocTable.get(self.h5s[i+1].filename+'-cyc')
				dat1=[]
				for j in xrange(len(dat)):
					if dat[j]!=3.14159265:
						dat1.append(dat[j])
					
				dat=dat1
				for j in xrange(len(dat)):
					dat[j]=str(int(dat[j]))
					for k in xrange(10-len(dat[j])):
						dat[j]='0'+dat[j]
				
				for j in xrange(len(dat)):
					self.cycles.append(dat[j])
				self.h5s[i+1].cycle=dat
				dat=preprocTable.get(self.h5s[i+1].filename+'-age')
				dat1=[]
				for j in xrange(len(dat)):
					if dat[j]!=3.14159265:
						dat1.append(dat[j])
				dat=dat1
				self.h5s[i+1].age=dat
				for j in xrange(len(dat)):
					self.ages.append(dat[j])
			try:
			    self.cycles = sorted(self.cycles, cmp=self.numeric_compare)
			except TypeError:
			    print "There was a problem sorting the cycles.  You may have problems later.  Please consider reloading(h5T) and trying again"
			
			try:
			    self.ages = sorted(self.ages, cmp=self.numeric_compare)        
			except TypeError:
			    None
		print self.h5sStarted	
		t2=t.time()
		print "Time: "+str(t2-self.t1)
		
	def numeric_compare(self,x, y):
		if int(x)>int(y):
			return 1
		elif int(x)==int(y):
			return 0
		else: # x<y
			return -1

	
	def add_data(self):
		if debug:
			print "add data"
		self.done.append(self.h5s[-1].filename)
		
		if len(self.done) == len(self.h5s)-1 or len(self.done) == len(self.h5s):
			self.all_done()
		
		
	#	Clears up memory upon file deletion
	def __del__(self):
		print 'File holder destruction event'
		self.terminate()
			
	# This function determines which cycle, which file, which storage mechanism (cattr or data) and returns it 
	def get(self, *args):
		if debug:
			print "get"
		#	This function takes in a variety of inputs
		#	option 1
		#	get(dataitem)
		#		fetches the dataitem for all cycles
		#	option 2
		#	get(cycle_list, dataitem)
		#		fetches the dataitem from the list of cycles
		#	option 3
		#	get(cycle_list, 'iso_massf', isotope)		isotope Must be in the form 'H-2'
		#		fetches the isotope data for a list of cycles

		#	Check out the inputs		
		if len(args) > 4:
			print 'Improper use of this function'
			return None 
		isotope_of_interest = []
		dat = []		
		cycle_list = []

		scale = 1
	#	print args
	#	print len(args)
		print'args', args, self.hattrs, self.cattrs
		
		if len(args) == 2:
			dataitem = args[0]
			if self.hattrs.count(dataitem) == 0:
				print 'a'
				cycle_list = []
				scale = int(args[1])
			else:
				print 'b'
				self.h5s[0] = mrT.File(self.h5s[0].filename,'r')
				dat = self.h5s[0].attrs.get(dataitem, None)
				return dat
		elif len(args) == 3:
			cycle_list = args[0]
			dataitem = args[1]
			scale = int(args[2])
		elif len(args) == 4:
			cycle_list = args[0]
			dataitem = args[1]
			isotope_of_interest = args[2]
			scale = int(args[3])		
		if dataitem=='yps' and dataitem not in self.dcols:
			print self.dcols
			dataitem='iso_massf'
		if dataitem=='iso_massf' and dataitem not in self.dcols:
			print self.dcols
			dataitem='yps'
			
	#	print dataitem

		#	if it is a cattr call, it will not have a cycle value
		if cycle_list == []:
			#	Open all the files.  
			for h5 in self.h5s:
				h5.h5 = mrT.File(h5.filename,'r')
		
			#	Decide which cycles are actually important 
			for x in range(len(self.cycles)/scale):
				cycle_list.append(self.cycles[scale*x])
			
			
			
			#	Fetch the data from each cycle
			for cyc in cycle_list:
				for h5 in self.h5s:
					if h5.cycle.count(cyc) or h5.cycle.count(str(cyc)):
						print 'hello'
						self.h5sStarted[self.h5s.index(h5)]=True
						h5.run()
						try:
							h5.wait()
						except AttributeError:
							print 'cheating again'
						temp = h5.fetch_data_one(dataitem,cyc)
						
						try:
							dat.append(h5.h5[self.cycle_header+str(cyc)].attrs.get(dataitem, None)[0])
						except TypeError:
							dat.append(h5.h5[self.cycle_header+str(cyc)].attrs.get(dataitem, None))
						except IndexError:
							print 'looking in the wrong place'		
			#	Close all the files when done
			for h5 in self.h5s:
				h5.h5.close()			
		else:
			c_len = len(str(self.h5s[0].cycle[0]))
			#print 'c_len', c_len
			if self.h5s[0].new:
				for z in xrange(len(cycle_list)):
					while len(str(cycle_list[z])) < c_len:
						cycle_list[z] = '0'+str(cycle_list[z])	

			#	We already have a list of cycles to check out
			for cyc in cycle_list:
				for h5 in self.h5s:
					if h5.cycle.count(int(cyc)) or h5.cycle.count(str(cyc)):
						
						try:
							if not self.h5sStarted[self.h5s.index(h5)]:
								print "not Sarted"
								self.h5sStarted[self.h5s.index(h5)]=True
								print self.h5sStarted
								h5.run()
								h5.wait()
								
								self.temp = h5.fetch_data_one(dataitem,cyc)
								
							else:
								#print "Sarted"
								self.temp = h5.fetch_data_one(dataitem,cyc)
								#print self.temp
								
						except :
							None
						#	Occasionally data comes out formatted in a funny way (arrays nested in arrays....)
						#	This strips the nested arrays until the actual data is found
						if dataitem != 'iso_massf' or dataitem != 'yps':
							#print 'hello'
							#print self.temp
							'''
							while np.ndim(self.temp) > 1:
								shape = np.shape(self.temp)
								self.temp = self.temp[0]
							'''	
							#print self.temp
						
						else:
							
							while np.ndim(self.temp) > 2:
								shape = np.shape(self.temp)
								self.temp = self.temp[0]
							
							while len(self.temp) < 2:
								self.temp = self.temp[0]
								
						
						try:
							dat.append(self.temp)
						except AttributeError:
							
							np.append(dat, self.temp)
									
									
		if len(dat) < 2 and (dataitem != 'iso_massf' or dataitem != 'yps'):
			try:
				dat = dat[0]
			except IndexError:
				None
			#	print dat
			#	print 'indexerror'
		print dat	
		return dat	
	# This function determines which cycle, which file, which storage mechanism (cattr or data) and returns it 
	def get1(self, *args):
		if debug:
			print "get1"
		#	This function takes in a variety of inputs
		#	option 1
		#	get(dataitem)
		#		fetches the dataitem for all cycles
		#	option 2
		#	get(cycle_list, dataitem)
		#		fetches the dataitem from the list of cycles
		#	option 3
		#	get(cycle_list, 'iso_massf', isotope)		isotope Must be in the form 'H-2'
		#		fetches the isotope data for a list of cycles

		#	Check out the inputs		
		if len(args) > 3:
			print 'Improper use of this function'
			return None 
		isotope_of_interest = []
		
		if len(args) == 1:
			dataitem = args[0]
	
			if self.hattrs.count(dataitem) == 0:
				cycle_list = self.cycles
			else:
				self.h5s[0] = mrT.File(self.h5s[0].filename,'r')
				dat = self.h5s[0].attrs.get(dataitem, None)
				return dat
		elif len(args) == 2:
			cycle_list = args[0]
			dataitem = args[1]
		elif len(args) == 3:
			cycle_list = args[0]
			dataitem = args[1]
			isotope_of_interest = args[2]		
		
			
		#	Just in case the user inputs integers 
		try:
			for x in xrange(len(cycle_list)):
				cycle_list[x] = str(cycle_list[x])
		except TypeError:
			cycle_list = [str(cycle_list)]


		try:
			if cycle_list.isdigit():
				cycle_list = [cycle_list]
				for cycle in cycle_list:
					if len(cycle) != len(self.cycles[0]):

						diff = len(self.cycles[0])-len(cycle)
						OO = ''
						while diff >=1:
							OO+='0'
	
						cycle = OO+cycle

		except AttributeError:
			if cycle_list[0].isdigit():
								
				for x in xrange(len(cycle_list)):
					if len(str(cycle_list[x])) != len(str(self.cycles[0])):
	
						diff = len(str(self.cycles[0]))-len(str(cycle_list[x]))
						
						OO = ''
						while diff >=1:
							OO+='0'
							diff-=1
					
						try:
							cycle_list[x] = OO+cycle_list[x]
						except TypeError:
							cycle_list[0] = OO+cycle_list[0]

		#	if it is a cattr call, it will not have a cycle value
		dat = []		
	#	for h5 in self.h5s:
	#		h5.h5 = mrT.File(h5.filename,'r')
	#'/rpod2/fherwig/tmp/tmp/'
		for cyc in cycle_list:
			for h5 in self.h5s:
				if h5.cycle.count(int(cyc)) or h5.cycle.count(str(cyc)):
					
					if not self.h5sStarted[self.h5s.index(h5)]:
						self.h5sStarted[self.h5s.index(h5)]=True
						
						h5.run()
						
						#try:
						h5.wait()
						#except:
						#	print 'failed thread:',  h5	
						temp = h5.fetch_data_one(dataitem,cyc)
						
					else:
						temp = h5.fetch_data_one(dataitem,cyc)
					#	Occasionally data comes out formatted in a funny way (arrays nested in arrays....)
					#	This strips the nested arrays until the actual data is found
					
					
					
					#print 'temp', temp
					

					
					#else:
						
					#	while np.ndim(temp) > 2:
					#		shape = np.shape(temp)
					#		temp = temp[0]
					#	
					#	while len(temp) < 2:
					#		temp = temp[0]
					if (dataitem == 'iso_massf' or dataitem == 'yps') and isotope_of_interest != []: #
						#	Figure out the index
						#print 'yps', dataitem
						index = 0
						for x, iso in enumerate(self.isotopes):
							print str(iso[0]+'-'+iso[1]), isotope_of_interest
							if str(iso[0]+'-'+iso[1]) == isotope_of_interest:
								index = x
								break
						#print 'iso_massf',temp
						temp = temp[:,index]
					#	Now add the information to the list we pass back
					elif (dataitem=='iso_massf' or dataitem=='yps'):
						print 'the right stuff', isotope_of_interest
					else:
						#print 'enter', dataitem
						while np.ndim(temp) > 1:
							shape = np.shape(temp)
							temp = temp[0]
						

					try:
						dat.append(temp)
						#print 'right append'
					except AttributeError:
						np.append(dat, temp)	
						#print 'bad append'
													
		if len(dat) < 2 and (dataitem != 'iso_massf'or dataitem != 'yps'):

			try:
				dat = dat[0]
			except IndexError:
				None
			except TypeError:
				None
		try:
			if len(dat) < 2 and isotope_of_interest != []:
				dat = dat[0]
		except TypeError:
			None	
		except IndexError:
			None
			
		
		
		print self.h5sStarted
		return dat	
		
	#	Determines the file-cycle match and gets the associated information
	def fetch_datas(self,dataitem1,dataitem2,cycle, scale):
		if debug:
			print "fetch_datas"
		dat = []
		if cycle == None:
			dat1 = []
			dat2 = []			
			index = 0		
			#	Loop through all the files and grab the appropriate data
			for x in xrange(len(self.h5s)):
				self.h5s[x].h5 = mrT.File(self.h5s[x].filename,'r')
			for y in xrange(len(self.cycles)/scale):
				for x in xrange(len(self.h5s)):
					try:
						node = self.h5s[x].h5[self.cycle_header+str(self.cycles[scale*y])].attrs
						dat1.append(node.get(dataitem1, None))
						dat2.append(node.get(dataitem2, None))	
						#print node
					except (IndexError,KeyError):
						1+1
						#print 'bad cycle', self.cycles[scale*y]
						#print self.h5sStarted														
				#for y in xrange(len(self.h5s[x].cycle)/scale):
				#	node = self.h5s[x].h5[self.cycle_header+str(self.h5s[x].cycle[scale*y])].attrs
				#	dat1.append(node.get(dataitem1, None))
				#	dat2.append(node.get(dataitem2, None))	
			for x in xrange(len(self.h5s)):	
				self.h5s[x].h5.close()
			try:
				while len(dat1) == 1:
					print 'trimmed dat1'
					#print dat1
					dat1 = dat1[0]
			except AttributeError:
				while dat1.shape[0] == 1:
					print 'trimmed dat1'
					#print dat1
					dat1 = dat1[0]

			try:
				while len(dat2) == 1:
					print 'timmed dat2'
				#	print dat2
					dat2 = dat2[0]
			except AttributeError:
				while dat2.shape[0] == 1:
					print 'timmed dat2'
				#	print dat2
					dat2 = dat2[0]					
			
			
			dat = [dat1,dat2]
		
		else:
			#	We have a user-defined list of cycles to check out
			try:
				for cyc in cycle:
					for h5 in self.h5s:
						if h5.cycle.count(int(cyc)):
							dat.append(h5.fetch_data(dataitem1,dataitem2,cyc))
						else:
							print 'not Matched'
			except TypeError:
				for h5 in self.h5s:
					if h5.cycle.count(cycle):
						dat.append(h5.fetch_data(dataitem1,dataitem2,cycle))
		return dat		
	
	def fetch_sp_kh_data(self, cycle_list, scale, ranger, element):
		if debug:
			print "fetch_sp_kh_data"
		t0 = t.time()
		for h5 in self.h5s:
			h5.h5 = mrT.File(h5.filename,'r')
	
		t2 = t.time()	
		
		Z = np.zeros([int(len(cycle_list)), 3000], float)
		zindex = 0

		for cyc in cycle_list:
			for h5 in self.h5s:
				if h5.cycle.count(cyc) or h5.cycle.count(str(cyc)):
					if element == -1:
						goal = "dcoeff"
					else:
						goal = "iso_massf"
					
					try:
						node = h5.h5[self.cycle_header+str(cyc)].__getitem__('SE_DATASET')
					except ValueError:
						node = h5.h5[self.cycle_header+str(cyc)]
					#temp1 = node.__getitem__("dcoeff")
					try:
						temp1 = node.__getitem__(goal)#node.col('convection_indicator')
					except:
						temp1 = node.__getitem__("yps")					
					temp2 = node.__getitem__("mass")#node.col('mass')	
					
					
					
					
					try:
						for x in xrange(len(temp2)):
							if float(temp1[x][element]) >= float(limit):
								try:
									Z[zindex][x] = temp2[x]
								except:
									print 'indexerror'
									break
						zindex += 1
					except IndexError:
						for x in xrange(len(temp2)):
					
							if float(temp1[x]) >= float(limit):
							
								try:
									Z[zindex][x] = temp2[x]
						
								except:
									print 'indexerror'
									break
						zindex += 1						
						
						
					#break
		for h5 in self.h5s:
			h5.h5.close()
		return Z
	
	
	def fetch_KH_data(self, scale, ranger, cycle_list):
		#cycle_list = []
		#age = []
		if debug:
			print "fetch_KH_data"
		t0 = t.time()
		for h5 in self.h5s:
			h5.h5 = mrT.File(h5.filename,'r')
		t1 = t.time()
		print 'Opening files took: ', str(t1-t0), scale
	
		t2 = t.time()
		print 'Buidling cycle List took: ', str(t2-t1)
		

		
		Z = np.zeros([int(len(cycle_list)), 3000], float)	# 	holds the convective information
			
		zindex = 0

		print 'starting check', len(cycle_list)
		
		for cyc in cycle_list:
			for h5 in self.h5s:
				if h5.cycle.count(cyc):
					#print  'checking cycle:', cyc
					try:
						try:
							node = h5.h5[self.cycle_header+str(cyc)].__getitem__('SE_DATASET')
						except ValueError:
							node = h5.h5[self.cycle_header+str(cyc)]
					except:
						print 'encountered a problem with cycle:',cyc
						break
						
					
					temp1 = node.__getitem__("convection_indicator")#node.col('convection_indicator')
					temp2 = node.__getitem__("mass")#node.col('mass')	
		
					
					maxi =  abs(max(temp1))
				#	print 'maxi', maxi
					
					#temp1 = abs(temp1/(temp1-1e-10))	
					#if any(temp1) > 1:
					for x in xrange(len(temp1)):
						if temp1[x] != 1:
							temp1[x] = 0
							
					
					
					Z[zindex][:len(temp2)] = temp1*temp2
					
					zindex += 1

							
						
		t3 = t.time()
		print 'Data fetch took: ', str(t3-t2)
		
		
		for h5 in self.h5s:
			h5.h5.close()	
			
		t4 = t.time()
		print 'file closing took: ', str(t4-t3)
	
		return Z#conv
	
		
	def fetch_sparse_yps(self,block_plot, scale,isotope_index,textEdit):
		if debug:
			print "fetch_sparse_yps"
		self.cycles.sort()
		data1 = []
		data2 = []
		all_mass = []
		cycle_list = []
		#scale = 4
		block_plot=True
		isotope_index = 0

		t0 = t.time()
		limits = [0.3645,0.04956,8.5e-4,2.4082e-3]
		
		for h5 in self.h5s:
			h5.h5 = mrT.File(h5.filename,'r')
		
		for x in range(len(self.cycles)/scale):
			cycle_list.append(self.cycles[scale*x])

		t1 = t.time()
		textEdit.append('Open Time: ' + str(t1-t0))
		textEdit.append('cycling...............')
		for cyc in cycle_list:
			for h5 in self.h5s:
				if h5.cycle.count(cyc):
					try:
						node = h5.h5['/cycle-'+str(cyc)].__getitem__('SE_DATASET')
		
						if block_plot:
							temp1 = []
							temp2 = []
							temp3 = []
							
							temp1 = node.__getitem__("convection_indicator")#node.col('convection_indicator')
							temp2 = node.__getitem__("mass")#node.col('mass')		
							#temp1 = np.dot(temp1,temp2)
							for x in range(len(temp1)):
								temp3.append(abs(temp2[x]*float(temp1[x])))

							data1.append(temp3)
							data2.append(temp2)
										
						
						else:
							try:							
								temp1 = node.__getitem__("iso_massf")#col('iso_massf')
							except:
								temp1 = node.__getitem__("yps")
							temp2 = node.__getitem__("mass")#col('mass')
							fake1 = []
							fake2 = []
							for x in range(len(temp1)/scale):
								fake1.append(temp1[scale*x][isotope_index])
								fake2.append(temp2[scale*x])			
							data.append(fake1)
							all_mass.append(fake2)
							del fake1, fake2, temp1, temp2
							
					except AttributeError:
						print 'Age does not exist'
		t2 = t.time()			
		textEdit.append( 'Cycle Time: ' + str(t2-t1))
		textEdit.append( 'closing a bunch of files')
		for h5 in self.h5s:
			h5.h5.close()
			textEdit.append(str( h5.h5))
		t3 = t.time()
		textEdit.append( 'Close TIme: ' + str(t3-t2))

		return [data2,data1]
		
	#	uses the index information to build list of isotopes from tables A,Z	
	def fetch_isotopes(self):
		isos = []
		try:
			for x in range(len(self.Tables[1])):
				isos.append([self.isotopes[int(self.Tables[1][x])], self.Tables[0][x]])
		except IndexError:
			None
		return isos
	
	
		
#	This wrapper class allows some automated activity when an h5file is initialized.
#	upon inmitialization the h5 file is opened and certain bits of data is read.  
#	This class also interacts with the h5fileholder class to access needed data.

class h5File(qc.QThread):
	h5 = None
	filename = None
	cycle = []
	age = []
	dcol = []
	hattr = []
	data = []
	skipped_nodes = 0
	ver = ''
	classname = ''
	cattr=[]
	Table = []
	isomeric_state = []
	new = True
	


	#	Initialize the class
	def __init__(self, filepath,deep_search, new):
		qc.QThread.__init__(self,None)
		print 'starting'
		#	Instantiate
		self.h5 = None
		self.filename = None
		self.cycle = []
		self.age = []
		self.dcol = []
		self.data = []
		self.skipped_nodes = 0
		self.ver = ''
		self.classname = ''
		self.hattr = []
		self.cattr=[]
		self.Table = []
		self.isomeric_state = []
		self.A = []
		self.Z = []
		self.new = True
		
		#	Build
		self.new = new
		self.filename = filepath
		self.deep_search = deep_search		
		
		self.filename = filepath

		if self.new:
			self.cycle_header = 'cycle'
		else:
			self.cycle_header = 'cycle-'
	
	def run(self):
	
		if self.deep_search:
			self.search_deep()
		else:
			self.search_shallow()
		try:
			self.h5.close()
		except:
			None#print 'error closing file: ', self.h5
		print 'done'
		#print 'Table: ', self.Table
	
	#	Fetches a single category of information
	def fetch_data_one(self,dataitem,cycle):
		#print 'fetching data one'
		self.h5 = mrT.File(self.filename,'r')
		
		try:
			data = self.h5.__getitem__(self.cycle_header+str(cycle)).__getitem__('SE_DATASET')[dataitem]
		except ValueError:
			try:
				data = self.h5.__getitem__(self.cycle_header+str(cycle)).attrs.get(dataitem, None)[0]
			except TypeError:
				data = self.h5.__getitem__(self.cycle_header+str(cycle))[dataitem]
				#print data

	
		try:
			while data.shape[0] < 2:
				data = data[0]
		except IndexError:
			None	
		
		self.h5.close()
		return data
	#	same as above, but for 2 sets of data		
	def fetch_data(self,dataitem1,dataitem2,cycle):
		self.h5 = mrT.File(self.filename,'r')
		#print 'cycle ',cycle
		
		try:
			dataset = self.h5.__getitem__(self.cycle_header+str(cycle)).__getitem__('SE_DATASET')
			
			data1 = dataset[dataitem1]
			data2 = dataset[dataitem2]
			
		
		except ValueError:
			
			dataset = self.h5.__getitem__(self.cycle_header+str(cycle)).attrs.values()
			data1 = dataset[self.cattr.index(dataitem1)][0]
			data2 = dataset[self.cattr.index(dataitem2)][0]

							
		
		self.h5.close()
		del dataset
		data = [data1,data2]
		return data	

	#	The typical search algirthm when a h5file class is initialized
	def search_shallow(self):
		self.h5 = mrT.File(self.filename,'r')        
		temp = self.h5.keys()
		
		for te in temp:
		    if te[0] == 'c':
			if te[5:].isdigit():
			    self.cycle.append(str(te[5:]))
			    
			    try:
				self.age.append(self.h5[te].attrs.get("age",None)[0])
			    except TypeError:
				self.age.append(self.h5[te].attrs.get("age",None))
			else:
			    self.cycle.append(str(te.split('-')[1]))
			    try:
				self.age.append(self.h5[te].attrs.get("age", None)[0])
			    except TypeError:
				self.age.append(self.h5[te].attrs.get("age",None))
				self.cycle.sort()        
				self.age.sort()

	def search_deep(self):
		self.h5 = mrT.File(self.filename,'r')        
		temp = self.h5.keys()
		
		#    Handles the change in cycle nomenclature
		self.new = True
		for te in temp:
		    if te.count('-'):
			self.new = False
			break
		
		if self.new:
		    self.cycle_header = 'cycle'
		    for te in temp:
			if te[0] == 'c':
			    if te[5:].isdigit():
				self.cycle.append(str(te[5:]))
				try:
				    self.age.append(self.h5[te].attrs.get("age",None)[0])
				except TypeError:
				    self.age.append(self.h5[te].attrs.get("age",None))
			    else:
				self.isomeric_state.append(self.h5[te]['data'])
			else:
	
			    obj = self.h5[te].__iter__()
			    if str(te).count('A'):
				holder = []
				for ob in obj:
				    holder.append(ob[0])
				self.Table.append(holder)
				self.A.append(holder)
			    elif str(te).count('Z'):
				holder = []
				for ob in obj:
				    holder.append(ob[0])
				self.Table.append(holder)
				self.Z.append(holder)    
			    else:
				holder = []
				for ob in obj:
				    holder.append(ob[0])
				self.Table.append(holder)
				self.isomeric_state.append(holder)    
			    
		    try:
			temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).__getitem__('SE_DATASET').dtype.__str__().split(',')
		    except ValueError:
			temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).dtype.__str__().split(',')
		
		else:        
		    self.cycle_header = 'cycle-'
		    for te in temp:
			try:
			    self.cycle.append(str(te.split('-')[1]))
	
			    try:
				self.age.append(self.h5[te].attrs.get("age", None)[0])
			    except TypeError:
				self.age.append(self.h5[te].attrs.get("age", None))
			except IndexError:
		    
			    obj =  self.h5[te].__iter__()
		    
			    if str(te).count('A'):
				holder = []
				for ob in obj:
				    holder.append(ob[0])
				self.Table.append(holder)
				self.A.append(holder)
			    elif str(te).count('Z'):
				holder = []
				for ob in obj:
				    holder.append(ob[0])
				self.Table.append(holder)
				self.Z.append(holder)    
			    else:
				holder = []
				for ob in obj:
				    holder.append(ob[0])
				self.Table.append(holder)
				self.isomeric_state.append(holder)                        
				
		    self.cycle.sort()    
		
		# This is kind of stupid, but I have not found a way to access this information directly.
		
		    try:
			temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).__getitem__('SE_DATASET').dtype.__str__().split(',')
		    except ValueError:
			temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).dtype.__str__().split(',')
	
		
		for tem in temp:
		    if tem.count('<') ==0:
			try:
			    self.dcol.append(tem.split('\'')[1])
			except IndexError:
			    None
	
		attrs = self.h5.attrs
		for at in attrs:
		    self.hattr.append(at)
		self.cattr = self.h5[self.cycle_header+str(self.cycle[0])].attrs.keys()
		
		table = []
		grp = self.h5[self.cycle_header+str(self.cycle[0])]
		for gr in grp:
		    try:
			table.append(float(gr[0]))
		    except ValueError:
			None
		
		self.h5.close()
		return None


        
        	
	"""	old search deep
	def search_deep(self):
		self.h5 = mrT.File(self.filename,'r')		
		temp = self.h5.keys()
		
		#	Handles the change in cycle nomenclature
		self.new = True
		for te in temp:
			#print te
			if te.count('-'):
				#print 'checked'
				self.new = False
				#print self.new

		attrs = self.h5.attrs
		for at in attrs:
			self.hattr.append(at)
		#print self.cycle
		
		
		if self.new:
			self.cycle_header = 'cycle'
			for te in temp:
				if te[0] == 'c':
					#print 'cycle', te
					if te[5:].isdigit():
						self.cycle.append(te[5:])
					else:
						self.isomeric_state.append(self.h5[te]['data'])
				else:
					#print 'table', te
					obj = self.h5[te].__iter__()
					holder = []
					for ob in obj:
						holder.append(ob[0])
					self.Table.append(holder)
					
			try:
				#print self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).__getitem__('SE_DATASET').dtype
				temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).__getitem__('SE_DATASET').dtype.__str__().split(',')
			except ValueError:
				temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).dtype.__str__().split(',')
				#print temp
		
		else:		
			self.cycle_header = 'cycle-'
			for te in temp:
				try:
					self.cycle.append(te.split('-')[1])

				except IndexError:

					obj =  self.h5[te].__iter__()
					holder = []
					for ob in obj:
						#print ob
						holder.append(ob[0])
						#print holder
					self.Table.append(holder)
		
		
		
		# This is kind of stupid, but I have not found a way to access this information directly.
		
			try:
				temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).__getitem__('SE_DATASET').dtype.__str__().split(',')
			except ValueError:
				temp =  self.h5.__getitem__(self.cycle_header+str(self.cycle[0])).dtype.__str__().split(',')
				#print temp
		#print self.cycle	
		
		for tem in temp:
			if tem.count('<') ==0:
				try:
					self.dcol.append(tem.split('\'')[1])
				except IndexError:
					None

		self.cattr = self.h5[self.cycle_header+str(self.cycle[0])].attrs.keys()
		#print self.cattr
		table = []
		grp = self.h5[self.cycle_header+str(self.cycle[0])]
		for gr in grp:
			try:
				table.append(float(gr[0]))
			except ValueError:
				None
	#	print self.cattr, self.age
		self.cycle.sort()	
		print 'end search deep'
		self.h5.close()
	"""
	
	def __del__(self):
		print 'Deleting H5File'
		self.terminate()
		print self.h5
		

#t0 = t.time()
	
#f = h5FileHolder(['M2.00Z0.010.0000001.out.h5','M2.00Z0.010.0001001.out.h5'], 'new_mesa_file/')

#datas = f.fetch_datas('logL','R_sol', f.cycles[0], 1)
#print datas
#t1 = t.time()
#print 'Total Time: ' + str(t1-t0)	
		
		
		
