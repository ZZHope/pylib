"""
Ascii_table.py

By Daniel Alexander Bertolino Conti
Fall 2010
If the user find any bugs or errors, please email 

This module will read in simple ascii tables and store the data within them

Assumptions:
Headers are always at the beginning of the file and start with a capital H
The next line after the header lines, is a line of column attribute names.
The column attribute names are separated by '  ' by default or whatever the 
	user dictates to the class.
All the data columns are of equal length.
Any file name that has 'trajectory' or 'Trajectory' in it, is assumed to be a 
	Trajectory type file

Here is an example run through.
bash-4.1$ python
Python 2.6.4 (r264:75706, Jun  4 2010, 18:20:16) 
[GCC 4.4.4 20100503 (Red Hat 4.4.4-2)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ppm import *
>>> p=ascii_table('c12cg.dat')
>>> p.hattrs
['1 12  6  1  1  1  0  0  0  1 13  7  0', '55   1.943', 'c12pg']
>>> a.dcols
['upper', 'lower', 'T9', 'ado(gs)', 'ado/CA88', 'tt/gs', 'ado/fit', 'CA88', 'fitted']
>>> a.get('upper')
[1.3400000000000001e-24, ... 1590.0]
>>>a.plot('T9', 'ado/CA88')
plots data
"""
from numpy import *
from data_plot import *
import matplotlib.pylab as pyl
import matplotlib.pyplot as pl
import os

class ascii_table(DataPlot):
	'''
	Data structure to read simple data tables and trajectory data tables	
	'''
	files = []
	sldir = ''
	hattrs=[]
	dcols=[]
	data ={}
	dataType=''
	headerLines=[]
	def __init__(self,fileName,sldir='.',sep='  ', dataType='normal',headers=[],dcols=[],data=[],read=True,headerLines=[]):
		'''
		Init method that reads in ascii type files and trajectory type files.
		By default this method reads ascii type files.  If the user wants
		a trajectory file read, either the file must have 'trajectory'
		in the filename or the user must set the dataType='trajectory'
		
		Input:
		sldir: Standard directory of fileName
		fileName: The name of the file we are looking at, or writeing to
		read: Boolean of weather this is reading or writing a file.
			Defaults to Reading
		sep: The seperator that seperates column attributes in FileName
		     Defaults to '  '
		dataType: What type of ascii table  you are reading and or writeing
			  The only two options currently are 'normal' and 'trajectory'
		Headers: A list of Header strings or if the file being writtin 
			 is  of type trajectory, this is a dictionary of header 
			 attributes and their associated values
		dcols: A list of data attributes
		data:  A list of lists (or of numpy arrays) of columns of data
		headerLines: Additional list of strings of header data, only used
			     in trajectory data Types.
		'''
		self.sldir=sldir
		self.files.append(fileName)
		self.dataType=dataType
		if 'trajectory' in fileName or 'Trajectory' in fileName:
			self.dataType='trajectory'
		
		if read:
			
			self.hattrs,self.data=self._readFile(sldir,fileName,sep)
			self.dcols=self.data.keys()
		else:
			a=self.write(fileName,headers,dcols,data,headerLines,sldir,sep)
			if a ==None:
				return None
			self.hattrs,self.data=self._readFile(sldir,fileName,sep)
		self.dcols=self.data.keys()

	def get(self, attri):
		'''
		Method that dynamically determines the type of attribute that is 
		passed into this method. Also it then returns that attribute's 
		associated data.
		Input:
		attri: The attribute we are looking for.
		'''
		isCol=False
		isHead=False
		
		if attri in self.dcols:
			isCol=True
		elif attri in self.hattrs:
			isHead=True
		else:
			print "That attribute does not exist in this File"
			print 'Returning None'
		
		if isCol:
			return self.getColData(attri)
		elif isHead:
			return hattrs	
	
	def getColData(self,attri):
		'''
		Method that returns column data
		Input:
		attri: The attribute we are looking for.
		'''
		return self.data[attri]
	
	def _readFile(self,sldir,fileName,sep):
		'''
		Private method that reads in the header and column data
		'''
		
		if sldir.endswith('/'):
			fileName = str(sldir)+str(fileName)
		else:
			fileName = str(sldir)+'/'+str(fileName)
		
		
		fileLines=[] #list of lines in the file
		header=[]    #list of Header lines
		dataCols=[]  #Dictionary of data column names
		data=[]	     #List of Data lists
		cols=[]      #List of column names
		
		f=open(fileName,'r')
		fileLines=f.readlines()
		i=0
		if self.dataType != 'trajectory':
			
			while i<len(fileLines):
				if fileLines[i].startswith('H'):
					tmp=fileLines[i].lstrip('H')
					header.append(tmp.strip())
				else:
					break
				i+=1
			
			cols=fileLines[i].split(sep)
			
			tmp=[]
			tmp1=[]
			for j in range(len(cols)):
					tmp1=cols[j].strip()
					if tmp1 !='':
						tmp.append(tmp1)
			cols=tmp		
			i+=1
		else:
			header={}
			while fileLines[i].startswith('#') or '=' in fileLines[i]:
				if fileLines[i].startswith('#') and cols==[]:
					cols=fileLines[i].strip('#')
					cols=cols.strip()
					cols=cols.split()
				elif fileLines[i].startswith('#'):
					tmp1=fileLines[i].strip('#')
					tmp1=tmp1.strip()
					self.headerLines.append(tmp1)
				elif not fileLines[i].startswith('#'):
					tmp=fileLines[i].split('=')
					tmp[0]=tmp[0].strip()
					tmp[1]=tmp[1].strip()
					if header=={}:
						header={str(tmp[0]):str(tmp[1])}
					else:
						header[str(tmp[0])]=str(tmp[1])
				i+=1
		while i<len(fileLines):
			tmp=fileLines[i].split(sep)
			for j in xrange(len(tmp)):
				tmp[j]=tmp[j].strip()
			data.append(tmp)
			i+=1
		tmp=[]
		tmp1=[]
		for j in range(len(data)):
			for k in range(len(data[j])):
				tmp1=data[j][k].strip()
				if tmp1 !='':
					tmp.append(tmp1)
			data[j]=tmp
			tmp=[]
		tmp=[]
		#print cols
		#print data
		for j in range(len(cols)):
			for k in range(len(data)):
				tmp.append(float(data[k][j]))
			tmp=array(tmp)	
						
			if j == 0:
				dataCols={cols[j]:tmp}
			else:
				dataCols[cols[j]]=tmp
			tmp=[]
			
		return header,dataCols

#Global methods

		
def write(fileName,headers,dcols,data,headerLines=[],sldir='.',sep='  ',trajectory=False):
		'''
		Method for writeing Ascii files.
		Note the attribute name at position i in dcols will be associated
		with the column data at index i in data.
		Also the number of data columns(in data) must equal the number
		of data attributes (in dcols)
		Also all the lengths of that columns must all be the same.
		Input:
		fileName: The file where this data will be written.
		Headers: A list of Header strings or if the file being writtin 
			 is  of type trajectory, This is a dictionary of header 
			 attributes and their associated values
		dcols: A list of data attributes
		data:  A list of lists (or of numpy arrays).
		headerLines: Additional list of strings of header data, only used
			     in trajectory data Types.
		sldir: Where this fill will be written.
		sep: What seperatesa the data column attributes
		trajectory: Boolean of if we are writeing a trajectory type file
		'''
		if sldir.endswith('/'):
			fileName = str(sldir)+str(fileName)
		else:
			fileName = str(sldir)+'/'+str(fileName)
		tmp=[] #temp variable
		lines=[]#list of the data lines
		lengthList=[]# list of the longest element (data or column name)
			     # in each column
			     
		if os.path.exists(fileName):
			print 'Warning this method will overwrite '+ fileName
			print 'Would you like to continue? (y)es or (n)no?'
			s = raw_input('--> ')
			if s=='Y' or s=='y' or s=='Yes' or s=='yes':
				print 'Yes selected'
				print 'Continuing as normal'
			else:
				print 'No Selected'
				print 'Returning None'
				return None
		
		if len(data)!=len(dcols):
			print 'The number of data columns does not equal the number of Data attributes'
			print 'returning none'
			return None
		if trajectory:
			
			keys=headers.keys()
			sep=' '
		for i in xrange(len(headers)):
			if not trajectory:
				tmp.append('H '+headers[i]+'\n')
			else:
				tmp.append(str(keys[i])+ ' = '+str(headers[keys[i]])+'\n')
		headers=tmp
		tmp=''
		
		for i in xrange(len(data)): #Line length stuff
			length=len(dcols[i])
			for j in xrange(len(data[i])):
				if len(str(data[i][j]))>length:
					length=len(str(data[i][j]))
			lengthList.append(length)
		
		tmp=''
		tmp1=''
		if trajectory:
			tmp='#'
		for i in xrange(len(dcols)):
			tmp1=dcols[i]
			if not trajectory:
				if len(dcols[i]) < lengthList[i]:
					j=lengthList[i]-len(dcols[i])
					for k in xrange(j):
						tmp1+=' '
				tmp+=sep+tmp1
			else:
				tmp+=' '+dcols[i]
		tmp+='\n'
		dcols=tmp
		tmp=''
		
		
		for i in xrange(len(data[0])):
			for j in xrange(len(data)):
				tmp1=str(data[j][i])
				if len(str(data[j][i])) < lengthList[j]:
					l=lengthList[j]-len(str(data[j][i]))
					for k in xrange(l):
						tmp1+=' '
				tmp+=sep+tmp1
			lines.append(tmp+'\n')
			tmp=''
			
		f=open(fileName,'w')
		if not trajectory:
			for i in xrange(len(headers)):
				f.write(headers[i])
			f.write(dcols)
		else:
			f.write(dcols)
			for i in xrange(len(headerLines)):
				f.write('# '+headerLines[i]+'\n')
			for i in xrange(len(headers)):
				f.write(headers[i])
		for i in xrange(len(lines)):
			f.write(lines[i])
		
		f.close()
		return None
