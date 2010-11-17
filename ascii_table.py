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
Data columns are seperated by spaces

Here is an example run through.
bash-4.1$ python
Python 2.6.4 (r264:75706, Jun  4 2010, 18:20:16) 
[GCC 4.4.4 20100503 (Red Hat 4.4.4-2)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from ppm import *
>>> p=AsciiTable('c12cg.dat')
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

class AsciiTable(DataPlot):
	'''
	Data structure to read simple data tables	
	'''
	files = []
	sldir = ''
	hattrs=[]
	dcols=[]
	data ={}
	def __init__(self,fileName,sldir='.',sep='  '):
		'''
		Init method
		Input:
		sldir: Standard directory of fileName
		fileName: The name of the file we are looking at
		sep: The seperator that seperates column attributes in FileName
		     Defaults to '  '
		'''
		self.sldir=sldir
		self.files.append(fileName)
		
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
		while i<len(fileLines):
			tmp=fileLines[i].split(' ')
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
				
						
			if j == 0:
				dataCols={cols[j]:tmp}
			else:
				dataCols[cols[j]]=tmp
			tmp=[]
			
		return header,dataCols
