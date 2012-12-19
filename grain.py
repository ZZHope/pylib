'''grain is a collection of routines to analyze presolar grain data (and can probably be extended to observational 
   data at a later stage). 

   This class load the current version of the presolar grain database for further processing. A private databse can
   be given as well as described later. Several routines (see below and NuGrid book) can be used to filter, plot, and
   retrieve grain data.
   The presolar grain database is supported by the group at Washington University, mainly Frank Gyngard. The database can
   be found at
   http://presolar.wustl.edu/PGD/Presolar_Grain_Database.html
   
   Important note: This script assumes that you have a full SVN tree checked out (or actually, that you have at least the 
   utils folder and the validation folder on the same level checked out.

   Usage of these tools:
   +++++++++++++++++++++

   
   +++
   For questions, bug reports, please contact trappitsch@uchicago.edu
   Reto Trappitsch for the NuGrid collaboration
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
from utils import *
from data_plot import *

class gdb():
   '''This class provides easy access to the presolar grain databse, as described in the header
   The databse is read in by default, however you can choose a private database and if you do so, if you want to use 
   the private database exclusively or together with the whole database.

   arguments:
   fname    -     filename to your private database, if not in the main tree structure with other databse, give full path
   gdbload  -     True or False: Do you want to load the grain database or not? Default is True
   '''

   def __init__(self,fname=None,gdbload=True):
      # grab data
      header_desc, header_data, desc, data = preprocessor(fname,gdbload)
      # make dictionary
      tmp = range(len(header_desc))
      descdict = dict(zip(header_desc,tmp))
      tmp = range(len(header_data))
      datadict = dict(zip(header_data,tmp))

      # make private instances w/ all the data
      self._header_desc = header_desc
      self._header_data = header_data
      self._desc = desc
      self._data = data
      self._descdict = descdict
      self._datadict = datadict
      # make the working data
      self.header_desc = header_desc
      self.header_data = header_data
      self.desc = desc
      self.data = data
      self.descdict = descdict
      self.datadict = datadict

   def __del__(self):
      print 'Presolar grain database available at: http://presolar.wustl.edu/PGD/Presolar_Grain_Database.html'

   def reset_filter(self):
      '''
      Resets the filter and goes back to initialized value
      '''
      self.header_desc = self._header_desc
      self.header_data = self._header_data
      self.desc = self._desc
      self.data = self._data
      self.descdict = self._descdict
      self.datadict = self._datadict


# subroutine that reads in data and splits into nice numpy arrays
def preprocessor(fname,gdbload):
   # path to validation folder 
   scriptpathtmp = __file__   # path of current script
   if len(scriptpathtmp.split('/')) == 1:   # if we are in the utils/pylib folder
      scriptpathtmp = os.path.abspath('.') + '/grain.py'
   svnpathtmp = ''   # current path where the full svn tree is
   for i in range(len(scriptpathtmp.split('/'))-3):   # go up to head folder where utils and validation folder is
      svnpathtmp += scriptpathtmp.split('/')[i] + '/'
   gdbdir = svnpathtmp + 'validation/grain_data/'   # grain data directory

   # Initialize private file if available
   if fname != None:
      data_pri = list()
      f_in = open(gdbdir + fname,'r')
      data_read = f_in.read().splitlines()
      f_in.close()
      for i in range(len(data_read)):
         if len(data_read[i]) > 1:
            data_pri.append(data_read[i].split('\t'))
      print 'Private file ' + fname + ' initialized.'

   
   # Initialize grain database
   if gdbload:
      # load the four files
      data_sic = list()
      data_gra = list()
      data_oxi = list()
      data_mis = list()
      # SiC
      f_in = open(gdbdir + 'SiC-All.txt','r')
      data_read = f_in.read().splitlines()
      f_in.close()
      for i in range(len(data_read)):
         if len(data_read[i]) > 1:
            data_sic.append(data_read[i].split('\t'))
      # Graphites
      f_in = open(gdbdir + 'graphite-All.txt','r')
      data_read = f_in.read().splitlines()
      f_in.close()
      for i in range(len(data_read)):
         if len(data_read[i]) > 1:
            data_gra.append(data_read[i].split('\t'))
      # Oxides and Silicates
      f_in = open(gdbdir + 'oxide-silicate-all.txt','r')
      data_read = f_in.read().splitlines()
      f_in.close()
      for i in range(len(data_read)):
         if len(data_read[i]) > 1:
            data_oxi.append(data_read[i].split('\t'))
      # Miscellaneous grains
      f_in = open(gdbdir + 'miscellaneous-SiN.txt','r')
      data_read = f_in.read().splitlines()
      f_in.close()
      for i in range(len(data_read)):
         if len(data_read[i]) > 1:
            data_mis.append(data_read[i].split('\t'))

      # now bring all files together into one database (if private is not the only file specified)
      header_data = list()   # header for data
      header_desc = list()   # header for description

      # SiC - first file 
      headswtch=True   # switch from description to data 
      for head in data_sic[0]:
         if len(head) != 0:
            if headswtch:   # we have a description header
               header_desc.append(head)
               if len(head) >=5:
                  if head[0:5] == 'Notes':
                     headswtch=False
            else:
               header_data.append(head)
      # Graphites
      headswtch=True   # switch from description to data 
      for head in data_gra[0]:
         if len(head) != 0:
            if headswtch:   # we have a description header
               writeswtch=True
               for i in range(len(header_desc)):   # check if entry already exists
                  if header_desc[i] == head:
                     writeswtch=False
               if writeswtch:
                  header_desc.append(head)
               # the data start after the column 'Notes' (always!)
               if len(head) >=5:
                  if head[0:5] == 'Notes':
                     headswtch=False
            else:
               writeswtch=True
               for i in range(len(header_data)):
                  if header_data[i] == head:
                     writeswtch=False
               if writeswtch:
                  header_data.append(head)
      # Oxides
      headswtch=True   # switch from description to data 
      for head in data_oxi[0]:
         if len(head) != 0:
            if headswtch:   # we have a description header
               writeswtch=True
               for i in range(len(header_desc)):   # check if entry already exists
                  if header_desc[i] == head:
                     writeswtch=False
               if writeswtch:
                  header_desc.append(head)
               # the data start after the column 'Notes' (always!)
               if len(head) >=5:
                  if head[0:5] == 'Notes':
                     headswtch=False
            else:
               writeswtch=True
               for i in range(len(header_data)):
                  if header_data[i] == head:
                     writeswtch=False
               if writeswtch:
                  header_data.append(head)
      # Miscellanious
      headswtch=True   # switch from description to data 
      for head in data_mis[0]:
         if len(head) != 0:
            if headswtch:   # we have a description header
               writeswtch=True
               for i in range(len(header_desc)):   # check if entry already exists
                  if header_desc[i] == head:
                     writeswtch=False
               if writeswtch:
                  header_desc.append(head)
               # the data start after the column 'Notes' (always!)
               if len(head) >=5:
                  if head[0:5] == 'Notes':
                     headswtch=False
            else:
               writeswtch=True
               for i in range(len(header_data)):
                  if header_data[i] == head:
                     writeswtch=False
               if writeswtch:
                  header_data.append(head)
      # Private
      if fname != None:
         headswtch=True   # switch from description to data 
         for head in data_pri[0]:
            if len(head) != 0:
               if headswtch:   # we have a description header
                  writeswtch=True
                  for i in range(len(header_desc)):   # check if entry already exists
                     if header_desc[i] == head:
                        writeswtch=False
                  if writeswtch:
                     header_desc.append(head)
                  # the data start after the column 'Notes' (always!)
                  if len(head) >=5:
                     if head[0:5] == 'Notes':
                        headswtch=False
               else:
                  writeswtch=True
                  for i in range(len(header_data)):
                     if header_data[i] == head:
                        writeswtch=False
                  if writeswtch:
                     header_data.append(head)
      
      # Prepare the data -> description and data, fill it into appropriate forms
      # total amount of data entries
      totdata = len(data_sic) + len(data_gra) + len(data_oxi) + len(data_mis) - 4
      if fname != None:
         totdata += len(data_pri) - 1
      totdata = int(totdata)

      # initialize the description list
      descr = np.zeros((totdata,len(header_desc)),dtype='|S1024')   # string with 1024 characters 

      # initialize data array
      data = np.zeros((totdata,len(header_data)))

      # fill description list
      
      # SiC
      totdata_i = 0  # counter
      for h in range(len(data_sic[0])):
         # description
         for d in range(len(header_desc)):
            if header_desc[d] == data_sic[0][h]:
               for dd in range(1,len(data_sic)):
                  descr[dd-1+totdata_i][d] = data_sic[dd][h]
               break
         # if not a description but data
         for e in range(len(header_data)):
            if header_data[e] == data_sic[0][h]:
               for ee in range(1,len(data_sic)):
                  if data_sic[ee][h] != '':
                     try:
                        data[ee-1][e] = float(data_sic[ee][h])
                     except ValueError:
                        #print data_sic[ee][h] 
                        None
               break
      # Graphite
      totdata_i = len(data_sic)-1  # counter
      for h in range(len(data_gra[0])):
         # description
         for d in range(len(header_desc)):
            if header_desc[d] == data_gra[0][h]:
               for dd in range(1,len(data_gra)):
                  descr[dd-1+totdata_i][d] = data_gra[dd][h]
               break
         # if not a description but data
         for e in range(len(header_data)):
            if header_data[e] == data_gra[0][h]:
               for ee in range(1,len(data_gra)):
                  if data_gra[ee][h] != '':
                     try:
                        data[ee-1+totdata_i][e] = float(data_gra[ee][h])
                     except ValueError:
                        None
                        #print data_gra[ee][h] 
               break
      # Oxides      
      totdata_i = len(data_sic) + len(data_gra) - 2 # counter
      for h in range(len(data_oxi[0])):
         # description
         for d in range(len(header_desc)):
            if header_desc[d] == data_oxi[0][h]:
               for dd in range(1,len(data_oxi)):
                  descr[dd-1+totdata_i][d] = data_oxi[dd][h]
               break
         # if not a description but data
         for e in range(len(header_data)):
            if header_data[e] == data_oxi[0][h]:
               for ee in range(1,len(data_oxi)):
                  if data_oxi[ee][h] != '':
                     try:
                        data[ee-1+totdata_i][e] = float(data_oxi[ee][h])
                     except ValueError:
                        #print data_oxi[ee][h] 
                        None
               break
      # Miscellaneous   
      totdata_i = len(data_sic) + len(data_gra) + len(data_oxi) - 3 # counter
      for h in range(len(data_mis[0])):
         # description
         for d in range(len(header_desc)):
            if header_desc[d] == data_mis[0][h]:
               for dd in range(1,len(data_mis)):
                  descr[dd-1+totdata_i][d] = data_mis[dd][h]
               break
         # if not a description but data
         for e in range(len(header_data)):
            if header_data[e] == data_mis[0][h]:
               for ee in range(1,len(data_mis)):
                  if data_mis[ee][h] != '':
                     try:
                        data[ee-1+totdata_i][e] = float(data_mis[ee][h])
                     except ValueError:
                        #print data_mis[ee][h] 
                        None
               break
      # Private
      if fname != None:
         totdata_i = len(data_sic) + len(data_gra) + len(data_oxi) + len(data_mis) - 4 # counter
         for h in range(len(data_pri[0])):
            # description
            for d in range(len(header_desc)):
               if header_desc[d] == data_pri[0][h]:
                  for dd in range(1,len(data_pri)):
                     descr[dd-1+totdata_i][d] = data_pri[dd][h]
                  break
            # if not a description but data
            for e in range(len(header_data)):
               if header_data[e] == data_pri[0][h]:
                  for ee in range(1,len(data_pri)):
                     if data_pri[ee][h] != '':
                        try:
                           data[ee-1+totdata_i][e] = float(data_pri[ee][h])
                        except ValueError:
                           #print data_pri[ee][h] 
                           None
                  break


      return header_desc, header_data, descr, data
