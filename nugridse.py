'''nugridse is a collection of plots of data in se-type h5 files

   sedir    - the directory on which the h5 files can be
   pattern - a string pattern that the file names have to
                    contain in order to be read

   

   usage:
   start by loading the module
   1 : import nugridse as mp

   you can get help:
   2 : help mp

   next, initiate a se class instance:   
   3 : pt=mp.se('path/to/dir/with/se-files')
   
   which would read all h5 files from the current directory; again, do
   help(mp.se) or mp.se? do get info for more options, as for example
   how to specify a match pattern to select only certain files in the
   directory; this may be useful if there are very many files and you
   want only every 10th; example:k
   3 : pt =mp.se('.','M1.65Z0.020.00') 

      note: initializing an instance with 120 files with 1000 packets
      each may take 3-4 minutes, which is too long for many people;
      therefore the initialization module will generate and write an
      index file; if such an index file is present initialization will
      take less then a second; the initialization module will
      automatically detect various situations in which the index file
      needs to be rewritten, for example if there is a new file.

   look at the data that is available in your instance
   4 : pt.sedir 
   5 : pt.sefiles

   The actual data is in pt.se and the following commands give you
   access to header and cycle attributes as well as the cycle profile
   data:
      pt.se.cattrs
      pt.se.hattrs
      pt.se.dcols

   Available cycles can be viewed with
      pt.se.cycles

   You can get any of the quantities via the 'get' method, which is
   relatively smart to give you things in various ways:
      pt.get('rho')

   would give you all rho vectors for all cycles, which is maybe more
   than you want, so try/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0012001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0013001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0014001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0015001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0016001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0017001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0018001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0019001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0020001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0021001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0022001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0023001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0024001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0025001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0026001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0027001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0028001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0029001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0030001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0031001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0032001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0033001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0034001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0035001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0036001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0037001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0038001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0039001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0040001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0041001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0042001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0043001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0044001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0045001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0046001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0047001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0048001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0049001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0050001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0051001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0052001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0053001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0054001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0055001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0056001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.0057001.se.h5
/rpod2/fherwig/M/work_operator0/Mop0_2.0MsunZ0.01/M2.00Z0.010.idx
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0000001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0001001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0002001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0003001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0004001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0005001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0006001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0007001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0008001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0009001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0010001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0011001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0012001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0013001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0014001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0015001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0016001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0017001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0018001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0019001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0020001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0021001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0022001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0023001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0024001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0025001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0026001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0027001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0028001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0029001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0030001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0031001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0032001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0033001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0034001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0035001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0036001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0037001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0038001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0039001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0040001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0041001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0042001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0043001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0044001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0045001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0046001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0047001.se.h5
/rpod2/fherwig/M/work_operator0_fullnet/Mop0_2.0MsunZ0.01/M2.00Z0.010.0048001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0000001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0001001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0002001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0003001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0004001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0005001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0006001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0007001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0008001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0009001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0010001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0011001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0012001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0013001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0014001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0015001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0016001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0017001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0018001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0019001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0020001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0021001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0022001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0023001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0024001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0025001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0026001.se.h5
/rpod2/fherwig/M/work_operator7/Mop7_2.0MsunZ0.01/M2.00Z0.010.0027001.se.h5

      pt.get(300,'rho')
   to get the rho vector for cycle 300, instead of one cycle you may
   also supply a cycle list for the first argument

   Do help(pt) (or whatever you instance is) for a full description of
   what else your instance has to offer, which includes methods to
   work with data

      note: a particularly nice feature is that you can work with
      instances of different types of data in a very flexible way, for
      example in lists:
         cases=[pt1,pt2,pt3]
         for this_case in cases:
            do something with this_case

   There are various plotting methods, including plotting abundance
   charts (abu_chart), isotopic abundance distributions and the
   generic 'plot' method that lets you plot every quantity against any
   other quantity that can possibly be plotted. Some of the methods
   (including the three just mentioned) are available via the super
   class data_plot, and are equally available in other python
   modules, as for example mesa.py or ppn.py. Several methods accept
   lists of cycles which implies to create a series of frames to be
   written to disk in order to make movies.

   6 : pt.plot('mass','rho',fname=3000)

   whereas pt.plot('mass','rho',fname=[3000,4000]) will produce two png files,
   one for each cycle, do help(m.plot) for a full list of options
   
   mppnp output allows to do plots with abundances:
   10: pt.iso_abund(23500)
   11: pt.abu_chart(23500)

   Also for iso_abund and abu_chart: if instead of a single
	cycle the user inputs a list of cycles, the method will then,
	instead of plotting them, will then save a .png for each
	cycle. Also if you just want a singular plot saved, the user
	can input their cycle, in a list like [0]. And that will save
	their plot.
'''

import h5T
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.pylab as pyl
from matplotlib.ticker import MultipleLocator
import os  
import re
import sys
import time
import glob
from utils import *
from data_plot import *

class se(DataPlot,Utils):
    ''' This class provides easy access to h5 files from the NuGrid
        project, along with some standard plots
        
    arguments: 
    sedir    - the directory on which the h5 files can be
    pattern - a string pattern that the file names have to
                    contain in order to be read
    rewrite - a boolean of if the user would like to rewrite the preprocessor
    	      file. Defaults to False 

    example: f=nu.plot_tools('.','260') reads all h5 files with the
                  string 260 in the name in the present directory
    '''

    sedir = ''
    sefiles = []
    se = []         # main data dictionary
    pattern=''

    def __init__(self, sedir='.', pattern='.h5', rewrite=False):
        
        slist = os.listdir(sedir)
        self.pattern=pattern
        expr = re.compile(pattern)
       
        sefiles=filter(expr.search,slist)
        
        self.se = []         # main data dictionary
        self.sedir=sedir
        self.sefiles=sefiles
        self._stable_names() # provides in addition to stable_el from 
                             # utils also just the stable element names
        
        self.se=h5T.Files(sedir,sefiles,rewrite=rewrite)
                
    def __del__(self):
        print 'Closing plot_tools'
        
    def get(self, cycle_list,dataitem=None,isotope=None, sparse=1):
    	'''
    	Simple function that simply calls h5T.py get method.
        There are three ways to call this function
            option 1
            get(dataitem)
                fetches the dataitem for all cycles, interperates the argument 
                cycle list, as Data Item
           option 2
           get(cycle_list, dataitem)
                fetches the dataitem from the list of cycles. If Dataitem is an 
                isotope, it then returns self.get(cycle_list,'iso_massf',dataitem)
           option 3
            get(cycle_list, 'iso_massf', isotope)
               isotope Must be in the form 'H-2'
                fetches the isotope data for a list of cycles
         Additional Input:
        	sparse - implements a sparsity factor on the fetched data
        '''
        return self.se.get(cycle_list,dataitem,isotope,sparse)


    def get_decayed(self, cycle_list, dataitem=None, isotope=None, sparse=1):
        '''
        This function gives back the fully decayed isotope. By default, it wants to
        beta-decays every isotope, however, then checks if this is okay or not. If not
        it actually does goes into a database (see routine below) to do the proper decay.
        This database might not be complete, so pay attention to the chart of the 
        nuclides.
        
        Standard input is as in regular get function:
        "There are three ways to call this function
            option 1
            get(dataitem)
                fetches the dataitem for all cycles, interperates the argument 
                cycle list, as Data Item
           option 2
           get(cycle_list, dataitem)
                fetches the dataitem from the list of cycles. If Dataitem is an 
                isotope, it then returns self.get(cycle_list,'iso_massf',dataitem)
           option 3
            get(cycle_list, 'iso_massf', isotope)
               isotope Must be in the form 'H-2'
                fetches the isotope data for a list of cycles
         Additional Input:
        	sparse - implements a sparsity factor on the fetched data"

         Additional features are:
         
         IN PROGRESS: CONTACT RETO: trappitsch@uchicago.edu
            - 
            - 

        '''
        return None

        
    def plot_prof_1(self,mod,species,xlim1,xlim2,ylim1,ylim2,symbol=None):

        ''' plot one species for cycle between xlim1 and xlim2
        If m
        species      - which species to plot
        mod          - model to plot
        xlim1, xlim2 - mass coordinate range 
        ylim1, ylim2 - mass fraction coordinate range
	symbol       - indicate which symbol you want to use, if required. '''
        DataPlot.plot_prof_1(self,species,mod,xlim1,xlim2,ylim1,ylim2,symbol)
        '''
        tot_mass=self.se.get(mod,'total_mass')    
        age=self.se.get(mod,'age')    
        mass=self.se.get(mod,'mass')    
        Xspecies=self.se.get(mod,'iso_massf',species)
        pyl.plot(mass,np.log10(Xspecies),'-',label=species)
        pyl.xlim(xlim1,xlim2)
        pyl.ylim(ylim1,ylim2)
        pyl.legend()

        pl.xlabel('$Mass$ $coordinate$', fontsize=20)
        pl.ylabel('$X_{i}$', fontsize=20)
        pl.title('Mass='+str(tot_mass)+', Time='+str(age)+' years, cycle='+str(mod))
        '''

    def plot_prof_2(self,mod,species,xlim1,xlim2):

        ''' plot one species for cycle between xlim1 and xlim2

        no log
        species      - which species to plot
        mod          - model to plot
        xlim1, xlim2 - mass coordinate range 
        '''

        mass=self.se.get(mod,'mass')    
        Xspecies=self.se.get(mod,'yps',species)
        pyl.plot(mass,Xspecies,'-',label=str(mod)+', '+species)
        pyl.xlim(xlim1,xlim2)
        pyl.legend()
    def write_deltatable(self,filename='default',decayed=True,dcycle=500,iniabufile='../../frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn'):
        '''
        This subroutine is to write out tables with delta values for cosmochemists to use in comparison with their data. No options are necessarily needed to load this routine, however, it might be useful to specify a filename.
        This file furthermore searches for thermal pulses, hence, it's only useful for non-explosive TP-AGB stars!
        filename: choose the filename you want to set here
        decayed:  Value if decayed massfractions of isotope should be taken or not
        dcycle:   difference between cycles to search for thermal pulses
        iniabufile: file with initial abundances. as a standard value, the GN93 file is used (in USEEPP folder). Important, input file has to be USEEPP conform -> see NuGrid book!
        '''
        print 'This is a preliminary version - contact Reto for more information!'
        # make filename
        filename = 'delta_outfile.txt'
        # which iso_massmf?
        if decayed:
            isomassmf = 'iso_massf_decay'
        else:
            isomassmf = 'iso_massf'
        # read in thermal pulse position and co_ratio (from private routine)
        tp_pos, co_return = self._tp_finder(dcycle)
        # define isotopes to read out and calculate delta values
        noele = 1   # there is at least one element available
        elelist = list()
        isolisttmp = list()
        isolist = list()
        isoabulist = list()
        isoabulisttmp = list()
        # load file w/ initial abundances
        inut = iniabu(iniabufile)
        # list of elements:
        itmp = inut.names[0][0:2].replace(' ','')
        elelist.append([itmp, inut.z[0]])
        isolisttmp.append(inut.a[0])
        isoabulisttmp.append(inut.abu[0])

        for i in range(1,len(inut.names)):
            if inut.names[i][0:2].replace(' ','') != itmp:
                # add to elelist
                itmp = inut.names[i][0:2].replace(' ','')
                elelist.append([itmp, inut.z[i]])       
                # append to isolist
                isolist.append(isolisttmp)
                isolisttmp = list()
                isolisttmp.append(inut.a[i])
                # append to isoabulist
                isoabulist.append(isoabulisttmp)
                isoabulisttmp = list()
                isoabulisttmp.append(inut.abu[i])
                # increase number of elements
                noele += 1
            else:
                # add to isolisttmp
                isolisttmp.append(inut.a[i])
                isoabulisttmp.append(inut.abu[i])
        isolist.append(isolisttmp)
        isoabulist.append(isoabulisttmp)
   
        ### CHOP TEMPORARILY ###
        elelist = elelist[5:len(elelist)-2]
        isolist = isolist[5:len(elelist)-2]
        isoabulist = isoabulist[5:len(elelist)-2]
        # make an index vector for at which position for each line the most abundand isotope is
        isoabuindex = zeros(len(isoabulist),dtype=int)   # as integer!
        for i in range(len(isoabulist)):
            tmp = 0
            for j in range(len(isoabulist[i])):
                if isoabulist[i][j] > tmp:
                    isoabuindex[i] = j
                    tmp = isoabulist[i][j]
         
        # and now - finally - make a list with all delta values that need to be calculated
        # 0: isotope name, 1: nominator isotope, 2: denominator isotope
        # also make solar system ratios list
        deltalist = list()
        solsysratio = list()
        for i in range(len(isolist)):
            if len(isolist[i]) > 1:
                for j in range(len(isolist[i])):
                    if j != isoabuindex[i]:
                        deltalist.append([elelist[i][0], isolist[i][j], isolist[i][isoabuindex[i]]])
                        solsysratio.append(isoabulist[i][j] / isoabulist[i][isoabuindex[i]])
        solsysratio = array(solsysratio)
        # Make an array to write all data into. first column, cycle number, second column, co ratio, then all delta values
        write_out = zeros((len(tp_pos), len(deltalist) + 2))

        # now add first three columns to write_out
        for i in range(len(write_out)):
            write_out[i][0] = tp_pos[i]
            write_out[i][1] = co_return[i]
   
        # header for delta values
        deltaheader = list()
   
        ### BIG DELTA LOOP ###
        for deli in range(len(deltalist)):
            ## read in the two isotopes of interest, iso1 and iso2 (iso1 as nominator)
            iso1 = zeros(len(tp_pos))
            iso2 = zeros(len(tp_pos))
            # isotope names:
            if len(deltalist[deli][0]) == 1:
                iso1name = deltalist[deli][0].upper() + '-' + str(int(deltalist[deli][1]))
                iso2name = deltalist[deli][0].upper() + '-' + str(int(deltalist[deli][2]))
            else:
                iso1name = deltalist[deli][0][0].upper() + deltalist[deli][0][1] + '-' + str(int(deltalist[deli][1]))
                iso2name = deltalist[deli][0][0].upper() + deltalist[deli][0][1] + '-' + str(int(deltalist[deli][2]))
            # readin loop
            for i in range(len(tp_pos)):
                iso1[i] = self.get(int(tp_pos[i]),isomassmf,iso1name)
                iso2[i] = self.get(int(tp_pos[i]),isomassmf,iso2name)
         
            ## make number fractions
            iso1 /= deltalist[deli][1]
            iso2 /= deltalist[deli][2]
      
            ## calculate deltavalue array
            deltavalues = (iso1/iso2 / solsysratio[deli] - 1.) * 1000.

            ## append values to write_out
            for i in range(len(deltavalues)):
                write_out[i][deli+2] = deltavalues[i]
         
            ## header
            deltaheader.append('del_' + deltalist[deli][0] + str(int(deltalist[deli][1])) + '_' + deltalist[deli][0] + str(int(deltalist[deli][2])))


        ### WRITE OUT DATAFILE ###
        outf = open(filename,'w')

        # write header
        outf.writelines('cycle_no,'.ljust(18) + 'c_o_ratio,'.ljust(18))
        for i in range(len(deltaheader)):
            if i != len(deltaheader)-1:
                outfstr = deltaheader[i] + ','
            else:
                outfstr = deltaheader[i]
            outf.writelines(outfstr.ljust(18))
        outf.writelines('\n')
   
        # write data
        for i in range(len(write_out)):
            for j in range(len(write_out[i])):
                if j != len(write_out[i])-1:
                    outfstr = format(write_out[i][j],'.7E') + ','
                else:
                    outfstr = format(write_out[i][j],'.7E')
                outf.writelines(outfstr.ljust(18))
            outf.writelines('\n')
   
        # close file
        outf.close()


    def plot_isoratios(self,xiso,yiso,graintype=None,tp_finding=False,deltax=True,deltay=True,logx=False,logy=False,title=None,legend=True,dcycle=500,errbar=True,iniabufile='iniab2.0E-02GN93.ppn'):
        ''' returns the ratios of model surface data for validation w/ grain data
        rt, 2012
        
        xiso:       give isotopes as ['Fe',57,'Fe',56]. x axis, yaxis
        yiso:       as x iso, but for y axis
        graintype:	Example: [['sic','M'],['oxides','1','2']], see utils, graindata_handler routine
        tp_finding:  Enables the finding of thermal pulses - especially for grain people important. Calls the private routine _tp_finder -> see below.
        deltax:     True for plotting delta values on x-axis. If False, plots ratios not normalized to solar
        deltay:     As deltax but for y axis
        logx:       Logarithm of x axis
        logy:       Logarithm of y axis
        title:      Title of the plot
        legend:     Legend in plot or not?
        dcycle:     Difference between cycles to take for thermal pulse searching, if searching is deactivated, dcycle describes how often cycles are sampled.
        errbar:     Do you want to plot error bars with the grain data? True or False
        iniabufile: Solar abundances that are used to normalize too. Lives in ../../frames/mppnp/USEEPP, can also be given absolute to avoid to put it in USEEPP relative to frame. Attention: we assume a standard svn tree (see relative path)!
        '''
        ### WORK ON PATH ###
        # define svn path form path where script runs, depending on standard input or not
        if len(iniabufile.split('/')) == 1 :   # means not an absolute path!
            scriptpathtmp = __file__
            if len(scriptpathtmp.split('/')) == 1:   # in folder where nugridse is
                scriptpathtmp = os.path.abspath('.') + '/nugridse.py'   # to get the current dir
            svnpathtmp = '/'
            for i in range(len(scriptpathtmp.split('/'))-3):   # -3 to go to folders up! 
                if scriptpathtmp.split('/')[i] != '':
                    svnpathtmp += scriptpathtmp.split('/')[i] + '/'
            iniabufile = svnpathtmp + 'frames/mppnp/USEEPP/' + iniabufile   # make absolute path for iniabufile
        # read in thermal pulse position and co_ratio (from private routine)
        print 'Need to get a lot of data - this might take a while the first time'
        if tp_finding: 
            tp_pos, co_return = self._tp_finder(dcycle) 
        else:
            # define equally spaced cycles to sample data, according to dcycle, we also call it tp_pos for now
            tp_pos = range(1000,int(self.se.cycles[len(self.se.cycles)-1])+dcycle,dcycle)   # start at cycle 1000, in case some weird things go on before
            # get c and o and form ratio of final vector
            c_nf = array(self.get(tp_pos,'elem_numf','C'))
            o_nf = array(self.get(tp_pos,'elem_numf','O'))
            co_return = c_nf / o_nf
        ## read in isotope data
        isox1 = array(self.get(tp_pos,xiso[0] + '-' + str(int(xiso[1]))))
        isox2 = array(self.get(tp_pos,xiso[2] + '-' + str(int(xiso[3]))))
        isoy1 = array(self.get(tp_pos,yiso[0] + '-' + str(int(yiso[1]))))
        isoy2 = array(self.get(tp_pos,yiso[2] + '-' + str(int(yiso[3]))))
        # ratios - number fractions!
        isox_ratio = (isox1 / xiso[1]) / (isox2 / xiso[3])
        isoy_ratio = (isoy1 / yiso[1]) / (isoy2 / yiso[3])
        # load initial abundance file if necessary
        inut = iniabu(iniabufile)
        # make delta values if necessary
        if deltax:
            deltax_solsys = inut.isoratio_init(xiso)
            ret_x = (isox_ratio / deltax_solsys - 1.) * 1000.
        else:
            ret_x = isox_ratio   # just the ratio
            deltax_solsys = None
        if deltay:
            deltay_solsys = inut.isoratio_init(yiso)
            ret_y = (isoy_ratio / deltay_solsys - 1.) * 1000.
        else:
            ret_y = isoy_ratio   # just the ratio
            deltay_solsys = None
        ### get graindata ###
        if graintype == None:
            graindata = None
        else:
            graindata = graindata_handler(xiso,isosy=yiso,graintype_in=graintype,deltax=deltax,deltay=deltay,iniabufile_in=iniabufile)
        ### now plot it using the appropriate too ###
        DataPlot.plot_ratios(self,graindata=graindata,misosx=ret_x,misosy=ret_y,solsysx=deltax_solsys,solsysy=deltay_solsys,m_co=co_return,misosxname=xiso,misosyname=yiso,deltax=deltax,deltay=deltay,logx=logx,logy=logy,title=title,legend=legend,iniabufile=iniabufile,errbar=errbar)
        
    def _tp_finder(self,dcycle):   # Private routine
        '''
        Routine to find thermal pulses in given star and returns an index vector that gives the cycle number in which the thermal pulse occure. The routine looks for the C/O ratio jumping up and up, so only useful in TP-AGB star. A vector is given back that indicates the position of the cycle that is at 95% of the thermal pulse (to make sure it's not in the next one and that most of the processing is done).
        The script also returns the co_ratio vector - the C/O ratio (number fraction) at the given thermal pulse.
        '''
        # read in c and o isotopes for all cycles, regarding deltacycle
        last_cycle = int(self.se.cycles[len(self.se.cycles)-1])
        cyc_tp = range(1,last_cycle + dcycle, dcycle)
        c_nf = array(self.get(cyc_tp,'elem_numf','C'))
        o_nf = array(self.get(cyc_tp,'elem_numf','O'))
        
        # search for thermal pulses
        co_ratio = c_nf / o_nf
        tp_guess     = 200   # this should be an upper limit!
        tp_guess_max = 200   # to through an error
        # guess variables, i is the actual break criterion, n a max counter
        gi = 0
        gn = 0
   
        while gi != 1 and gn < 10000:
            tp_ind = list()
            i = 0
            while i < len(co_ratio)-2:
                gcompar= 1. / (dcycle*tp_guess*100.)
                slope1 = (co_ratio[i+1]-co_ratio[i])/(dcycle)
                slope2 = (co_ratio[i+2]-co_ratio[i+1])/dcycle
                if slope1 > gcompar and slope2 < gcompar and co_ratio[i+1] > co_ratio[i]:
                    tp_ind.append(i+1)
                    i += 3   # jump three cycles to avoid defining a single cycle twice!
                else:
                    i += 1
    
            if abs(len(tp_ind) - tp_guess) < tp_guess/2:   # gotta be within factor two of guess
                gi = 1
            else:
                gn += 1
                tp_guess /= 2
        # check w/ maximum of thermal pulses allowed
        if len(tp_ind) > tp_guess_max:
            print 'Problem detected with number of pulses'
        # create thermal pulse vector
        tp_startf = zeros(len(tp_ind))   # found start
        for i in range(len(tp_startf)):
            tp_startf[i] = cyc_tp[tp_ind[i]]   
        # read out isotopic composition at 95% of the thermal pulse and the initial of the star
        # set up thermal pulse positions
        tp_limits = zeros(len(tp_startf)+1)
        for i in range(len(tp_startf)):
            tp_limits[i] = tp_startf[i]
        tp_limits[len(tp_limits)-1] = int(self.se.cycles[len(self.se.cycles)-1])
        # thermal pulse position (where to read the isotope ratio)
        tp_pos = list()
        for i in range(len(tp_startf)):
            tp_pos.append(int(tp_limits[i] + 0.95 * (tp_limits[i+1] - tp_limits[i]))) 
        # create co_ret vector to return c/o ratio vector
        co_return = zeros(len(tp_pos))
        for i in range(len(tp_pos)):
            co_return[i] = co_ratio[tp_ind[i]]
        # return the two vectors
        return tp_pos,co_return


    def plot4iso_exp(self,isotope_list,shift=0,graintype=None,deltax=False,deltay=False,logx=False,logy=False,addiso=None,weighting=None,co_toggle='c',pl_title=None,modlegend=None,errbar=True,plt_show=True,plt_symb='o',plt_col='b',plt_markersize=None,plt_linewidth=None,plt_sparse=10,plt_massrange=False,iniabufile='iniab2.0E-02GN93.ppn'):
        '''
        This subroutine plots 4 isotope plots for explosive stars, assuming one model initialized in nugridse
        Description here on how it is done. It plots the model along w/ presolar grain data if wanted.
        Grain data are plotted along with the models (if wanted)
        Change once grainclass is done

        Input:
        - isotope_list:   List of isotopes, format ['C-12','C-13','N-14','N-15']. First pair: nominator, denominator of x-axis, second pair: y-axis
        - shift:          By how much do you want to shift the models back. If zero, the last cycle is taken
        - graintype:      Graintype, e.g., 'all', 'sic', ['sic','oxides'], [['sic','N']], ... see graindata_handler in utils.py
        - deltax, deltay: Plot delta values on  x and y axes?
        - logx, logy:     Logarithmic axes?
        - addiso:         Add an isotope. Format ['C-12', 0.5 ,'N-12'] to add N12 to C12 and multiply it w/ a factor of 0.5. Multiple isotopes can be added, the factor is optional and does not have to be given. Isotopes can be added to other isotopes as well, i.e., [['C-12','N-12'],['C-13','N-13']]
        - weighting:      None -> plot every profile separately, 'zone' -> average each zone, 'all' -> average all selected zones
        - co_toggle:      Select 'c' for selecting zones with C/O >= 1, select 'o' for C/O <= 1, 'a' takes the whole star
        - pl_title:       Plot title
        - modlegend:      Legend for model, if None the script uses the folder the data are in as a title
        - errbar:         Do you want errorbars for presolar grain data? True or False
        - plt_show:       Show the plot or not. Set to False if you're overprinting several models.
        - plt_symb:       Defines the symbol to plot the model data with
        - plt_col:        Defines the color of the plotted model symbol
        - plt_markersize: Array that defines all the markersized. If you use it, make sure there are enough sizes given!
        - plt_linewidth:  Linewidth for everything, as markersize and symbols. Make sure you have enough values
        - plt_sparse:     Every so many datapoint is plotted for model data
        - plt_massrange:  Plot mass of shell with first and last datapoint of each zone.
        - iniabufile:     File where the initial abundances are stored. Give absolute path or relative to USEEPP!
        '''

        ### WORK ON PATH ###
        # define svn path form path where script runs, depending on standard input or not
        if len(iniabufile.split('/')) == 1 :   # means not an absolute path!
            scriptpathtmp = __file__
            if len(scriptpathtmp.split('/')) == 1:   # in folder where nugridse is
                scriptpathtmp = os.path.abspath('.') + '/nugridse.py'   # to get the current dir
            svnpathtmp = '/'
            for i in range(len(scriptpathtmp.split('/'))-3):   # -3 to go to folders up! 
                if scriptpathtmp.split('/')[i] != '':
                    svnpathtmp += scriptpathtmp.split('/')[i] + '/'
            iniabufile = svnpathtmp + 'frames/mppnp/USEEPP/' + iniabufile   # make absolute path for iniabufile
        ### END PATH WORK ###
        
        ### Input compatibility ###
        co_toggle = co_toggle.lower()

        # cycle
        cyc_no = self.se.cycles[len(self.se.cycles)-1-shift]

        mass =   array(self.get(cyc_no,'mass'))
        c_elem = array(self.get(cyc_no,'C-12'))
        c_elem +=array(self.get(cyc_no,'C-13'))
        o_elem = array(self.get(cyc_no,'O-16'))
        o_elem +=array(self.get(cyc_no,'O-17'))
        o_elem +=array(self.get(cyc_no,'O-18'))
        co_ratio = c_elem / o_elem * (16. / 12.)
        # get isotopes
        isotope_profile = []
        for i in range(4):
            isotope_profile.append(array(self.get(cyc_no,isotope_list[i])))

        # add radioactive isotopes (if given)
        if addiso != None:
            if type(addiso[0] == list):   # then list of lists
                for i in range(len(addiso)):
                    for j in range(4):
                        if isotope_list[j] == addiso[i][0]:
                            multiplicator_addiso = 1.
                            try:
                                multiplicator_addiso = float(addiso[i][1])
                                starter = 2
                            except ValueError:
                                starter = 1
                            for k in range(starter,len(addiso[i])):
                                isotope_profile[j] += array(self.get(cyc_no,addiso[i][k])) * multiplicator_addiso
            else:
                for j in range(4):
                    if isotope_list[j] == addiso[0]:
                        multiplicator_addiso = 1.
                        try:
                            multiplicator_addiso = float(addiso[1])
                            starter = 2
                        except ValueError:
                                starter = 1
                        for k in range(starter,len(addiso)):
                            isotope_profile[j] += array(self.get(cyc_no,addiso[k])) * multiplicator_addiso

        # search for carbon / oxygen rich layers
        crich = []   # alternating start stop values. if odd number, then surface is c-rich, but add stop number
        dumb = True
        if co_toggle != 'a':
            for i in range(len(co_ratio)):
                if co_toggle == 'c':   # carbon rich
                    if dumb:
                        if co_ratio[i] >= 1:
                            crich.append(i)
                            dumb = False
                            continue
                    else: 
                        if co_ratio[i] < 1:
                            crich.append(i)
                            dumb = True
                elif co_toggle == 'o':   # oxygen richt
                    if dumb:
                        if co_ratio[i] <= 1:
                            crich.append(i)
                            dumb = False
                            continue
                    else: 
                        if co_ratio[i] > 1:
                            crich.append(i)
                            dumb = True
                else:
                    print 'Select your enrichment!'
                    return None
        else:   # take whole star
            print 'Using all profiles to mix'
            crich.append(0)
            crich.append(len(co_toggle))

        if len(crich)%2 == 1:
            crich.append(len(co_ratio)-1)

        if len(crich) == 0:
            print 'Star did not get rich in C or O, depending on what you specified'
            return None
        
        # make isotope_profile into array and transpose
        isotope_profile = array(isotope_profile).transpose()

        
        # Ask user which zones to use
        if co_toggle != 'a':
            if co_toggle == 'c':
                print '\n\nI have found the following carbon rich zones:\n'
            elif co_toggle == 'o':
                print '\n\nI have found the following oxygen rich zones:\n'
            
            mass_tmp = zeros((len(crich)))
            for i in range(len(crich)):
                mass_tmp[i] = mass[crich[i]]

            j = 1
            for i in range(len(crich)/2):
                print 'Mass range (' + str(j) + '):\t' + str(mass_tmp[2*i]) + ' - ' + str(mass_tmp[2*i+1]) 
                j += 1
            print '\n'
            usr_zones = input('Please select which mass range you want to use. Select 0 for all zones. Otherwise give one zone or a list of zones separated by comma (e.g.: 1, 2, 4): ')

            crich_dumb = crich
            if usr_zones == 0:
                print 'I continue w/ all zones then'
            elif type(usr_zones) == int:   # only one zone selected
                tmp = int(usr_zones)-1
                crich = crich_dumb[2*tmp:2*tmp+2]
            else:
                crich = []
                for i in range(len(usr_zones)):
                    tmp = int(usr_zones[i])-1
                    crich.append(crich_dumb[2*tmp])
                    crich.append(crich_dumb[2*tmp + 1])

        # weight profiles according to weighting factor using the selected crich
        # define isos_to_use variable for later
        if weighting == None:
            isos_to_use = []
            for i in range(len(crich)/2):
                isos_dumb = []
                n = crich[2*i]
                while n <= crich[2*i+1]:
                    isos_dumb.append(isotope_profile[n])
                    n += 1
                isos_to_use.append(array(isos_dumb))


        elif weighting.lower() == 'zone' or weighting.lower() == 'zones':
            # make array w/ mass weigted isotope ratio (4) for all mass zones
            isotope_profile_cweight = zeros((len(crich)/2,4))
            mass_tot = []
            for i in range(len(isotope_profile_cweight)):   # 2*i is start, 2*i+1 is stop value
                if crich[2*i] == 0:
                    print 'C- / O-rich in first shell (core).'
                else:
                    dumb = crich[2*i + 1]
                    j = crich[2*i]
                    mass_tmp = 0
                    while j <= dumb:
                        mass_shell = mass[j] - mass[j-1]
                        mass_tmp += mass_shell
                        for k in range(4):
                            isotope_profile_cweight[i][k] += isotope_profile[j][k]*mass_shell
                        j += 1
                    mass_tot.append(mass_tmp)
            for i in range(len(isotope_profile_cweight)):
                for j in range(4):
                    isotope_profile_cweight[i][j] /= mass_tot[i]
            isos_to_use = [array(isotope_profile_cweight)]

        elif weighting.lower() == 'all':   # average all zones by mass
            isos_tmp = zeros((1, len(isotope_profile[0])))
            for i in range(len(isotope_profile)-1):   # neglect surface effects
                for j in range(len(isos_tmp[0])):
                    mass_shell = mass[i+1] - mass[i]
                    isos_tmp[0][j] += isotope_profile[i][j]*mass_shell
            # weight all
            isos_tmp /= sum(mass)
            isos_to_use = [isos_tmp]

        # legend for model
        if modlegend == None:
            modlegend = self.se.filename.split('/')
            modlegend = modlegend[len(modlegend)-1]   # depth of folder

        # change to isotope numbers from mass!
        for i in range(len(isos_to_use)):
            for j in range(len(isos_to_use[i])):
                for k in range(len(isos_to_use[i][j])):
                    isos_to_use[i][j][k] /= float(isotope_list[k].split('-')[1])

        # do the ratios and stuff
        inut = iniabu(iniabufile)
        ratiox = []
        ratioy = []
        for i in range(len(isos_to_use)):
            ratiox_dumb = []
            ratioy_dumb = []
            for j in range(len(isos_to_use[i])):
                ratiox_dumb.append(isos_to_use[i][j][0] / isos_to_use[i][j][1])
                ratioy_dumb.append(isos_to_use[i][j][2] / isos_to_use[i][j][3])
            ratiox.append(array(ratiox_dumb))
            ratioy.append(array(ratioy_dumb))
        
        ratiox_solsys = inut.isoratio_init(isotope_list[0:2])
        ratioy_solsys = inut.isoratio_init(isotope_list[2:4])
        if deltax:
            ratiox = (ratiox / ratiox_solsys - 1.) * 1000.
        if deltay:
            ratioy = (ratioy / ratioy_solsys - 1.) * 1000.
        
        # create massrange array if necessary
        if plt_massrange:
            plt_massrange = []
            for i in range(len(ratiox)):
                plt_massrange.append([ratiox[i][0], ratioy[i][0], mass[crich[2*i]], ratiox[i][len(ratiox[i])-1], ratioy[i][len(ratioy[i])-1], mass[crich[2*i+1]]])   # start: x-ratio, y-ratio, mass label, then same for stop

        ### get graindata ###
        if graintype == None:
            graindata = None
        else:
            graindata = graindata_handler(isotope_list[0:2],isosy=isotope_list[2:4],graintype_in=graintype,deltax=deltax,deltay=deltay,iniabufile_in=iniabufile)

        ### send to data_plot.py -> plot_ratios ###
        DataPlot.plot_ratios(self,ratiox,ratioy,solsysx=ratiox_solsys,solsysy=ratioy_solsys,graindata=graindata,m_co=None,misosxname=isotope_list[0:2],misosyname=isotope_list[2:4],deltax=deltax,deltay=deltay,logx=logx,logy=logy,title=pl_title,legend=True,iniabufile=iniabufile,modlegend=modlegend,calling_routine='4iso_exp',plt_symb=plt_symb,plt_col=plt_col,plt_modms=plt_markersize,plt_modlw=plt_linewidth,plt_sparse=plt_sparse,plt_show=plt_show,plt_mrng=plt_massrange,errbar=errbar)


    def plot4(self,num):
        self.plot_prof_1(num,'H-1',0.,5.,-5,0.)
        self.plot_prof_1(num,'He-4',0.,5.,-5,0.)
        self.plot_prof_1(num,'C-12',0.,5.,-5,0.)
        self.plot_prof_1(num,'O-16',0.,5.,-5,0.)
        pyl.legend(loc=3)

    def plot4_nolog(self,num):
        self.plot_prof_2(num,'H-1',0.,5.)
        self.plot_prof_2(num,'He-4',0.,5.)
        self.plot_prof_2(num,'C-12',0.,5.)
        self.plot_prof_2(num,'O-16',0.,5.)
        pyl.legend(loc=3)

    def plot_prof_sparse(self,mod,species,xlim1,xlim2,ylim1,ylim2,sparse,symbol):

        ''' plot one species for cycle between xlim1 and xlim2

        species      - which species to plot
        mod          - model to plot
        xlim1, xlim2 - mass coordinate range 
        ylim1, ylim2 - mass fraction coordinate range 
        sparse       - sparsity factor for points 
        symbol (str) - which symbol you want to use ? '''
        mass=self.se.get(mod,'mass')    
        Xspecies=self.se.get(mod,'yps',species)
        pyl.plot(mass[0:len(mass):sparse],np.log10(Xspecies[0:len(Xspecies):sparse]),symbol)
        pyl.xlim(xlim1,xlim2)
        pyl.ylim(ylim1,ylim2)
        pyl.legend()

    def trajectory(self,ini,end,delta,mass_coo,age_in_sec=False):
        
        ''' create a trajectory out of a stellar model 
        
        ini       - initial model
        end       - final model 
        delta     - sparsity factor of the frames
        mass_coo  - mass coordinate for the traj
        age_in_sec (default False) set to True if age in se file is in seconds (like in MESA)

        Output:
        returns radius_at_mass_coo, density_at_mass_coo, temperature_at_mass_coo, age_all
        
        plus writes a file with the trajectory information to be used with ppn
        warning: remove the old trajectory, if you have any for the same mass coordinate.
        you are appending data, not overwriting.

        Update: this method works for output types with indexes going
        from the outside in (MESA) or the other way around. Also the
        requested quantities are linearly interpolated in the mass
        shell.

        '''
        
        f = open('traj_'+str(mass_coo)+'.dat','a')
        radius_at_mass_coo=[]
        density_at_mass_coo=[]
        temperature_at_mass_coo=[]
        age_all=[]
        for step in range(ini,end+1,delta):
                age=self.se.get(step,'age')
                if age_in_sec:
                    age /= constants.one_year
                mass=self.se.get(step,'mass')  
                temperature=self.se.get(step,'temperature')
                rho=self.se.get(step,'rho')
                radius=self.se.get(step,'radius')
                my_things=[temperature,rho,radius]

                if mass[0]>mass[len(mass)-1]:
                    zone_above=where(mass>mass_coo)[0][-1]
                    zone_below=zone_above+1
                else:
                    zone_above=where(mass>mass_coo)[0][0]
                    zone_below=zone_above-1
                   
                if mass[zone_below]>mass[zone_above]:
                    sys.exit("ERROR: finding of zone index confused")
                all_things_interplt=[]
                for thing in my_things:
                    thing_interplt=thing[zone_below]+(mass_coo-mass[zone_below])* \
                        (thing[zone_above]-thing[zone_below])/(mass[zone_above]-mass[zone_below])
                    all_things_interplt.append(thing_interplt)
                this_temperature,this_rho,this_radius=all_things_interplt

                string = str(step)+'  '+str(age)+'  '+str(this_temperature)+'  '+str(this_rho) 
                f.write(string+"\n")
                radius_at_mass_coo.append(this_radius)
                density_at_mass_coo.append(this_rho)
                temperature_at_mass_coo.append(this_temperature)
                age_all.append(age)
        f.close()
        return radius_at_mass_coo, density_at_mass_coo, temperature_at_mass_coo, age_all


    def abund_at_masscoorinate(self,ini,mass_coo):

        ''' create a file with distribution at a given mass coord, and at a given time step. 
        This for instance may be used as intial distribution for function trajectory, to
        reproduce local conditions in ppn.   

        ini       - model number
        mass_coo  - mass coordinate for the traj
  
        warning: remove the old abundance file, if you have any for the same mass coordinate.
        you are appending data, not overwriting.
        '''
        
        f = open('massf_'+str(mass_coo)+'.dat','a')

        age=self.se.get(ini,'age')
        mass=self.se.get(ini,'mass')  
        temperature=self.se.get(ini,'temperature')
        rho=self.se.get(ini,'rho')
        for i in range(len(mass)):
          if mass_coo == mass[i]:
             mass_coo_new = mass[i]
             zone = int(i) 
          elif mass_coo > mass[i]:
            	try:
                  	dum = mass[i+1]   
                        if mass_coo <= mass[i+1]:
                            mass_coo_new = mass[i+1]
                            zone = int(i+1) 
                except IndexError:
                        mass_coo_new = mass[i]
                        zone = int(i)
        string = 'C Model number | time (yr)  |  Temperature (GK)  | density (cm^-3) ' 
        f.write(string+"\n")
        string = 'C '+str(ini)+'  '+str(age)+'  '+str(temperature[zone])+'  '+str(rho[zone]) 
        f.write(string+"\n")
        string = 'C Mass coordinate that is really used' 
        f.write(string+"\n")
        string = 'C '+str(mass_coo_new)
        f.write(string+"\n")
        abund_string = self.se.dcols[5]
        print abund_string
        # the for loop below maybe optimized, I believe 
        # defining before isotopes and abundances. Marco 13 Jannuary 2011
        for i in range(len(self.se.isotopes)):
           string = 'D '+str(self.se.isotopes[i].split('-')[0]).upper()+'   '+str(self.se.isotopes[i].split('-')[1]).upper()+'    '+str(self.se.get(ini,abund_string)[zone][i])
           f.write(string+"\n")
        f.close()
         
         


    def kip(self,cycle_end,mix_thresh,xaxis,sparse):
        '''
        This function uses a threshold diffusion coefficient, above which the
        the shell is considered to be convective, to plot a Kippenhahn diagram.

            mix_thresh: the threshold diffusion coefficient
            xaxis:      age, cycle, log_age or log_time_left
            sparse:     sparsity factor when plotting from cyclelist

        ex: pt=mp.se('/ngpod1/swj/see/mppnp_out/scratch_data/M25.0Z1e-02','.h5')
            pt.kip(10000,'log_time_left',100)
        '''

        original_cyclelist = self.se.cycles
        cyclelist = original_cyclelist[0:cycle_end:sparse]

        xx = self.se.ages[:cycle_end:sparse]
        totalmass = []
        m_ini = float(self.se.get('mini'))


        fig = pl.figure(1)
        ax = pl.subplot(1,1,1)
        fsize = 12

        def getlims(d_coeff,massco):
            '''This function returns the convective boundaries for a cycle,
            given the cycle's dcoeff and massco columns, taking into account
            whether surface or centre are at the top'''
            plotlims = []
            if massco[0] > massco[-1]:
                for j in range(-1,-len(d_coeff)-1,-1):
                    if j == -1:
                        if d_coeff[j] >= mix_thresh:
                            plotlims.append(massco[j])
                        else:
                            pass
                    elif (d_coeff[j]-mix_thresh)*(d_coeff[j+1]-mix_thresh) < 0:
                        plotlims.append(massco[j])
                    if j == -len(d_coeff):
                        if d_coeff[j] >= mix_thresh:
                            plotlims.append(massco[j])
                return plotlims       
            else:
                for j in range(len(d_coeff)):
                    if j == 0:
                        if d_coeff[j] >= mix_thresh:
                            plotlims.append(massco[j])
                        else:
                            pass
                    elif (d_coeff[j]-mix_thresh)*(d_coeff[j-1]-mix_thresh) < 0:
                        plotlims.append(massco[j])
                    if j == len(d_coeff)-1:
                        if d_coeff[j] >= mix_thresh:
                            plotlims.append(massco[j])
                return plotlims

            

        if xaxis == 'age':
            ax.set_xlabel('Age [yr]',fontsize=fsize)
        elif xaxis == 'cycle':
            xx = cyclelist
            ax.set_xlabel('Cycle',fontsize=fsize)
        elif xaxis == 'log_age':
            for i in range(len(xx)):
                xx[i] = np.log10(xx[i])
            ax.set_xlabel('log$_{10}$(age) [yr]',fontsize=fsize)
        elif xaxis == 'log_time_left':
            for i in range(len(xx)):
                xx[i] = np.log10(max(xx)-xx[i])
            xx[-2] = xx[-3]-abs(xx[-4]-xx[-3])
            xx[-1] = xx[-2]-abs(xx[-3]-xx[-2])
            ax.set_xlabel('log$_{10}$(time until collapse) [yr]',fontsize=fsize)

        #centre-surface flag:
        flag = False

        if self.se.get(cyclelist[1],'mass')[0] > self.se.get(cyclelist[1],'mass')[-1]:
            flag = True

        for i in range(len(cyclelist)):
            if flag == True:
                totalmass.append(self.se.get(cyclelist[i],'mass')[0])
            else:
                totalmass.append(self.se.get(cyclelist[i],'mass')[-1])
            percent = int(i*100/len(cyclelist))
            sys.stdout.flush()
            sys.stdout.write("\rcreating color map " + "...%d%%" % percent)
            d_coeff = self.se.get(cyclelist[i],'dcoeff')
            massco = self.se.get(cyclelist[i],'mass')
            plotlims = getlims(d_coeff,massco)
            for k in range(0,len(plotlims),2):
                ax.axvline(xx[i],ymin=plotlims[k]/m_ini,ymax=plotlims[k+1]/m_ini,color='b',linewidth=0.5)
        

        ax.plot(xx, totalmass, color='black', linewidth=1)
        if xaxis == 'log_time_left':
            ax.axis([xx[0],xx[-1],0.,m_ini])
        else:
            ax.axis([min(xx),max(xx),0.,m_ini])
        ax.set_ylabel('Mass [$M_{\odot}$]',fontsize=fsize)

        pl.show()

    def kip_cont(self,modstart,modstop,dcoeff_thresh=1.e12,xres=1000,ylims=[0.,0.],xlims=[0.,0.],yres=2000,ixaxis='log_time_left',outfile='',landscape_plot=False,codev='KEP'):
        '''
        This function creates a Kippenhahn diagram as a contour plot of the
        se output using a diffusion coefficient threshold.
        ---------------------------------------------------
        ex: pt.kip_cont(0,-1,xres=1000,yres=1000,ixaxis='log_time_left')
        unfortunately, although the implicit do loops are faster, I cannot
        implement a % bar, but this example call should not take more than 1 minute.
        ---------------------------------------------------
        modstart, modstop first and last cycle numbers to plot
        dcoeff_thresh     diffusion coefficient threshold, above which assume
                          convection is present
        xres,yres         x and y resolution
        xlims, ylims      x and y plot limits
        ixaxis            'age', 'model_number' or 'log_time_left'
        outfile           name of output file including extension, which will
                          be saved in 300 dpi
        landscape_plot    for a landscape plot
        codev             'KEP', 'MES' or 'GNV'
        '''

        original_cyclelist = np.array([int(self.se.cycles[i]) for i in range(len(self.se.cycles))])
        nmodels=len(self.se.cycles[modstart:modstop])
        print 'nmodels           = ', nmodels
        sparse = int(max(1,nmodels/xres))
        print sparse
        cyclelist = original_cyclelist[modstart:modstop:sparse]
        cyclelist=cyclelist[0:-1]
        print cyclelist
        xx = np.array([self.se.ages[int(cyclelist[i])] for i in range(len(cyclelist))])
        age_unit=self.get('age_unit')
        if codev == 'KEP':
            oneyear = 365*24*3600
        elif codev == 'MES':
            oneyear=self.get('one_year')
        elif codev == 'GNV':
            oneyear=1.e0


        m_min=ylims[0]
        m_max=ylims[1]
        if m_min==0.:
            m_min=0.
            ylims[0] = m_min
        if m_max==0.:
            m_max=float(self.se.get('mini'))
            ylims[1] = m_max


        dy = (m_max-m_min)/float(yres)
        print 'dy = ', dy
        y = np.arange(m_min, m_max, dy)
        print y

        Z = np.zeros([len(y),len(xx)],float)

        print cyclelist

#        print 'getting mass_matrix...'
#        mass_matrix=self.se.get([cyclelist[i] for i in range(len(cyclelist))],'mass')
#        print 'getting dcoeff_matrix...'
#        dcoeff_matrix=self.se.get(cyclelist,'dcoeff')
#        print 'getting total_mass_vector...'
#        total_mass_vector=self.se.get(cyclelist,'total_mass')

        print 'calculating convection matrix... '
        conv_i_vec_matrix=np.array([np.interp(y,self.se.get(cyclelist[i],'mass')[::-1],self.se.get(cyclelist[i],'dcoeff')[::-1],0.,self.se.get(cyclelist[i],'mass')[::-1][-1]) for i in range(len(cyclelist))])
#        conv_i_vec_matrix=[]
#        for i in range(len(cyclelist)):
#            conv_i_vec_matrix.append(np.interp(y,self.se.get(cyclelist[i],'mass')[::-1],self.se.get#(cyclelist[i],'dcoeff')[::-1],0.,self.se.get(cyclelist[i],'mass')[::-1][-1]))
#            percent = int(i*100/len(cyclelist))
#            sys.stdout.flush()
#            sys.stdout.write("\rcreating color map " + "...%d%%" % percent)

        print 'getting stellar mass...'
        if codev == 'MES':
            mtot = np.array([self.se.get(cyclelist[i],'total_mass')/1.9891e33 for i in range(len(cyclelist))])
        else:
            mtot = np.array([self.se.get(cyclelist[i],'total_mass') for i in range(len(cyclelist))])
        print mtot
#        conv_i_vec_matrix=np.array([np.interp(y,mass_matrix[i][::-1],dcoeff_matrix[i][::-1],0.,mass_matrix[i][::-1][-1]) for i in range(len(cyclelist))])
#        print 'getting stellar mass...'
#        mtot = np.array([total_mass_vector[i]/1.9891e33 for i in range(len(cyclelist))])

        print np.shape(conv_i_vec_matrix)
        print np.shape(np.transpose(conv_i_vec_matrix))
        print np.shape(Z)
        print np.shape(xx)
        print np.shape(y)

        print conv_i_vec_matrix

       # print xx[0]
       # # Set up x-axis
       # if xax is 'logtimerev': 
       #     xx=np.log10(np.max(xx)+t_eps-xx)                
       #     ax.set_xlabel('$\log(t_{final} - t)$  $\mathrm{[yr]}$',fontsize=fsize)
       # elif xax is 'time':
       #     ax.set_xlabel('$t$ $\mathrm{[yr]}$',fontsize=fsize)
       #     if tzero is True:
       #         print 'zero time is '+str(xx[0])
       #         xx=xx-xx[0]
       #         ax.set_xlabel('$t - t_0$ $\mathrm{[yr]}$',fontsize=fsize)

        #########################################################################
        # PLOT:
        fsize=20

	params = {'axes.labelsize':  fsize,
	  'text.fontsize':   fsize,
	  'legend.fontsize': fsize,
	  'xtick.labelsize': fsize*0.8,
	  'ytick.labelsize': fsize*0.8,
	  'text.usetex': False}
	pl.rcParams.update(params)
        fig = pl.figure(1)
        ax = pl.axes()
        if landscape_plot == True:
		fig.set_size_inches(9,4)
		fsize=20
	        pl.gcf().subplots_adjust(bottom=0.15)
	        pl.gcf().subplots_adjust(right=0.85)

        ax.set_ylabel('$\mathrm{Mass }(M_\odot)$',fontsize=fsize)

        xxx=[]
        age_array = self.se.ages
        if codev == 'KEP':
            original_ages = np.array(self.se.ages[0:modstop])
            for i in range(1,len(original_ages)):
                if (original_ages[i]-original_ages[i-1]) < 0.:
                    age_at_restart_idx = i-1
                    age_at_restart = original_ages[i-1]
                    print 'age restart found at cycle = '+str(age_at_restart_idx)+', age = '+str(age_at_restart)
                    KEPLER = True
                    break

            for i in range(age_at_restart_idx+1,len(original_ages)):
                original_ages[i] = original_ages[i] + age_at_restart
            age_array = original_ages

	if ixaxis == 'log_time_left':
	# log of time left until core collapse
	    gage= np.array(age_array)/oneyear
	    lage=np.zeros(len(gage))
	    agemin = max(abs(gage[-1]-gage[-2])/5.,1.e-10)
	    for i in np.arange(len(gage)):
	        if gage[-1]-gage[i]>agemin:
	            lage[i]=np.log10(gage[-1]-gage[i]+agemin)
	        else :
	            lage[i]=np.log10(agemin)
	    xxx = lage[modstart:modstop:sparse]
            xxx = xxx[0:-1]
	    print 'plot versus time left'
	    ax.set_xlabel('$\mathrm{log}_{10}(t^*) \, \mathrm{(yr)}$',fontsize=fsize)
            if xlims == [0.,0.]:
                xlims = [xxx[0],xxx[-1]]
	elif ixaxis =='model_number':
	    xxx= cyclelist
	    print 'plot versus model number'
	    ax.set_xlabel('Model number',fontsize=fsize)
            if xlims == [0.,0.]:
                xlims = [cyclelist[0],cyclelist[-1]]
	elif ixaxis =='age':
	    xxx = age_array/oneyear/1.e6
	    print 'plot versus age'
	    ax.set_xlabel('Age [Myr]',fontsize=fsize)
            if xlims == [0.,0.]:
                xlims = [xxx[0],xxx[-1]]

        cmapMIX=mpl.colors.ListedColormap(['w','#8B8386']) # rose grey
        #cmap = mpl.colors.ListedColormap(['w','b'])

        print 'plotting contours'
        print len(xxx),len(y)
        ax.autoscale(False)
        low_level = dcoeff_thresh - dcoeff_thresh/10
        high_level= dcoeff_thresh + dcoeff_thresh/10
        ax.contourf(xxx,y,np.transpose(conv_i_vec_matrix), cmap=cmapMIX, alpha=0.6,levels=[low_level,1.e99])
        ax.contour(xxx,y,np.transpose(conv_i_vec_matrix), cmap=cmapMIX,levels=[low_level,high_level])
        ax.plot(xxx,mtot,color='k')
        #for i in range(len(xxx)):
        #  for j in range(len(y)):
        #    pl.plot(xxx[i],conv_i_vec_matrix[i,j],'bo',alpha=0.5)
        ax.axis([xlims[0],xlims[-1],ylims[0],ylims[1]])
        if outfile != '':
            fig.savefig(outfile,dpi=300)
        pl.show()

    def kip_cont2(self,sparse,cycle_start=0,cycle_end=0,plot=['dcoeff'],thresholds=[1.0E+12],xax='log_time_left',alphas=[1.0],yllim=0.,yulim=0.,y_res=2000,xllim=0.,xulim=0.,age='seconds',sparse_intrinsic=20, engen=False,netnuc_name='eps_nuc',engenalpha=0.6,outfile='plot.pdf',annotation='',KEPLER=False):
        '''
        !! EXPERIMENTAL FEATURE !!
        This function creates a Kippenhahn diagram as a contour plot of the
        .se.h5 or .out.h5 files using any continuous variable (columns in the
        hdf5 cycle data.
        Multiple columns may be plotted, their name indicated in the list "plot",
        and theirthresholds in the list "thresholds".
        
        Currently, this is only designed to take one threshold for each variable
        but future versions will hopefully be able to plot with multiple
        thresholds, however you may circumvent this issue by repeating the variable
        in "plots" and entering a second threshold for it in "thresholds".

        Arguments:
        --------------------
        sparse:              x-axis (timestep) sparsity;
                             true sparsity = sparse*sparse_intrinsic
                             Try 100 or 500 for .se.h5 and 20 for .out.h5 files
                             for preliminaxllimry plots
        sparse_intrinsic:    sparsity of timesteps in the data provided (usually
                             20 for .out.h5 
                             files and 1 for .se.h5 files
        cycle_start:         cycle from which you wish to plot (defaults to 0)
        cycle_end:           maximum cycle that you wish to plot (if 0, = last
                             cycle available).
        plot:                1-D array containing the variables to be plotted (as 
                             strings, e.g. plots=['dcoeff','C-13']. I recommend
                             always plotting 'dcoeff' as plots[0])
        thresholds:          1-D array containing the thresholds corresponding to
                             the variables in "plots". The threshold for 'dcoeff'
                             defaults to 1.0E+12
        xax:                 x-axis quantity; either 'log_time_left' or 'cycles'
        alphas:              array containing the opacity (0 to 1) of the contour
                             for each variable.
        yllim:               lower plot limit for y-axis (mass co-ordinate)
        yulim:               lower plot limit for y-axis (mass co-ordinate)
        y_res:               y-axis resolution. Defaults to 2000 but increasing
                             to as much as 10000
                             does not significantly affect the plotting time.
        age:                 either 'years' or 'seconds', depending on the data
        engen:               boolean indicating whether the user would like to
                             plot Kippenhahn of convective zones and energy
                             generation. If True, please still include
                             plots=['dcoeff'] and thresholds=[1.0E+12'] in your call.
                             This will require the data to have an 'eps_nuc' column,
                             so the plot is only working for .se.h5 files from MESA
                             in the current se library. This is the most recent
                             addition, so probably the most buggy. The plot script
                             will automatically calculate and assign multiple
                             thresholds according to the model.
        engenalpha:          opacity of the energy generation contours.
        netnuc_name:         the name of the column containing (eps_nuc-eps_neu).
                             If you do not have available eps_neu, then you can give
                             netnuc_name="eps_nuc" to just plot energy generation.
        '''

        # Organize cycles and ages:
        original_cyclelist = self.se.cycles
        if cycle_end==0.:
            cycle_end = original_cyclelist[-1]
        cycle_end = int(cycle_end)/sparse_intrinsic - 1
        if cycle_start==0:
            pass
        else:
            cycle_start = int(cycle_start)/sparse_intrinsic - 1
        cyclelist = original_cyclelist[cycle_start:cycle_end:sparse]
        # fix for KEPLER restart counting at O burning:
        original_ages= self.se.ages
        age_at_restart = 0.
        age_at_restart_idx = 0
        if KEPLER == True:
            for i in range(1,len(original_ages)):
                if (original_ages[i]-original_ages[i-1]) < 0.:
                    age_at_restart_idx = i-1
                    age_at_restart = original_ages[i-1]
                    print 'age restart found at cycle = '+str(age_at_restart_idx)+', age = '+str(age_at_restart)
                    KEPLER = True
                    break

            for i in range(age_at_restart_idx+1,len(original_ages)):
                original_ages[i] = original_ages[i] + age_at_restart


        # Figure:
        fig = pl.figure()
        ax = pl.axes()
	params = {'axes.labelsize':  15,
          'text.fontsize':   15,
          'legend.fontsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'text.usetex': False}
        fsize=18
	pl.rcParams.update(params)
        # X-axis:
        if xax == 'log_time_left':
            if KEPLER == True:
                xxtmp = original_ages[cycle_start:cycle_end:sparse]
            else:
                xxtmp = self.se.ages[cycle_start:cycle_end:sparse]
            if age == 'years':
                if KEPLER == True:
                    pass
                else:
                    xxtmp = self.se.ages[cycle_start:cycle_end:sparse]
            elif age == 'seconds':
                if KEPLER == True:
                    for i in range(len(cyclelist)):
                        xxtmp[i] = original_ages[cycle_start:cycle_end:sparse][i]/31558149.984
                else:
                    for i in range(len(cyclelist)):
                        xxtmp[i] = self.se.ages[cycle_start:cycle_end:sparse][i]/31558149.984
        if xax == 'cycles':
            xx = cyclelist
            xxtmp = cyclelist
        # Set up x-axis according to whether ages are in years or seconds and
        # re-write as log(time left). The last entry will always be -inf in this 
        # way so we calculate it by extrapolating the anti-penultimate and
        # penultimate entries.
        # Modified from the GENEC gdic (R. Hirschi)
        if xax == 'log_time_left':
            xx=np.zeros(len(xxtmp))
            agemin = max(abs(xxtmp[-1]-xxtmp[-2])/5.,1.e-10)
            for i in np.arange(len(xxtmp)):
                 if xxtmp[-1]-xxtmp[i]>agemin:
                     xx[i]=np.log10(xxtmp[-1]-xxtmp[i]+agemin)
                 else :
                     xx[i]=np.log10(agemin)
            ax.set_xlabel('$log_\mathrm{10}(t_\mathrm{end}-t)\;[\mathrm{yr}]$',fontsize=fsize-1)
        if xax == 'cycles':
            ax.set_xlabel('$\mathrm{CYCLE}$',fontsize=fsize-1)
        # Y-axis limits and resolution:
        totalmass = []
        try:
            m_ini = float(self.se.get('mini'))
        except:
            mini=m.se.get(0,'total_mass')
            mini=mini/constants.mass_sun
            print 'getting Mini from 1st cycle'

        if yulim==0.:
            yulim = m_ini 
        dy = m_ini/float(y_res)
        vlinetol = 1.0E-8
        # Set up (y-axis) vector and a 3-D (hist) array to store all of the
        # contours.
        y = np.arange(0., m_ini, dy)
        if engen == True:
            Z = np.zeros([len(y),len(xxtmp),len(plot)+2],float)
        else:
            Z = np.zeros([len(y),len(xxtmp),len(plot)],float)

        # Define function extracting the contour boundaries which will be
        # called for every cycle in cyclelist, for every variable to be plotted
        # along with its corresponding threshold(s).
        def getlims(variable_array,thresh,massco_array):
            '''This function returns the variable boundaries (in mass) for a cycle,
            given the cycle's variable and mass columns, ensuring that the boundaries
            are ordered centre to surface (as some .se.h5 files are the opposite).'''
            plotlims = []
            # Just a fix for some get problem I was having:
            if len(massco_array) == 2:
                massco_array = massco_array[0]
            if len(variable_array) == 2:
                variable_array = variable_array[0]
            if massco_array[0] > massco_array[-1]:
                for j in range(-1,-len(variable_array)-1,-1):
                    if j == -1:
                        if variable_array[j] >= thresh:
                            plotlims.append(massco_array[j])
                        else:
                            pass
                    elif (variable_array[j]-thresh)*(variable_array[j+1]-thresh) < 0:
                        plotlims.append(massco_array[j])
                    if j == -len(variable_array):
                        if variable_array[j] >= thresh:
                            plotlims.append(massco_array[j])
                return plotlims       
            else:
                for j in range(len(variable_array)):
                    if j == 0:
                        if variable_array[j] >= thresh:
                            plotlims.append(massco_array[j])
                        else:
                            pass
                    elif (variable_array[j]-thresh)*(variable_array[j-1]-thresh) < 0:
                        plotlims.append(massco_array[j])
                    if j == len(variable_array)-1:
                        if variable_array[j] >= thresh:
                            plotlims.append(massco_array[j])
                return plotlims
        # Flag preventing plotting any other variables on an energy generation
        # Kippenhahn plot:
        if engen == True:
            plot = ['dcoeff']
        # This loop gets the mass co-ordinate array and the variable arrays,
        # calls to get the boundaries in order, and populates the contour array.
        ypscoeff = [-1,-1,-1] # this should have same length as plot - quick fix for yps.
        total_massco = []
        for i in range(len(cyclelist)):
#            print 'CYCLE: ', cyclelist[i]
            massco = self.se.get(cyclelist[i],'mass')
            total_massco.append(max(massco))
            plotlimits=[]
            for j in range(len(plot)):
                if plot[j][1] == '-' or plot[j][2] == '-':
                    # Assume file has yps, not iso_massf
                    ypsthere = True
                    try: variables = self.se.get(cyclelist[i],'yps')
                    # If this is not the case, do the usual se calls for iso_massf
                    except KeyError:
                        variables = self.se.get(cyclelist[i],plot[j])
                        ypsthere = False
                    # If yps is there, ask which indices correspond to the 
                    # elements that are to be plotted, one by one.
                    if ypsthere == True:
                        if ypscoeff[j] == -1:
                            ypscoeff[j] = int(raw_input("What integer is your element "+str(plot[j])+" in the 'yps' array? "))
                        else:
                            pass
                        variables = self.se.get(cyclelist[i],'yps')[:,ypscoeff[j]]
                else:
                    variables = self.se.get(cyclelist[i],plot[j])
                plotlims = getlims(variables,thresholds[j],massco)
                plotlimits.append(plotlims)
            percent = int(i*100/len(cyclelist))
            sys.stdout.flush()
            sys.stdout.write("\rcreating color map " + "...%d%%" % percent)
            for g in range(len(plot)):
                for k in range(0,len(plotlimits[g]),2):
                    llimit = plotlimits[g][k]
                    ulimit = plotlimits[g][k+1]
                    #if xx[i] >= 0:
                    for f in range(y_res):
                        if llimit<=y[f] and ulimit>y[f]:
                            Z[f,i,g]=1.
                    #else:
                     #   ax.axvline(xx[i],ymin=llimit/m_ini,ymax=ulimit/m_ini,color='#8B8386',alpha=alphas[0],linewidth=0.5)
        # This function determines the adjacent two mass cells to a point in the
        # y-vector (which contains mass co-ordinates centre to surface, split into
        # y_res chunks), returning their index in the mass co-ordinate vector for
        # that timestep/cycle.
        def find_nearest(array,value):
            '''Returns [lower,upper] indexes locating adjacent mass cells (in the
               massco vector) around y-value (one of y_res points equally spaced
               between centre and surface).'''
            idx=(np.abs(array-value)).argmin()
            lims=np.zeros([2],int)
            if idx == len(array)-1: # SJONES post-mod
                lims[0] = idx - 1   # SJONES post-mod
                lims[1] = idx       # SJONES post-mod
                return lims
            if array[idx] < value:
                if array[idx]-array[idx+1] < 0.:
                    lims[0] = idx
                    lims[1] = idx-1
                    return lims
                else:
                    lims[0] = idx
                    lims[1] = idx+1
                    return lims
            elif array[idx] > value:
                if array[idx]-array[idx+1] < 0.:
                    lims[0] = idx+1
                    lims[1] = idx
                    return lims
                else:
                    lims[0] = idx-1
                    lims[1] = idx
                    return lims
        # This flag enebles the loop below it to populate the contour array for
        # energy generation. It does not take threshold arguments, as the array
        # contains the log of the energy generation rather than "above" or "below".
        # Because of this, contour boundaries are automatically calculated
        # according to the max energy generation in the model.
        dummy_engen=[]
        engen_signs = []
        if engen == True:
        # Requires eps_nuc array in the data. Produces energy generation contour
        # by linearly interpolating eps_nuc between mass co-ordinates according
        # to the y-resolution:
            for i in range(len(cyclelist)):
#                print 'CYCLE: ', cyclelist[i]
                max_energy_gen = 0.
                min_energy_gen = 0.
                massco = self.se.get(cyclelist[i],'mass')
                if len(massco) <= 10:
                    massco=massco[0]
                dummy_engen = self.se.get(cyclelist[i],netnuc_name)
                if len(dummy_engen) <= 10:
                    dummy_engen = dummy_engen[0]
                for f in range(len(dummy_engen)):
                    # make all values absolute, but note in engen_signs which were negative:
                    if dummy_engen[f] == 0.:
                        engen_signs.append(1.)
                    else:
                        engen_signs.append(dummy_engen[f]/abs(dummy_engen[f]))
                    if abs(engen_signs[f]) != 1.:
                        print 'engen sign not +/- 1!!'
                        print 'engen_signs['+str(f)+'] = ',engen_signs[f]
                        print 'dummy_engen[f] = ', dummy_engen[f]
                        sys.exit()
                    dummy_engen[f] = abs(dummy_engen[f])
                log_epsnuc = np.log10(dummy_engen)
                # now insert the correct signs again:
                for f in range(len(log_epsnuc)):
                    log_epsnuc[f] = log_epsnuc[f]*engen_signs[f]
                #for f in range(len(log_epsnuc)):
                    #if str(log_epsnuc[f]) == 'nan': log_epsnuc[f] = 0.
#                print log_epsnuc
                percent = int(i*100/len(cyclelist))
                sys.stdout.flush()
                sys.stdout.write("\rcreating color map " + "...%d%%" % percent)
                for j in range(len(y)):
                    if j == len(y)-1:
                        energy_here = 0.
                    elif j == 0:
                        energy_here = log_epsnuc[-1]
                    elif y[j] > max(massco):
                        energy_here = 0.
                    else:
                        lims = find_nearest(massco,y[j])
                        frac = (y[j]-massco[lims[0]])/(massco[lims[1]]-massco[lims[0]])
                        energy_here = frac*(log_epsnuc[lims[1]]-log_epsnuc[lims[0]]) + log_epsnuc[lims[0]]
                    if energy_here > max_energy_gen:
                        max_energy_gen = energy_here
                    if energy_here < min_energy_gen:
                        min_energy_gen = energy_here
                    if abs(max_energy_gen) > 100.:
                        print y[j]
                        print engen_signs[f], log_epsnuc[f], frac, lims[0], lims[1], massco[lims[0]], massco[lims[1]]
                        print (massco[lims[1]]-massco[lims[0]]), (y[j]-massco[lims[0]])
                        print max_energy_gen
                        print 'exit due to energy generation > 100'
                        sys.exit()
#                    print energy_here
#                    print max_energy_gen
#                    if energy_here >= 0.:
                    #Z[j,i,1] = 10**energy_here #SJONES comment
                    if energy_here < 0.:
                        energy_here = 0.
                    Z[j,i,1] = energy_here
#                    if energy_here < 0.:
#                        Z[j,i,2] = 10**energy_here


        # Here we define the colourmap for the energy generation and an array
        # containing a list of colours in which to plot each variable (in the
        # order that the variables appear in "plots") iso_colours is obsolete
        # but was for when we tried plotting isotopes with just their boundary
        # lines as opposed to shading (for clarity). Colourmaps of these choices
        # are written to cmap (array).
        engen_cmap=mpl.cm.get_cmap('Blues')
        engen_cmap.set_under(color='w',alpha=engenalpha)

        enloss_cmap=mpl.cm.get_cmap('Reds')
        colours = ['#8B8386','m','g','b']
        iso_colours = ['b','r','y']
        cmap = []
        for i in range(len(plot)):
            cmap.append(mpl.colors.ListedColormap(['w',colours[i]]))

        print 'plotting contours'

        if xllim==0. and xulim==0.:
            ax.axis([float(xx[0]),float(xx[-1]),yllim,yulim])
        else:
            ax.axis([xllim,xulim,yllim,yulim])

        # Plot all of the contours. Levels indicates to only plot the shaded
        # regions and not plot the white regions, so that they are essentially
        # transparent. If engen=True, then the energy generation levels
        # (boundary values) are calculated (in dex) from 2 to the maximum in
        # steps of 2.
        for i in range(len(plot)):
            ax.contourf(xx,y,Z[:,:,i],levels=[0.5,1.5],colors=colours[i], alpha=alphas[i])
        if engen == True:
            #ceiling = int(max_energy_gen+1)
            #floor = int(min_energy_gen+1)
            #cburn = ax.contourf(xx,y,Z[:,:,1],cmap=engen_cmap,locator=mpl.ticker.LogLocator(),alpha=engenalpha) # SJONES comment
            cburn = ax.contourf(xx,y,Z[:,:,1],cmap=engen_cmap,alpha=engenalpha,levels=range(5,32,5))
            cbarburn = pl.colorbar(cburn)
#            if min_energy_gen != 0:
#                closs = ax.contourf(xx,y,Z[:,:,2],cmap=enloss_cmap,locator=mpl.ticker.LogLocator(),alpha=engenalpha)
#                cbarloss = pl.colorbar(closs)

#            cbarburn.set_label('$\epsilon_\mathrm{nuc}-\epsilon_{\\nu} \; (\mathrm{erg\,g}^{-1}\mathrm{\,s}^{-1})$, > 0',fontsize=fsize)
            cbarburn.set_label('$log_{10}(\epsilon_\mathrm{nuc}) \; (\mathrm{erg\,g}^{-1}\mathrm{\,s}^{-1})$, > 0',fontsize=fsize)


        pl.plot(xx,total_massco,color='k')
        pl.text(0.9,0.9,annotation,horizontalalignment='right',transform = ax.transAxes,fontsize=fsize)
        pl.ylabel('$\mathrm{Mass}\;[M_\odot]$',fontsize=fsize-1)
        pl.savefig(outfile)
        print outfile+' is done.'
        pl.show()


    def abup_se_plot(mod,species):

        ''' plot species from one ABUPP file and the se file
        
        You must use this function in the directory where the ABP
        files are and an ABUP file for model mod must exist.

        mod - model to plot, you need to have an ABUPP file for that model'''

# Marco, you have already implemented finding headers and columns in
# ABUP files. You may want to transplant that into here?
        species='C-12'

        filename = 'ABUPP%07d0000.DAT' % mod
        print filename
        mass,c12=np.loadtxt(filename,skiprows=4,usecols=[1,18],unpack=True)
        c12_se=self.se.get(mod,'iso_massf','C-12')
        mass_se=self.se.get(mod,'mass')

        pyl.plot(mass,c12)
        pyl.plot(mass_se,c12_se,'o',label='cycle '+str(mod))
        pyl.legend()
        

    def read_iso_abund_marco(self, mass_range, cycle):
        ''' plot the abundance of all the chemical species
        inputs:
        mass_range - a 1x2 array containing the lower and upper mass range.  
                     if None, it will plot over the entire range
        cycle       - a string/integer of the cycle of interest.
        stable     - a boolean of whether to filter out the unstables.
            
        '''

	import utils as u

        masses = []
        #    Check the inputs
        #if not self.se.cycles.count(str(cycle)):
        #    print 'You entered an cycle that doesn\'t exist in this dataset:', cycle
        #    print 'I will try and correct your format.'
        #    cyc_len = len(self.se.cycles[-1])
        
        #    print cyc_len, len(str(cycle))
        #
        #    while len(str(cycle)) < cyc_len:
        #        cycle = '0'+str(cycle)
        #        print cycle
        
        #    if not self.se.cycles.count(str(cycle)):
        #        print 'I was unable to correct your cycle.  Please check that it exists in your dataset.'

        masses = self.se.get(cycle,'mass')
        if mass_range == None:
            print 'Using default mass range'
            mass_range = [min(masses),max(masses)]    
	# what this was for??? Marco        
	#masses.sort()
        #mass_range.sort()

        
        
        print 'Using The following conditions:'
        print '\tmass_range:', mass_range[0], mass_range[1]
        print '\tcycle:', cycle

        isotope_names = self.se.isotopes
	u.convert_specie_naming_from_h5_to_ppn(isotope_names)
	names_ppn_world = u.spe
	number_names_ppn_world = u.n_array	
	u.define_zip_index_for_species(names_ppn_world,number_names_ppn_world)
                        
        # from here below I read the abundance.
        
        #name_specie_in_file=self.se.dcols[5]
	# I am using directly 'iso_massf' only because somehow m20 explosive do not have dcols....
	name_specie_in_file='iso_massf'
	abunds=self.se.get(cycle,name_specie_in_file)

        global used_masses
        used_masses = []
        self.mass_frac = []
        for i in range(len(masses)):
            if mass_range[0] <=  masses[i]  and mass_range[1] >=  masses[i] :
                used_masses.append(masses[i])
            	self.mass_frac.append(abunds[i])


    def decay(self,mass_frac):

        '''this module simply calculate abundances of isotopes after decay
        it requires before being used a call to read_iso_abund_marco and stable_species.
        see e.g., function interface'''

	import utils as u

        global decayed_multi_d
        decayed_multi_d=[]
	#print len(mass_frac)
	#print len(decay_raw)
        for iii in range(len(mass_frac)):
            jj=-1
            decayed=[]
            for i in range(len(u.decay_raw)):
                 if u.jdum[i] > 0.5:
                    jj=jj+1
                    dummy=0.
                    for j in range(len(u.decay_raw[i])):             
                        try:    
                            dum_str = u.decay_raw[i][j]    
                            dummy = dummy + float(self.mass_frac[iii][u.cl[dum_str.lower().capitalize()]])
			    #print cl[dum_str.lower().capitalize()]		
                            #print dum_str, mass_frac[iii][cl[dum_str.capitalize()]]
                        except KeyError:
                            None            
                            #print 'I am not in the network:',decay_raw[i][j]
                        except IndexError:
                            None
                            #print 'I am not read',cl[decay_raw[i][j].lower().capitalize()],decay_raw[i][j]    
                    decayed.append(dummy) 
            decayed_multi_d.append(decayed)    
	#print 'I am here'
        #print decayed_multi_d[0][back_ind['CU 63']]
        #print mass_frac[0][cl[('CU 63').lower().capitalize()]],spe[cl[('CU 63').lower().capitalize()]]
	#print back_ind        
        #print decayed_multi_d[0][back_ind['ZR 90']]
        #print mass_frac[0][cl[('zr 90').capitalize()]],spe[cl[('zr 90').capitalize()]]
        #print decayed_multi_d[0][back_ind['ZR 91']]
        #print mass_frac[0][cl[('zr 91').capitalize()]],spe[cl[('zr 91').capitalize()]]
        #print decayed_multi_d[0][back_ind['ZR 92']]
        #print mass_frac[0][cl[('zr 92').capitalize()]],spe[cl[('zr 92').capitalize()]]
        #print decayed_multi_d[0][back_ind['ZR 94']]
        #print mass_frac[0][cl[('zr 94').capitalize()]],spe[cl[('zr 94').capitalize()]]
        #print decayed_multi_d[0][back_ind['ZR 96']]
        #print mass_frac[0][cl[('zr 96').capitalize()]],spe[cl[('zr 96').capitalize()]]
        #print decayed_multi_d[0][back_ind['CS133']]
        #print mass_frac[0][cl[('cs133').capitalize()]],spe[cl[('cs133').capitalize()]]
        #print spe,len(u.spe)
        #print cl,len(cl)    
    
    
    def burnstage(self, **keyw):
        """
        This function calculates the presence of burning stages and outputs
        the ages when key isotopes are depleted and uses them to calculate burning
        lifetimes.
    
        burnstage()
        
        A list is returned containing the following information:
        
        [burn_cycles, burn_ages, burn_abun, burn_type, burn_lifetime]
    
        Cycles contain the cycle numbers for the various points where the abundance
        is abun.  The age of the star at each point and the type of burning is
        indicated by those arrays.  The lifetimes are calculated by 
    
        The following keywords can also be used:
    
        | Keyword argument | Default Value:
        ------------------------------------------------------------------------------
         abund               "iso_massf"
         isoa                "A"
         isoz                "Z"
         mass                "mass"
         cycle               "cycle"
         cyclefin             0
        
        All arguments change the name of fields used to read data from HDF5 files,
        other than cyclefin.  Cyclefin is the last timestep to use when reading 
        files
        """

        if ("isoa" in keyw) == False:
            keyw["isoa"] = "A"
        if ("isoz" in keyw) == False:
            keyw["isoz"] = "Z"
        if ("mass" in keyw) == False:
            keyw["mass"] = "mass"
        if ("age" in keyw) == False:
            keyw["age"] = "age"
        if ("abund" in keyw) == False:
            keyw["abund"] = "iso_massf"
        if ("cycle" in keyw) == False:
            keyw["cycle"] = "cycle"
        if ("cyclefin" in keyw) == False:
            cyclefin = 1.e99
        else:
            cyclefin = keyw["cyclefin"]

        burn_cycles = []
        burn_ages = []
        burn_abun = []
        burn_type = []
    
        burn_lifetime = []
    
        firstfile = True
    
        hemax, cmax, nemax, omax = 0, 0, 0, 0
        
            
        cycles_list = self.se.cycles
        xm_list     = self.se.get(keyw["mass"])
        age_list    = self.se.get(keyw["age"])
        abund_list  = self.se.get(keyw["abund"])
        tables = self.se.Tables()
        isoa_list   = tables[0]
        isoz_list   = tables[1]
    
        # Check the order of the input data
        xm_init = xm_list[0]
        centre = 0
        if isinstance(xm_init, list) == True:
            if xm_init[0] > xm_init[1]:
                # mass is descending with shell number and the centre of the star
                # is the last shell
                centre = -1
        
        # Determine yps indices for certain abundances
        for i in range(len(isoa_list)):
            try:
                A = int(isoa_list[i])
                Z = int(isoz_list[i])
            except TypeError:
                A = int(isoa_list[i][0])
                Z = int(isoz_list[i][0])                
            
            if A == 1 and Z == 1:
                h1index = i
            elif A == 4 and Z == 2:
                he4index = i
            elif A == 12 and Z == 6:
                c12index = i
            elif A == 16 and Z == 8:
                o16index = i
            elif A == 20 and Z == 10:
                ne20index = i
            elif A == 28 and Z == 14:
                si28index = i
 
        if firstfile == True:
            hmax = abund_list[0][centre][h1index]
            firstfile = False
            
        # Try and determine the location of a convective core using the central and
        # next from central shells
        for i in range(1, len(cycles_list)-1):
            if cycles_list[i] > cyclefin and cyclefin != 0:
                pair = False
                age1 = -1
                for i in range(len(burn_type)):
                    if 'start' in burn_type[i] and pair == False:
                        age1 = burn_ages[i]
                        pair = True
                    elif 'end' in burn_type[i] and pair == True:
                        age2 = burn_ages[i]
                        pair = False
                        if age1 != -1:
                            burn_lifetime.append(age2 - age1)
                            age1 = -1
                            
                return [burn_cycles, burn_ages, burn_abun, burn_type,
                burn_lifetime]
            
            # H-burning
            hcen  = abund_list[i][centre][h1index]
            hcennext = abund_list[i+1][centre][h1index]
            if hcen >1.e-10:
                if hcennext < hmax-0.003 and hcen >= hmax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(hcen)
                    burn_type.append('H_start')
        
                if hcennext < 1.e-1 and hcen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('H')
                    
                if hcennext < 1.e-2 and hcen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('H')

                if hcennext < 1.e-3 and hcen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('H')
                    
                if hcennext < 1.e-4 and hcen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('H')

                if hcennext < 1.e-5 and hcen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('H_end')
                    hemax = abund_list[i][centre][he4index]

                if hcennext < 1.e-6 and hcen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('H')
                    
                if hcennext < 1.e-9 and hcen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('H')
                    
            # He-burning
            hecen  = abund_list[i][centre][he4index]
            hecennext = abund_list[i+1][centre][he4index]
            if hcen < 1.e-5 and hecen > 1.e-10:
                if hecennext < hemax-0.003 and hecen >= hemax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(hecen)
                    burn_type.append('He_start')
        
                if hecennext < 1.e-1 and hecen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('He')
                    
                if hecennext < 1.e-2 and hecen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('He')

                if hecennext < 1.e-3 and hecen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('He')
                    
                if hecennext < 1.e-4 and hecen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('He')

                if hecennext < 1.e-5 and hecen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('He_end')
                    cmax = abund_list[i][centre][c12index]

                if hecennext < 1.e-6 and hecen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('He')
                    
                if hecennext < 1.e-9 and hecen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('He')
    
            # C-burning
            ccen  = abund_list[i][centre][c12index]
            ccennext = abund_list[i+1][centre][c12index]
            if hcen < 1.e-5 and hecen < 1.e-5 and ccen > 1.e-10:
                if ccennext < cmax-0.003 and ccen >= cmax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(ccen)
                    burn_type.append('C_start')
        
                if ccennext < 1.e-1 and ccen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('C')
                    
                if ccennext < 1.e-2 and ccen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('C')

                if ccennext < 1.e-3 and ccen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('C')
                    
                if ccennext < 1.e-4 and ccen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('C')

                if ccennext < 1.e-5 and ccen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('C_end')
                    nemax = abund_list[i][centre][ne20index]

                if ccennext < 1.e-6 and ccen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('C')
                    
                if ccennext < 1.e-9 and ccen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('C')

            # Ne-burning
            necen  = abund_list[i][centre][ne20index]
            necennext = abund_list[i+1][centre][ne20index]
            if hcen < 1.e-5 and hecen < 1.e-5 and ccen < 1.e-3 and necen > 1.e-10:
                if necennext < nemax-0.003 and necen >= nemax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(necen)
                    burn_type.append('Ne_start')
        
                if necennext < 1.e-1 and necen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('Ne')
                    
                if necennext < 1.e-2 and necen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('Ne')
                    

                if necennext < 1.e-3 and necen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('Ne_end')
                    omax = abund_list[i][centre][o16index]
                    
                if necennext < 1.e-4 and necen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('Ne')

                if necennext < 1.e-5 and necen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('Ne')

                if necennext < 1.e-6 and necen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('Ne')

                if necennext < 1.e-9 and necen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('Ne')
                    
            # O-burning
            ocen  = abund_list[i][centre][o16index]
            ocennext = abund_list[i+1][centre][o16index]
            if hcen < 1.e-5 and hecen < 1.e-5 and ccen < 1.e-3 and ocen > 1.e-10:
                if ocennext < omax-0.003 and ocen >= omax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(ccen)
                    burn_type.append('O_start')
        
                if ocennext < 1.e-1 and ocen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('O')
                    
                if ocennext < 1.e-2 and ocen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('O')

                if ocennext < 1.e-3 and ocen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('O')
                    
                if ocennext < 1.e-4 and ocen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('O')

                if ocennext < 1.e-5 and ocen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('O_end')

                if ocennext < 1.e-6 and ocen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('O')
                    
                if ocennext < 1.e-9 and ocen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('O')
                        
        
        pair = False
        age1 = -1
        for i in range(len(burn_type)):
            if 'start' in burn_type[i] and pair == False:
                age1 = burn_ages[i]
                pair = True
            elif 'end' in burn_type[i] and pair == True:
                age2 = burn_ages[i]
                pair = False
                if age1 != -1:
                    burn_lifetime.append(age2 - age1)
                    age1 = -1
                    
        return [burn_cycles, burn_ages, burn_abun, burn_type, burn_lifetime]
    
    def burnstage_upgrade(self, **keyw):
        """
        This function calculates the presence of burning stages and outputs
        the ages when key isotopes are depleted and uses them to calculate burning
        lifetimes.
    
        burnstage()
        
        A list is returned containing the following information:
        
        [burn_cycles, burn_ages, burn_abun, burn_type, burn_lifetime]
    
        Cycles contain the cycle numbers for the various points where the abundance
        is abun.  The age of the star at each point and the type of burning is
        indicated by those arrays.  The lifetimes are calculated by 
    
        The following keywords can also be used:
    
        | Keyword argument | Default Value:
        ------------------------------------------------------------------------------
         abund               "iso_massf"
         isoa                "A"
         isoz                "Z"
         mass                "mass"
         cycle               "cycle"
         cyclefin             0
        
        All arguments change the name of fields used to read data from HDF5 files,
        other than cyclefin.  Cyclefin is the last timestep to use when reading 
        files
        """

        if ("isoa" in keyw) == False:
            keyw["isoa"] = "A"
        if ("isoz" in keyw) == False:
            keyw["isoz"] = "Z"
        if ("mass" in keyw) == False:
            keyw["mass"] = "mass"
        if ("age" in keyw) == False:
            keyw["age"] = "age"
        if ("abund" in keyw) == False:
            keyw["abund"] = "iso_massf"
        if ("cycle" in keyw) == False:
            keyw["cycle"] = "cycle"
        if ("cyclefin" in keyw) == False:
            cyclefin = 1.e99
        else:
            cyclefin = keyw["cyclefin"]

        burn_cycles = []
        burn_ages = []
        burn_abun = []
        burn_type = []
    
        burn_lifetime = []
    
        firstfile = True
    
        hemax, cmax, nemax, omax = 0., 0., 0., 0.
	hburn_logic = True
	heburn_logic = True
	cburn_logic = True
	neburn_logic = True
	oburn_logic = True		

	hburn_start_logic = False
	heburn_start_logic = False
	cburn_start_logic = False
	neburn_start_logic = False
	oburn_start_logic = False		        
            
        #cycles_list = self.se.cycles
	cyc = self.se.cycles
	sparsity_factor = int(1)	
	cycles_list = range(int(cyc[0]),int(cyc[len(cyc)-1]),((int(cyc[1])-int(cyc[0])))*sparsity_factor)
        age_list    = self.se.get(keyw["age"])
	# I want to read only the isotopes I need to identify the burning stages.
	all_isos=self.se.isotopes
	#list_all_isos=all_isos.tolist()
	useful_species = species_list("burn_stages")
	useful_index  = []
	for iso in useful_species:
		#useful_index.append(all_isos.index(iso))
		useful_index.append(useful_species.index(iso))

        specie_index={}
        for a,b in zip(useful_species,useful_index):
            specie_index[a] = b  

        # Check the order of the input data
        xm_init = self.se.get(0,keyw["mass"])
        central_zone = 0
	external_zone = -1
        if isinstance(xm_init, list) == True:
            if xm_init[0] > xm_init[1]:
                # mass is descending with shell number and the centre of the star
                # is the last shell
                central_zone = -1
		external_zone = 0

	# central zone first	
	zone = 0
	xm_cyc=[]
	xm_list=[]
	for i in cycles_list:
		xm_cyc  = self.se.get(i,keyw["mass"])[central_zone]
		xm_list.append(xm_cyc)

	
	abund_list = []
	for i in cycles_list:
		abund_tmp = []	
		for iso in useful_species:
			abund_cyc = self.se.get(i,keyw["abund"],iso)[central_zone]
			abund_tmp.append(abund_cyc)
		abund_list.append(abund_tmp)	

        if firstfile == True:
            hsurf = self.se.get(0,keyw["abund"],'H-1')[external_zone]
            #hesurf = self.se.get(0,keyw["abund"],'He-4')[external_zone]
            firstfile = False
            
        # Try and determine the location of a convective core using the central and
        # next from central shells
        for i in range(1, len(cycles_list)-1):
            if cycles_list[i] > cyclefin and cyclefin != 0:
                pair = False
                age1 = -1
                for i in range(len(burn_type)):
                    if 'start' in burn_type[i] and pair == False:
                        age1 = burn_ages[i]
                        pair = True
                    elif 'end' in burn_type[i] and pair == True:
                        age2 = burn_ages[i]
                        pair = False
                        if age1 != -1:
                            burn_lifetime.append(age2 - age1)
                            age1 = -1
                            
                return [burn_cycles, burn_ages, burn_abun, burn_type,
                burn_lifetime]
            	
		print 'passa 3'

            # H-burning
	    if hburn_logic:
 	           hcen  = abund_list[i][specie_index['H-1']]
                   hcennext = abund_list[i+1][specie_index['H-1']]
            	   if hcen >1.e-10:
                	if hcennext < hsurf-0.003 and hcen >= hsurf-0.003:
                    		burn_cycles.append(cycles_list[i])
                    		burn_ages.append(age_list[i])
                    		burn_abun.append(hcen)
                    		burn_type.append('H_start')
				hburn_start_logic = True
        
                	if hcennext < 1.e-1 and hcen >= 1.e-1:
                    		burn_cycles.append(cycles_list[i])
                    		burn_ages.append(age_list[i])
                    		burn_abun.append(1.e-1)
                    		burn_type.append('H')
                    
                	if hcennext < 1.e-2 and hcen >= 1.e-2:
                    		burn_cycles.append(cycles_list[i])
                    		burn_ages.append(age_list[i])
                    		burn_abun.append(1.e-2)
                    		burn_type.append('H')

                	if hcennext < 1.e-3 and hcen >= 1.e-3:
                    		burn_cycles.append(cycles_list[i])
                    		burn_ages.append(age_list[i])
                    		burn_abun.append(1.e-3)
                   		burn_type.append('H')
                    
                	if hcennext < 1.e-4 and hcen >= 1.e-4:
                    		burn_cycles.append(cycles_list[i])
                    		burn_ages.append(age_list[i])
                    		burn_abun.append(1.e-4)
                    		burn_type.append('H')

                	if hcennext < 1.e-5 and hcen >= 1.e-5:
                    		burn_cycles.append(cycles_list[i])
                    		burn_ages.append(age_list[i])
                    		burn_abun.append(1.e-5)
                    		burn_type.append('H_end')
                    		hemax = abund_list[i][specie_index['He-4']]
				if hburn_start_logic:
					hburn_logic == False

                	if hcennext < 1.e-6 and hcen >= 1.e-6:
                    		burn_cycles.append(cycles_list[i])
                    		burn_ages.append(age_list[i])
                    		burn_abun.append(1.e-6)
                    		burn_type.append('H')
                    
                	#if hcennext < 1.e-9 and hcen >= 1.e-9:
                    	#	burn_cycles.append(cycles_list[i])
                    	#	burn_ages.append(age_list[i])
                    	#	burn_abun.append(1.e-9)
                    	#	burn_type.append('H')
                    
            # He-burning
            hecen  = abund_list[i][specie_index['He-4']]
            hecennext = abund_list[i+1][specie_index['He-4']]
            if hcen < 1.e-5 and hecen > 1.e-10:
                if hecennext < hemax-0.003 and hecen >= hemax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(hecen)
                    burn_type.append('He_start')
        
                if hecennext < 1.e-1 and hecen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('He')
                    
                if hecennext < 1.e-2 and hecen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('He')

                if hecennext < 1.e-3 and hecen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('He')
                    
                if hecennext < 1.e-4 and hecen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('He')

                if hecennext < 1.e-5 and hecen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('He_end')
                    cmax = abund_list[i][specie_index['C-12']]

                if hecennext < 1.e-6 and hecen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('He')
                    
                if hecennext < 1.e-9 and hecen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('He')
    
            # C-burning
            ccen  = abund_list[i][specie_index['C-12']]
            ccennext = abund_list[i+1][specie_index['C-12']]
            if hcen < 1.e-5 and hecen < 1.e-5 and ccen > 1.e-10:
                if ccennext < cmax-0.003 and ccen >= cmax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(ccen)
                    burn_type.append('C_start')
        
                if ccennext < 1.e-1 and ccen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('C')
                    
                if ccennext < 1.e-2 and ccen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('C')

                if ccennext < 1.e-3 and ccen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('C')
                    
                if ccennext < 1.e-4 and ccen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('C')

                if ccennext < 1.e-5 and ccen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('C_end')
                    nemax = abund_list[i][specie_index['Ne-20']]

                if ccennext < 1.e-6 and ccen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('C')
                    
                if ccennext < 1.e-9 and ccen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('C')

            # Ne-burning
            necen  = abund_list[i][specie_index['Ne-20']]
            necennext = abund_list[i+1][specie_index['Ne-20']]
            if hcen < 1.e-5 and hecen < 1.e-5 and ccen < 1.e-3 and necen > 1.e-10:
                if necennext < nemax-0.003 and necen >= nemax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(necen)
                    burn_type.append('Ne_start')
        
                if necennext < 1.e-1 and necen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('Ne')
                    
                if necennext < 1.e-2 and necen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('Ne')
                    

                if necennext < 1.e-3 and necen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('Ne_end')
                    omax = abund_list[i][specie_index['O-16']]
                    
                if necennext < 1.e-4 and necen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('Ne')

                if necennext < 1.e-5 and necen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('Ne')

                if necennext < 1.e-6 and necen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('Ne')

                if necennext < 1.e-9 and necen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('Ne')
                    
            # O-burning
            ocen  = abund_list[i][specie_index['O-16']]
            ocennext = abund_list[i+1][specie_index['O-16']]
            if hcen < 1.e-5 and hecen < 1.e-5 and ccen < 1.e-3 and ocen > 1.e-10:
                if ocennext < omax-0.003 and ocen >= omax-0.003:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(ccen)
                    burn_type.append('O_start')
        
                if ocennext < 1.e-1 and ocen >= 1.e-1:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-1)
                    burn_type.append('O')
                    
                if ocennext < 1.e-2 and ocen >= 1.e-2:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-2)
                    burn_type.append('O')

                if ocennext < 1.e-3 and ocen >= 1.e-3:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-3)
                    burn_type.append('O')
                    
                if ocennext < 1.e-4 and ocen >= 1.e-4:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-4)
                    burn_type.append('O')

                if ocennext < 1.e-5 and ocen >= 1.e-5:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-5)
                    burn_type.append('O_end')

                if ocennext < 1.e-6 and ocen >= 1.e-6:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-6)
                    burn_type.append('O')
                    
                if ocennext < 1.e-9 and ocen >= 1.e-9:
                    burn_cycles.append(cycles_list[i])
                    burn_ages.append(age_list[i])
                    burn_abun.append(1.e-9)
                    burn_type.append('O')
                        
        
	print 'passa 4'

        pair = False
        age1 = -1
        for i in range(len(burn_type)):
            if 'start' in burn_type[i] and pair == False:
                age1 = burn_ages[i]
                pair = True
            elif 'end' in burn_type[i] and pair == True:
                age2 = burn_ages[i]
                pair = False
                if age1 != -1:
                    burn_lifetime.append(age2 - age1)
                    age1 = -1
                    
        return [burn_cycles, burn_ages, burn_abun, burn_type, burn_lifetime]


    
    def cores(self, incycle, **keyw):
        """        
        This function uses the abundances as a function of time to return
        core masses.  These core masses are:
        
        M_alpha  : mass of the helium core
        M_CO     : mass of the carbon-oxygen core
        M_si     : mass of the silicon core
        M_fe     : mass of the iron core
        M_rem    : mass of the remnant star post-supernova
    
        i)  cores(cycle)
    
        where cycle is the cycle to choose where to take the core masses
        
        A list is returned containing the following information:
        
        [cores, core_type, core_info]
    
        cores contains the core masses for the types found in core_type.  Core_info
        stores a variety of information about the abundances at certain mass coordinates.
    
        The following keywords can also be used:
    
        | Keyword argument | Default Value:
        ------------------------------------------------------------------------------
         abund               "iso_massf"
         isoa                "A"
         isom                "isomeric_state"
         isoz                "Z"
         mass                "mass"
         core_opt             0:  Use this for stellar evolution output
                              1:  Use this for MPPNP output
        
        All arguments change the name of fields used to read data from HDF5 files,
        other than option.  Option controls which scheme to use for calculating the
        core masses:
            core_opt = 0:  This uses alpha-isotopes to calculate the silicon and iron
                           cores.
                           Si_core = Si28+S32+Ar36+Ca40+Ti44
                           Fe_core = Cr48+Fe52+Ni56
            core_opt = 1:  This uses all elements in a network with proton number
                           Z=23 as a boundary.
                           Si_core = Sum of all isotopes with 14 <= Z <= 22
                           Fe_core = Sum of all isotopes with 23 <= Z <= 28
        """

        def infomod(core_opt, *inp):
            """
            This calls the correct infomod function depending on the value of core_opt
            """
            if core_opt == 0:
                return infomod1(inp[0], inp[1], inp[2])
            elif core_opt == 1:
                return infomod2(inp[0], inp[1], inp[2], inp[3], inp[4])
                
        def infomod1(shell, yps, isoindex):
            """
            Function for defining data to print into a string.
            This is used for the case of core_opt = 0
            """
            
            xh = yps[shell][isoindex[0]]
            xhe = yps[shell][isoindex[1]]
            xc = yps[shell][isoindex[2]]
            xo = yps[shell][isoindex[3]]
            xsi = yps[shell][isoindex[4]]
            xs = yps[shell][isoindex[5]]
            xa = yps[shell][isoindex[6]]
            xca = yps[shell][isoindex[7]]
            xti = yps[shell][isoindex[8]]
            xcr = yps[shell][isoindex[9]]
            xfe = yps[shell][isoindex[10]]
            xni = yps[shell][isoindex[11]]
            
            xsicore = xsi + xs + xa + xca + xti
            xfecore = xcr + xfe + xni
            
            return ' h=' + "%12.4e"%(xh) + ' he=' + "%12.4e"%(xhe) + \
            ' c=' + "%12.4e"%(xc) + ' o=' + "%12.4e"%(xo) + \
            ' si=' + "%12.4e"%(xsicore) + ' fe=' + "%12.4e"%(xfecore)
    
        def infomod2(shell, yps, isoindex, xsicore, xfecore):
            """
            Function for defining data to print into a string.
            This is used for the case of core_opt = 1
            """
            xsicore = 0.
            xsicoren = 0.
            xfecore = 0.
            xfecoren = 0.
            
            xh = yps[shell][isoindex[0]]
            xhe = yps[shell][isoindex[1]]
            xc = yps[shell][isoindex[2]]
            xo = yps[shell][isoindex[3]]
            
            for j in iso_si:
                xsicore += yps[i][j]
                xsicoren += yps[i+iter][j]
            for j in iso_fe:
                xfecore += yps[i][j]
                xfecoren += yps[i+iter][j]
            
            return 'shellnb = ' + str(shell) + ' h=' + "%12.4e"%(xh) + \
            ' he=' + "%12.4e"%(xhe) + ' c=' + "%12.4e"%(xc) + \
            ' o=' + "%12.4e"%(xo) + ' si=' + "%12.4e"%(xsicore) + \
            ' fe=' + "%12.4e"%(xfecore)

        if ("isoa" in keyw) == False:
            keyw["isoa"] = "A"
        if ("isoz" in keyw) == False:
            keyw["isoz"] = "Z"
        if ("mass" in keyw) == False:
            keyw["mass"] = "mass"
        if ("abund" in keyw) == False:
            keyw["abund"] = "iso_massf"
        if ("coreopt" in keyw) == False:
            core_opt = 1
        else:
            core_opt = keyw["coreopt"]
            
        cores = [0,0,0,0,0,0]
        core_type = ["","","","","",""]
        core_info = []
        iso_si = []
        iso_fe = []
        first = True
                        
        xm_list     = self.se.get(incycle, keyw["mass"])
        abund_list  = self.se.get(incycle, keyw["abund"])
        tables = self.se.Tables
        isoa_list   = tables[0]
        isoz_list   = tables[1]
                        
        # Check the order of the input data
        xmass = xm_list
        yps = abund_list
        centre = 0
        surface = len(xmass)-1
        iter = -1
        if isinstance(xmass, list) == True:
            if xmass[0] > xmass[1]:
                # mass is descending with shell number and the centre of the star
                # is the last shell
                centre = len(xmass)-1
                surface = 0
                iter = 1
            
        mco01 = 0.
        na    = 0
        nco   = 0

        h1index = -1
        he4index = -1
        c12index = -1
        o16index = -1
        si28index = -1
        s32index = -1
        a36index = -1
        ca40index = -1
        ti44index = -1
        cr48index = -1
        fe52index = -1
        ni56index = -1
        
        # Determine yps indices for certain abundances
        for i in range(len(isoa_list)):
            try:
                A = int(isoa_list[i])
                Z = int(isoz_list[i])
            except TypeError:
                A = int(isoa_list[i][0])
                Z = int(isoz_list[i][0])
            
            if A == 1 and Z == 1:
                h1index = i
            elif A == 4 and Z == 2:
                he4index = i
            elif A == 12 and Z == 6:
                c12index = i
            elif A == 16 and Z == 8:
                o16index = i
            
            if core_opt == 0:
                if A == 28 and Z == 14:
                    si28index = i
                elif A == 32 and Z == 16:
                    s32index = i
                elif A == 36 and Z == 18:
                    a36index = i
                elif A == 40 and Z == 20:
                    ca40index = i
                elif A == 44 and Z == 22:
                    ti44index = i
                elif A == 48 and Z == 24:
                    cr48index = i
                elif A == 52 and Z == 26:
                    fe52index = i
                elif A == 56 and Z == 38:
                    ni56index = i

        if h1index == -1 or he4index == -1 or c12index == -1 or o16index == -1:
            print "A key isotope(s) is not found in network!  Please check for the \
            presence of H1, He4, C12 and O16."
            os.sys.exit()
        
        # Si_core = Si28+S32+Ar36+Ca40+Ti44
        # Fe_core = Cr48+Fe52+Ni56  
        if core_opt == 0:
            if si28index == -1 or s32index == -1 or a36index == -1 or ca40index == -1 \
            or ti44index == -1 or cr48index == -1 or fe52index == -1 or ni56index == -1:
                print "Key isotopes for measuring the core mass with core_opt = 0 \
                are missing.  Setting core_opt = 1."
                core_opt = 1

        # Si_core = Sum of all isotopes with 14 <= Z <= 22
        # Fe_core = Sum of all isotopes with 23 <= Z <= 28
        if core_opt == 1:
            for i in range(len(isoa_list)):
                try:
                    A = int(isoa_list[i])
                    Z = int(isoz_list[i])
                except TypeError:
                    A = int(isoa_list[i][0])
                    Z = int(isoz_list[i][0])
                if Z >= 14 and Z <= 22:
                    iso_si.append(i)
                if Z >= 23 and Z <= 28:
                    iso_fe.append(i)

        isoindex = [h1index, he4index, c12index, o16index]
        if core_opt == 0:
            isoindex.extend([si28index, s32index, a36index, ca40index,
            ti44index, cr48index, fe52index, ni56index])
        
        if first == True:
            first = False
            cores[0] = xmass[surface]
            core_type[0] = "Final"
            core_info.append(infomod(core_opt, surface, yps, isoindex, iso_si, iso_fe))
        
        # Iterate over shells to determine the core masses
        for i in np.arange(surface, centre, iter):
            
            xsicore = 0.
            xfecore = 0.
            xsicoren = 0.
            xfecoren = 0.
            
            # Abundances of key isotopes at this shell (and the next one)
            xh = yps[i][h1index]
            xhe = yps[i][he4index]
            xc = yps[i][c12index]
            xo = yps[i][o16index]

            xhn = yps[i+iter][h1index]
            xhen = yps[i+iter][he4index]
            xcn = yps[i+iter][c12index]
            xon = yps[i+iter][o16index]
            
            xcocore = xc + xo
            xcocoren = xcn + xon
            
            if core_opt == 0:
                xsi = yps[i][si28index]
                xs = yps[i][s32index]
                xa = yps[i][a36index]
                xca = yps[i][ca40index]
                xti = yps[i][ti44index]
                xcr = yps[i][cr48index]
                xfe = yps[i][fe52index]
                xni = yps[i][ni56index]
                
                xsin = yps[i+iter][si28index]
                xsn = yps[i+iter][s32index]
                xan = yps[i+iter][a36index]
                xcan = yps[i+iter][ca40index]
                xtin = yps[i+iter][ti44index]
                xcrn = yps[i+iter][cr48index]
                xfen = yps[i+iter][fe52index]
                xnin = yps[i+iter][ni56index]
                
                xsicore = xsi + xs + xa + xca + xti
                xfecore = xcr + xfe + xni
                
                xsicoren = xsin + xsn + xan + xcan + xtin
                xfecoren = xcrn + xfen + xnin
                
            elif core_opt == 1:
                for j in iso_si:
                    xsicore += yps[i][j]
                    xsicoren += yps[i+iter][j]
                for j in iso_fe:
                    xfecore += yps[i][j]
                    xfecoren += yps[i+iter][j]

            # M_alpha:
            if xh >= 1.e-9:
                if xhen > 0.75 and xhe <= 0.75:
                    if cores[1] == 0:
                        cores[1] = xmass[i]
                        core_type[1] = 'Malpha 75%'
                        core_info.append('Malpha 75% '+str(xmass[i-iter])+ \
                        infomod(core_opt, i-iter, yps, isoindex, iso_si, iso_fe))
                        core_info.append('Malpha 75% '+str(xmass[i])+ \
                        infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                        core_info.append('Malpha 75% '+str(xmass[i+iter])+ \
                        infomod(core_opt, i+iter, yps, isoindex, iso_si, iso_fe))
                if xhen > 0.5 and xhe <= 0.5:
                    core_info.append('Malpha 50% '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                if xhn < 1.e-2 and xh >= 1.e-2:
                    na=i
                    core_info.append('Malpha 1.e-2 '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                if xhn < 1.e-4 and xh >= 1.e-4:
                    core_info.append('Malpha 1.e-4 '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                if xhn < 1.e-5 and xh >= 1.e-5:
                    core_info.append('Malpha 1.e-5 '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
        
            # M_CO:
            if xhe >= 1.e-9 and xh < 1.e-9:
                if xcocoren > 0.5 and xcocore <= 0.5:
                    core_info.append('Mco 50% '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                if xcocoren > 0.75 and xcocore <= 0.75:
                    core_info.append('Mco 75% '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                if xhen < 1.e-2 and xhe >= 1.e-2:
                    if cores[2] == 0:
                        mco01=xmass[i]
                        nco=i
                        cores[2] = xmass[i]
                        core_type[2] = 'Mco'
                        core_info.append('Mco 1.e-2 '+str(xmass[i-iter])+ \
                        infomod(core_opt, i-iter, yps, isoindex, iso_si, iso_fe))
                        core_info.append('Mco 1.e-2 '+str(xmass[i])+ \
                        infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                        core_info.append('Mco 1.e-2 '+str(xmass[i+iter])+ \
                        infomod(core_opt, i+iter, yps, isoindex, iso_si, iso_fe))
                if xhen < 1.e-4 and xhe >= 1.e-4:
                    core_info.append('Mco 1.e-4 '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                if xhen < 1.e-5 and xhe >= 1.e-5:
                    core_info.append('Mco 1.e-5 '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))

            # M_Si:
            if xcocore >= 1.e-7:
                if xsicoren > 0.5 and xsicore <= 0.5:
                    if cores[3] == 0:
                        core_info.append('Msi 50% '+str(xmass[i-iter])+ \
                        infomod(core_opt, i-iter, yps, isoindex, iso_si, iso_fe))
                        core_info.append('Msi 50% '+str(xmass[i])+ \
                        infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                        core_info.append('Msi 50% '+str(xmass[i+iter])+ \
                        infomod(core_opt, i+iter, yps, isoindex, iso_si, iso_fe))
                        cores[3] = xmass[i]
                        core_type[3] = 'Msi 50%'
                if xsicoren > 0.75 and xsicore <= 0.75:
                    core_info.append('Msi 75% '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                if xcocoren > 1.e-2 and xcocore <= 1.e-2:
                    core_info.append('Msi 1.e-2 '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
            
            # M_fe:
            if xfecoren > 0.5 and xfecore <= 0.5:
                if cores[4] == 0:
                    core_info.append('Mfe 50% '+str(xmass[i-iter])+ \
                    infomod(core_opt, i-iter, yps, isoindex, iso_si, iso_fe))
                    core_info.append('Mfe 50% '+str(xmass[i])+ \
                    infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                    core_info.append('Mfe 50% '+str(xmass[i+iter])+ \
                    infomod(core_opt, i+iter, yps, isoindex, iso_si, iso_fe))
                    cores[4] = xmass[i]
                    core_type[4] = 'Mfe 50%'
            if xfecoren > 0.75 and xfecore <= 0.75:
                core_info.append('Mfe 75% '+str(xmass[i])+ \
                infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
            if xsicoren > 1.e-2 and xsicore <= 1.e-2:
                core_info.append('Mfe 1.e-2 '+str(xmass[i])+ \
                infomod(core_opt, i, yps, isoindex, iso_si, iso_fe))
                    
        mcof=mco01
        i=na
        while i < nco:
            i+=1
            xco = yps[i][c12index] + yps[i][o16index]
            xcop = yps[i-iter][c12index] + yps[i-iter][o16index]
            mcof+=(xco+xcop)/2.*(xmass[i-iter]-xmass[i])
        cores[2] = mcof
        core_type[2] = 'Mco int'
    
        # M_rem vs M_CO relationship
        xmco = np.array([0.72,1.01,1.38,2.10,2.97,4.94,7.20,13.88,24.75,38.12,58.73])
        xmre = np.array([0.72,1.01,1.23,1.42,1.65,2.15,2.72,4.3,7.6,11.4,17.1])
    
        if 'Mco' in core_type:
            for i in np.arange(len(xmco)):
                if cores[2] < xmco[i]:
                    xco = yps[i][c12index] + yps[i][o16index]
                    xcop = yps[i-iter][c12index] + yps[i-iter][o16index]
                    cores[5] = xmre[i-1]+(xmre[i]-xmre[i-1])*(cores[j]-xmco[i-1])/(xmco[i]-xmco[i-1])
                    core_type[5] = 'Mrem'
                    break
                
        return [cores, core_type, core_info]    
    
    def presnyields(self, *cycles, **keyw):
        """
        This function calculates the presupernova yields of a full structure profile
        from a remnant mass, mrem, to the surface.  The command is:
    
        presnyields(cycle1, cycle2, mrem = 0.0)
    
        where cycle2 is optional.  cycle1 is the cycle to perform the
        presupernova yields calculations.  If cycle2 is also specified, the
        yields are outputted using `initial' abundances from cycle2,
        otherwise the ejected masses are outputted.
        
        mrem is specified using a keyword argument and tells the program
        where to begin integrating.
                
        The following keywords can be used:
    
        | Keyword argument | Default Value:
        ------------------------------------------------------------------------------
         abund               "iso_massf"
         xm                  "mass"
         mrem                0
    
        abund and xm are used when the variables within the input file differ in
        their names.  The default values are set to the output typically found in an MPPNP
        output file.  For example, if the table for the abundances is called "abundances"
        instead of the default value, use abund = "abundances" as a keyword argument.
        """
        abund_list  = []
        xm_list  = []
        
        if ("xm" in keyw) == False:
            keyw["xm"] = "mass"
        if ("abund" in keyw) == False:
            keyw["abund"] = "iso_massf"
        if ("mrem" in keyw) == False:
            mrem = 0.
        else:
            mrem = keyw["mrem"]
    
        # Only two cycles are required in this program.
        # Any more will be ignored
        cylen = len(cycles)
        if cylen > 2:
            cylen = 2
        for i in range(cylen):
            cycle = cycles[i]
            abund_list.append(self.se.get(cycle, keyw['abund']))
            xm_list.append(self.se.get(cycle, keyw['xm']))
                
        isoX = abund_list[0]
        xm = xm_list[0]
        if cylen == 2:
            isoXini = abund_list[1] # initial abundances
        
        niso = len(isoX[0,:])
        nshells = len(xm)
        
        X_i = np.zeros([niso], float)
        ybound = np.zeros([niso], float)
        xarray = np.zeros([nshells+1], float)
        yarray = np.zeros([nshells+1], float)
    
        # This part determines the index of the mass coordinate array which
        # is closest to the specified remnant mass.  This is used in the next
        # part, which is used to interpolate the abundances at the remnant mass
        for k in range(nshells):
            k1 = k
            if mrem<=xm[k]:
                break
                
        # This part is the interpolation of the isotopes found at the remnant
        # mass.  
        for i in range(niso):
            if k1>=1:
                if isoX[k1-1,i]!=0.0:
                    m=(isoX[k1,i]-isoX[k1-1,i])/(xm[k1]-xm[k1-1])
                    ybound[i] = isoX[k1-1,i] +m*(mrem-xm[k1-1])
                else:
                    ybound[i]=1.e-99
            if k1==0:
                if isoX[k1,i]!=0.0:
                    ybound[i]=isoX[k1,i]
                else:
                    ybound[i]=1.e-99
        
        # This part merges the interpolated data and the existing arrays into
        # the arrays xarray and yarray.  Once this is done, the summation is
        # made. 
        xarray[0] = mrem
        xarray[1:nshells-k1+1] = xm[k1:nshells]
        for i in range(niso):
            yarray[0] = ybound[i]
            for j in range(nshells-k1):
                if isoX[k1+j,i] != 0.0:
                    yarray[j+1] = isoX[k1+j,i]
                else:
                    yarray[j+1] = 1.e-99
        
            if cylen == 1:
                # Calculate the ejected masses
                for j in range(nshells-k1):
                    X_i[i] = X_i[i] + ((0.5*(yarray[j+1] + yarray[j])) * \
                    (xarray[j+1] - xarray[j]))
            elif cylen == 2:
                # Calculate the SN yield.
                for j in range(nshells-k1):
                    X_i[i] = X_i[i] + ((0.5*(yarray[j+1] + yarray[j]) - isoXini[-1,i]) * \
                    (xarray[j+1] - xarray[j]))

        
        return X_i
    
    def windyields(self, ini, end, delta, **keyw):
        """
        This function returns the wind yields and ejected masses.  The command
        is:

        X_i, E_i = data.windyields(ini, end, delta)

        where `ini' is the starting cycle, `end' is the finishing cycle and `delta' is the
        cycle interval.  `data' is a plot_tools object.  The function returns a list of
        wind yields, X_i, and a list of ejected masses, E_i, in whichever mass units were
        used (usually solar masses).
 
        The following keywords cand also be used:

        | Keyword argument | Default Value:
        ------------------------------------------------------------------------------
        abund               "iso_massf"
        tmass               "mass"
        cycle               "cycle"
  
        The keyword arguments are used when the variables within the input file differ in
        name from their default values typically found in an MPPNP output file.  If the
        data table differs in name, use these keywords.  For example, if the table for the
        abundances is called "abundances" instead of "iso_massf", then use
        abund = "abundances" as a keyword argument.
        """

        if ("tmass" in keyw) == False:
            keyw["tmass"] = "mass"
        if ("abund" in keyw) == False:
            keyw["abund"] = "iso_massf"
        if ("cycle" in keyw) == False:
            keyw["cycle"] = "cycle"        

        print "Windyields() initialised.  Reading files..."

        ypsinit = [] 
        niso = 0
        X_i = []
        E_i = []
        totalmass = []
        ypssurf  = []
        cycles = []
        first = True

        # The following statements copy global functions into local memory,
        # which is called faster, speeding up the code slightly
        wc = self.windcalc
        cycleret = self.se.cycles
        retrieve = self.se.get
        capp = cycles.extend
        tapp = totalmass.extend
        yapp = ypssurf.extend

        # Retrieve the data from the files
        for i in xrange(ini,end+1,delta):
            step = int(i)
            capp([int(cycleret[i-ini])])
            tapp([retrieve(step,keyw["tmass"])])
            yapp([retrieve(step,keyw["abund"])])
        
        print "Reading complete.  Calculating yields and ejected masses..."

        nsteps = len(cycles)-1
        niso = len(ypssurf[0])
        X_i = np.zeros([niso], float)
        E_i = np.zeros([niso], float)
        # Call the windyields calculator
        X_i, E_i = wc(first, totalmass, nsteps, niso, ypssurf, \
        ypsinit, X_i, E_i, cycles)
        return X_i, E_i


    def windcalc(self, first, totalmass, nsteps, niso, ypssurf, ypsinit, \
    X_i, E_i, cycles):
        # This function calculates the windyields and ejected masses as called from
        # windyields().  It uses a summation version of the formulae used in Hirschi
        # et al. 2005, "Yields of rotating stars at solar metallicity".
 
        # If it is the first file, the arrays need to be created and the initial
        # abundances set

        if first == True:
            X_i = np.zeros([niso], float)
            E_i = np.zeros([niso], float)
            ypsinit = ypssurf[0]
            for m in xrange(niso):
                for n in xrange(nsteps):
                    X_i[m] = X_i[m] + ((totalmass[n] - totalmass[n+1]) * \
                    (0.5 * (ypssurf[n][m] + ypssurf[n+1][m]) - ypsinit[m]))
                    E_i[m] = E_i[m] + ((totalmass[n] - totalmass[n+1]) * \
                    (0.5 * (ypssurf[n][m] + ypssurf[n+1][m])))
        
        else:
            for m in range(niso):
                for n in xrange(nsteps):
                    X_i[m] = X_i[m] + ((totalmass[n] - totalmass[n+1]) * \
                    (0.5 * (ypssurf[n][m] + ypssurf[n+1][m]) - ypsinit[m]))
                    E_i[m] = E_i[m] + ((totalmass[n] - totalmass[n+1]) * \
                    (0.5 * (ypssurf[n][m] + ypssurf[n+1][m])))

        return X_i, E_i


    def average_iso_abund_marco(self,mass_range,cycle,stable,i_decay):
    	''' Interface to average over mass_range. 
    	directory     -  location of h5 file to plot. Needed for plot_tools
    	name_h5_file  -  name of h5 file. Needed for plot_tools
    	mass_range    - required to plot data in a certain mass range. Needed for read_iso_abund_marco
    	cycle         - which cycle from the h5 file?. Needed for read_iso_abund_marco
    	stable        - logic if want to plot only stable or not.  
    	i_decay       - if = 1 I plot not decayed, if = 2 I plot decayed. Make sense only if stable is true.'''


    	import utils as u

    	if not stable and i_decay == 2:
        	print 'ERROR: choose i_decay = 1'  
        	return
    
    	#data=mp.se(directory,name_h5_file)
    	self.read_iso_abund_marco(mass_range,cycle)
    	#print spe
    	if i_decay == 2:
        	u.stable_specie()
        	self.decay(self.mass_frac)
    
    	# here I am calculating average mass fraction for all isotopes in given mass range, and then
    	# if needed calculating average over decayed.
    	# warning: mass_range is bigger than used_masses range, by definition. Should I use it?
    	print 'average over used_masses range, not over original mass_range'
    	print used_masses[0],used_masses[len(used_masses)-1],'instead of',mass_range[0],mass_range[1]
    
    	global average_mass_frac
    	average_mass_frac = []
    
    	if len(used_masses) >= 2:
     		dm_tot = abs(used_masses[len(used_masses)-1]-used_masses[0])
        	for j in range(len(u.spe)-1):
            		temp = 0.
            		for i in range(len(used_masses)-1):
	        		dm_i = abs(used_masses[i+1]-used_masses[i])
            			temp = float(self.mass_frac[i][j]*dm_i/dm_tot) + temp
            		average_mass_frac.append(temp)
        	#print average_mass_frac
    	elif  len(used_masses) == 1:
        	print 'case with 1 mass zone only, not implemented yet'
    
    
        
    	somma = 0.
    	somma = sum(average_mass_frac)
   	print 'departure from 1 of sum of average_mass_frac=',abs(1. - somma)
    
    	# not let's do it over decayed also, if i_decay = 2
    	if i_decay == 2:
        	global average_mass_frac_decay
        	average_mass_frac_decay = []
		dm_tot = abs(used_masses[len(used_masses)-1]-used_masses[0])
		#
		#print len(decayed_multi_d[0]),decayed_multi_d[0]        
		for j in range(len(u.back_ind)):
        		temp = 0.
           		for i in range(len(used_masses)-1):
		        	dm_i = abs(used_masses[i+1]-used_masses[i])
	            		temp = float(decayed_multi_d[i][j]*dm_i/dm_tot) + temp
            		average_mass_frac_decay.append(temp)

        	somma = 0.
		somma = sum(average_mass_frac_decay)
        	print 'departure from 1 of sum of average_mass_frac_decay=',abs(1. - somma)

    	# now I have the average abundances. We can do the plot.
    
    def get_elem_names(self):
    	''' returns for one cycle an element name dictionary.'''

    	import utils as u

	# provide library for Z versus element names, and Z for elements        
	#element_name = self.se.elements
	element_name = self.elements_names
	u.give_zip_element_z_and_names(element_name)
	self.z_of_element_name = u.index_z_for_elements      
	# notice that z is not the element index! Z and indexing of elements can be different
	

    def get_abundance_iso_decay(self,cycle):
    	''' returns the decayed stable isotopes.
	If average is True, than is peforming an average over the given mass range.'''

    	import utils as u

	masses_for_this_cycle = self.se.get(cycle,'mass')
    	self.read_iso_abund_marco([min(masses_for_this_cycle),max(masses_for_this_cycle)],cycle)

        u.stable_specie()
        self.decay(self.mass_frac)

	self.index_for_all_species = u.cl
	self.index_for_stable_species = u.back_ind

	self.decayed_stable_isotopes_per_cycle = decayed_multi_d	

	# from here read solar abundances
	solar_factor = 2.
	u.solar('iniab1.0E-02.ppn_GN93',solar_factor)

	self.stable_isotope_identifier=u.jjdum
	self.stable_isotope_list=u.stable

	self.isotopic_production_factors=[]
	for i in range(len(masses_for_this_cycle)):
		pf_dum=[]
		jj=0
    		for j in range(len(self.stable_isotope_identifier)):
        		if self.stable_isotope_identifier[j] == 1:
            			pf_dum.append(float(self.mass_frac[i][self.index_for_all_species[self.stable_isotope_list
[jj].capitalize()]]/u.solar_abundance[self.stable_isotope_list[jj].lower()]))
				jj=jj+1
        		#elif self.stable_isotope_identifier[j] == 0:
            		#	pf_dum.append(float(0.))
		self.isotopic_production_factors.append(pf_dum)	

	self.isotopic_production_factors_decayed=[]
	for i in range(len(masses_for_this_cycle)):
		pf_dum_d=[]
		jj=0
    		for j in range(len(self.stable_isotope_identifier)):
        		if self.stable_isotope_identifier[j] == 1:
				pf_dum_d.append(float(self.decayed_stable_isotopes_per_cycle[i][self.index_for_stable_species[self.stable_isotope_list
[jj].upper()]]/u.solar_abundance[self.stable_isotope_list[jj].lower()]))
				jj=jj+1
	        self.isotopic_production_factors_decayed.append(pf_dum_d)	



    def get_abundance_elem(self,cycle):
    	''' returns the undecayed element profile (all elements that are in elem_names).'''

    	import utils as u

	masses_for_this_cycle = self.se.get(cycle,'mass')
    	self.read_iso_abund_marco([min(masses_for_this_cycle),max(masses_for_this_cycle)],cycle)

        u.stable_specie()
        self.decay(self.mass_frac)

	# provide library for Z versus element names, and Z for elements        	
	element_name = self.se.elements
	u.give_zip_element_z_and_names(element_name)
	# from here read solar abundances
	solar_factor = 2.
	u.solar('iniab1.0E-02.ppn_GN93',solar_factor)

	self.stable_isotope_identifier=u.jjdum
	self.stable_isotope_list=u.stable
		

	self.element_abundance_not_decayed=[]
	self.element_abundance_decayed =[]
	self.element_production_factors=[]
	self.element_production_factors_decayed=[]
	for i in range(len(masses_for_this_cycle)):
		mass_fractions_array_decayed = decayed_multi_d[i]
		mass_fractions_array_not_decayed = self.mass_frac[i]
		u.element_abund_marco(2,self.stable_isotope_list,self.stable_isotope_identifier,mass_fractions_array_not_decayed,mass_fractions_array_decayed)
		self.element_abundance_not_decayed.append(u.elem_abund)
		self.element_abundance_decayed.append(u.elem_abund_decayed)
		self.element_production_factors.append(u.elem_prod_fac)
		self.element_production_factors_decayed.append(u.elem_prod_fac_decayed)
		
	#self.decayed_stable_isotopes_per_cycle = decayed_multi_d	
	

def obsolete_plot_iso_abund_marco(directory,name_h5_file,mass_range,cycle,logic_stable,i_decay,file_solar,solar_factor):
    ''' Interface to plot average over mass_range. 
    directory     -  location of h5 file to plot. Needed for plot_tools
    name_h5_file  -  name of h5 file. Needed for plot_tools
    mass_range    - required to plot data in a certain mass range. Needed for read_iso_abund_marco
    cycle         - which cycle from the h5 file?. Needed for read_iso_abund_marco
    logic_stable  - logic if want to plot only stable or not.  
    i_decay       - if = 1 I plot not decayed, if = 2 I plot decayed. Make sense only if stable is true
    file_solar        - file where to take solar abundances
    solar_factor      - float to correct initial abundances to solar, e.g. for Z=0.01 and AG89 solar_factor = 2.'''


    # provide library for Z versus element names, and Z for elements 	
    u.give_zip_element_z_and_names()
    # solar abundances are read here
    u.solar(file_solar,solar_factor)
    # from here I have average abundances in mass_range to plot
    average_iso_abund_marco(mass_range,cycle,logic_stable,i_decay)
    
    fig = pl.figure()            # Figure object
    ax = fig.add_subplot(1,1,1)     # Axes object: one row, one column, first plot (one plot!)
    # Tick marks
    xminorlocator = MultipleLocator(1)
    xmajorlocator = MultipleLocator(10)
    ax.xaxis.set_major_locator(xmajorlocator)
    ax.xaxis.set_minor_locator(xminorlocator)
    yminorlocator = MultipleLocator(0.1)
    ymajorlocator = MultipleLocator(1)
    ax.yaxis.set_major_locator(ymajorlocator)
    ax.yaxis.set_minor_locator(yminorlocator)

    ax.set_yscale('log')
 
    if not logic_stable:
        for i in range(len(u.spe)):
            pl.plot(amass_int[u.cl[spe[i]]],average_mass_frac[u.cl[spe[i]]],'ko')

        pl.xlabel('$Mass$ $number$', fontsize=20)
        pl.ylabel('$X_{i}$', fontsize=20)

        pl.ylim(1.0e-10,10.)
        pl.xlim(55,110)
    
    elif logic_stable:
    # plot stable
        #for i in range(len(stable)):
        #    pl.plot(amass_int[cl[stable[i].capitalize()]],average_mass_frac[cl[stable[i].capitalize()]]/u.solar_abundance[stable[i].lower()],'ko')
        
        if i_decay == 2:
            for j in range(len(stable)):
                    #print cl[stable[j].capitalize()],stable[j].capitalize(),amass_int[cl[stable[j].capitalize()]]
                    pl.plot(amass_int[u.cl[stable[j].capitalize()]],u.mass_fractions_array_decayed[u.back_ind[stable[j]]]/u.solar_abundance[stable[j].lower()],'Dk')
    
        for i in range(len(stable)):
            for j in range(len(stable)): 
                if stable[i][:2] == stable[j][:2]:
                    if stable[i] == stable[j-1]:
                        adum  =[amass_int[u.cl[stable[i].capitalize()]],amass_int[u.cl[stable[j].capitalize()]]]
                        mfdum =[float(average_mass_frac[u.cl[stable[i].capitalize()]])/float(u.solar_abundance[stable[i].lower()]),float(average_mass_frac[u.cl[stable[j].capitalize()]])/float(u.solar_abundance[stable[j].lower()])]
                        mfddum=[float(u.average_mass_frac_decay[u.back_ind[stable[i]]])/float(u.solar_abundance[stable[i].lower()]),float(u.average_mass_frac_decay[u.back_ind[stable[j]]])/float(u.solar_abundance[stable[j].lower()])]
                        #pl.plot(adum,mfdum,'k-')
			# I had to add this try/except...why? I guess is someone related to H2, that I spotted that was wrong in stable_raw...
			# should deal without this. Have to be solved when I have time Marco (June 7 2011)
                    	if i_decay == 2:
				try:
                        		pl.plot(adum,mfddum,'k-')
				except UnboundLocalError:
					continue	  

    pl.xlabel('$Mass$ $number$', fontsize=20)
    pl.ylabel('$X_{i}/X_{sun}$', fontsize=20)

    pl.ylim(1.0e-3,5000.)
    pl.xlim(55,210)
    
    
    
    pl.grid()
    pl.show() 
      


                         

def obsolete_plot_el_abund_marco(directory,name_h5_file,mass_range,cycle,logic_stable,i_decay,file_solar,solar_factor,symbol='ko'):
    ''' Interface to plot elements abundances averaged over mass_range. 
    directory     -  location of h5 file to plot. Needed for plot_tools
    name_h5_file  -  name of h5 file. Needed for plot_tools
    mass_range    - required to plot data in a certain mass range. Needed for read_iso_abund_marco
    cycle         - which cycle from the h5 file?. Needed for read_iso_abund_marco
    logic_stable  - logic if want to plot only stable or not.  
    i_decay       - if = 1 I plot not decayed, if = 2 I plot decayed. Make sense only if stable is true
    file_solar        - file where to take solar abundances
    solar_factor      - float to correct initial abundances to solar, e.g. for Z=0.01 and AG89 solar_factor = 2.'''

    # provide library for Z versus element names, and Z for elements 	
    u.give_zip_element_z_and_names()
    # solar abundances are read here
    u.solar(file_solar,solar_factor)
    # from here I have average abundances in mass_range to plot
    average_iso_abund_marco(mass_range,cycle,logic_stable,i_decay)
    # element abundances are calculated here
    mass_fractions_array_decayed = average_mass_frac_decay 	
    mass_fractions_array_not_decayed = average_mass_frac 	
    u.element_abund_marco(i_decay,stable,jjdum,mass_fractions_array_not_decayed,mass_fractions_array_decayed)
    
    
    fig = pl.figure()            # Figure object
    ax = fig.add_subplot(1,1,1)     # Axes object: one row, one column, first plot (one plot!)
    # Tick marks
    xminorlocator = MultipleLocator(1)
    xmajorlocator = MultipleLocator(10)
    ax.xaxis.set_major_locator(xmajorlocator)
    ax.xaxis.set_minor_locator(xminorlocator)
    yminorlocator = MultipleLocator(0.1)
    ymajorlocator = MultipleLocator(1)
    ax.yaxis.set_major_locator(ymajorlocator)
    ax.yaxis.set_minor_locator(yminorlocator)

    ax.set_yscale('log')
 
    if not logic_stable:
        for i in range(u.z_bismuth):
            pl.plot(z_for_elem[i],elem_prod_fac[i],symbol,markersize=10.)

        pl.xlabel('$Atomic$ $number$', fontsize=20)
        pl.ylabel('$X_{i}/X_{sun}$', fontsize=20)

        pl.ylim(1.0e-2,1000.)
        pl.xlim(0,95)
    
    elif logic_stable:
        for i in range(u.z_bismuth):
            if index_stable[i] == 1:
		continue
                #pl.plot(z_for_elem[i],elem_prod_fac[i],'ko')
        if i_decay == 2:
            for i in range(u.z_bismuth):
                if index_stable[i] == 1:
                    pl.plot(z_for_elem[i],elem_prod_fac_decayed[i],symbol,markersize=10.)

        pl.xlabel('$Atomic$ $number$', fontsize=20)
        pl.ylabel('$X_{i}/X_{sun}$', fontsize=20)

        pl.ylim(1.0e-2,1000.)
        pl.xlim(0,95)
    
    
    
    pl.grid()
    pl.show() 
    
    
