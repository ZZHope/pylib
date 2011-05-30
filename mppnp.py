'''mppnp is a collection of plots of data in se-type h5 files

   sedir    - the directory on which the h5 files can be
   pattern - a string pattern that the file names have to
                    contain in order to be read

   

   usage:
   start by loading the module
   1 : import mppnp as mp

   you can get help:
   2 : help mp

   next, initiate a se class instance:   
   3 : pt=mp.se('path/to/dir/with/se-files')
   
   which would read all h5 files from the current directory; again, do
   help(mp.se) or mp.se? do get info for more options, as for example
   how to specify a match pattern to select only certain files in the
   directory; this may be useful if there are very many files and you
   want only every 10th; example:
   3 : pt =mp.se('.','M1.65Z0.020.00') 

      note: initialising an instance with 120 files with 1000 packets
      each may take 3-4 minutes, which is too long for many people;
      therefore the initialisation module will generate and write an
      index file; if such an index file is present initialisation will
      take less then a second; the initialisation module will
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

   Available cylces can be viewed with
      pt.se.cycles

   You can get any of the quanties via the 'get' method, which is
   relatively smart to give you things in various ways:
      pt.get('rho')

   would give you all rho vectors for all cycles, which is maybe more
   than you want, so try
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
            do someting with this_case

   There are various plotting methodes, including plotting abundance
   charts (abu_chart), isotopic abundance distributions and the
   generic 'plot' method that lets you plot every quantity against any
   other quantity that can possibly be plotted. Some of the methods
   (including the three just mentioned) are available via the super
   class data_plot, and are equally available in other python
   modulues, as for example mesa.py or ppn.py. Several methods accept
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

    def trajectory(self,ini,end,delta,mass_coo):
        
        ''' create a trajectory out of a stellar model 
        
        ini       - initial model
        end       - final model 
        delta     - sparsity factor of the frames
        mass_coo  - mass coordinate for the traj
        
        warning: remove the old trajectory, if you have any for the same mass coordinate.
        you are appending data, not overwriting.'''
        
        f = open('traj_'+str(mass_coo)+'.dat','a')
        for step in range(ini,end+1,delta):
                age=self.se.get(step,'age')
                mass=self.se.get(step,'mass')  
                temperature=self.se.get(step,'temperature')
                rho=self.se.get(step,'rho')

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
                                        

                string = str(step)+'  '+str(age)+'  '+str(temperature[zone])+'  '+str(rho[zone]) 
                f.write(string+"\n")

        f.close()

    def abund_at_masscoorinate(self,ini,end,delta,mass_coo):

        ''' create a trajectory out of a stellar model 

        ini       - initial model
        end       - final model 
        delta     - sparsity factor of the frames
        mass_coo  - mass coordinate for the traj
  
        warning: remove the old trajectory, if you have any for the same mass coordinate.
        you are appendind data, not overwriting.
        '''
        
        f = open('traj_'+str(mass_coo)+'.dat','a')
        for step in range(ini,end+1,delta):
            age=self.se.get(step,'age')
            mass=self.se.get(step,'mass')  
            temperature=self.se.get(step,'temperature')
            rho=self.se.get(step,'rho')
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
            string = str(step)+'  '+str(age)+'  '+str(temperature[zone])+'  '+str(rho[zone]) 
            f.write(string+"\n")

        f.close()
         


    def kip(self, mix_thresh,xaxis,sparse):
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
        cyclelist = original_cyclelist[0:len(original_cyclelist):sparse]

        xx = self.se.get(cyclelist,'age')
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
            xx[-1] = xx[-2]*1.001
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
            d_coeff = self.se.get(cyclelist[i],'dcoeff')
            massco = self.se.get(cyclelist[i],'mass')
            plotlims = getlims(d_coeff,massco)
            for k in range(0,len(plotlims),2):
                ax.axvline(xx[i],ymin=plotlims[k]/m_ini,ymax=plotlims[k+1]/m_ini,color='b',linewidth=0.5)
        

        ax.plot(xx, totalmass, color='black', linewidth=1)
        if xaxis == 'log_time_left':
            ax.axis([max(xx),min(xx),0.,m_ini])
        else:
            ax.axis([min(xx),max(xx),0.,m_ini])
        ax.set_ylabel('Mass [$M_{\odot}$]',fontsize=fsize)

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
        
    '''
    def iso_abund(self, mass_range, cycle, stable):
        plot the abundance of all the chemical species
        inputs:
            mass_range - a 1x2 array containing the lower and upper mass range.  
                                     if None, it will plot over the entire range
            cycle       - a string/integer of the cycle of interest.
            stable     - a boolean of whether to filter out the unstables.
                
        
        elem_list = []
        elem_index = []
        isotope_to_plot = self.se.isotopes 
        abunds = self.se.get(cycle,'iso_massf')
        masses = []
        masses = self.se.get(cycle,'mass')
        if mass_range == None:
            print 'Using default mass range'
            mass_range = [min(masses),max(masses)]    
        masses.sort()
        mass_range.sort()
        
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
        
        
        
        print 'Using The following conditions:'
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
        #print elem_index
        #print '\n'
        if stable: #removing unstable elements
            
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

       
        while len(abunds) == 1:
            abunds = abunds[0]

        abund_plot = []
        mass_num = []
        
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
                                        print abunds[l][elem_index[index[j]][k]],abs(masses[l+1]-masses[l])
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
        
        
        
        
        #temp3 = []
        mass_num = []
        for j in xrange(len(index)):
            temp = []
            #temp2 = []    
            try:
                
                for k in xrange(len(elem_index[index[j]])):
                
                    temp.append(float(isotope_to_plot[elem_index[index[j]][k]].split('-')[1]))

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
        #print len(index)
        
        
        #print abund_plot,mass_num
        for j in xrange(len(index)):        #Loop through the elements of interest
            #    Process the line
            #print 'processing line'
            for l in xrange(len(abund_plot[j])):
                if abund_plot[j][l] == 0:
                    abund_plot[j][l] = 1e-20
                    
                    
            
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
                coordinates = [mass_num[j][abund_plot[j].index(max(abund_plot[j]))],max(abund_plot[j])]     
                pl.text(coordinates[0],coordinates[1], elem_list[index[j]])
                    #print self.parent.elem_list[self.index[j]],coordinates[0],coordinates[1]
            except ValueError:
                None
                #print 'Empty var:  ', abund_plot[j]
            except IndexError:
                None
                #print 'out of bounds: ', len(abund_plot), j
                    
            try:
                pl.semilogy(mass_num[j],abund_plot[j],'bo')
            except OverflowError:
                None
                #print 'div by zero', len(mass_num[j]), len(abund_plot[j])
            except IndexError:
                None
                #print 'out of bounds II: ', len(mass_num),len(abund_plot),j
        title = str('Abundance of Isotopes over range %4.2f' %mass_range[0]) + str('-%4.2f' %mass_range[1]) +\
                str(' for cycle %d' %int(cycle))
        
        pl.ylim([1e-13,10])
        pl.title(title)
        pl.xlabel('Mass Number')
        pl.ylabel('Relative Abundance')
        pl.grid()
        pl.show()
        return
        '''

    def read_iso_abund_marco(self, mass_range, cycle):
        ''' plot the abundance of all the chemical species
        inputs:
        mass_range - a 1x2 array containing the lower and upper mass range.  
                     if None, it will plot over the entire range
        cycle       - a string/integer of the cycle of interest.
        stable     - a boolean of whether to filter out the unstables.
            
        '''

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
        masses.sort()
        mass_range.sort()

        
        
        print 'Using The following conditions:'
        print '\tmass_range:', mass_range[0], mass_range[1]
        print '\tcycle:', cycle

        
        spe_rude1 = []
        spe_rude2 = []
        spe_rude3 = []
        for i in range(len(self.se.isotopes)):
            spe_rude1.append(self.se.isotopes[i].split('-')[0])
            spe_rude2.append(self.se.isotopes[i].split('-')[1])
        # spe_rude1 is elem name and spe_rude2 is mass number.
        #print spe_rude1,spe_rude2
        k = 0
        for i in range(len(spe_rude1)):
            try: 
            	if int(spe_rude2[i]) < 10:
            	    spe_rude3.append(str(spe_rude1[i][0:2])+str('  ')+str(spe_rude2[i][0:3]))    
            	elif int(spe_rude2[i]) >= 10 and int(spe_rude2[i]) < 100 :    
            	    spe_rude3.append(str(spe_rude1[i][0:2])+str(' ')+str(spe_rude2[i][0:3]))    
            	elif int(spe_rude2[i]) >= 100 :    
            	    spe_rude3.append(str(spe_rude1[i][0:2])+str(spe_rude2[i][0:3]))    
            except ValueError:
                k = k+1
                None

        global spe
        spe = []
        n_array = []
        for i in range(len(spe_rude3)):
            if len(str(spe_rude1[i])) == 1:
                spe.append(str(spe_rude3[i][0:1])+str(' ')+str(spe_rude3[i][1:4]))        
            else:
                spe.append(spe_rude3[i]) 
	    n_array.append(i)
        if spe[0]=='Ne  1':
		spe[0] = 'N   1' 
                                          
        # spe_rude is the isotope name, in agreement with what we use in ppn, etc.
        # need to do this to can use other functions without changing them drastically.

                 
        
        # here I skip isomers...
        global amass_int
        amass_int=np.zeros(len(spe_rude2)) 
        for i in range(len(spe_rude2)-k):   
            amass_int[i]=int(spe_rude2[i])
            #print amass_int


        # here I have to create an array for the atomic numbers.
        # I need to this when I calculate and plot element abundances 
        
        global znum_int
        znum_int=np.zeros(len(spe)) 
        
        for i in range(len(spe)):
            if str(spe[i][0:2]) == 'H ':
                znum_int[i] = 1        
            elif str(spe[i][0:2]) == 'He':
                znum_int[i] = 2
            elif str(spe[i][0:2]) == 'Li':
                znum_int[i] = 3
            elif str(spe[i][0:2]) == 'Be':
                znum_int[i] = 4
            elif str(spe[i][0:2]) == 'B ':
                znum_int[i] = 5
            elif str(spe[i][0:2]) == 'C ':
                znum_int[i] = 6
            elif str(spe[i][0:2]) == 'N ':
                znum_int[i] = 7
            elif str(spe[i][0:2]) == 'O ':
                znum_int[i] = 8
            elif str(spe[i][0:2]) == 'F ':
                znum_int[i] = 9
            elif str(spe[i][0:2]) == 'Ne':
                znum_int[i] = 10
            elif str(spe[i][0:2]) == 'Na':
                znum_int[i] = 11
            elif str(spe[i][0:2]) == 'Mg':
                znum_int[i] = 12
            elif str(spe[i][0:2]) == 'Al':
                znum_int[i] = 13
            elif str(spe[i][0:2]) == 'Si':
                znum_int[i] = 14
            elif str(spe[i][0:2]) == 'P ':
                znum_int[i] = 15
            elif str(spe[i][0:2]) == 'S ':
                znum_int[i] = 16
            elif str(spe[i][0:2]) == 'Cl':
                znum_int[i] = 17
            elif str(spe[i][0:2]) == 'Ar':
                znum_int[i] = 18
            elif str(spe[i][0:2]) == 'K ':
                znum_int[i] = 19
            elif str(spe[i][0:2]) == 'Ca':
                znum_int[i] = 20
            elif str(spe[i][0:2]) == 'Sc':
                znum_int[i] = 21
            elif str(spe[i][0:2]) == 'Ti':
                znum_int[i] = 22
            elif str(spe[i][0:2]) == 'V ':
                znum_int[i] = 23
            elif str(spe[i][0:2]) == 'Cr':
                znum_int[i] = 24
            elif str(spe[i][0:2]) == 'Mn':
                znum_int[i] = 25
            elif str(spe[i][0:2]) == 'Fe':
                znum_int[i] = 26
            elif str(spe[i][0:2]) == 'Co':
                znum_int[i] = 27
            elif str(spe[i][0:2]) == 'Ni':
                znum_int[i] = 28
            elif str(spe[i][0:2]) == 'Cu':
                znum_int[i] = 29
            elif str(spe[i][0:2]) == 'Zn':
                znum_int[i] = 30
            elif str(spe[i][0:2]) == 'Ga':
                znum_int[i] = 31
            elif str(spe[i][0:2]) == 'Ge':
                znum_int[i] = 32
            elif str(spe[i][0:2]) == 'As':
                znum_int[i] = 33
            elif str(spe[i][0:2]) == 'Se':
                znum_int[i] = 34
            elif str(spe[i][0:2]) == 'Br':
                znum_int[i] = 35
            elif str(spe[i][0:2]) == 'Kr':
                znum_int[i] = 36
            elif str(spe[i][0:2]) == 'Rb':
                znum_int[i] = 37
            elif str(spe[i][0:2]) == 'Sr':
                znum_int[i] = 38
            elif str(spe[i][0:2]) == 'Y ':
                znum_int[i] = 39
            elif str(spe[i][0:2]) == 'Zr':
                znum_int[i] = 40
            elif str(spe[i][0:2]) == 'Nb':
                znum_int[i] = 41
            elif str(spe[i][0:2]) == 'Mo':
                znum_int[i] = 42
            elif str(spe[i][0:2]) == 'Tc':
                znum_int[i] = 43
            elif str(spe[i][0:2]) == 'Ru':
                znum_int[i] = 44
            elif str(spe[i][0:2]) == 'Rh':
                znum_int[i] = 45
            elif str(spe[i][0:2]) == 'Pd':
                znum_int[i] = 46
            elif str(spe[i][0:2]) == 'Ag':
                znum_int[i] = 47
            elif str(spe[i][0:2]) == 'Cd':
                znum_int[i] = 48
            elif str(spe[i][0:2]) == 'In':
                znum_int[i] = 49
            elif str(spe[i][0:2]) == 'Sn':
                znum_int[i] = 50
            elif str(spe[i][0:2]) == 'Sb':
                znum_int[i] = 51
            elif str(spe[i][0:2]) == 'Te':
                znum_int[i] = 52
            elif str(spe[i][0:2]) == 'I ':
                znum_int[i] = 53
            elif str(spe[i][0:2]) == 'Xe':
                znum_int[i] = 54
            elif str(spe[i][0:2]) == 'Cs':
                znum_int[i] = 55
            elif str(spe[i][0:2]) == 'Ba':
                znum_int[i] = 56
            elif str(spe[i][0:2]) == 'La':
                znum_int[i] = 57
            elif str(spe[i][0:2]) == 'Ce':
                znum_int[i] = 58
            elif str(spe[i][0:2]) == 'Pr':
                znum_int[i] = 59
            elif str(spe[i][0:2]) == 'Nd':
                znum_int[i] = 60
            elif str(spe[i][0:2]) == 'Pm':
                znum_int[i] = 61
            elif str(spe[i][0:2]) == 'Sm':
                znum_int[i] = 62
            elif str(spe[i][0:2]) == 'Eu':
                znum_int[i] = 63
            elif str(spe[i][0:2]) == 'Gd':
                znum_int[i] = 64
            elif str(spe[i][0:2]) == 'Tb':
                znum_int[i] = 65
            elif str(spe[i][0:2]) == 'Dy':
                znum_int[i] = 66
            elif str(spe[i][0:2]) == 'Ho':
                znum_int[i] = 67
            elif str(spe[i][0:2]) == 'Er':
                znum_int[i] = 68
            elif str(spe[i][0:2]) == 'Tm':
                znum_int[i] = 69
            elif str(spe[i][0:2]) == 'Yb':
                znum_int[i] = 70
            elif str(spe[i][0:2]) == 'Lu':
                znum_int[i] = 71
            elif str(spe[i][0:2]) == 'Hf':
                znum_int[i] = 72
            elif str(spe[i][0:2]) == 'Ta':
                znum_int[i] = 73
            elif str(spe[i][0:2]) == 'W ':
                znum_int[i] = 74
            elif str(spe[i][0:2]) == 'Re':
                znum_int[i] = 75
            elif str(spe[i][0:2]) == 'Os':
                znum_int[i] = 76
            elif str(spe[i][0:2]) == 'Ir':
                znum_int[i] = 77
            elif str(spe[i][0:2]) == 'Pt':
                znum_int[i] = 78
            elif str(spe[i][0:2]) == 'Au':
                znum_int[i] = 79
            elif str(spe[i][0:2]) == 'Hg':
                znum_int[i] = 80
            elif str(spe[i][0:2]) == 'Tl':
                znum_int[i] = 81
            elif str(spe[i][0:2]) == 'Pb':
                znum_int[i] = 82
            elif str(spe[i][0:2]) == 'Bi':
                znum_int[i] = 83
            elif str(spe[i][0:2]) == 'Po':
                znum_int[i] = 84
            elif str(spe[i][0:2]) == 'At':
                znum_int[i] = 85

        # here the index to connect name and atomic numbers.
        global i_znum
        i_znum = {}    
        for a,b in zip(spe,znum_int):
            i_znum[a]=b
                        
        #connect name to column number
        global cl
        cl={}
        for a,b in zip(spe,n_array):
            cl[a] = b  

        # from here below I read the abundance.
        
        abunds = []
        name_specie_in_file=self.se.dcols[5]

        for i in range(len(spe)):
            #print spe[i],self.se.isotopes[i]
            abunds.append(self.se.get(cycle,name_specie_in_file,self.se.isotopes[i]))
            #print abunds[0][0],abunds[0][len(masses)-1]
            # abunds[i] is giving abundances for isotope i. I want mass_frac[i],
        # the isotopic distribution for mass zone i
        global used_masses
        used_masses = []
        global mass_frac
        mass_frac = []
        for i in range(len(masses)):
            if mass_range[0] <=  masses[i]  and mass_range[1] >=  masses[i] :
                used_masses.append(masses[i])
            	temp = []
            	for j in range (len(spe)):
                	temp.append(abunds[j][i])
            		mass_frac.append(temp)    
        #print mass_frac[len(masses)-1][cl['H   1']]        
        #print len(mass_frac)


    def decay(self):

        '''this module simply calculate abundances of isotopes after decay
        it requires before being used a call to read_iso_abund_marco and stable_species.
        see e.g., function interface'''

        global decayed_multi_d
        decayed_multi_d=[]
	#print len(mass_frac)
	#print len(decay_raw)
        for iii in range(len(mass_frac)):
            jj=-1
            decayed=[]
            for i in range(len(decay_raw)):
                 if jdum[i] > 0.5:
                    jj=jj+1
                    dummy=0.
                    for j in range(len(decay_raw[i])):             
                        try:    
                            dum_str = decay_raw[i][j]    
                            dummy = dummy + float(mass_frac[iii][cl[dum_str.lower().capitalize()]])
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
        #print spe,len(spe)
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
 
        The following keywords can also be used:

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

# this belongs into a superclass (check with Marco?)
def solar(filename_solar):
    ''' read solar abundances from filename_solar'''

    f0=open(filename_solar)
    sol=f0.readlines()
    f0.close 
    sol[0].split("         ")

    # Now read in the whole file and create a hashed array:
    global names_sol
    names_sol=[]
    global z_sol
    z_sol=[]    
    yps=np.zeros(len(sol))
    mass_number=np.zeros(len(sol))
    for i in range(len(sol)):
        z_sol.append(int(sol[i][1:3]))
        names_sol.extend([sol[i].split("         ")[0][4:]])
        yps[i]=float(sol[i].split("         ")[1])
        mass_number[i]=int(names_sol[i][2:5])
    #  convert 'h   1' in prot, not needed any more??
    #names_sol[0] = 'prot '
    
    
    # now zip them together:
    global solar_abundance
    solar_abundance={}
    for a,b in zip(names_sol,yps):
        solar_abundance[a] = b


# this belongs into a superclass (check with Marco?)
def stable_specie():

    ''' provide the list of stable species, and decay path feeding stables '''


    #from numpy import *
    #from define_input import *


    stable_raw=[]
    stable_raw = ['H   1', 'H    2',\
    'HE  3', 'HE  4',\
    'LI  6', 'LI  7',\
    'BE  9',\
    'B  10', 'B  11',\
    'C  12', 'C  13',\
    'N  14', 'N  15',\
    'O  16', 'O  17', 'O  18',\
    'F  19',\
    'NE 20', 'NE 21', 'NE 22',\
    'NA 23',\
    'MG 24', 'MG 25', 'MG 26',\
    'AL 27',\
    'SI 28', 'SI 29', 'SI 30',\
    'P  31',\
    'S  32', 'S  33', 'S  34', 'S  36',\
    'CL 35', 'CL 37',\
    'AR 36', 'AR 38', 'AR 40',\
    'K  39', 'K  40', 'K  41',\
    'CA 40', 'CA 42', 'CA 43', 'CA 44', 'CA 46', 'CA 48',\
    'SC 45',\
    'TI 46', 'TI 47', 'TI 48', 'TI 49', 'TI 50',\
    'V  50', 'V  51',\
    'CR 50', 'CR 52', 'CR 53', 'CR 54',\
    'MN 55',\
    'FE 54', 'FE 56', 'FE 57', 'FE 58',\
    'CO 59',\
    'NI 58', 'NI 60', 'NI 61', 'NI 62', 'NI 64',\
    'CU 63', 'CU 65',\
    'ZN 64', 'ZN 66', 'ZN 67', 'ZN 68', 'ZN 70',\
    'GA 69', 'GA 71',\
    'GE 70', 'GE 72', 'GE 73', 'GE 74', 'GE 76',\
    'AS 75',\
    'SE 74', 'SE 76', 'SE 77', 'SE 78', 'SE 80', 'SE 82',\
    'BR 79', 'BR 81',\
    'KR 78', 'KR 80', 'KR 82', 'KR 83', 'KR 84', 'KR 86',\
    'RB 85', 'RB 87',\
    'SR 84', 'SR 86', 'SR 87', 'SR 88',\
    'Y  89',\
    'ZR 90', 'ZR 91', 'ZR 92', 'ZR 94', 'ZR 96',\
    'NB 93',\
    'MO 92', 'MO 94', 'MO 95', 'MO 96', 'MO 97', 'MO 98', 'MO100',\
    'RU 96', 'RU 98', 'RU 99', 'RU100', 'RU101', 'RU102', 'RU104',\
    'RH103',\
    'PD102', 'PD104', 'PD105', 'PD106', 'PD108', 'PD110',\
    'AG107', 'AG109',\
    'CD106', 'CD108', 'CD110', 'CD111', 'CD112', 'CD113', 'CD114', 'CD116',\
    'IN113',  'IN115',\
    'SN112', 'SN114', 'SN115', 'SN116', 'SN117', 'SN118', 'SN119', 'SN120', 'SN122', 'SN124',\
    'SB121', 'SB123',\
    'TE120', 'TE122', 'TE123', 'TE124', 'TE125', 'TE126', 'TE128', 'TE130',\
    'I 127',\
    'XE124', 'XE126', 'XE128', 'XE129', 'XE130', 'XE131', 'XE132', 'XE134', 'XE136',\
    'CS133',\
    'BA130', 'BA132', 'BA134', 'BA135', 'BA136', 'BA137', 'BA138',\
    'LA138', 'LA139',\
    'CE136', 'CE138', 'CE140', 'CE142',\
    'PR141',\
    'ND142', 'ND143', 'ND144', 'ND145', 'ND146', 'ND148', 'ND150',\
    'SM144', 'SM147', 'SM148', 'SM149', 'SM150', 'SM152', 'SM154',\
    'EU151', 'EU153',\
    'GD152', 'GD154', 'GD155', 'GD156', 'GD157', 'GD158', 'GD160',\
    'TB159',\
    'DY156', 'DY158', 'DY160', 'DY161', 'DY162', 'DY163', 'DY164',\
    'HO165',\
    'ER162', 'ER164', 'ER166', 'ER167', 'ER168', 'ER170',\
    'TM169',\
    'YB168', 'YB170', 'YB171', 'YB172', 'YB173', 'YB174', 'YB176',\
    'LU175', 'LU176',\
    'HF174', 'HF176', 'HF177', 'HF178', 'HF179', 'HF180',\
    'TA180', 'TA181',\
    'W 180', 'W 182', 'W 183', 'W 184', 'W 186',\
    'RE185', 'RE187',\
    'OS184', 'OS186', 'OS187', 'OS188', 'OS189', 'OS190', 'OS192',\
    'IR191', 'IR193',\
    'PT190', 'PT192', 'PT194', 'PT195', 'PT196', 'PT198',\
    'AU197',\
    'HG196', 'HG198', 'HG199', 'HG200', 'HG201', 'HG202', 'HG204',\
    'TL203', 'TL205',\
    'PB204', 'PB206', 'PB207', 'PB208',\
    'BI209']

    jj=-1
    global count_size_stable
    count_size_stable=[]
    global stable
    stable=[]
    global jdum
    jdum=np.zeros(len(stable_raw))
    global jjdum
    jjdum=np.zeros(len(spe))
    for i in range(len(stable_raw)):
       dum_str = stable_raw[i]
       for j in range(len(spe)):
               if stable_raw[i].capitalize() == spe[j]:
                   stable.append(stable_raw[i]) 
                   jdum[i]=1
		   jjdum[j]=1
                   jj=jj+1
                   count_size_stable.append(int(jj))
    #print stable
    # back_ind is an index to go back, to use the order of stable
    # useful for example for decayed yields.
    global back_ind
    back_ind={}     
    for a,b in zip(stable,count_size_stable):
         back_ind[a]=b       
    #print 'in stable:',back_ind['SE 74']   
    # definition of decay paths
    global decay_raw
    decay_raw=[]
    decay_raw=[['H   1'],\
    ['H   2'],\
    ['HE  3'],\
    ['HE  4','B   8'],\
    ['LI  6'],\
    ['LI  7','BE  7'],\
    ['BE  9'],\
    ['B  10','BE 10'],\
    ['B  11','C  11','BE 11'],\
    ['C  12'],\
    ['C  13','N  13','O  13'],\
    ['N  14','C  14','O  14'],\
    ['N  15','C  15','O  15','F  15'],\
    ['O  16'],\
    ['O  17','F  17'],\
    ['O  18','F  18','NE 18'],\
    ['F  19','O  19','NE 19'],\
    ['NE 20','F  20','NA 20'],\
    ['NE 21','F  21','NA 21'],\
    ['NE 22','NA 22','MG 22','F  22'],\
    ['NA 23','MG 23','NE 23'],\
    ['MG 24','AL 24','NA 24','NE 24',],\
    ['MG 25','NA 25','AL 25'],\
    ['MG 26','SI 26','AL 26','NA 26','AL*26'],\
    ['AL 27','SI 27','MG 27'],\
    ['SI 28','AL 28','MG 28'],\
    ['SI 29','P  29','AL 29','MG 29'],\
    ['SI 30','S  30','P  30','AL 30','MG 30'],\
    ['P  31','S  31','SI 31'],\
    ['S  32','P  32','SI 32'],\
    ['S  33','CL 33','P  33','SI 33'],\
    ['S  34','CL 34','P  34'],\
    ['S  36','P  36'],\
    ['CL 35','S  35','AR 35'],\
    ['CL 37','S  37','AR 37'],\
    ['AR 36','CL 36'],\
    ['AR 38','CL 38'],\
    ['AR 40','CL 40'],\
    ['K  39','AR 39'],\
    ['K  40'],\
    ['K  41','AR 41','CA 41'],\
    ['CA 40'],\
    ['CA 42','K  42','AR 42'],\
    ['CA 43','K  43','AR 43'],\
    ['CA 44','K  44','AR 44'],\
    ['CA 46','K  46'],\
    ['CA 48','K  48'],\
    ['SC 45','CA 45','K  45'],\
    ['TI 46','SC 46'],\
    ['TI 47','SC 47','CA 47'],\
    ['TI 48','SC 48'],\
    ['TI 49','SC 49','CA 49'],\
    ['TI 50','SC 50'],\
    ['V  50'],\
    ['V  51','CR 51','TI 51'],\
    ['CR 50'],\
    ['CR 52','MN 52','V  52','TI 52'],\
    ['CR 53','MN 53','V  53'],\
    ['CR 54','MN 54','V  54'],\
    ['MN 55','FE 55','CR 55'],\
    ['FE 54'],\
    ['FE 56','NI 56','CO 56','MN 56','CR 56'],\
    ['FE 57','NI 57','CO 57','MN 57'],\
    ['FE 58','CO 58'],\
    ['CO 59','FE 59'],\
    ['NI 58'],\
    ['NI 60','CO 60'],\
    ['NI 61','CO 61','FE 61','CU 61'],\
    ['NI 62','CO 62','CU 62'],\
    ['NI 64','CU 64'],\
    ['CU 63','NI 63','ZN 63'],\
    ['CU 65','NI 65','ZN 65'],\
    ['ZN 64','CU 64'],\
    ['ZN 66','CU 66','NI 66','GA 66','GE 66'],\
    ['ZN 67','CU 67','GA 67'],\
    ['ZN 68','GA 68','GE 68','CU 68'],\
    ['ZN 70'],\
    ['GA 69','ZN 69','GE 69'],\
    ['GA 71','ZN 71','GE 71'],\
    ['GE 70','GA 70'],\
    ['GE 72','GA 72','ZN 72','AS 72'],\
    ['GE 73','GA 73','AS 73'],\
    ['GE 74','GA 74'],\
    ['GE 76'],\
    ['AS 75','GE 75','GA 75','SE 75'],\
    ['SE 74','BR 74','KR 74'],\
    ['SE 76','AS 76','BR 76'],\
    ['SE 77','AS 77','GE 77','BR 77'],\
    ['SE 78','AS 78','GE 78','BR 75'],\
    ['SE 80'],\
    ['SE 82'],\
    ['BR 79','SE 79','AS 79','KR 79'],\
    ['BR 81','SE 81','KR 81'],\
    ['KR 78','RB 78','SR 78'],\
    ['KR 80','BR 80'],\
    ['KR 82','BR 82','RB 82'],\
    ['KR 83','BR 83','RB 83'],\
    ['KR 84','BR 84','RB 84'],\
    ['KR 86'],\
    ['RB 85','KR 85','SR 85','KR*85'],\
    ['RB 87','KR 87'],\
    ['SR 84','Y  84','ZR 84'],\
    ['SR 86','RB 86'],\
    ['SR 87','Y  87'],\
    ['SR 88','RB 88','KR 88','Y  88','ZR 88'],\
    ['Y  89','SR 89','RB 89','KR 89','ZR 89'],\
    ['ZR 90','Y  90','SR 90','NB 90'],\
    ['ZR 91','Y  91','SR 91','NB 91'],\
    ['ZR 92','Y  92','SR 92','NB 92'],\
    ['ZR 94','Y  94'],\
    ['ZR 96'],\
    ['NB 93','ZR 93','MO 93'],\
    ['MO 92','TC 92','RU 92'],\
    ['MO 94','NB 94','TC 94','RU 94'],\
    ['MO 95','NB 95','ZR 95','Y  95'],\
    ['MO 96','NB 96','TC 96'],\
    ['MO 97','NB 97','ZR 97','TC 97'],\
    ['MO 98','NB 98'],\
    ['MO100'],\
    ['RU 96','RH 96','PD 96'],\
    ['RU 98','TC 98','RH 98','PD 98'],\
    ['RU 99','TC 99','MO 99'],\
    ['RU100','TC100','RH100'],\
    ['RU101','TC101','RH101','MO101'],\
    ['RU102','MO102','TC102'],\
    ['RU104','TC104'],\
    ['RH103','RU103','TC103','PD103'],\
    ['PD102','AG102','CD102'],\
    ['PD104','RH104'],\
    ['PD105','RH105','RU105','AG105','AG105'],\
    ['PD106','RH106','RU106','AG106'],\
    ['PD108'],\
    ['PD110'],\
    ['AG107','PD107','CD107'],\
    ['AG109','PD109','CD109'],\
    ['CD106','IN106','SN106'],\
    ['CD108','AG108','IN108','SN108'],\
    ['CD110','AG110'],\
    ['CD111','AG111','PD111','IN111'],\
    ['CD112','AG112','PD112'],\
    ['CD113','AG113'],\
    ['CD114'],\
    ['CD116'],\
    ['IN113','SN113','SB113','TE113'],\
    ['IN115','CD115'],\
    ['SN112','IN112','SB112','TE112'],\
    ['SN114','IN114','SB114','TE114'],\
    ['SN115','SB115','TE115','I 115'],\
    ['SN116','IN116'],\
    ['SN117','IN117','SB117','CD117'],\
    ['SN118','SB118'],\
    ['SN119','SB119'],\
    ['SN120','SB120'],\
    ['SN122'],\
    ['SN124'],\
    ['SB121','SN121','TE121'],\
    ['SB123','SN123'],\
    ['TE120','I 120','XE120'],\
    ['TE122','SB122'],\
    ['TE123','I 123'],\
    ['TE124','SB124','I 124'],\
    ['TE125','SB125','I 125','SN125'],\
    ['TE126','SB126','SN126'],\
    ['TE128','SB128','SN128','I 128'],\
    ['TE130'],\
    ['I 127','TE127','XE127'],\
    ['XE124','CS124','BA124'],\
    ['XE126','CS126','BA126'],\
    ['XE128','I 128'],\
    ['XE129','I 129','TE129','SB129'],\
    ['XE130','I 130'],\
    ['XE131','I 131','TE131','CS131'],\
    ['XE132','I 132','TE132','CS132'],\
    ['XE134','I 134'],\
    ['XE136'],\
    ['CS133','XE133','I 133','BA133'],\
    ['BA130','LA130','CE130','PR130'],\
    ['BA132','LA132','CE132','PR132'],\
    ['BA134','CS134'],\
    ['BA135','CS135','XE135','I 135'],\
    ['BA136','CS136'],\
    ['BA137','CS137','XE137','LA137'],\
    ['BA138','CS138','XE138'],\
    ['LA138'],\
    ['LA139','BA139','CS139','CE139'],\
    ['CE136','PR136','ND136'],\
    ['CE138','PR138','ND138','PM138','SM138'],\
    ['CE140','LA140','BA140'],\
    ['CE142','LA142','BA142'],\
    ['PR141','CE141','LA141','BA141'],\
    ['ND142','PR142'],\
    ['ND143','PR143','CE143','LA143','PM143'],\
    ['ND144','PR144','CE144','LA144','PM144'],\
    ['ND145','PR145','CE145','PM145'],\
    ['ND146','PR146','CE146'],\
    ['ND148'],\
    ['ND150'],\
    ['SM144','EU144','GD144'],\
    ['SM147','PM147','ND147'],\
    ['SM148','PM148'],\
    ['SM149','PM149','ND149'],\
    ['SM150','PM150'],\
    ['SM152','PM152','ND152','EU152'],\
    ['SM154','PM154'],\
    ['EU151','SM151','PM151','ND151'],\
    ['EU153','SM153','PM153','GD153'],\
    ['GD152','EU152','TB152','DY152','HO152'],\
    ['GD154','EU154'],\
    ['GD155','EU155','SM155','TB155'],\
    ['GD156','EU156','SM156'],\
    ['GD157','EU157','SM157','TB157'],\
    ['GD158','EU158','SM158'],\
    ['GD160'],\
    ['TB159','GD159','EU159'],\
    ['DY156','HO156','ER156','TM156','YB156'],\
    ['DY158','HO158','ER158','TM158','YB158'],\
    ['DY160','TB160'],\
    ['DY161','TB161','GD161'],\
    ['DY162','TB162','GD162'],\
    ['DY163','TB163','HO163'],\
    ['DY164','TB164'],\
    ['HO165','DY165','ER165'],\
    ['ER162','TM162','YB162','LU162'],\
    ['ER164','TM164','YB164','LU164'],\
    ['ER166','HO166','DY166'],\
    ['ER167','HO167','DY167'],\
    ['ER168','HO168','DY168','TB161','TM168'],\
    ['ER170'],\
    ['TM169','ER169','HO169','YB169'],\
    ['YB168','LU168','HF168','TA168','W 168'],\
    ['YB170','TM170'],\
    ['YB171','TM171','ER171','LU171'],\
    ['YB172','TM172','ER172','LU172'],\
    ['YB173','TM173','ER173','LU173'],\
    ['YB174','TM174','ER174','LU174'],\
    ['YB176'],\
    ['LU175','YB175','TM175','HF175'],\
    ['LU176'],\
    ['HF174','TA174','W 174','RE174','OS174'],\
    ['HF176'],\
    ['HF177','LU177','YB177'],\
    ['HF178','LU178','YB178'],\
    ['HF179','LU179','TA179'],\
    ['HF180','LU180'],\
    ['TA180'],\
    ['TA181','HF181','W 181'],\
    ['W 180','RE180','OS180','IR180','PT180'],\
    ['W 182','TA182','HF182'],\
    ['W 183','TA183','HF183','RE183'],\
    ['W 184','TA184','HF184','RE184'],\
    ['W 186','TA186'],\
    ['RE185','W 185','TA185'],\
    ['RE187','W 187'],\
    ['OS184','IR184','PT184','AU184'],\
    ['OS186','RE186'],\
    ['OS187'],\
    ['OS188','RE188','W 188','IR188'],\
    ['OS189','RE189','W 189','IR189'],\
    ['OS190','RE190','W 190','IR190'],\
    ['OS192'],\
    ['IR191','OS191','RE191','PT191'],\
    ['IR193','OS193','PT193'],\
    ['PT190','AU190','HG190','TL190','PB190'],\
    ['PT192','IR192'],\
    ['PT194','IR194','OS194'],\
    ['PT195','IR195','OS195','AU195'],\
    ['PT196','IR196','AU196'],\
    ['PT198'],\
    ['AU197','PT197','HG197'],\
    ['HG196','TL196','PB196'],\
    ['HG198','AU198'],\
    ['HG199','AU199','PT199'],\
    ['HG200','AU200','PT200'],\
    ['HG201','AU201','PT201'],\
    ['HG202','AU202','PT202'],\
    ['HG204'],\
    ['TL203','HG203','PB203'],\
    ['TL205','HG205','PB205'],\
    ['PB204','TL204'],\
    ['PB206','TL206','HG206','PO210'],\
    ['PB207','TL207','HG207','BI207'],\
    ['PB208','TL208','HG208','BI208'],\
    ['BI209','PB209','TL209']]
    #print decay_raw

# is this redundant with  def iso_abund(self, mass_range, cycle, stable) ???
def average_iso_abund_marco(directory,name_h5_file,mass_range,cycle,stable,i_decay):
    ''' Interface to average over mass_range. 
    directory     -  location of h5 file to plot. Needed for plot_tools
    name_h5_file  -  name of h5 file. Needed for plot_tools
    mass_range    - required to plot data in a certain mass range. Needed for read_iso_abund_marco
    cycle         - which cycle from the h5 file?. Needed for read_iso_abund_marco
    stable        - logic if want to plot only stable or not.  
    i_decay       - if = 1 I plot not decayed, if = 2 I plot decayed. Make sense only if stable is true.'''


    #import nuh5p 
    import mppnp as mp 	    

    if not stable and i_decay == 2:
        print 'ERROR: choose i_decay = 1'  
        return
    
    data=mp.se(directory,name_h5_file)
    data.read_iso_abund_marco(mass_range,cycle)
    #print spe
    if i_decay == 2:
        mp.stable_specie()
        data.decay()


    
    # here I am calculating average mass fraction for all isotopes in given mass range, and then
    # if needed calculating average over decayed.
    # warning: mass_range is bigger than used_masses range, by definition. Should I use it?
    print 'average over used_masses range, not over original mass_range'
    print used_masses[0],used_masses[len(used_masses)-1],'instead of',mass_range[0],mass_range[1]
    
    global average_mass_frac
    average_mass_frac = []
    
    if len(used_masses) >= 2:
        dm_tot = abs(used_masses[len(used_masses)-1]-used_masses[0])
        for j in range(len(spe)):
            temp = 0.
            for i in range(len(used_masses)-1):
	        dm_i = abs(used_masses[i+1]-used_masses[i])
            	temp = float(mass_frac[i][j]*dm_i/dm_tot) + temp
            	#temp = float(mass_frac[i][j]*abs(used_masses[i+1]-used_masses[i])) + temp
            	#temp = temp/abs(used_masses[len(used_masses)-1]-used_masses[0])    
            average_mass_frac.append(temp)
        #print average_mass_frac
    elif  len(used_masses) == 1:
        print 'case with 1 mass zone only, not implemented yet'
    
    
        
    somma = 0.
    for i in range(len(spe)):
        somma = float(average_mass_frac[i]) + somma
    print 'departure from 1 of sum of average_mass_frac=',abs(1. - somma)
    
    # not let's do it over decayed also, if i_decay = 2
    if i_decay == 2:
        global average_mass_frac_decay
        average_mass_frac_decay = []
	dm_tot = abs(used_masses[len(used_masses)-1]-used_masses[0])
	#
	#print len(decayed_multi_d[0]),decayed_multi_d[0]        
	for j in range(len(back_ind)):
        	temp = 0.
           	for i in range(len(used_masses)-1):
		        dm_i = abs(used_masses[i+1]-used_masses[i])
	            	temp = float(decayed_multi_d[i][j]*dm_i/dm_tot) + temp
           		#temp = float(decayed_multi_d[i][j]*abs(used_masses[i+1]-used_masses[i])) + temp
            		#temp = temp/abs(used_masses[len(used_masses)-1]-used_masses[0])    
            	average_mass_frac_decay.append(temp)

        somma = 0.
        for i in range(len(back_ind)):
            somma = float(average_mass_frac[i]) + somma
        print 'departure from 1 of sum of average_mass_frac_decay=',abs(1. - somma)

    # now I have the average abundances. We can do the plot.
    

def plot_iso_abund_marco(directory,name_h5_file,mass_range,cycle,logic_stable,i_decay,file_solar,solar_factor):
    ''' Interface to plot average over mass_range. 
    directory     -  location of h5 file to plot. Needed for plot_tools
    name_h5_file  -  name of h5 file. Needed for plot_tools
    mass_range    - required to plot data in a certain mass range. Needed for read_iso_abund_marco
    cycle         - which cycle from the h5 file?. Needed for read_iso_abund_marco
    logic_stable  - logic if want to plot only stable or not.  
    i_decay       - if = 1 I plot not decayed, if = 2 I plot decayed. Make sense only if stable is true
    file_solar        - file where to take solar abundances
    solar_factor      - float to correct initial abundances to solar, e.g. for Z=0.01 and AG89 solar_factor = 2.'''


    # solar abundances are read here
    solar(file_solar)
    # from here I have average abundances in mass_range to plot
    average_iso_abund_marco(directory,name_h5_file,mass_range,cycle,logic_stable,i_decay)
    
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
        for i in range(len(spe)):
            pl.plot(amass_int[cl[spe[i]]],average_mass_frac[cl[spe[i]]],'ko')

        pl.xlabel('$Mass$ $number$', fontsize=20)
        pl.ylabel('$X_{i}$', fontsize=20)

        pl.ylim(1.0e-10,10.)
        pl.xlim(55,110)
    
    elif logic_stable:
    # plot stable
        for i in range(len(stable)):
            pl.plot(amass_int[cl[stable[i].capitalize()]],average_mass_frac[cl[stable[i].capitalize()]]/solar_abundance[stable[i].lower()]/solar_factor,'ko')
        
        if i_decay == 2:
            for j in range(len(stable)):
                    #print cl[stable[j].capitalize()],stable[j].capitalize(),amass_int[cl[stable[j].capitalize()]]
                    pl.plot(amass_int[cl[stable[j].capitalize()]],average_mass_frac_decay[back_ind[stable[j]]]/solar_abundance[stable[j].lower()]/solar_factor,'Dg')
    
        for i in range(len(stable)):
            for j in range(len(stable)): 
                if stable[i][:2] == stable[j][:2]:
                    if stable[i] == stable[j-1]:
                        adum  =[amass_int[cl[stable[i].capitalize()]],amass_int[cl[stable[j].capitalize()]]]
                        mfdum =[float(average_mass_frac[cl[stable[i].capitalize()]])/float(solar_abundance[stable[i].lower()]*solar_factor),float(average_mass_frac[cl[stable[j].capitalize()]])/float(solar_abundance[stable[j].lower()]*solar_factor)]
                        mfddum=[float(average_mass_frac_decay[back_ind[stable[i]]])/float(solar_abundance[stable[i].lower()]*solar_factor),float(average_mass_frac_decay[back_ind[stable[j]]])/float(solar_abundance[stable[j].lower()]*solar_factor)]
                        pl.plot(adum,mfdum,'k-')
                    if i_decay == 2:
                              pl.plot(adum,mfddum,'g-')  

    pl.xlabel('$Mass$ $number$', fontsize=20)
    pl.ylabel('$X_{i}/X_{sun}$', fontsize=20)

    pl.ylim(1.0e-2,1000.)
    pl.xlim(55,210)
    
    
    
    pl.grid()
    pl.show() 
      


def element_abund_marco(i_decay,solar_factor):
    ''' Here elements abundances, solar element abundances and production factors for elements are calculated'''


    # this way is done in a really simple way. May be done better for sure, in a couple of loops.
    # I keep this, since I have only to copy over old script. Falk will probably redo it.
    

    global z_bismuth
    z_bismuth = 83

    global z_for_elem
    z_for_elem = []
    global index_stable
    index_stable = []
    global elem_abund
    elem_abund = np.zeros(z_bismuth)
    global elem_abund_decayed
    elem_abund_decayed = np.zeros(z_bismuth)
    global solar_elem_abund
    solar_elem_abund = np.zeros(z_bismuth)
    global elem_prod_fac
    elem_prod_fac = np.zeros(z_bismuth)
    global elem_prod_fac_decayed
    elem_prod_fac_decayed = np.zeros(z_bismuth)
    
    i_for_stable = 1
    i_for_unstable = 0
    for i in range(z_bismuth):
        z_for_elem.append(int(i+1))
    	# the only elements below bismuth with no stable isotopes are Tc and Pm
    	if i+1 == 43 or i+1 == 61:
        	index_stable.append(i_for_unstable) 
    	else:
        	index_stable.append(i_for_stable)
        
    # notice that elem_abund include all contribution, both from stables and unstables in
    # that moment.
    for i in range(z_bismuth):
        dummy = 0.
        for j in range(len(spe)):
            if znum_int[j] == i+1 and jjdum[j] > 0.5:
                dummy = dummy + float(average_mass_frac[j])
    	elem_abund[i] = dummy

    for i in range(z_bismuth):
        dummy = 0.
        for j in range(len(solar_abundance)):
            if z_sol[j] == i+1:
                dummy = dummy + float(solar_abundance[names_sol[j]])
    	solar_elem_abund[i] = dummy


    for i in range(z_bismuth):
        if index_stable[i] == 1:
            elem_prod_fac[i] = float(elem_abund[i]/solar_elem_abund[i]/solar_factor)
        elif index_stable[i] == 0:
            elem_prod_fac[i] = 0.    


    if i_decay == 2:
        for i in range(z_bismuth):
            dummy = 0.
            for j in range(len(average_mass_frac_decay)):
                if znum_int[cl[stable[j].capitalize()]] == i+1:
                    #print znum_int[cl[stable[j].capitalize()]],cl[stable[j].capitalize()],stable[j]
                    dummy = dummy + float(average_mass_frac_decay[j])
       	    elem_abund_decayed[i] = dummy


        for i in range(z_bismuth):
            if index_stable[i] == 1:
                elem_prod_fac_decayed[i] = float(elem_abund_decayed[i]/solar_elem_abund[i]/solar_factor)
            elif index_stable[i] == 0:
                elem_prod_fac_decayed[i] = 0.    


                         

def plot_el_abund_marco(directory,name_h5_file,mass_range,cycle,logic_stable,i_decay,file_solar,solar_factor,symbol='ko'):
    ''' Interface to plot elements abundances averaged over mass_range. 
    directory     -  location of h5 file to plot. Needed for plot_tools
    name_h5_file  -  name of h5 file. Needed for plot_tools
    mass_range    - required to plot data in a certain mass range. Needed for read_iso_abund_marco
    cycle         - which cycle from the h5 file?. Needed for read_iso_abund_marco
    logic_stable  - logic if want to plot only stable or not.  
    i_decay       - if = 1 I plot not decayed, if = 2 I plot decayed. Make sense only if stable is true
    file_solar        - file where to take solar abundances
    solar_factor      - float to correct initial abundances to solar, e.g. for Z=0.01 and AG89 solar_factor = 2.'''

    # solar abundances are read here
    solar(file_solar)
    # from here I have average abundances in mass_range to plot
    average_iso_abund_marco(directory,name_h5_file,mass_range,cycle,logic_stable,i_decay)
    # element abundances are calculated here
    element_abund_marco(i_decay,solar_factor)
    
    
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
        for i in range(z_bismuth):
            pl.plot(z_for_elem[i],elem_prod_fac[i],symbol,markersize=10.)

        pl.xlabel('$Atomic$ $number$', fontsize=20)
        pl.ylabel('$X_{i}/X_{sun}$', fontsize=20)

        pl.ylim(1.0e-2,1000.)
        pl.xlim(0,95)
    
    elif logic_stable:
        for i in range(z_bismuth):
            if index_stable[i] == 1:
		continue
                #pl.plot(z_for_elem[i],elem_prod_fac[i],'ko')
        if i_decay == 2:
            for i in range(z_bismuth):
                if index_stable[i] == 1:
                    pl.plot(z_for_elem[i],elem_prod_fac_decayed[i],symbol,markersize=10.)

        pl.xlabel('$Atomic$ $number$', fontsize=20)
        pl.ylabel('$X_{i}/X_{sun}$', fontsize=20)

        pl.ylim(1.0e-2,1000.)
        pl.xlim(0,95)
    
    
    
    pl.grid()
    pl.show() 
    
    
