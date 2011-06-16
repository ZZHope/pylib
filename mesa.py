''' MESA output data loading and plotting

    Falk Herwig for the MESA collaboration (v0.1, 23JUN2010)

    mesa.py provides tools to get MESA stellar evolution data output
    into your favourite python session. In the LOGS directory MESA
    outputs two types of files: star.log is a time evolution output,
    printing one line per so many cycles (e.g. each cycle) of all
    sorts of things. lognnn.data files are profile data files. nnn is
    the number of log.data files that is translated into model cycles
    in the profiles.index file.

    MESA allows users to freely define what should go into these two
    types of outputs, which means that column numbers can and do
    change. mesa.py reads in both types of files and present them (as
    well as any header attributes) as arrays that can be referenced by
    the actual column name as defined in the header section of the
    files. mesa.py then defines a (hopefully growing) set of standard
    plots that make use of the data just obtained.

    mesa.py is organised as a module that can be imported into any
    python or ipython session. It is related to nmuh5.py which is a
    similar module to deal with 'se' output, used by the NuGrid
    collaboration. mesa.py does not need se libraries. The 'se' output
    files that can be written with MESA can be read and processed with
    the nuh5.py tool.

    mesa.py is providing two class objects, profile and star_log. The
    first makes profile data available, the second reads and plots the
    star.log file. Note that several instances of these can be
    initiated within one session and data from different instances
    (i.e. models, tracks etc) can be overplotted.

    Here is how a simple session could look like that is plotting an
    HRD (I prefer to load ipython with matplotlib and numpy support
    via the alias
    alias mpython='ipython -pylab -p numpy -editor emacsclient')

        vortex$ mpython
        Python 2.5.2 (r252:60911, Feb 22 2008, 07:57:53) 
        Type "copyright", "credits" or "license" for more information.

        IPython 0.9.1 -- An enhanced Interactive Python.
        ?         -> Introduction and overview of IPython's features.
        %quickref -> Quick reference.
        help      -> Python's own help system.
        object?   -> Details about 'object'. ?object also works, ?? prints more.

        IPython profile: numpy

          Welcome to pylab, a matplotlib-based Python environment.
          For more information, type 'help(pylab)'.

        In [1]: import mesa as ms

        In [2]: help ms
        ------> help(ms)

        In [4]: s=ms.star_log('.')

        In [5]: s.hrd()
    
     In order to find out what header attributes and columns are
     available in star.log use:
     
        In [6]: s.header_attr
        Out[6]: 
        {'burn_min1': 50.0,
         'burn_min2': 1000.0,
         'c12_boundary_limit': 0.0001,
         'h1_boundary_limit': 0.0001,
         'he4_boundary_limit': 0.0001,
         'initial_mass': 2.0,
         'initial_z': 0.01}

        In [7]: s.cols
        Out[7]: 
        {'center_c12': 38,
         'center_h1': 36,
         'center_he4': 37,
          ...  
          
    In order to read the profile data from the first log.data file in
    profiles.index, and then get the mass and temperature out and
    finally plot them try:

        In [9]: a1=ms.mesa_profile('LOGS',1)
        100 in profiles.index file..
        The 1. log.data file is 44
        reading ./log44.data ...

        In [10]: T=a1.get('temperature')

        In [11]: mass=a1.get('mmid')

        In [12]: plot(mass,T)
        Out[12]: [<matplotlib.lines.Line2D object at 0x8456ed0>]

    Or, you could have had it easier in the following way:
        In [13]: a1.plot('mass','na23',logY=True)
    where the superclass plot method interprets data column headers
    correctly and does all the work for you.

    Of course, a1.cols etc are available here as well and many other
    things. E.g. a.model contains an array with all the models for
    which log.data are available. You may initiate a profile object
    with a model number:

        In [14]: a2=ms.mesa_profile('.',55000,num_type='model')
        100 in profiles.index file ...
        reading ./log87.data ...

'''
from data_plot import *
import numpy as np
import matplotlib.pylab as pyl
import matplotlib.pyplot as pl
import os



class mesa_profile(DataPlot):
    ''' read profiles.index and prepare reading MESA profile files

    starts with reading profiles.index and creates hash array
    log.data can then be accessed via prof_plot
    
    '''

    sldir = ''
    
    def __init__(self,sldir,num,num_type='log_i',prof_ind_name='profiles.index',log_prefix='log',data_suffix='.data'):
        '''read a log.data profile file

        input:
        sldir       directory path of LOGS

        num         by default this is the i. log.data available
                    (e.g. num=1 is the 1. available profile file),
                    however if you give 
        num_type    as 'log_num' then num will be interpreted as the
                    log.data number log_num (log_num is the number
                    that appears in the file names of type
                    log23.data), or try 'model' to get the prfile
                    log.data file for model (or cycle number) used by
                    the stellar evolution code
        prof_ind_name    use this optional argument if the profiles.index 
                         file hasn an alternative name, for example, do 
                         superpro=ms.profile('LOGS',1,prof_ind_name='super.prof') 
        log_prefix, data_suffix are optional arguments that allow you to change
                    the defaults for the log.data profile files. '''

        self.prof_ind_name = prof_ind_name
        self.sldir         = sldir

        if num_type is 'model':
            self.profiles_index()
            try:
                log_num=self.log_ind[num]
            except KeyError:
                print 'There is no log.data file for this model'
                return
        elif num_type is 'log_i':
            log_num=self.log_file_ind(num)
            if log_num == -1:
                print "Could not find a log.data file with that number"
                return
        elif num_type is 'log_num':
            log_num = num
        else:
            print 'unknown num_type'
            return

        filename=self.sldir+'/'+log_prefix+str(log_num)+data_suffix
        
        print 'reading '+filename+' ...'
        header_attr = read_mesafile(filename,only='header_attr')
        num_zones=int(header_attr['num_zones'])
        header_attr,cols,data = read_mesafile(filename,data_rows=num_zones,only='all')

        self.cols        = cols
        self.header_attr = header_attr
        self.data        = data


    def __del__(self):
        print 'Closing profile tool ...'

    def profiles_index(self):
        ''' read profiles.index and make hash array

        log_ind     hash array that returns log.data file number from model number
        model       the models for which log.data is available'''

        prof_ind_name = self.prof_ind_name 

        f = open(self.sldir+'/'+prof_ind_name,'r')
        line = f.readline()
        numlines=int(line.split()[0])
        print str(numlines)+' in profiles.index file ...'

        model=[]
        log_file_num=[]
        for line in f:
            model.append(int(line.split()[0]))
            log_file_num.append(int(line.split()[2]))

        log_ind={}    # log.data number from model
        for a,b in zip(model,log_file_num):
            log_ind[a] = b
            
        self.log_ind=log_ind
        self.model=model

# let's start with functions that aquire data

    def log_file_ind(self,inum):
        ''' information about available log.data files
        
        inmu       attempt to get number of inum's log.data file
        inum_max   max number of log.data files available'''
        
        self.profiles_index()
        if inum <= 0:
            print "Smallest argument is 1"
            return

        inum_max = len(self.log_ind)
        inum -= 1
        
        if inum > inum_max:
            print 'There are only '+str(inum_max)+' log.data file available.'
            log_data_number = -1
            return log_data_number
        else:
            log_data_number=self.log_ind[self.model[inum]]
            print 'The '+str(inum+1)+'. log.data file is '+ \
                  str(log_data_number)
            return log_data_number

    def get(self,str_name):
        ''' return a column of data with the name str_name
        
        str_name is the name of the column as printed in the
        lognnn.data file; get the available columns from self.cols
        (where you replace self with the name of your instance)'''

        column_array = self.data[:,self.cols[str_name]-1].astype('float')
        return column_array



        
class star_log(DataPlot):
    ''' read star.log MESA output and plot various things, including
    HRD, Kippenhahn etc
    
    sldir              - which LOGS directory
    slname='star.log'  - optional argument if star.log file has alternative name,
    clean_starlog=True - request new cleaning of star.log, makes star.logsa which 
                         is the file that is actually read and plotted
    use like this: another=ms.star_log('LOGS',slname='anothername')
    '''

    sldir  = ''
    slname = ''
    header_attr = []
    cols = [] 
    
    def __init__(self,sldir,slname='star.log',clean_starlog=True):
        self.sldir  = sldir
        self.slname = slname
        self.clean_starlog  = clean_starlog

        if not os.path.exists(sldir+'/'+slname):
            print 'error: no star.log file found in '+sldir
        else:
            self.read_starlog()

    def __del__(self):
        print 'Closing star_log tool ...'

# let's start with functions that aquire data
    def read_starlog(self):
        ''' read star.log file again'''

        sldir   = self.sldir
        slname  = self.slname
        slaname = slname+'sa'
        if self.clean_starlog and os.path.exists(sldir+'/'+slaname):
            os.remove(sldir+'/'+slaname)
            
        if not os.path.exists(sldir+'/'+slaname):
            print 'No star.logsa file found, create new one from star.log.'
            cleanstarlog(sldir+'/'+slname)
        else:
            print 'Using old star.logsa file ...'
            
        cmd=os.popen('wc '+sldir+'/'+slaname)    
        cmd_out=cmd.readline()
        cnum_cycles=cmd_out.split()[0]
        num_cycles=int(cnum_cycles) - 6

        filename=sldir+'/'+slaname

        header_attr,cols,data = read_mesafile(filename,data_rows=num_cycles)

        self.cols        = cols
        self.header_attr = header_attr
        self.data        = data
        
    def get(self,str_name):
        ''' return a column of data with the name str_name
        
        str_name is the name of the column as printed in star.log
        get the available columns from self.cols (where you replace
        self with the name of your instance'''

        column_array = self.data[:,self.cols[str_name]-1].astype('float')
        return column_array
        
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
    
    def kippenhahn_CO(self,num_frame,xax,t0_model=0,title='Kippenhahn diagram',\
                       tp_agb=0.):
		''' Kippenhahn plot as a function of time or model with CO ratio
		
		num_frame    number of frame to plot this plot into 
                xax          string that is either model or time to
                             indicate what is to be used on the x-axis

                t0_model     model for the zero point in time, for AGB
                             plots this would be usually the model of
                             the 1st TP, which can be found with the
                             Kippenhahn plot 
                title        figure title

                tp_agb       if >= 0 then 
                             ylim=[h1_min*1.-tp_agb/100 : h1_max*1.+tp_agb/100] 
                             with h1_min, h1_max the min and max H-free 
                             core mass coordinate
                '''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction.'+\
			  ' needs to be "time" or "model"'
		
                t0_mod=xaxisarray[t0_model]
	    
		h1_boundary_mass  = self.get('h1_boundary_mass')
		he4_boundary_mass = self.get('he4_boundary_mass')
		star_mass         = self.get('star_mass')
		mx1_bot           = self.get('mx1_bot')*star_mass
		mx1_top           = self.get('mx1_top')*star_mass
		mx2_bot           = self.get('mx2_bot')*star_mass
		mx2_top           = self.get('mx2_top')*star_mass
		surface_c12       = self.get('surface_c12')
		surface_o16       = self.get('surface_o16')
	
		COratio=(surface_c12*4.)/(surface_o16*3.)
	
		pyl.plot(xaxisarray[t0_model:]-t0_mod,COratio[t0_model:],'-k',label='CO ratio')
		pyl.ylabel('C/O ratio')
		pyl.legend(loc=4)

		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')
	
		pyl.twinx()
		pyl.plot(xaxisarray[t0_model:]-t0_mod,h1_boundary_mass[t0_model:],label='h1_boundary_mass')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,he4_boundary_mass[t0_model:],label='he4_boundary_mass')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_bot[t0_model:],',r',label='conv bound')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_top[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_bot[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_top[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,star_mass[t0_model:],label='star_mass')
		pyl.ylabel('mass coordinate')
		pyl.legend(loc=2)
                if tp_agb > 0.:
                    h1_min = min(h1_boundary_mass[t0_model:])
                    h1_max = max(h1_boundary_mass[t0_model:])
                    h1_min = h1_min*(1.-tp_agb/100.)
                    h1_max = h1_max*(1.+tp_agb/100.)
                    print 'setting ylim to zoom in on H-burning:',h1_min,h1_max                    
                    pyl.ylim(h1_min,h1_max)

    def kippenhahn(self,num_frame,xax,t0_model=0,title='Kippenhahn diagram',\
                       tp_agb=0.):
		''' Kippenhahn plot as a function of time or model
		
		num_frame    number of frame to plot this plot into 
                xax          string that is either model or time to
                             indicate what is to be used on the x-axis

                t0_model     model for the zero point in time, for AGB
                             plots this would be usually the model of
                             the 1st TP, which can be found with the
                             Kippenhahn plot 
                title        figure title

                tp_agb       if >= 0 then 
                             ylim=[h1_min*1.-tp_agb/100 : h1_max*1.+tp_agb/100] 
                             with h1_min, h1_max the min and max H-free 
                             core mass coordinate
                '''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 'kippenhahn_error: invalid string for x-axis selction.'+\
			  ' needs to be "time" or "model"'
		
                t0_mod=xaxisarray[t0_model]
	    
		h1_boundary_mass  = self.get('h1_boundary_mass')
		he4_boundary_mass = self.get('he4_boundary_mass')
		star_mass         = self.get('star_mass')
		mx1_bot           = self.get('mx1_bot')*star_mass
		mx1_top           = self.get('mx1_top')*star_mass
		mx2_bot           = self.get('mx2_bot')*star_mass
		mx2_top           = self.get('mx2_top')*star_mass
	

		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')
	
		pyl.plot(xaxisarray[t0_model:]-t0_mod,h1_boundary_mass[t0_model:],label='h1_boundary_mass')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,he4_boundary_mass[t0_model:],label='he4_boundary_mass')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_bot[t0_model:],',r',label='conv bound')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx1_top[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_bot[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,mx2_top[t0_model:],',r')
		pyl.plot(xaxisarray[t0_model:]-t0_mod,star_mass[t0_model:],label='star_mass')
		pyl.ylabel('mass coordinate')
		pyl.legend(loc=2)
                if tp_agb > 0.:
                    h1_min = min(h1_boundary_mass[t0_model:])
                    h1_max = max(h1_boundary_mass[t0_model:])
                    h1_min = h1_min*(1.-tp_agb/100.)
                    h1_max = h1_max*(1.+tp_agb/100.)
                    print 'setting ylim to zoom in on H-burning:',h1_min,h1_max                    
                    pyl.ylim(h1_min,h1_max)

    def t_surfabu(self,num_frame,xax,t0_model=0,title='surface abundance'):
		''' t_surfabu plots surface abundance evolution as a function of time
		
		num_frame    number of frame to plot this plot into 
                xax          string that is either model or time to
                             indicate what is to be used on the x-axis

                t0_model     model for the zero point in time, for AGB
                             plots this would be usually the model of
                             the 1st TP, which can be found with the
                             Kippenhahn plot 
                title        figure title
                '''
	
		pyl.figure(num_frame)
		
		if xax == 'time':
		    xaxisarray = self.get('star_age')
		elif xax == 'model':
		    xaxisarray = self.get('model_number')
		else:
		    print 't-surfabu error: invalid string for x-axis selction.'+\
			  ' needs to be "time" or "model"'
	    
		star_mass         = self.get('star_mass')
		surface_c12       = self.get('surface_c12')
		surface_c13       = self.get('surface_c13')
		surface_n14       = self.get('surface_n14')
		surface_o16       = self.get('surface_o16')                
                
                target_n14 = -3.5

	
		COratio=(surface_c12*4.)/(surface_o16*3.)
                t0_mod=xaxisarray[t0_model]
                log10_c12=np.log10(surface_c12[t0_model:])

                eps=1.e-3
                logxaxisarray=log10(max(xaxisarray[t0_model:])+eps-xaxisarray[t0_model:])

		pyl.plot(logxaxisarray,log10_c12,\
                             label='$^{12}\mathrm{C}$')
		pyl.plot(logxaxisarray,np.log10(surface_c13[t0_model:]),\
                             label='$^{13}\mathrm{C}$')
		pyl.plot(logxaxisarray,np.log10(surface_n14[t0_model:]),\
                             label='$^{14}\mathrm{N}$')
		pyl.plot(logxaxisarray,np.log10(surface_o16[t0_model:]),\
                             label='$^{16}\mathrm{O}$')
#                pyl.plot([min(xaxisarray[t0_model:]-t0_mod),max(xaxisarray[t0_model:]-t0_mod)],[target_n14,target_n14])

		pyl.ylabel('mass fraction $\log X$')
		pyl.legend(loc=2)

		if xax == 'time':
		    pyl.xlabel('t / yrs')
		elif xax == 'model':
		    pyl.xlabel('model number')
	

		pyl.twinx()
		pyl.plot(xaxisarray[t0_model:]-t0_mod,COratio[t0_model:],'-k',label='CO ratio')
		pyl.ylabel('C/O ratio')
		pyl.legend(loc=4)
                pyl.title(title)

# ... end t_surfabu

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
# below are some utilities that the user typically never calls directly

def read_mesafile(filename,data_rows=0,only='all'):
    ''' private routine that is not directly called by the user
    '''
    f=open(filename,'r')
    vv=[]
    v=[]
    lines = [] 
    line  = ''
    for i in range(0,6):
        line = f.readline()
        lines.extend([line])
    
    hval  = lines[2].split()
    hlist = lines[1].split()
    header_attr = {}
    for a,b in zip(hlist,hval):
        header_attr[a] = float(b)  
    if only is 'header_attr':
        return header_attr

    cols    = {}
    colnum  = lines[4].split()
    colname = lines[5].split()
    for a,b in zip(colname,colnum):
        cols[a] = int(b)
            
    data = []
    for i in range(data_rows):
        line = f.readline()
        v=line.split()
        try: 
            vv=np.array(v,dtype='float64')
        except ValueError:
            for item in v:
                if item.__contains__('.') and not item.__contains__('E'):
                    v[v.index(item)]='0'
        data.append(vv)
    f.close()
    a=np.array(data) 
    data = []
    return header_attr, cols, a


def cleanstarlog(file_in):
    ''' cleaning star.log, e.g. to take care of repetitive restarts
    
    private, should not be called by user directly

    file_in     typically the filename of the mesa output star.log file,
                creates a clean file called star.logsa

    (thanks to Raphael for providing this tool)            
    '''

    file_out=file_in+'sa'
    f = open(file_in)
    lignes = f.readlines()
    f.close()

    nb    = np.array([],dtype=int)   # model number
    nb    = np.concatenate((nb    ,[  int(lignes[len(lignes)-1].split()[ 0])])) 
    nbremove = np.array([],dtype=int)   # model number
    i=-1

    for i in np.arange(len(lignes)-1,0,-1):
        line = lignes[i-1]
        
        if i > 6 and line != "" :
            if int(line.split()[ 0])>=nb[-1]:
                nbremove = np.concatenate((nbremove,[i-1])) 
            else:
                nb = np.concatenate((nb    ,[  int(line.split()[ 0])])) 
    i=-1
    for j in nbremove:
        lignes.remove(lignes[j])
 
    fout=file(file_out,'w')
    for j in np.arange(len(lignes)):
        fout.write(lignes[j])
    fout.close()
 
