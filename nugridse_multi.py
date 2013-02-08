'''nugridse_multi allows to create plots of of data from multiple run dirs with se-type h5 files

   sedir    - the directory on which the h5 files can be
   pattern - a string pattern that the file names have to
                    contain in order to be read

   

   usage:
   start by loading the module
   1 : import nugridse_multi as mp_multi

   you can get help:
   2 : help mp_multi

   next, read the different run dirs in 'path/to/dir/with/run/dirs':   
   3 : pt=mp_multi.read_run_dirs('path/to/dir/with/run/dirs')
   
   which would read all h5 files from the run directories; again, do
   help(mp.se) or mp.se? 

    ........
'''


import matplotlib.pyplot as plt
import numpy as np
import os
from nugridse import *

##########################
#### Experimental class  ####
##########################


class read_run_dirs(se):

	'''
		This class allows access to multiple runs in the path dir.
		Instead of defining one path, an array of paths , mult_dir can be
		defined. In this case each path in multi_dir represents a run.
		
	'''

	def __init__(self,rundir='.',multi_dir=[]):

		slist = os.listdir(rundir)


		pattern='.h5'
		expr = re.compile(pattern)
		expr = re.compile(pattern)
		self.runs_H5_out=[]
		self.runs_H5_surf=[]
		self.run_dirs_name=[]	
		for element in slist:
			run_path=rundir+"/"+element
			if os.path.isdir(run_path+"/H5_out") and os.path.isdir(run_path+"/H5_surf"):
				sefiles = os.listdir(run_path+"/H5_out")
				if (filter(expr.search,sefiles)) <1:
					print "Warning: No hdf5 files found in "+run_path+"/H5_out"	
				sefiles = os.listdir(run_path+"/H5_surf")
				if (filter(expr.search,sefiles)) <1:
					print "Warning: No hdf5 files found in "+run_path+"/H5_surf"				
				self.runs_H5_surf.append(rundir+"/"+element+"/H5_surf")
				self.runs_H5_out.append(rundir+"/"+element+"/H5_out") 
				self.run_dirs_name.append(element)
				print "Read "+element		
			
			
			
	def weighted_yields_massgrid(self,runs=[],isotopes=[],cycles=[],title=''):
		'''
			Plots behaviour star mass dependent yields of different isotopes - specify dirs at the beginning
			Beware  of different z
			The isotopes should be being produced because of the logscale...
			At the moment 4 isotopes supported

			runs : If array is empty function uses all available directories.					
				
			" -Final wind yields - isotopic composition folded with IMF" 
			dirs=['M1.650Z0.0001', 'M2.000Z0.0001','M3.000Z0.0001','M5.000Z0.0001']
			dirs=['M1.650Z0.0010', 'M2.000Z0.0010','M3.000Z0.0010','M5.000Z0.0010']
			dirs=['M1.650Z0.0060', 'M2.000Z0.0060','M3.000Z0.0060','M5.000Z0.0060']		
		
			cycles=[[0,80000,1000],[0,80000,1000],[0,80000,1000],[0,120000,1000]]	
		
			pap.weighted_yields_massgrid(runs=['M1.650Z0.0001', 'M2.000Z0.0001','M3.000Z0.0001','M5.000Z0.0001'],isotopes=["He-4","C-12","O-16"],cycles=[[0,10000,1000],[0,10000,1000],[0,10000,1000],[0,10000,1000]])
	
		'''

		HDF5_surf=[]
		if len(runs) ==0:
			HDF5_surf=self.runs_H5_surf
		else:
			for i in range(len(self.run_dirs_name)):
				if self.run_dirs_name[i] in runs:
					HDF5_surf.append(self.runs_H5_surf[i])

		legend_k=0
		sefiles=[]
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		

		color=['r','b','g','k']
		marker_type=['o','D','s','p']
		line_style=['--','-','-.',':']
		fig=plt.figure()
		plt.rcParams.update({'font.size': 16})
		plt.rc('xtick', labelsize=16) 
		plt.rc('ytick', labelsize=16) 
		ax = fig.add_subplot(1,1,1)
		ax.set_yscale('log')
	
	
		for j in range(len(isotopes)):
			legend_k=0
			yields=[]
			star_mass_array=[]	
			for i in range(len(HDF5_surf)):
				#if (legend_k ==0):
				star_mass=sefiles[i].get("mini")[0]				
				star_mass_array.append(star_mass)
				#below way not very efficient
				iso_yield_folded,iso_yield_unfolded =self.weighted_yields(sefiles[i],cycles[i][0],cycles[i][1],cycles[i][2],isotopes,star_mass,False,legend=isotopes[j],color=color[j],title="",plot_fig=False)
						
				yields.append(iso_yield_folded[j])				
				#	legend_k=1
				#	print yields			
					#plt.plot(star_mass,yields[j],marker='*',markersize=8,mfc=color[j],linestyle='None',label=isotopes[j])			
				#else:
					#yields.append(weighted_yields(sefiles[i],cycles[i][0],cycles[i][1],cycles[i][2],isotopes,star_mass,False,color=color[j],title="",plot_fig=False))
			order_indices=[i[0] for i in sorted(enumerate(star_mass_array), key=lambda x:x[1])]
			star_mass_array = [ star_mass_array[i] for i in order_indices]		
			yields = [ yields[i] for i in order_indices]		
			plt.plot(star_mass_array,yields,marker=marker_type[j],markersize=10,mfc=color[j],linestyle=line_style[j],label=isotopes[j])				
					#plt.plot(star_mass,yields[j],marker='*',markersize=8,mfc=color[j],linestyle='None')
		plt.legend()				
		plt.xlabel("M/M$_{\odot}$",fontsize=20)
		plt.minorticks_on()
		plt.ylabel("Weighted stellar yields",fontsize=20)
		plt.title(title)		
		plt.xlim(0,max(star_mass_array)+2)
	


	

	def multi_surface_plots(self,runs=[],cycles=[[3470,50000,1000],[4740,60000,1000]],x_range=[-1.8,0.1],y_range=[]):
		'''
		Plot multiple runs
		cycle must be 2d array
		
	
		# creating prefix variable, we have agreed that this is an 11 character
		# name with the initial mass and metallicity encoded: M5.00Z0.010 
		
		
		'''	
	
	
		line_style=['--','-','-.',':']
		color=['r','b','g','k']

		sefiles=[]
		legend=[]
		HDF5_surf=[]
		if len(runs) ==0:
			HDF5_surf=self.runs_H5_surf
			runs=self.run_dirs_name
		else:
			for i in range(len(self.run_dirs_name)):
				if self.run_dirs_name[i] in runs:
					HDF5_surf.append(self.runs_H5_surf[i])

		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		

		for i in range(len(HDF5_surf)):							
			mass=sefiles[i].get("mini")[0]	
			z=sefiles[i].get("zini")[0]	
			legend.append(str(mass)+"M$_{\odot}$ Z= "+str(z))					
			self.surface_plots(HDF5_surf[i],cycles[i],legend[i],color[i],line_style[i],title="")
	
		#HDF5_dirs=["/nfs/rpod3/critter/Results/test2_template_mppnp/M1.650Z0.0001","/nfs/rpod3/critter/Results/test2_template_mppnp/M5.000Z0.0001"],	
			
	
		plt.figure(2)	
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()		
		plt.plot([-1000,1000],[0,0],"-",c="k",linewidth=2)
		plt.plot([0,0],[-1000,1000],"-",c="k",linewidth=2)		
		if xmax<0:
			xmax=0.1
		if xmin >0:
			xmin=-0.1	
		if ymax<0:
			ymax=0.1	
		if ymin >0:
			ymin=-0.1		
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)					
		plt.figure(3)	
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()	
		plt.plot([-1000,1000],[0,0],"-",c="k",linewidth=2)
		plt.plot([0,0],[-1000,1000],"-",c="k",linewidth=2)	
		if xmax<0:
			xmax=0.1
		if xmin >0:
			xmin=-0.1	
		if ymax<0:
			ymax=0.1	
		if ymin >0:
			ymin=-0.1
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)			
			
		plt.figure(4)	
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()	
		plt.plot([-1000,1000],[0,0],"-",c="k",linewidth=2)
		plt.plot([0,0],[-1000,1000],"-",c="k",linewidth=2)	
		if xmax<0:
			xmax=0.1
		if xmin >0:
			xmin=-0.1	
		if ymax<0:
			ymax=0.1	
		if ymin >0:
			ymin=-0.1
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)		
	
		plt.figure(5)	
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()	
		plt.plot([-1000,1000],[0,0],"-",c="k",linewidth=2)
		plt.plot([0,0],[-1000,1000],"-",c="k",linewidth=2)	
		if xmax<0:
			xmax=0.1
		if xmin >0:
			xmin=-0.1	
		if ymax<0:
			ymax=0.1	
		if ymin >0:
			ymin=-0.1
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)		
	
		if len(x_range)>0:
			plt.figure(1);plt.xlim(x_range[0],x_range[1])
			plt.figure(2);plt.xlim(x_range[0],x_range[1])
			plt.figure(3);plt.xlim(x_range[0],x_range[1])
			plt.figure(4);plt.xlim(x_range[0],x_range[1])
			plt.figure(5);plt.xlim(x_range[0],x_range[1])

	#if len(y_range)>0:
	#	figure(1);plt.ylim(y_range[])
	#	figure(2);plt.ylim(y_range)
	#	figure(2);plt.ylim(y_range)

	
	def surface_plots(self,HDF5surf_dir,cycles,legend="",color="r",line_style="-",title="",hs=["Ba","La","Nd","Sm"],ls=["Sr","Y","Zr"]): 
	
		'''surface hdf5
	
			ls: Sr,Y,Zr  - elements
			hs: Ba,La,Nd,Sm -elements
		
			isotopes,elements,HDF5_dir,cycles,legend,title,x_range=[],y_range=[]):
	
		'''	
		
		
		
		sefiles=se(HDF5surf_dir)

		t0_model=cycles[0]


		#2 choices to plot
		#modelrange=[10,30,10]
		#t0_model=0
		t0_time=sefiles.se.get(t0_model,"age")
		starage=[]
		ls_abu=[]
		hs_abu=[]
		fe_abu=[]
		rb_abu=[]####n-density + efficiency of N22(a,n)
		sr_abu=[]
		quot_abu_mg_26_24=[]
		quot_abu_mg_25_24=[]
		quot_abu_gd_152_154=[]
		quot_abu_zr_96_94=[]		
		hs_isotopes=[]
		ls_isotopes=[]
		rb_isotopes=[]
		fe_isotopes=[]
		sr_isotopes=[]
		mg_isotopes=['Mg-24','Mg-25','Mg-26']
		gd_isotopes=['Gd-152','Gd-154']
		zr_isotopes=['Zr-96','Zr-94']
		for iso in sefiles.se.isotopes:
			for i in range(len(ls)):
				if ls[i] in iso:
					ls_isotopes.append(iso)
			for i in range(len(hs)):
				if hs[i] in iso:
					hs_isotopes.append(iso)
			if "Rb" in iso:
					rb_isotopes.append(iso)
			if "Fe" in iso:
					fe_isotopes.append(iso)
			if "Sr" in iso:
					sr_isotopes.append(iso)		
		print hs_isotopes,ls_isotopes,rb_isotopes,fe_isotopes,sr_isotopes
		for cycle in np.arange(cycles[0],cycles[1],cycles[2]):
			#print cycle
			starage.append(sefiles.se.get(cycle,"age")-t0_time)
			ls_abu_1=0
			hs_abu_1=0
			rb_abu_1=0
			fe_abu_1=0
			sr_abu_1=0			

			for iso in sefiles.se.isotopes:
				if iso in ls_isotopes:
					ls_abu_1+=sefiles.get(cycle,"iso_massf",iso)
				if iso in hs_isotopes:
					hs_abu_1+=sefiles.get(cycle,"iso_massf",iso)		
				if iso in rb_isotopes:
					rb_abu_1+=sefiles.get(cycle,"iso_massf",iso)
				if iso in sr_isotopes:
					sr_abu_1+=sefiles.get(cycle,"iso_massf",iso)
				if iso in fe_isotopes:
					fe_abu_1+=sefiles.get(cycle,"iso_massf",iso)
			#create hs and ls as mean values of the elements		
			hs_abu.append(hs_abu_1/len(hs))
			ls_abu.append(ls_abu_1/len(ls))
			fe_abu.append(fe_abu_1)
			rb_abu.append(rb_abu_1)
			sr_abu.append(sr_abu_1)			
			quot_abu_mg_26_24.append( sefiles.get(cycle,"iso_massf",mg_isotopes[2])/sefiles.get(cycle,"iso_massf",mg_isotopes[0]) )
			quot_abu_mg_25_24.append( sefiles.get(cycle,"iso_massf",mg_isotopes[1])/sefiles.get(cycle,"iso_massf",mg_isotopes[0] ))
			quot_abu_gd_152_154.append( sefiles.get(cycle,"iso_massf",gd_isotopes[0])/sefiles.get(cycle,"iso_massf",gd_isotopes[1])) 
			quot_abu_zr_96_94.append( sefiles.get(cycle,"iso_massf",zr_isotopes[1])/sefiles.get(cycle,"iso_massf",zr_isotopes[0]) )
	
		print quot_abu_mg_26_24
		print quot_abu_mg_25_24
		print quot_abu_gd_152_154
		print quot_abu_zr_96_94	
		
		#Plotting different figures, beware of square bracket notation
		plt.figure(1)
		plt.plot(starage,np.log10(np.array(hs_abu)/np.array(ls_abu)),marker='*',markersize=8,mfc=color,linestyle=line_style,label=legend)
		plt.xlabel("Star age")
		plt.ylabel("[hs/ls]")
		plt.legend()
		plt.title(title)
		plt.draw()
		#Paper 1 plot
		plt.figure(2);plt.plot(np.log10(np.array(hs_abu)/np.array(ls_abu)),np.log10(np.array(ls_abu)/np.array(fe_abu)),c=color,linestyle=line_style,label=legend)  #,marker='*',markersize=8,
		plt.xlabel("[hs/ls]");plt.ylabel("[ls/Fe]")
		plt.legend()
		plt.title(title)
		plt.draw()	
	
		#Paper 1 plot
		plt.figure(3);plt.plot(np.log10(np.array(hs_abu)/np.array(ls_abu)),np.log10(np.array(rb_abu)/np.array(sr_abu)),c=color,linestyle=line_style,label=legend)
		plt.xlabel("[hs/ls]")
		plt.ylabel("[Rb/Sr]")
		plt.legend()
		plt.title(title)
		plt.draw()

		#Paper 1 plot
		plt.figure(4);plt.plot(quot_abu_mg_25_24,quot_abu_mg_26_24,c=color,linestyle=line_style,label=legend)
		plt.xlabel("$^{25}$Mg/$^{24}$Mg")
		plt.ylabel("$^{26}$Mg/$^{24}$Mg")
		plt.legend()
		plt.title(title)
		plt.draw()
	
		#Paper 1 plot
		plt.figure(5);plt.plot(quot_abu_zr_96_94,quot_abu_gd_152_154,c=color,linestyle=line_style,label=legend)
		plt.xlabel("$^{96}$Zr/$^{94}$Zr")
		plt.ylabel("$^{152}$Gd/$^{154}$Gd")
		plt.legend()
		plt.title(title)
		plt.draw()
		
		return [starage,hs_abu,ls_abu,fe_abu,rb_abu,sr_abu]






	#####################Internal use

	def weighted_yields(self,sefiles,cyclestart=11000,cycleend=12000,sparse=100,isotopes=["H-1","H-2","H-3","He-4"],star_mass=1.65,label=True,legend="",color="r",title=" -Final wind yields - isotopic composition folded with IMF",plot_fig=True):
		import utils as u
		import re	
		import nugridse as mp
		'''
				Uses H5_surf
	  	    This function returns the wind yields and ejected masses for stable species.  The command
     		   is
		'''
	
	
		X_i, E_i =sefiles.windyields(cyclestart, cycleend, sparse)
		###calculate isotopes
		#isotopes=['Fe-56','Co-60','Ni-61','Cu-65','Zn-67',]
		#isotopes['Ga-71','Ge-73','As-75','Se-77','Br-82','Kr-84','Rb-87','Sr-88','Y-89','Zr-92','Nb-94','Zr-96','Mo-96','Ba-137','Ba-138','La-139','Ce-140','Nd-142','Sm-144','Pb-206']					
		#206
		plotiso=[]
		mass=[]
		name=[]
		isotope_names = sefiles.se.isotopes
		sefiles.mass_frac = [X_i] #double array important for decay
	
		#convert isotope name scheme
		u.convert_specie_naming_from_h5_to_ppn(isotope_names) #-define u.spe
		names_ppn_world = u.spe 
		number_names_ppn_world = u.n_array
		u.define_zip_index_for_species(names_ppn_world,number_names_ppn_world) #--define u.cl!	
	
		#let species decay
		u.stable_specie()
		sefiles.decay([X_i]) #double array is important for decay()
		X_i_decayed=decayed_multi_d[0]
		
		#convert specified isotopes in correct name scheme
		u.convert_specie_naming_from_h5_to_ppn(isotopes) #-define u.spe
		isotopes_ppn_world = u.spe
	
		for i in range(len(isotopes_ppn_world)):
			isotopes_ppn_world[i]=isotopes_ppn_world[i].upper()
		
		for j in range(len(u.stable)):
			if u.stable[j] in isotopes_ppn_world:
				plotiso.append(X_i_decayed[j])
				mass.append(int(re.findall(r'\d+',u.stable[j])[0]))
				name.append(u.stable[j])
		###Folding and plotting###			
		#Salpeter mass function,Kroupa 2001:(lazyness..)
		#f(M)=M**-2.3
		f_m=(star_mass**(-2.3))
		plotiso_IMF=[]
		for i in range(len(plotiso)):
			plotiso_IMF.append(f_m*plotiso[i])
		if plot_fig == True:	
			plt.plot(mass,plotiso_IMF,marker='*',markersize=8,mfc=color,linestyle='None',label=legend)
			plt.xlabel("Mass number A")
			plt.ylabel("Abundance")
			plt.title(title)	
			if label ==True:
				k=-1
				for j in range(len(plotiso)):
					k+=1
					if isotopes_ppn_world[j] != name[j]:
						k+=1
						print "Input isotope >>"+isotopes_ppn_world[j]+"<< not stable and will be removed"
					plt.annotate(name[j], xytext = (0, 10),textcoords = 'offset points' ,xy=(mass[j], plotiso_IMF[j]))
			mina=mass[0]-4 	
			if mina<0:
				mina=0
			plt.xlim(mina,mass[-1]+4)
			plt.legend()
		return plotiso_IMF, plotiso #folded and unfolded - both stable
		
	
	
	def is_stable(self,isotope):
		stable_file="/astro/critter/Documents/Master/Simulations/pythonscripts/stable.txt"
		stable_list = np.loadtxt(stable_file, dtype='string')
		isotope="".join(isotope.split())
		for i in range(len(stable_list)):
		#print "|"+isotope+stable_list[i][0]+"|"	
			if isotope == stable_list[i][0]:
				return stable_list[i][1]
				break

	
	
	
	
