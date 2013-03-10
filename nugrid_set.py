''' Tool for analysis of mesa and mppnp data of paper2 calculations

    v0.1, 15JAN2013: Christian Ritter


	##########################
	#### Still experimental tool  ####
	##########################



	This tools contains two classes which can deal with set data
	of MESA and MPPNP runs. The class mesa_set allows the analysis
	of multiple mesa runs and the class mppnp_set the analysis
	of multiple MPPNP dirs.
	
	Both classes allow during their initialisation the set of 
	one dir path - containing the set calculations or the set of 
	multiple run dirs.


	Python 2.7.3 (default, Jul 24 2012, 10:05:38) 
	Type "copyright", "credits" or "license" for more information.

	IPython 0.12 -- An enhanced Interactive Python.
	?         -> Introduction and overview of IPython's features.
	%quickref -> Quick reference.
	help      -> Python's own help system.
	object?   -> Details about 'object', use 'object??' for extra details.

	IPython profile: numpy

	Welcome to pylab, a matplotlib-based Python environment [backend: GTKAgg].
	For more information, type 'help(pylab)'.

	#for working with mesa set - e.g. from inside the set directory

	In [1]: import nugrid_set as set
	In [2]: setresults=set.mesa_set(".")
	In [3]: setresults.multi_DUP()

	e.g for calculation of DUP parameter evolution of all set runs.


	#for working with mppnp set 

        In [2]: setresults=set.mppnp_set(".")
        In [3]: setresults.weighted_yields_massgrid()

	#for calculation of yields of the set 

'''
##need to be changed, matplotlib import
import matplotlib.pyplot as plt
from matplotlib.pyplot import * 
import numpy as np
import os
from nugridse import *
from mesa import *
from paper2_analysis import *
import glob

import utils
symbs=utils.symbol_list('lines1')



class mppnp_set(se):

	'''
		This class allows access to multiple runs in the path dir.
		Instead of defining one path, an array of paths , mult_dir can be
		defined. In this case each path in multi_dir represents a run.
		

                rundir - dir with run directory
                multi_dir - specify dirs at different location
                extra_label - set labels of runs, else run dir name will be chosen      

                for paper2 agb analysis:
                mppnpset=set.mppnp_set(multi_dir=["M1.650Z0.0001","M2.000Z0.0001/","M3.000Z0.0001","M5.000Z0.0001/"])
                mppnpset=set.mppnp_set(multi_dir=["M1.650Z0.0010","M2.000Z0.0010/","M3.000Z0.0010","M5.000Z0.0010/"])     
                mppnpset=set.mppnp_set(multi_dir=["M1.650Z0.0060","M2.000Z0.0060/","M3.000Z0.0060","M5.000Z0.0060/"])


        
                e.g.:
        
		
                        
                        or
                                

                Using the set_plot* functions you can instantly plot diagrams
                for all runs defined in MPPNP set


	


	'''

	def __init__(self,rundir='.',multi_dir=[],extra_label=[]):

                if len(multi_dir)==0:
			slist = os.listdir(rundir)
		else:
                        slist=multi_dir

		pwd=os.getcwd()
		pattern='.h5'
		expr = re.compile(pattern)
		expr = re.compile(pattern)
		self.runs_H5_out=[]
		self.runs_H5_surf=[]
		self.run_dirs_name=[]
		self.extra_label=[]
                runs_H5_out=[]
                runs_H5_surf=[]
                run_dirs_name=[]
		extra_label_1=[]
		i=0
		for element in slist:
                        
			if len(multi_dir)==0:
                                run_path=pwd+"/"+rundir+"/"+element
                        else:
                                if multi_dir[i][0] == "/":
                                        run_path=multi_dir[i]
                                else:
                                        run_path=pwd+"/"+multi_dir[i]

			if os.path.isdir(run_path+"/H5_out") and os.path.isdir(run_path+"/H5_surf"):
				sefiles = os.listdir(run_path+"/H5_out")
				if (filter(expr.search,sefiles)) <1:
					print "Warning: No hdf5 out files found in "+run_path+"/H5_out"
				else:
					runs_H5_out.append(run_path+"/H5_out")	
				sefiles = os.listdir(run_path+"/H5_surf")
				if (filter(expr.search,sefiles)) <1:
					print "Warning: No hdf5 surf files found in "+run_path+"/H5_surf"
				else:
					runs_H5_surf.append(run_path+"/H5_surf")				
				if len(extra_label)>0:
					extra_label_1.append(extra_label[i])
				else:
					extra_label_1.append(element)
				run_dirs_name.append(element)
				print "Read "+run_path		
			i+=1
		###order lists###
		sorted_indices=sorted(range(len(run_dirs_name)),key=lambda k: run_dirs_name[k])	
		print sorted_indices
		print run_dirs_name
		k=0
		for i in sorted_indices:
			print i
			self.run_dirs_name.append(run_dirs_name[i])
			self.extra_label.append(extra_label_1[i])
			self.runs_H5_surf.append(runs_H5_surf[i])
			self.runs_H5_out.append(runs_H5_out[i])
			k+=1
		print self.run_dirs_name
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
		if (len(runs) == 0):
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
		marker_type=['o','p','s','D']
		line_style=['--','-','-.',':']
		line_width=12*[5,3,4,6]
		fig=plt.figure()
		plt.rcParams.update({'font.size': 16})
		plt.rc('xtick', labelsize=16) 
		plt.rc('ytick', labelsize=16) 
		ax = fig.add_subplot(1,1,1)
		ax.set_yscale('log')
		z_index_files=[]
		z_values=[]
		j=-1
		for i in range(len(HDF5_surf)):
			j+=1
			star_z=sefiles[i].get("zini")[0]
			if star_z not in z_values:
				z_values.append(star_z)
				z_index_files.append([])
			z_index_files[ z_values.index(star_z)].append(i)
		print "index files",z_index_files		

		color_iso=-1
		yields_1=[]
		t=0
		for j in range(len(isotopes)):
			yields=[]
			star_mass_array=[]
			color_iso+=1
			legend_k=-1
			t=0	
			for w in range(len(z_index_files)):
				star_mass_array=[]
				yields=[]
				legend_k+=1
				for k in z_index_files[w]: 
					#if (legend_k ==0):
					star_mass=sefiles[k].get("mini")[0]
					star_z=sefiles[k].get("zini")[0]				
					star_mass_array.append(star_mass)
					#below way not very efficient
					if cycles[k][1]==-1:
						endcycle=sefiles[k].get("model_number")[-1] - 5000
					else:
						endcycle=cycles[k][1]
					if j==0:
						iso_name,iso_yield_folded,iso_yield_unfolded,iso_prod_name,iso_prod =self.weighted_yields_production_factor(sefiles[k],cycles[k][0],endcycle,cycles[k][2],isotopes,star_mass,False,legend=isotopes[j],color="k",title="")
						yields_1.append(iso_yield_folded)
					yields.append(yield_1[t][j])
					t+=1				
				#	legend_k=1
				#	print yields			
					#plt.plot(star_mass,yields[j],marker='*',markersize=8,mfc=color[j],linestyle='None',label=isotopes[j])			
				#else:
					#yields.append(weighted_yields(sefiles[i],cycles[i][0],cycles[i][1],cycles[i][2],isotopes,star_mass,False,color=color[j],title="",plot_fig=False))
				order_indices=[i[0] for i in sorted(enumerate(star_mass_array), key=lambda x:x[1])]
				star_mass_array = [ star_mass_array[i] for i in order_indices]		
				yields = [ yields[i] for i in order_indices]		
				#print "plotting graph",len(star_mass)
				plt.plot(star_mass_array,yields,marker=marker_type[legend_k],markersize=10,mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[legend_k],label=iso_name[j]+" , "+str(star_z)+"$Z_{\odot}$"  )				
					#plt.plot(star_mass,yields[j],marker='*',markersize=8,mfc=color[j],linestyle='None')
		plt.legend()				
		plt.xlabel("M/M$_{\odot}$",fontsize=20)
		plt.minorticks_on()
		plt.ylabel("Weighted stellar yields",fontsize=20)
		plt.title(title)		
		plt.xlim(0,max(star_mass_array)+2)



	def set_pocket(self,runs=[],cycle=[45532,47566],mass_cell=[0.641981,0.641981],isotopes=[],x_charge=False,mass_number_range=[],decayed=False,iso_label=False,title="PDCZ",colors=["red","blue","green","black"],yax_log=True,filename_norm="iniab2.0E-02GN93.ppn",elem_norm=''):

		'''
			plots isotopic composition for different runs, each with specific cycle number cycle[i] and specific mass cell mass_cell[i].
			Initial abundace file for normalization can be chosen as filename_norm, but MUST be in run directory

			decayed - plot only decayed isotopes
			mass_number_range - array of min mass number and max mass number, isotopes inbetween will be plottet, if set, isotope array will be ignored
			iso_label - set labels of isotopes True or False
			title - 
			colors -
			filename_norm - not yet included, normalization with initial abundance file
			....
			e.g.
			
			setse.set_pocket(cycle=[32840,29390],mass_cell=[0.6300,0.6076],isotopes=["He-4","C-12","O-16","Ne-22"],mass_number_range=[80,140],decayed=True)
			

	
		'''
	
		import nugridse as mp	
                sefiles=[]
                legend=[]
                HDF5_out=[]
                extra_label=[]

		if len(runs) ==0:
			HDF5_out=self.runs_H5_out
			runs=self.run_dirs_name
			extra_label=self.extra_label
		else:
			for i in range(len(self.run_dirs_name)):
				if self.run_dirs_name[i] in runs:
					HDF5_out.append(self.runs_H5_out[i])
                                        extra_label.append(self.extra_label[i])
		for i in range(len(HDF5_out)):
			reload(mp)
			file_norm=HDF5_out[i][:-6]+filename_norm
			sefiles=se(HDF5_out[i])
                        mass=sefiles.get("mini")[0]
                        z=sefiles.get("zini")[0]
                        legend=str(mass)+"M$_{\odot}$ Z= "+str(z)+", "+extra_label[i]
			if len(isotopes)==0:
                                self.plot_abu_atmasscoord(mp,sefiles,cycle[i],mass_cell[i],["He-4","C-12","O-16","Ne-22"],x_charge,mass_number_range,decayed,file_norm,elem_norm,yax_log,label=iso_label,legend=legend,color=colors[i],title=title)

				#plt.figure(111);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=["He-4","C-12","O-16","Ne-22"],label=True,legend=legend,color=colors[i],title=title)
				#plt.figure(222);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=['Fe-56','Co-60','Ni-61','Cu-65','Zn-67','Ga-71','Ge-73','As-75','Se-77','Br-82','Kr-84','Rb-87','Sr-88','Y-89','Zr-92','Nb-94','Zr-96','Mo-96','Ba-137','Ba-138','La-139','Ce-140','Nd-142','Sm-144','Pb-206'],label=True,legend=legend,color=colors[i],title=title)
			else:
				#plt.figure(111);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=isotopes,label=True,legend=legend,color=colors[i],title=title)
				self.plot_abu_atmasscoord(mp,sefiles,cycle[i],mass_cell[i],isotopes,x_charge,mass_number_range,decayed,file_norm,elem_norm,yax_log,label=iso_label,legend=legend,color=colors[i],title=title)


        def set_surface_abundances(self,runs=[],cycle=[-1,-1],isotopes=[],x_charge=True,mass_number_range=[],decayed=False,iso_label=False,title="",colors=["red","blue","green","black"],yax_log=True,filename_norm="iniab2.0E-02GN93.ppn",elem_norm=''):

                '''
                        plots isotopic composition for different runs, each with specific cycle number cycle[i] and specific mass cell mass_cell[i].
                        Initial abundace file for normalization can be chosen as filename_norm, but MUST be in run directory

                        decayed - plot only decayed isotopes
                        mass_number_range - array of min mass number and max mass number, isotopes inbetween will be plottet, if set, isotope array will be ignored
                        iso_label - set labels of isotopes True or False
                        title - 
                        colors -
                        filename_norm - not yet included, normalization with initial abundance file
                        ....
                        e.g.
                        
                        setse.set_pocket(cycle=[32840,29390],mass_cell=[0.6300,0.6076],isotopes=["He-4","C-12","O-16","Ne-22"],mass_number_range=[80,140],decayed=True)
                        

        
                '''


                import nugridse as mp
                sefiles=[]
                legend=[]
                HDF5_out=[]
                extra_label=[]


                if len(runs) ==0:
                        HDF5_out=self.runs_H5_out
                        runs=self.run_dirs_name
                        extra_label=self.extra_label
                else:
                        for i in range(len(self.run_dirs_name)):
                                if self.run_dirs_name[i] in runs:
                                        HDF5_out.append(self.runs_H5_out[i])
                                        extra_label.append(self.extra_label[i])
                for i in range(len(HDF5_out)):
                        reload(mp)
                        file_norm=HDF5_out[i][:-6]+filename_norm
                        sefiles=se(HDF5_out[i])
			if cycle[i] ==-1:
				cycle_1=sefiles.se.cycles[-1]
			else:
				cycle_1=cycle[i]
			masscell=sefiles.get(cycle_1,"mass")[0]
                        mass=sefiles.get("mini")[0]
                        z=sefiles.get("zini")[0]
                        legend=str(mass)+"M$_{\odot}$ Z= "+str(z)+", "+extra_label[i]
                        if len(isotopes)==0:
                                self.plot_abu_atmasscoord(mp,sefiles,cycle_1,masscell,isotopes,x_charge,mass_number_range,decayed,file_norm,elem_norm,yax_log,label=iso_label,legend=legend,color=colors[i],title=title)

                                #plt.figure(111);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=["He-4","C-12","O-16","Ne-22"],label=True,legend=legend,color=colors[i],title=title)
                                #plt.figure(222);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=['Fe-56','Co-60','Ni-61','Cu-65','Zn-67','Ga-71','Ge-73','As-75','Se-77','Br-82','Kr-84','Rb-87','Sr-88','Y-89','Zr-92','Nb-94','Zr-96','Mo-96','Ba-137','Ba-138','La-139','Ce-140','Nd-142','Sm-144','Pb-206'],label=True,legend=legend,color=colors[i],title=title)
                        else:
                                #plt.figure(111);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=isotopes,label=True,legend=legend,color=colors[i],title=title)
                                self.plot_abu_atmasscoord(mp,sefiles,cycle_1,masscell,isotopes,x_charge,mass_number_range,decayed,file_norm,elem_norm,yax_log,label=iso_label,legend=legend,color=colors[i],title=title)


	def plot_saga_data(self,fig=7,label_plot=["C-rich RGB","CEMP RGB","EMP RGB"],path="/nfs/rpod3/critter/DPG/saga_platform"):
		
		figure_data=plt.figure(fig)
		pwd=os.getcwd()
		os.chdir(path)
		execfile(path+"/"+"read_files_images.py")
		print pwd
		os.chdir(pwd)
		print "done"	
		


	def set_surface_plots(self,runs=[],cycles=[[10000,1000],[10000,1000]],t0_model=[],decayed=True,mesarunpath="",x_range=[-1.8,0.1],y_range=[]):
		'''
		Plot surface abundance evolution of multiple runs

		plots [hs/ls] vs star age
		plots  [ls/Fe] vs [hs/ls] 
		[Rb/Sr] [hs/ls]
		plots Mg25/Mg24 vs Mg26/Mg24
          	plots plt.xlabel("$^{96}$Zr/$^{94}$Zr") plt.ylabel("$^{152}$Gd/$^{154}$Gd")
		plots [ls/Fe] vs star age
		cycles - 2d array, containing the start and stop cycle for each run and the sparsity
		t0_model - the starting model to plot, for the time dependent plot, t0_model will be..
		x_range - define range of 
	
		!!Check/c setse.run_dirs_name array for order of runs for setting the cycles parameter!!	
	
		# creating prefix variable, we have agreed that this is an 11 character
		# name with the initial mass and metallicity encoded: M5.00Z0.010 
		
		e.g. 
		
		set_surface_plots(self,runs=[],cycles=[[3113,50000,1000],[3031,60000,1000]],t0_model=[3113,3031],x_range=[],y_range=[]	
		test2_template: [3620, 5145, 3030, 4757, 11997, 5288, 9470, 5123, 2996, 3312, 5828, 4749]

		['./M3.000Z0.0060/LOGS',
		 './M5.000Z0.0010/LOGS',
		 './M2.000Z0.0001/LOGS',
		 './M3.000Z0.0001/LOGS',
		 './M1.650Z0.0060/LOGS',
		 './M3.000Z0.0010/LOGS',
		 './M1.650Z0.0010/LOGS',
		 './M5.000Z0.0060/LOGS',
		 './M2.000Z0.0010/LOGS',
		 './M2.000Z0.0060/LOGS',
	 	'./M1.650Z0.0001/LOGS',
		 './M5.000Z0.0001/LOGS
	
		'''	
		

                color=['r','b','k','g']	#'r','r','r','b','b','k','g','k','k','g','b','g']
                marker_type=['o','o','o','p','p','p','D','D','D','s','s','s']
		line_style=3*['--','-','-.',':']
		#IDEA PASS PARAMETE BETWEEN BOTH CLASSES
	        #if len(runs)==0:
                #	t0_model=self.set_find_first_TP()
		#cycles=[]
		#for i in range(len(cyc)):
		#	cycles.append([t0_model[i],cyc[i][0],cyc[i][1]])
		
		print "Do you want a TP-dependece on time axis?"
		print "Read corresponding MESA runs?"
		#answer=raw_input("Read corresponding MPPNP runs?")
		#if answer ==True:
		if len(mesarunpath) >0:
		#mesa_path=raw_input("Read corresponding MPPNP runs?")
			dir_1=mesarunpath
			multi_dir=[]
			for i in range(len(self.run_dirs_name)):
				multi_dir.append(mesarunpath+"/"+self.run_dirs_name[i])
			mesaset=mesa_set(multi_dir=multi_dir)
			peak_lum_model_array_1,h1_mass_min_DUP_model_array = mesaset.multi_DUP(t0_model=[],plot_fig=False)
			t0_model=[]
			peak_lum_model_array=[]
			cycles=[]
			for i in range(len(peak_lum_model_array_1)):
				cycles.append(map(int,np.array(peak_lum_model_array_1[i])))
				t0_model.append(int(peak_lum_model_array_1[i][0]))
		print "Cycles: ",cycles

		sefiles=[]
		legend=[]
		HDF5_surf=[]
		extra_label=[]
		if len(runs) ==0:
			HDF5_surf=self.runs_H5_surf
			runs=self.run_dirs_name
			extra_label=self.extra_label
		else:
			for i in range(len(self.run_dirs_name)):
				if self.run_dirs_name[i] in runs:
					HDF5_surf.append(self.runs_H5_surf[i])
					extra_label.append(self.extra_label[i])
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		

		for i in range(len(HDF5_surf)):							
			mass=sefiles[i].get("mini")[0]	
			z=sefiles[i].get("zini")[0]
			if z <0.001:
				color_1=color[0]
			elif z <0.006:
				color_1=color[1]
			else:
				color_1=color[2]	
			legend=str(mass)+"M$_{\odot}$ Z= "+str(z)+", "+extra_label[i]
			self.surface_plots(HDF5_surf[i],cycles[i],t0_model[i],decayed,legend,marker_type[i],color_1,line_style[i],title="")
	
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


		plt.figure(6)



	def set_triple_isotope_plot(self,runs=[],xiso=['Zr',96,'Zr',94],yiso=['Zr',92,'Zr',94],graintype=['sic'],C_star_only=True,dcycle=100,USEEPP_path='USEEPP',sens_runs=True):
		'''
        		xiso:       give isotopes as ['Fe',57,'Fe',56]. x axis, yaxis
	


		'''

		import os
		pwd=os.getcwd()
		iniabufiles=[]			
		
                sefiles=[]
                legend=[]
                HDF5_surf=[]
                extra_label=[]
                if len(runs) ==0:
                        HDF5_surf=self.runs_H5_surf
                        runs=self.run_dirs_name
                        extra_label=self.extra_label
                else:
                        for i in range(len(self.run_dirs_name)):
                                if self.run_dirs_name[i] in runs:
                                        HDF5_surf.append(self.runs_H5_surf[i])
                                        extra_label.append(self.extra_label[i])
		if sens_runs==True:
			import glob
			import os
			color=['r','r','b','b','k','g','k','g','k','b','r','g','g','b']
			marker_type=['o','D','s','p','o','D','s','p','p','o','p','o','s','p']
			line_style=['--','-','-.',':']
			first_plot=True
			for i in range(len(HDF5_surf)):
				sefiles=se(HDF5_surf[i])
                        	mass=sefiles.get("mini")[0]
                        	z=sefiles.get("zini")[0]
	                        new_notation='{:.1E}'.format(float(z))
        	                iniabufile="iniab"+new_notation+"GN93.ppn"     #pwd+"/"+USEEPP_path+'/'+"iniab"+new_notation+"GN93.ppn"
                        	legend=str(mass)+"M$_{\odot}$, "+str(z)+"Z$_{\odot}$" #+", "+extra_label[i]
				if first_plot==True:
					sefiles.plot_isoratios(xiso=xiso,yiso=yiso,graintype=graintype,tp_finding=False,C_star_only=C_star_only,deltax=True,deltay=True,logx=False,logy=False,title=None,legend=True,extra_legend=legend,dcycle=dcycle,errbar=True,iniabufile=iniabufile,plt_symb=marker_type[i],plt_col=color[i])		
					first_plot=False			
				else:
					print ""				
					sefiles.plot_isoratios(xiso=xiso,yiso=yiso,tp_finding=False,C_star_only=C_star_only,deltax=True,deltay=True,logx=False,logy=False,title=None,legend=True,extra_legend=legend,dcycle=dcycle,errbar=True,iniabufile=iniabufile,plt_symb=marker_type[i],plt_col=color[i])		


		#plot_isoratios(self,xiso,yiso,graintype=None,tp_finding=False,deltax=True,deltay=True,logx=False,logy=False,title=None,legend=True,dcycle=500,errbar=True,iniabufile='iniab2.0E-02GN93.ppn'):





	def write_prod_fact_stellarwinds(self,cyclerange=[5000,8000,1000],isotopes=["H-1","H-2","H-3","He-4"], label="fig:testlabel",caption="",table_header=[],file_1="test",table_size="normal",ascii_1=True,latex=True):

		'''
			
		'''
		production_factors=[]
		for i in range(len(self.runs_H5_surf)):
                	sefiles=se(self.runs_H5_surf[i])
			if len(cyclerange)==0:
				cyclestart=sefiles.se.cycles[0]
				cycleend=sefiles.se.cycles[-1]
				cyclesparse=500
			else:
				cyclestart=cyclerange[0]
				cycleend=cyclerange[1]
				cyclesparse=cyclerange[2]
			mass=sefiles.get("mini")[0]
		    	iso_folded,yield_folded,yield_unfolded,iso_production_factor_name,production_factor=self.weighted_yields_production_factor(sefiles=sefiles,cyclestart=cyclestart,cycleend=cycleend,sparse=cyclesparse,isotopes=isotopes,star_mass=mass,label=True,legend="",color="r",title="")

			production_factors.append(production_factor)
		cols=[iso_production_factor_name]+production_factors
		if latex==True:
			self.write_latex_table(header_file="Isotopic production factor stellar winds",cols=cols,label=label,caption=caption,table_header=table_header,latexfile=file_1,table_size=table_size)
		if ascii_1==True:
			self.write_ascii_table(header_file="Isotopic production factor stellar winds",cols=cols,txtfile=file_1)
	



	def write_latex_table(self,row_type="specie",header_file="Table element",cols=[[],[]],label="fig:testlabel",caption="The element abun",table_header=[],latexfile="test.tex",table_size="normal"):


		'''

			cols - 2d array containing rows and cols as shown below
                        header - header of table in latex
			label - label of the array, for reference in latex
			caption - latex-typical description of plot
			header_latexfile - header in latexfile 
			table_size - if "small" , set latex table size to small, else else normall
			
			writes tables of the form:

			cols[0][0]	cols[1][0]	... 
			cols[0][1]	cols[1][1]	...
			..		...	
			..		...		...


		'''

		if latexfile[-4:] != ".tex":
			latexfile=latexfile+".tex"
		file=open(latexfile,"w")
		#write header
		file.write("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"+"\n")
		file.write("%	"+header_file+"\n")
                file.write("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"+"\n")
		file.write("\n")
		file.write("\n")
		file.write("\\begin{table}"+"\n")
		file.write("\\begin{center}"+"\n")
		if table_size=="small":
			file.write("\\tabletypesize{\\small}")
		file.write("\\caption{"+caption+"}"+"\n")
		file.write("\\begin{tabular}{"+("c"*(1+len(self.runs_H5_out)))+"}"+"\n")
		file.write("\\hline"+"\n")
		#write content
		if len(table_header)==0:
			for i in range(len(self.runs_H5_out)):
                                sefiles=se(self.runs_H5_out[i])
                                mass=sefiles.get("mini")[0]
                                z=sefiles.get("zini")[0]
                                new_notation='{:.1E}'.format(float(z))
                                legend=str(mass)+"M$_{\odot}$ Z= "+str(z)
				table_header.append(legend)
		#write table header
		header_write=row_type
		for i in range(0,len(table_header)):
			header_write+=" & "+table_header[i]
		file.write(header_write+"\\\\"+"\n")	
		file.write("\\hline"+"\n")
		#write columns
		for i in range(len(cols[0])):
			line=cols[0][i]
			for j in range(1,len(cols)):
				a=cols[j][i]
				line+= " & "+str(round(float(str(a).split("e")[0]),3))+"E"+str(a).split("e")[1]
			file.write(line+"\\\\"+"\n")			

		#file end
		file.write("\\hline"+"\n")
		file.write("\\noalign{\\smallskip}"+"\n")
		file.write("\\hline"+"\n")
		file.write("\\end{tabular}"+"\n")
		file.write("\\label{"+label+"}"+"\n")
		file.write("\\end{center}"+"\n")
		file.write("\\end{table}"+"\n")
		
		

        def write_ascii_table(self,row_type="specie",header_file="Table element",cols=[[],[]],table_header=[],txtfile="test.txt",table_size="normal"):


                '''

                        cols - 2d array containing rows and cols as shown below
                        header - header of table in latex
                        label - label of the array, for reference in latex
                        caption - latex-typical description of plot
                        header_latexfile - header in latexfile 
                        table_size - if "small" , set latex table size to small, else else normall
                        
                        writes tables of the form:

                        cols[0][0]      cols[1][0]      ... 
                        cols[0][1]      cols[1][1]      ...
                        ..              ...     
                        ..              ...             ...


                '''
                if txtfile[-4:] != ".txt":
                        txtfile=txtfile+".txt"


                file_1=open(txtfile,"w")
                #write header
                file_1.write("#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"+"\n")
                file_1.write("#   "+header_file+"\n")
                file_1.write("#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"+"\n")
                file_1.write("\n")
                file_1.write("\n")
                #write content
                if len(table_header)==0:
                        for i in range(len(self.runs_H5_out)):
                                sefiles=se(self.runs_H5_out[i])
                                mass=sefiles.get("mini")[0]
                                z=sefiles.get("zini")[0]
                                new_notation='{:.1E}'.format(float(z))
                                legend=str(mass)+"Msun Z= "+str(z)
                                table_header.append(legend)
                #write table header
                header_write=row_type
                for i in range(0,len(table_header)):
                        header_write+="  "+table_header[i]
                file_1.write(header_write+"\n")
                file_1.write("\n")
                #write columns
                for i in range(len(cols[0])):
                        line=cols[0][i]
                        for j in range(1,len(cols)):
                                a=cols[j][i]
                                line+= "	"+str(round(float(str(a).split("e")[0]),3))+"E"+str(a).split("e")[1]
                        file_1.write(line+"\n")




	#########################################
        #####Internal or single run functions for MPPNP ####
	#########################################



	def weighted_yields_production_factor(self,sefiles,cyclestart=11000,cycleend=12000,sparse=100,isotopes=["H-1","H-2","H-3","He-4"],star_mass=1.65,label=True,legend="",color="r",title=" -Final wind yields - isotopic composition folded with IMF"):
		'''
				Uses H5_surf
		      This function returns the isotope name with corresponding wind yields, folded, unfolded. And isotope name with corresponding production factors of stable isotopes
        	is
		'''
	
                import utils as u
                import re
                import nugridse as mp
	
		X_i, E_i =sefiles.windyields(cyclestart, cycleend, sparse)
		print "#################",X_i[:12],sefiles.se.isotopes[0:12],"+#################"
		###calculate isotopes
		#isotopes=['Fe-56','Co-60','Ni-61','Cu-65','Zn-67',]
		#isotopes['Ga-71','Ge-73','As-75','Se-77','Br-82','Kr-84','Rb-87','Sr-88','Y-89','Zr-92','Nb-94','Zr-96','Mo-96','Ba-137','Ba-138','La-139','Ce-140','Nd-142','Sm-144','Pb-206']					
		#206
		iso_abu=[]
		production_factors=[]
		mass=[]
		name=[]
		initial_abu=[]
		isotope_names = sefiles.se.isotopes
		sefiles.mass_frac = [X_i] #double array important for decay
	
		##get initial abundance for specific isotopes
		for i in range(len(isotopes)):
			initial_abu.append(sefiles.get(0,isotopes[i]))
			
		#convert isotope name scheme
		u.convert_specie_naming_from_h5_to_ppn(isotope_names) #-define u.spe
		names_ppn_world = u.spe 
		number_names_ppn_world = u.n_array
		u.define_zip_index_for_species(names_ppn_world,number_names_ppn_world) #--define u.cl!
		
		#let species decay
		u.stable_specie()
		sefiles.decay([X_i]) #double array is important for decay()
		X_i_decayed=mp.decayed_multi_d[0]
		
		#convert specified isotopes in correct name scheme
		u.convert_specie_naming_from_h5_to_ppn(isotopes) #-define u.spe
		specific_isotopes_ppn_world = u.spe
		
		for i in range(len(specific_isotopes_ppn_world)):
			specific_isotopes_ppn_world[i]=specific_isotopes_ppn_world[i].upper()
		
		isotope_name_prod_factor=[]
		isotope_name_yield=[]
		for j in range(len(u.stable)):
			if u.stable[j] in specific_isotopes_ppn_world:
				idx=specific_isotopes_ppn_world.index(u.stable[j])
				print u.stable[j]
				print specific_isotopes_ppn_world
				print idx
				iso_abu.append(X_i_decayed[j])
				isotope_name_yield.append(isotopes[idx])
				print X_i_decayed[j],initial_abu[idx]
				if initial_abu[idx] <10e-20:
					print "isotope "+u.stable[j]+" not in initial abundance, and production factor wont be calculated"
				else:
					production_factors.append(X_i_decayed[j]/initial_abu[idx])
					isotope_name_prod_factor.append(isotopes[idx])
				mass.append(int(re.findall(r'\d+',u.stable[j])[0]))
				name.append(u.stable[j])
		###Folding and plotting###			
		#Salpeter mass function,Kroupa 2001:(lazyness..)
		#f(M)=M**-2.3
		f_m=(star_mass**(-2.3))
		iso_abu_IMF=[]
		for i in range(len(iso_abu)):
			iso_abu_IMF.append(f_m*iso_abu[i])
		print "RETURN: ",isotope_name_yield,iso_abu_IMF	
		return isotope_name_yield,iso_abu_IMF,iso_abu,isotope_name_prod_factor, production_factors #isotope name, folded, unfolded yields, and isotope name for production factor and production factor - all stable

 
	def weighted_yields(self,sefiles,cyclestart=11000,cycleend=12000,sparse=100,isotopes=["H-1","H-2","H-3","He-4"],star_mass=1.65,label=True,legend="",color="r",title=" -Final wind yields - isotopic composition folded with IMF",plot_fig=True):
		'''
				Uses H5_surf
	  	    This function returns the wind yields and ejected masses for stable species.  The command
     		   is
		'''
	

                import utils as u
                import re
                import nugridse as mp

	
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
		iso_IMF=[]
		for i in range(len(plotiso)):
			iso_IMF.append(f_m*plotiso[i])
		if plot_fig == True:	
			plt.plot(mass,iso_IMF,marker='*',markersize=8,mfc=color,linestyle='None',label=legend)
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
					plt.annotate(name[j], xytext = (0, 10),textcoords = 'offset points' ,xy=(mass[j], iso_IMF[j]))
			mina=mass[0]-4 	
			if mina<0:
				mina=0
			plt.xlim(mina,mass[-1]+4)
			plt.legend()
		return iso_IMF, plotiso #folded and unfolded - both stable
	


	def plot_abu_atmasscoord(self,mp,sefiles,cycle,masscell,isotopes_input,x_charge,mass_number_range,decayed,filename_solar_norm,elem_norm,yax_log,label,legend,color,title):
		import utils as u

		'''
			if mass_number_range greater 0, ignore isotopes input
	
		'''
		
		####decide which isotopes to process
		isotopes_stable=[]	
		isotopes=[]
		#if mass number range is needed
		if len(mass_number_range)>0:
			#get all isotopes in mass number range
			if mass_number_range[0] == 0:
				mass_number_min=sefiles.se.isotopes[0]
			else:
				mass_number_min=mass_number_range[0]
			if mass_number_range[1] == 0:
				mass_number_max=sefiles.se.isotopes[-1]
			else:
				mass_number_max=mass_number_range[1]	
			for i in range(len(sefiles.se.isotopes)):
				iso=sefiles.se.isotopes[i]
				#print iso.split("-")[1]
				#print mass_number_min 
				#print mass_number_max
				iso_number= iso.split("-")[1]
				iso_number= re.findall ('[\d ]+',iso_number)[0]
				if ( (int(iso_number)>=mass_number_min) and (int(iso_number)<=mass_number_max) ):
					if (is_stable(iso) == 't'):
						isotopes_stable.append(iso)
					isotopes.append(iso)
		else:
			isotopes=isotopes_input
			for i in range(len(sefiles.se.isotopes)):
				iso=sefiles.se.isotopes[i]
                                #print iso.split("-")[1]
                                #print mass_number_min 
                                #print mass_number_max
                                iso_number= iso.split("-")[1]
                                iso_number= re.findall ('[\d ]+',iso_number)[0]
                                if (is_stable(iso) == 't'):
                                	isotopes_stable.append(iso)
		
		print isotopes

		#############get isotopic distribution
		u.convert_specie_naming_from_h5_to_ppn(sefiles.se.isotopes)
		masses_for_this_cycle = sefiles.se.get(cycle,'mass') #<- find mass index you want
		idx=0	
		#print sefiles.isotopic_production_factors
		idx=np.abs(np.array(masses_for_this_cycle)-masscell).argmin()
		print "closest value: ", masses_for_this_cycle[idx]	
		plotiso=[]
		mass=[]
		plotiso_norm=[]
		iso_abu=[]
		iniabu_array=[i]
		name=[]
		
		mass=sefiles.se.A
		charge=sefiles.se.Z
		#sefiles.stable_isotope_list[isotope] <> get isotopes you want
		sefiles.iso_index=[]
		iniabu_norm=0	
		for j in range(len(sefiles.se.isotopes)):
			if sefiles.se.isotopes[j] in isotopes:
				name.append(sefiles.se.isotopes[j])
				if decayed == False:
					iso_abu.append(sefiles.get(cycle,sefiles.se.isotopes[j])[idx])
				sefiles.iso_index.append(j)
			if len(elem_norm)==1:
				if elem_norm in sefiles.se.isotopes[j]:
					iniabu_norm+=sefiles.get(cycle,sefiles.se.isotopes[j])[0]
					print "iniabu_norm....",iniabu_norm

		if decayed ==True:
               		########isotope decay, new abundances#######
			#isomers are excluded for now, in convert_specie naming...last ones 5
			#to decay whole network use all isotopes, change iso_abu
			iso_abu=sefiles.get(cycle,'iso_massf')[idx]
                	isotope_names = sefiles.se.isotopes
                	sefiles.mass_frac = [iso_abu] #double array important for decay

                	u.convert_specie_naming_from_h5_to_ppn(isotope_names) #-define u.spe
                	names_ppn_world_sefiles = u.spe
                	number_names_ppn_world = u.n_array
                	u.define_zip_index_for_species(names_ppn_world_sefiles,number_names_ppn_world) #--define u.cl!

                	#let species decay
                	u.stable_specie()
			print name
			print u.stable
			print "lentest:",len(u.stable),len(mp.decayed_multi_d[0])
			#print iso_abu
                	sefiles.decay([iso_abu]) #double array is important for decay()
                	iso_abu_decayed=mp.decayed_multi_d[0]
			print mp.decayed_multi_d


		#convert specified isotopes in correct name scheme
		u.convert_specie_naming_from_h5_to_ppn(isotopes) #-define u.spe
		names_ppn_world_chosen_isotopes = u.spe
	
		for i in range(len(names_ppn_world_chosen_isotopes)):
			names_ppn_world_chosen_isotopes[i]=names_ppn_world_chosen_isotopes[i].upper()
		
		#for normalization
		#u.solar(filename_solar_norm,solar_factor=1)


		mass_specific_isotopes=[]
		charge_specific_isotopes=[]
		label_specific_isotopes=[]
		###########Calculate production factor
		if decayed == True:
			for j in range(len(u.stable)):
				if u.stable[j] in names_ppn_world_chosen_isotopes:
					#if sefiles isotope in both name schemes have same order, than following is correct
					iso_idx=names_ppn_world_sefiles.index(u.stable[j].capitalize())
					#print u.stable[j].capitalize(), sefiles.se.isotopes[iso_idx]
					if len(elem_norm)==0:
						iniabu_norm=sefiles.get(0,sefiles.se.isotopes[iso_idx])[0]
						print "this is not "
                                        print iniabu_norm
					if iniabu_norm <1e-20:
                                                	print "Warning: Isotope "+sefiles.se.isotopes[i]+" not in initial abundance file, no plotting"
                                        else: 
						#norm=iniabu_array[iso_idx]		#u.solar_abundance[isotopes_ppn_world]
						plotiso.append(iso_abu_decayed[j]/iniabu_norm)
						mass_specific_isotopes.append(mass[iso_idx])
						label_specific_isotopes.append(sefiles.se.isotopes[iso_idx])
		else:			
			j=0
			for i in range(len(sefiles.se.isotopes)):
				if sefiles.se.isotopes[i] in isotopes:
					print sefiles.se.isotopes[i],isotopes
					if len(elem_norm)==0:
						print "count"+str(i)
						iniabu_norm=1
						print "print info:",sefiles.get(0,sefiles.se.isotopes[i])[0]
						print "aaaa",iniabu_norm
					print "aaaaaaaaaa",iniabu_norm
					#if len(iniabu_norm)>1:
					#	iniabu_norm=iniabu_norm[0]
					if iniabu_norm <1e-20:
						print "Warning: Isotope "+sefiles.se.isotopes[i]+" not in initial abundance file, no plotting"	
					else:	
						plotiso.append(iso_abu[j]/iniabu_norm)
						mass_specific_isotopes.append(mass[i])
						charge_specific_isotopes.append(charge[i])
						label_specific_isotopes.append(sefiles.se.isotopes[i])
					j+=1
		
		####Plotting


		fig=figure(0)
		if yax_log==True:
                	ax = fig.add_subplot(1,1,1)
                        ax.set_yscale('log')

		if x_charge==False:
			plt.plot(mass_specific_isotopes,plotiso,'*',label=legend,marker='*',markersize=8,mfc=color,linestyle='None')

	                if label ==True:
        	                for j in range(len(plotiso)):
                	                plt.annotate(label_specific_isotopes[j], xytext = (0, 10),textcoords = 'offset points' ,xy=(mass_specific_isotopes[j], plotiso[j]))


		else:
			plt.plot(charge_specific_isotopes,plotiso,'*',label=legend,marker='*',markersize=8,mfc=color,linestyle='None')
			if label ==True:
				for j in range(len(plotiso)):
					plt.annotate(label_specific_isotopes[j], xytext = (0, 10),textcoords = 'offset points' ,xy=(charge_specific_isotopes[j], plotiso[j]))

                plt.rcParams.update({'font.size': 16})
                plt.rc('xtick', labelsize=16)
                plt.rc('ytick', labelsize=16)
		plt.minorticks_on()
		plt.xlim(mass_specific_isotopes[0],mass_specific_isotopes[-1])
		plt.legend()
		plt.title(title)
 		plt.ylabel("$X_{i}/X_{ini}$",fontsize=20)
		if x_charge==False:
			plt.xlabel("Mass Number (A)",fontsize=20)
		else:
			plt.xlabel("Charge (Z)",fontsize=20)


	def set_ba(self):

		for i in range(len(self.runs_H5_surf)):
			sefiles=se(self.runs_H5_surf[i])
			y=[]
			mass=sefiles.get("mini")[0]
                        z=sefiles.get("zini")[0]
                        legend=str(mass)+"M$_{\odot}$ Z= "+str(z)
			for j in arange(0,len(sefiles.se.cycles)-1000,500):
    				y.append(sefiles.get(j,"elem_massf","Ba"))
			plt.plot(arange(0,len(sefiles.se.cycles)-1000,500),y,label=legend)



	
	def surface_plots(self,HDF5surf_dir,cycles,t0_model,decayed=True,legend="",marker_type='o',color="r",line_style="-",title="",hs=["Ba","La","Nd","Sm"],ls=["Sr","Y","Zr"]): 
	
		'''surface hdf5
	
			ls: Sr,Y,Zr  - elements
			hs: Ba,La,Nd,Sm -elements
		
			isotopes,elements,HDF5_dir,cycles,legend,title,x_range=[],y_range=[]):
		if cycle longer than 3 assume special sequence	
	
		plots [hs/ls] vs star age
	plots  [ls/Fe] vs [hs/ls] 
		[Rb/Sr] [hs/ls]
plots Mg25/Mg24 vs Mg26/Mg24
          plt.xlabel("$^{96}$Zr/$^{94}$Zr")
                plt.ylabel("$^{152}$Gd/$^{154}$Gd")
	'''	
		
		import utils as u
		import nugridse as mp	
                if len(cycles)>3:
                        cycle_range=cycles
			t0_model=cycle_range[0]
                else:
                        cycle_range=np.arange(cycles[0],cycles[1],cycles[2])
		
		sefiles=se(HDF5surf_dir)

		#if len(t0_model)==0:
	#		t0_model=cycles[0]


		#2 choices to plot
		#modelrange=[10,30,10]
		#t0_model=0
		t0_time=sefiles.se.get(t0_model,"age")
		t_begin_model=sefiles.se.get(0,"model_number")
		#t0_model= t0_model-t_begin_model
		starage=[]
		ls_abu=[]
		hs_abu=[]
		fe_abu=[]
		rb_abu=[]####n-density + efficiency of N22(a,n)
		sr_abu=[]
		ba_abu=[]
		eu_abu=[]
		ls_abu_ini=[]
		hs_abu_ini=[]
		ba_abu_ini=[]
		eu_abu_ini=[]
		fe_abu_ini=[]
		quot_abu_mg_26_24=[]
		quot_abu_mg_25_24=[]
		quot_abu_gd_152_154=[]
		quot_abu_zr_96_94=[]		
		hs_isotopes=[]
		ls_isotopes=[]
		rb_isotopes=[]
		fe_isotopes=[]
		sr_isotopes=[]
		ba_isotopes=[]
		eu_isotopes=[]
		mg_isotopes=['Mg-24','Mg-25','Mg-26']
		gd_isotopes=['Gd-152','Gd-154']
		zr_isotopes=['Zr-94','Zr-96']
		for iso in sefiles.se.isotopes:
			for i in range(len(ls)):
				if ls[i] in iso and (is_stable(iso) == 't'):
					ls_isotopes.append(iso)
			for i in range(len(hs)):
				if hs[i] in iso and (is_stable(iso) == 't'):
					hs_isotopes.append(iso)
			if "Rb" in iso:
					rb_isotopes.append(iso)
			if "Fe" in iso and (is_stable(iso) == 't'):
					fe_isotopes.append(iso)
			if "Sr" in iso and (is_stable(iso) == 't'):
					sr_isotopes.append(iso)
                        if "Ba" in iso and (is_stable(iso) == 't'):
                                        ba_isotopes.append(iso) 
                        if "Eu" in iso and (is_stable(iso) == 't'):
                                        eu_isotopes.append(iso) 	
		print "0000000000000000000"
		#if decayed ==True:
              #  	u.convert_specie_naming_from_h5_to_ppn(ls_isotopes) #-define u.spe
                 #       ls_isotopes = u.spe
	          #      for i in range(len(ls_isotopes)):
        	   #             ls_isotopes[i]=ls_isotopes[i].upper()
                    #    u.convert_specie_naming_from_h5_to_ppn(hs_isotopes) #-define u.spe
                     #   hs_isotopes = u.spe
                      #  for i in range(len(hs_isotopes)):
                       #         hs_isotopes[i]=hs_isotopes[i].upper()
            #           	u.convert_specie_naming_from_h5_to_ppn(rb_isotopes) #-define u.spe
             #           rb_isotopes = u.spe
              #          for i in range(len(rb_isotopes)):
               ##                 rb_isotopes[i]=rb_isotopes[i].upper()
		#	u.convert_specie_naming_from_h5_to_ppn(sr_isotopes) #-define u.spe
                 #       sr_isotopes = u.spe
                  #      for i in range(len(sr_isotopes)):
                   #             sr_isotopes[i]=sr_isotopes[i].upper()
                    #    u.convert_specie_naming_from_h5_to_ppn(ba_isotopes) #-define u.spe
                     #   ba_isotopes = u.spe
       #                 for i in range(len(ba_isotopes)):
        #                        ba_isotopes[i]=ba_isotopes[i].upper()
         #               u.convert_specie_naming_from_h5_to_ppn(eu_isotopes) #-define u.spe
         #               eu_isotopes = u.spe
          #              for i in range(len(eu_isotopes)):
           #                     eu_isotopes[i]=eu_isotopes[i].upper()
            #            u.convert_specie_naming_from_h5_to_ppn(fe_isotopes) #-define u.spe
             #           fe_isotopes = u.spe
              #          for i in range(len(fe_isotopes)):
               #                 fe_isotopes[i]=fe_isotopes[i].upper()
#
		print hs_isotopes,ls_isotopes,rb_isotopes,fe_isotopes,sr_isotopes,ba_isotopes,eu_isotopes
		#ini_abu=sefiles.get(0,"iso_massf_decay")
		a=u.iniabu('iniab2.0E-02GN93.ppn')#'../../frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn')
		for cycle in cycle_range:
			reload(mp)
			#print cycle
			starage.append(sefiles.se.get(cycle,"age")-t0_time)
			ls_abu_1=0
			hs_abu_1=0
			rb_abu_1=0
			fe_abu_1=0
			sr_abu_1=0
			ba_abu_1=0
			eu_abu_1=0
			ls_abu_ini_1=0
			hs_abu_ini_1=0
			ba_abu_ini_1=0
			eu_abu_ini_1=0
			fe_abu_ini_1=0
			print "11111111111111"	
			isotopes_complete=sefiles.se.isotopes
			if decayed ==True:
				isotope_abundances=sefiles.get(cycle,"iso_massf_decay")
			else:
				isotope_abundances=sefiles.get(cycle,"iso_massf")
                                ########isotope decay, new abundances#######
                                #isomers are excluded for now, in convert_specie naming...last ones 5
                                #to decay whole network use all isotopes, change iso_abu

	                #if decayed ==True:
			#	print "111111111111111233124324"
                        #	isotope_names = sefiles.se.isotopes
                        #	sefiles.mass_frac = [isotope_abundances] #double array important for decay
                        #	u.convert_specie_naming_from_h5_to_ppn(isotope_names) #-define u.spe
                        #	names_ppn_world_sefiles = u.spe
                        #	number_names_ppn_world = u.n_array
                        #	u.define_zip_index_for_species(names_ppn_world_sefiles,number_names_ppn_world) #--define u.cl!
                        #	#let species decay
                        #	u.stable_specie()
                        #	sefiles.decay([isotope_abundances]) #double array is important for decay()
                        #	iso_abu_decayed=mp.decayed_multi_d
			#	isotopes_complete=u.stable
			#	isotope_abundances=iso_abu_decayed[0]
			#	print "  "
			#	print isotopes_complete
			#	print len(isotope_abundances)		

			u.convert_specie_naming_from_h5_to_ppn(isotopes_complete)
			names_ppn_world=u.spe
			for i in range(len(names_ppn_world)):
				print "##############",len(isotopes_complete),len(names_ppn_world)
				names=names_ppn_world[i].lower()
				if isotopes_complete[i] in ls_isotopes:
					ls_abu_1+=(sefiles.get(cycle,"iso_massf_decay",isotopes_complete[i]))
					ls_abu_ini_1+=a.habu[names]
				if isotopes_complete[i] in hs_isotopes:
					hs_abu_1+=sefiles.get(cycle,"iso_massf_decay",isotopes_complete[i])
					hs_abu_ini_1+=a.habu[names]
				if isotopes_complete[i] in rb_isotopes:
					rb_abu_1+=sefiles.get(cycle,"iso_massf_decay",isotopes_complete[i])
				if isotopes_complete[i] in sr_isotopes:
					sr_abu_1+=sefiles.get(cycle,"iso_massf_decay",isotopes_complete[i])
                                if isotopes_complete[i] in ba_isotopes:
					#name=
                                        ba_abu_1+=sefiles.get(cycle,"iso_massf_decay",isotopes_complete[i])
					ba_abu_ini_1+=a.habu[names]
                                if isotopes_complete[i] in eu_isotopes:

                                        eu_abu_1+=sefiles.get(cycle,"iso_massf_decay",isotopes_complete[i])
					eu_abu_ini_1+=a.habu[names]
				if isotopes_complete[i] in fe_isotopes:
					#name=
					print isotopes_complete[i],sefiles.get(cycle,"iso_massf_decay",isotopes_complete[i])
					fe_abu_1+=sefiles.get(cycle,"iso_massf_decay",isotopes_complete[i])
					fe_abu_ini_1+=a.habu[names]
			#eu_abu_1=sefiles.get(cycle,"elem_massf_decay","Eu")
			print fe_abu_1 	
			print "fe: ",sefiles.get(cycle,"elem_massf","Fe")
			print "fe: ",sefiles.get(cycle,"elem_massf_decay","Fe")
			print fe_abu_1
			print ba_abu_1
			print "ba :",sefiles.get(cycle,"elem_massf_decay","Ba")	
			print ba_abu_1
			#create hs and ls as mean values of the elements		
			hs_abu.append(hs_abu_1/len(hs))
			ls_abu.append(ls_abu_1/len(ls))
			fe_abu.append(fe_abu_1)
			rb_abu.append(rb_abu_1)
			sr_abu.append(sr_abu_1)	
			ba_abu.append(ba_abu_1)
			eu_abu.append(eu_abu_1)	
			hs_abu_ini.append(hs_abu_ini_1)
			ls_abu_ini.append(ls_abu_ini_1)
			ba_abu_ini.append(ba_abu_ini_1)
			eu_abu_ini.append(eu_abu_ini_1)
			fe_abu_ini.append(fe_abu_ini_1)
			#quot_abu_mg_26_24.append( sefiles.get(cycle,"iso_massf",mg_isotopes[2])/sefiles.get(cycle,"iso_massf",mg_isotopes[0]) )
			#quot_abu_mg_25_24.append( sefiles.get(cycle,"iso_massf",mg_isotopes[1])/sefiles.get(cycle,"iso_massf",mg_isotopes[0] ))
			#quot_abu_gd_152_154.append( sefiles.get(cycle,"iso_massf",gd_isotopes[0])/sefiles.get(cycle,"iso_massf",gd_isotopes[1])) 
			#quot_abu_zr_96_94.append( sefiles.get(cycle,"iso_massf",zr_isotopes[1])/sefiles.get(cycle,"iso_massf",zr_isotopes[0]) )
	
		print "hsabu:" 
		print hs_abu
		print "lsabu"
		print ls_abu
		print "feabu"
		print fe_abu
		print "rbabu"
		print rb_abu
		print "srabu"
		print sr_abu
		print "babu"
		print ba_abu
		print "euabu"
		print eu_abu
		print "starage"
		print starage
		
		#Plotting different figures, beware of square bracket notation
		plt.figure(1)
		plt.plot(np.array(starage),np.log10(np.array(hs_abu)/np.array(ls_abu)),marker='*',markersize=8,mfc=color,linestyle=line_style,label=legend)
		plt.xlabel("$star age [t-t_0]$")
		plt.ylabel("[hs/ls]")
		plt.legend()
		plt.title(title)
		plt.draw()
		#Paper 1 plot
		print "plot2"
		plt.figure(2);plt.plot(np.log10(np.array(hs_abu)/np.array(ls_abu)),np.log10(np.array(ls_abu)/np.array(fe_abu)),c=color,linestyle=line_style,label=legend)  #,marker='*',markersize=8,
		plt.xlabel("[hs/ls]");plt.ylabel("[ls/Fe]")
		plt.legend()
		plt.title(title)
		plt.draw()	
		print "print3"
		#Paper 1 plot
		plt.figure(3);plt.plot(np.log10(np.array(hs_abu)/np.array(ls_abu)),np.log10(np.array(rb_abu)/np.array(sr_abu)),c=color,linestyle=line_style,label=legend)
		plt.xlabel("[hs/ls]")
		plt.ylabel("[Rb/Sr]")
		plt.legend()
		plt.title(title)
		plt.draw()
		print "plot4"
		#Paper 1 plot
		#plt.figure(4);plt.plot(quot_abu_mg_25_24,quot_abu_mg_26_24,c=color,linestyle=line_style,label=legend)
		plt.xlabel("$^{25}$Mg/$^{24}$Mg")
		plt.ylabel("$^{26}$Mg/$^{24}$Mg")
		plt.legend()
		plt.title(title)
		plt.draw()
	
		#Paper 1 plot
		#plt.figure(5);plt.plot(quot_abu_zr_96_94,quot_abu_gd_152_154,c=color,linestyle=line_style,label=legend)
		plt.xlabel("$^{96}$Zr/$^{94}$Zr")
		plt.ylabel("$^{152}$Gd/$^{154}$Gd")
		plt.legend()
		plt.title(title)
		plt.draw()
		#new
                plt.figure(6)
                plt.plot(np.array(starage),np.log10( np.array(hs_abu)*np.array(ls_abu_ini)/(np.array(ls_abu)*np.array(hs_abu_ini))),marker=marker_type,markersize=12,mfc=color,linestyle=line_style,label=legend)
                plt.xlabel("$t-t_0 [yr]$",fontsize=20)
                plt.ylabel("[hs/ls]",fontsize=20)
                plt.legend()
                plt.title(title)
                plt.draw()
                plt.rcParams.update({'font.size': 16})
                plt.rc('xtick', labelsize=16)
                plt.rc('ytick', labelsize=16)
                plt.minorticks_on()

                plt.rcParams.update({'font.size': 16})
                plt.rc('xtick', labelsize=16)
                plt.rc('ytick', labelsize=16)
                plt.minorticks_on()


                #new
                plt.figure(7)
                plt.plot( np.log10( np.array(eu_abu)*np.array(fe_abu_ini) /(np.array(fe_abu)*np.array(eu_abu_ini)  )), np.log10( np.array(ba_abu)*np.array(fe_abu_ini)/(np.array(fe_abu)*np.array(ba_abu_ini))),marker=marker_type,markersize=12,mfc=color,linestyle=line_style,label=legend)
                plt.xlabel("[Eu/Fe]",fontsize=20)
                plt.ylabel("[Ba/Fe]",fontsize=20)
                plt.legend()
                plt.title(title)
                plt.draw()
                plt.rcParams.update({'font.size': 16})
                plt.rc('xtick', labelsize=16)
                plt.rc('ytick', labelsize=16)
                plt.minorticks_on()




		
		return [starage,hs_abu,ls_abu,fe_abu,rb_abu,sr_abu]


	def pocket_composition(self,mp,sefiles,cycle,massbot,masstop,isotopes,label,legend,color,title):  ###hdf5out
        	#mass_range    - required to plot data in a certain mass range. Needed for read_iso_abund_marco
        	#cycle         - which cycle from the h5 file?. Needed for read_iso_abund_marco
        	#stable        - logic if want to plot only stable or not.  
        	#i_decay       - if = 1 I plot not decayed, if = 2 I plot decayed. Make sense only if stable is true.'''
		mass_range=[massbot,masstop]
		sefiles.average_iso_abund_marco(mass_range,cycle,stable=False,i_decay=1)
		mass=[]
		plotiso=[]
		startyields=[]
		plotiso_massfrac=[]
		for i in range(len(isotopes)):
			startyields.append(sefiles.get(sefiles.se.cycles[0],isotopes[i])[0])
		for j in range(len(sefiles.se.isotopes)):
			if sefiles.se.isotopes[j] in isotopes:
				plotiso.append(mp.average_mass_frac[j])
				mass.append(sefiles.se.A[j])
		for i in range(len(isotopes)):
			plotiso_massfrac.append(plotiso[i]/startyields[i])
		plt.plot(mass,plotiso_massfrac,marker='*',markersize=8,mfc=color,linestyle='None',label=legend)
		#plt.plot(mass,plotiso_massfrac,marker='*')
		if label ==True:
			for j in range(len(isotopes)):
				plt.annotate(isotopes[j], xytext = (0, 10),textcoords = 'offset points' ,xy=(mass[j], plotiso_massfrac[j]))
		mass=np.array(mass)
		plt.xlim(mass.min()-4,mass.max()+4)
		plt.legend()
		plt.title(title)
		plt.xlabel("mass number")
		plt.ylabel("Isotopic production factors")


########################################################################################
############################ START OF MESA SPECIFIC CONTENT#############################
########################################################################################








class mesa_set(history_data):

	'''
		This class allows access to multiple MESA runs in the path dir.
		Instead of defining one path, an array of paths , mult_dir can be
		defined. In this case each path in multi_dir represents a run.

		rundir - dir with run directory
		multi_dir - specify dirs at different location
		extra_label - set labels of runs, else run dir name will be chosen	
	
		e.g.:
	
			setres=set.mesa_set(".",extra_label=["double-f","reference"])		
			
			or
				
			setres=set.mesa_set(multi_dir=["M3.00Z2.0e-02_1_re_double_f","M3.00Z2.0e-02_re_referencerun"],extra_label=["double-f","reference]])

		for paper2 agb analysis:
		mesaset=set.mesa_set(multi_dir=["M1.650Z0.0001","M2.000Z0.0001/","M3.000Z0.0001","M5.000Z0.0001/"])
		mesaset=set.mesa_set(multi_dir=["M1.650Z0.0010","M2.000Z0.0010/","M3.000Z0.0010","M5.000Z0.0010/"])	
		mesaset=set.mesa_set(multi_dir=["M1.650Z0.0060","M2.000Z0.0060/","M3.000Z0.0060","M5.000Z0.0060/"])

		Using the set_plot* functions you can instantly plot diagrams
		for all runs defined in mesa set

                	setres.set_plot_kipp_CO()

	'''

	def __init__(self,rundir='.',multi_dir=[],extra_label=[]):

		if len(multi_dir)==0:
			slist = os.listdir(rundir)
		else:
			slist = multi_dir
		self.runs_H5=[]
		self.run_LOGS=[]
		self.run_historydata=[]
		self.run_label=[]
		i=0
		for element in slist:
			if len(multi_dir)==0:
				run_path=rundir+"/"+element
			else:
				if multi_dir[i][0] == "/":
					run_path=multi_dir[i]
				else:
					run_path=os.getcwd()+"/"+multi_dir[i]	
			if os.path.isdir(run_path):
				run_path_1=glob.glob(run_path+"/*/*.h5")
				if len(run_path_1)>0:
					h5_dir=run_path_1[0].split("/")[-2]		
					self.runs_H5.append(run_path+"/"+h5_dir)
				else:
					print "Warning: h5 files are not available in "+run_path+"/"+h5_dir
				if (len(glob.glob(run_path+"/*/*.data"))>0) or (len(glob.glob(run_path+"/*/*.log"))>0):
					self.run_LOGS.append(run_path+"/LOGS")
					print "Read "+run_path
					self.run_historydata.append(history_data(run_path+"/LOGS"))
					if len(extra_label)>0:
						self.run_label.append(self.create_label(self.run_historydata[-1],extra_label[i]))
					else:
						self.run_label.append(self.create_label(self.run_historydata[-1],element))
					i+=1
				else:
					#if len(multi_dir)>=0:
					print "Error: not history.data or star.log file found in "+run_path		

	
	def multi_DUP(self,dirs=[],path=path,t0_model=[],h_core_mass=False,plot_fig=True):
		'''
			z1e-2:1.65-5:[13690,3120,3163,5306]
			z2e-2:1.65-5:[16033,6214,3388,5368]
			/rpod3/fherwig/SEE/data/set1.2/see_wind
			/rpod3/fherwig/SEE/data/set1.1/see_wind
			dirs=["M1.65Z2.0e-02/LOGS","M2.00Z2.0e-02/LOGS","M3.00Z2.0e-02/LOGS","M5.00Z2.0e-02/LOGS"]
			dirs=["M1.65Z1.0e-02/LOGS","M2.00Z1.0e-02/LOGS","M3.00Z1.0e-02/LOGS","M5.00Z1.0e-02/LOGS"]
			note: M3_1e-2 shows decrease in core mass in the end
		
			pap.multi_DUP(dirs=["M1.65Z1.0e-02","M2.00Z1.0e-02","M3.00Z1.0e-02","M5.00Z1.0e-02"],path="/rpod3/fherwig/SEE/data/set1.1/see_wind",t0_model=[13690,3120,3163,5306])
			pap.multi_DUP(dirs=["M1.65Z2.0e-02","M2.00Z2.0e-02","M3.00Z2.0e-02","M5.00Z2.0e-02"],path="/rpod3/fherwig/SEE/data/set1.2/see_wind",t0_model=[16033,6214,3388,5368])
		
		'''
		if len(dirs)==0:
			t0_model=self.set_find_first_TP()
		historydata=[]
		if (len(dirs)) == 0:
			dirs=self.run_LOGS
			historydata=self.run_historydata
		else:
			for i in range(len(dirs)):
				historydata.append(history_data(dirs[i]+"/LOGS"))		
		color=['r','b','g','k','g','b','r','k','r','g','b','k']
		marker_type=['D','D','s','s','o','D','s','p','o','D','s','o']
		#line_style=['--','-','-.',':']
		peak_lum_model_array=[]
		h1_mass_min_DUP_model_array=[]
		for i in range(len(dirs)):
			#color,marker_type)
			peak_lum_model,h1_mass_min_DUP_model = historydata[i].find_TP_attributes(0,t0_model[i],color[i],marker_type[i],h_core_mass,no_fig=plot_fig)			
			peak_lum_model_array.append(peak_lum_model)
			h1_mass_min_DUP_model.append(h1_mass_min_DUP_model)

		return peak_lum_model_array,h1_mass_min_DUP_model_array

				
	###the following methods allow are parts of Falks vis3.py file, thanks Falk

	def set_plot_hrd(self,end_model=[]):
		'''
			end_model - array, control how far in models a run is plottet, if -1 till end

		'''
		m=self.run_historydata
		figure(1)
    		i=0
    		for case in m:
			t1_model=-1
			if end_model[i] != -1:
				t1_model=end_model[i]
			t0_model=case.get("model_number")[0]
        		logTeff=case.get('log_Teff')[:(t1_model-t0_model)]
        		logL=case.get('log_L')[:(t1_model-t0_model)]
        		plot(logTeff,logL,symbs[i],label=self.run_label[i])
        		i += 1
		case.xlimrev()
                ax = plt.gca()
                plt.rcParams.update({'font.size': 16})
                plt.rc('xtick', labelsize=16)
                plt.rc('ytick', labelsize=16)
    		legend(loc=4)
   		xlabel('log Teff',fontsize=18)
    		ylabel('log L',fontsize=18)
 
	def set_plot_t_mass(self):
    				
		m=self.run_historydata
		figure(2)
    		i=0
    		for case in m:
        		star_age=case.get('star_age')
       			h1_boundary_mass=case.get('h1_boundary_mass')
        		star_mass=case.get('star_mass')
        		plot(star_age[noffset:],h1_boundary_mass[noffset:],symbs[i],label=self.run_label[i])
        		plot(star_age[noffset:],star_mass[noffset:],symbs[i],label=self.run_label[i])
        	i += 1
    		legend(loc=6)
    		xlabel('time / yr')
    		ylabel('star mass, core mass')


	def set_plot_model_mass(self):
    		m=self.run_historydata
		figure(3)
    		i=0
    		for case in m:
        		model_number=case.get('model_number')
        		h1_boundary_mass=case.get('h1_boundary_mass')
      			star_mass=case.get('star_mass')
			plot(model_number[noffset:],h1_boundary_mass[noffset:],symbs[i],label=self.run_label[i])
        		plot(model_number[noffset:],star_mass[noffset:],symbs[i],label='')
        		i += 1
    		legend(loc=1)
    		xlabel('model number')
    		ylabel('star mass, core mass')

	def set_plot_CO(self):
		m=self.run_historydata    
		figure(4)
		i=0
	    	for case in m:
	       		star_age=case.get('star_age')
			model_number=case.get('model_number')
        		C=case.get('surface_c12')
        		O=case.get('surface_o16')
			CO=C*4./(O*3.)
        		plot(model_number[noffset:],CO[noffset:],symbs[i],label=self.run_label[i])
			i += 1
	    	legend(loc=2)
		xlabel('model number')
	    	ylabel('C/O ratio')

	def set_plot_CO_mass(self,end_model=[]):
    		m=self.run_historydata
		figure(4)
    		i=0
    		for case in m:
			t1_model=-1
                        if end_model[i] != -1:
                                t1_model=end_model[i]
                        t0_model=case.get("model_number")[0]
        		star_age=case.get('star_age')[:(t1_model-t0_model)]
        		model_number=case.get('model_number')[:(t1_model-t0_model)]
        		star_mass=case.get('star_mass')[:(t1_model-t0_model)]
        		C=case.get('surface_c12')[:(t1_model-t0_model)]
			O=case.get('surface_o16')[:(t1_model-t0_model)]
        		CO=C*4./(O*3.)
			plot(star_mass,CO,symbs[i],label=self.run_label[i])
        		i += 1
   		legend(loc=2)
    		ax = plt.gca()
                plt.rcParams.update({'font.size': 16})
                plt.rc('xtick', labelsize=16)
                plt.rc('ytick', labelsize=16)
		ax.invert_xaxis()
    		plt.xlabel('star mass $[M_{\odot}]$',fontsize=18)
		plt.ylabel('C/O Verhaeltnis', fontsize=18)

        def set_plot_CO_mH(self):
                m=self.run_historydata
                figure(4)
                i=0
                for case in m:
                        star_age=case.get('star_age')
                        model_number=case.get('model_number')
                        h1_mass=case.get('h1_boundary_mass')
			star_mass=case.get('star_mass')
			C=case.get('surface_c12')
                        O=case.get('surface_o16')
                        CO=C*4./(O*3.)
                        plot(star_mass[noffset:],CO[noffset:],symbs[i],label=self.run_label[i])
			plot(star_mass[noffset:],h1_mass[noffset:],symbs[i],label=self.run_label[i])
                        i += 1
                legend(loc=2)
                ax = plt.gca()
                ax.invert_xaxis()
                xlabel('star mass [M$_{H}$]')
                ylabel('C/O ratio')



	def set_plot_vsurf(self):
    		m=self.run_historydata
		figure(55)
   		i=0
    		for case in m:
			case.plot('model_number','v_div_csound_surf',legend=self.run_label[i],shape=symbs[i])
        	i += 1
    		legend(loc=2)
    		xlabel('model number')
    		ylabel('v_div_csound_surf')
    		if xlim_mod_min >= 0:
        		xlim(xlim_mod_min,xlim_mod_max)
    		ylim(-10.,10.)

	def set_plot_mdot(self,xtime=True):

                if xtime==True:
                        t0_model=self.set_find_first_TP()

		else:
			t0_model=len(self.run_historydata)*[0]

    		m=self.run_historydata
		figure(6)
   		i=0
    		for case in m:
			if xtime==True:

	                        t0_time=case.get('star_age')[t0_model[i]]
				print t0_time
				star_mass=case.get('star_mass')[t0_model[i]]
        	                star_age=case.get('star_age')[t0_model[i]:]  -t0_time
				mdot=case.get('log_abs_mdot')[t0_model[i]:]
				plt.plot(star_age,mdot,symbs[i],label=self.run_label[i])
			else:
        			case.plot('model_number','log_abs_mdot',legend=self.run_label[i],shape=symbs[i])
       			i += 1
    		legend(loc=2)
		if xtime==True:
			plt.xlabel('star age')
		else:
			plt.xlabel('model number')
    		ylabel('log_abs_mdot')
    		#ylim(-7,-3.5)

	def set_plot_R(self):
		m=self.run_historydata
	    	figure(7)
    		i=0
    		for case in m:
        		case.plot('star_age','log_R',legend=self.run_label[i],shape=symbs[i])
        	i += 1
    		legend(loc=2)
    		xlabel('model number')
    		ylabel('log_R')
    		if xlim_mod_min >= 0:
        		xlim(xlim_mod_min,xlim_mod_max)
	
	def set_plot_kipp(self,xtime=True,startfirstTP=True):
    		m=self.run_historydata
		i=0
                if startfirstTP==True:
                        t0_model=self.set_find_first_TP()

                else:
                        t0_model=len(self.run_historydata)*[0]
		
		if xtime==True:
    			for case in m:
        			case.kippenhahn(i,'time',t0_model=t0_model[i],c12_bm=False)
        			title(self.run_label[i])
        			i += 1
		else:
                	for case in m:
                                case.kippenhahn(i,'model',t0_model=t0_model[i],c12_bm=False)
                                title(self.run_label[i])
                                i += 1



	def set_plot_kipp_CO(self, startfirstTP=False):
		if startfirstTP==True:
			t0_model=self.set_find_first_TP()
		else:
			t0_model=len(self.run_historydata)*[0]
		m=self.run_historydata
    		i=10
		j=0
    		for case in m:
        		case.kippenhahn_CO(i,'model',t0_model=t0_model[j])
        		title(self.run_label[i-10])
        		i += 1
			j+=1


	def set_plot_surfabu(self):
    		m=self.run_historydata
		i=40
    		for case in m:
    			case.t_surfabu(i,'model')
			title(self.run_label[i-40])
        		legend(loc=4)
        		i += 1

	def set_plot_lumi(self):              
		m=self.run_historydata
                i=40
                for case in m:
			figure(i)
                        case.t_lumi(i,'model')
                        title(self.run_label[i-40])
                        legend(loc=4)
                        i += 1

	

	def set_thing(self,thing=['log_LH','log_LHe']):
		'''
			Specify what to plot
		'''
    		m=self.run_historydata
		i=70
   		for case in m:
        		figure(i)
			j = 0
			for thing in things:
            			case.plot('model_number',thing,legend=thing,shape=symbs[j])
            			j += 1
        		title(self.run_label[i-70])
        		i += 1
	
	def set_find_first_TP(self,dirs=[]):


                historydata=[]
                if (len(dirs)) == 0:
                        dirs=self.run_LOGS
                        historydata=self.run_historydata
                else:
                        for i in range(len(dirs)):
                                historydata.append(history_data(dirs[i]+"/LOGS"))
		t0_models=[]
		for j in range(len(dirs)):
			h1_boundary_mass  = historydata[j].get('h1_boundary_mass')
                	he4_boundary_mass = historydata[j].get('he4_boundary_mass')
                	star_mass         = historydata[j].get('star_mass')
                	#mx1_bot           = historyda.get('mx1_bot')*star_mass
                	#mx1_top           = historydata.get('mx1_top')*star_mass
                	mx2_bot           = historydata[j].get('mx2_bot')*star_mass
                	#mx2_top           = historydata.get('mx2_top')*star_mass
			he_lumi		  = historydata[j].get('log_LHe')
			h_lumi 		  = historydata[j].get('log_LH')
			#model_number            = historydata[j].get('model_number')
			lum_array=[]
			activate=False
			models=[]
			for i in range(len(h1_boundary_mass)):
				if (h1_boundary_mass[i]-he4_boundary_mass[i] <0.1) and (he4_boundary_mass[i]>0.2):
					if (mx2_bot[i]>he4_boundary_mass[i]) and (he_lumi[i]>h_lumi[i]):
						activate=True
						lum_array.append(he_lumi[i])
						models.append(i)
					if (activate == True) and (he_lumi[i]<h_lumi[i]):
						break	
			t0_models.append(models[np.argmax(lum_array)]  )
		return t0_models		

	#####################Internal use

	def create_label(self,historydata,extra_label):
		mass=historydata.header_attr['initial_mass']		
		z=historydata.header_attr['initial_z']
		if extra_label==" ":
			return str(mass)+"$M_{\odot}$, z="+str(z)
		else:
                        return str(mass)+"$M_{\odot}$, z="+str(z)+" , "+extra_label
			
def is_stable(isotope):
                
	isotopes={"N-1":"t",
	"H-1":"t","H-2":"t","He-3":"t","He-4":"t","Be-7":"f",
	"B-8":"f","Li-7":"t","C-11":"f","B-11":"t","C-12":"t",
	"C-13":"t","N-13":"f","N-14":"t","C-14":"f","N-15":"t",
	"O-16":"t","O-17":"t","O-18":"t","F-17":"f","F-18":"f",
	"F-19":"t","Ne-20":"t","Ne-21":"t","Ne-22":"t","Na-22":"f",
	"Na-23":"f","Mg-23":"f","Mg-24":"t","Mg-25":"t","Mg-26":"t",
	"Al-26":"f","Al-27":"t","Si-27":"f","Si-28":"t","Si-29":"t",
	"Si-30":"t","P-31":"t","S-31":"f","Be-8":"f","O-14":"f",
	"O-15":"f","Na-21":"f","Al-25":"f","P-29":"f","P-30":"f",
	"Pb-206":"t","Pb-207":"t","Bi-211":"f","Po-210":"f","F-20":"f",
	"Ne-19":"f","Na-24":"f","Mg-27":"f","Mg-28":"f","Al-28":"f",
	"Al-29":"f","Si-31":"f","Si-32":"f","P-32":"f","P-33":"f",
	"P-34":"f","P-35":"f","S-32":"t","S-33":"t","S-34":"t",
	"S-35":"f","S-36":"t","S-37":"f","S-38":"f","Cl-34":"f",
	"Cl-35":"t","Cl-36":"l","Cl-37":"t","Cl-38":"f","Cl-39":"f",
	"Cl-40":"f","Ar-35":"f","Ar-36":"t","Ar-37":"f","Ar-38":"t",
	"Ar-39":"f","Ar-40":"t","Ar-41":"f","Ar-42":"f","Ar-43":"f",
	"Ar-44":"f","K-38":"f","K-39":"t","K-40":"t","K-41":"t",
	"K-42":"f","K-43":"f","K-44":"f","K-45":"f","K-46":"f",
	"Ca-39":"f","Ca-40":"t","Ca-41":"f","Ca-42":"t","Ca-43":"t",
	"Ca-44":"t","Ca-45":"f","Ca-46":"t","Ca-47":"f","Ca-48":"t",
	"Ca-49":"f","Sc-43":"f","Sc-44":"f","Sc-45":"t","Sc-46":"f",
	"Sc-47":"f","Sc-48":"f","Sc-49":"f","Sc-50":"f","Ti-44":"f",
	"Ti-45":"f","Ti-46":"t","Ti-47":"t","Ti-48":"t","Ti-49":"t",
	"Ti-50":"t","Ti-51":"f","Ti-52":"f","V-47":"f","V-48":"f",
	"V-49":"f","V-50":"t","V-51":"t","V-52":"f","V-53":"f",
	"Cr-48":"f","Cr-49":"f","Cr-50":"t","Cr-51":"f","Cr-52":"t",
	"Cr-53":"t","Cr-54":"t","Cr-55":"f","Cr-56":"f","Mn-51":"f",
	"Mn-52":"f","Mn-53":"f","Mn-54":"f","Mn-55":"t","Mn-56":"f",
	"Mn-57":"f","Fe-52":"f","Fe-53":"f","Fe-54":"t","Fe-55":"f",
	"Fe-56":"t","Fe-57":"t","Fe-58":"t","Fe-59":"f","Fe-60":"f",
	"Fe-61":"f","Co-55":"f","Co-56":"f","Co-57":"f","Co-58":"f",
	"Co-59":"t","Co-60":"f","Co-61":"f","Co-62":"f","Co-63":"f",
	"Ni-56":"f","Ni-57":"f","Ni-58":"t","Ni-59":"l","Ni-60":"t",
	"Ni-61":"t","Ni-62":"t","Ni-63":"l","Ni-64":"t","Ni-65":"f",
	"Ni-66":"f","Ni-67":"f","Ni-68":"f","Cu-60":"f","Cu-61":"f",
	"Cu-62":"f","Cu-63":"t","Cu-64":"f","Cu-65":"t","Cu-66":"f",
	"Cu-67":"f","Cu-68":"f","Cu-69":"f","Cu-70":"f","Cu-71":"f",
	"Zn-62":"f","Zn-63":"f","Zn-64":"t","Zn-65":"f","Zn-66":"t",
	"Zn-67":"t","Zn-68":"t","Zn-69":"f","Zn-70":"t","Zn-71":"f",
	"Zn-72":"f","Zn-73":"f","Zn-74":"f","Ga-65":"f","Ga-66":"f",
	"Ga-67":"f","Ga-68":"f","Ga-69":"t","Ga-70":"f","Ga-71":"t",
	"Ga-72":"f","Ga-73":"f","Ga-74":"f","Ga-75":"f","Ge-66":"f",
	"Ge-67":"f","Ge-68":"f","Ge-69":"f","Ge-70":"t","Ge-71":"f",
	"Ge-72":"t","Ge-73":"t","Ge-74":"t","Ge-75":"f","Ge-76":"t",
	"Ge-77":"f","Ge-78":"f","As-69":"f","As-70":"f","As-71":"f",
	"As-72":"f","As-73":"f","As-74":"f","As-75":"t","As-76":"f",
	"As-77":"f","As-78":"f","As-79":"f","As-80":"f","As-81":"f",
	"Se-72":"f","Se-73":"f","Se-74":"t","Se-75":"f","Se-76":"t",
	"Se-77":"t","Se-78":"t","Se-79":"l","Se-80":"t","Se-81":"f",
	"Se-82":"t","Se-83":"f","Se-84":"f","Br-74":"f","Br-75":"f",
	"Br-76":"f","Br-77":"f","Br-78":"f","Br-79":"t","Br-80":"f",
	"Br-81":"t","Br-82":"f","Br-83":"f","Br-84":"f","Br-85":"f",
	"Br-86":"f","Br-87":"f","Kr-76":"f","Kr-77":"f","Kr-78":"t",
	"Kr-79":"f","Kr-80":"t","Kr-81":"f","Kr-82":"t","Kr-83":"t",
	"Kr-84":"t","Kr-85":"l","Kr-86":"t","Kr-87":"f","Kr-88":"f",
	"Kr-89":"f","Kr-90":"f","Rb-79":"f","Rb-80":"f","Rb-81":"f",
	"Rb-82":"f","Rb-83":"f","Rb-84":"f","Rb-85":"t","Rb-86":"f",
	"Rb-87":"t","Rb-88":"f","Rb-89":"f","Rb-90":"f","Rb-91":"f",
	"Sr-80":"f","Sr-81":"f","Sr-82":"f","Sr-83":"f","Sr-84":"t",
	"Sr-85":"f","Sr-86":"t","Sr-87":"t","Sr-88":"t","Sr-89":"f",
	"Sr-90":"f","Sr-91":"f","Sr-92":"f","Sr-93":"f","Sr-94":"f",
	"Y-85":"f","Y-86":"f","Y-87":"f","Y-88":"f","Y-89":"t",
	"Y-90":"f","Y-91":"f","Y-92":"f","Y-93":"f","Y-94":"f",
	"Y-95":"f","Y-96":"f","Zr-86":"f","Zr-87":"f","Zr-88":"f",
	"Zr-89":"f","Zr-90":"t","Zr-91":"t","Zr-92":"t","Zr-93":"l",
	"Zr-94":"t","Zr-95":"f","Zr-96":"t","Zr-97":"f","Zr-98":"f",
	"Nb-89":"f","Nb-90":"f","Nb-91":"f","Nb-92":"f","Nb-93":"t",
	"Nb-94":"f","Nb-95":"f","Nb-96":"f","Nb-97":"f","Nb-98":"f",
	"Nb-99":"f","Mo-90":"f","Mo-91":"f","Mo-92":"t","Mo-93":"f",
	"Mo-94":"t","Mo-95":"t","Mo-96":"t","Mo-97":"t","Mo-98":"t",
	"Mo-99":"f","Mo-100":"t","Mo-101":"f","Mo-102":"f","Tc-93":"f",
	"Tc-94":"f","Tc-95":"f","Tc-96":"f","Tc-97":"f","Tc-98":"l",
	"Tc-99":"f","Tc-100":"f","Tc-101":"f","Tc-102":"f","Tc-103":"f",
	"Tc-104":"f","Tc-105":"f","Ru-94":"f","Ru-95":"f","Ru-96":"t",
	"Ru-97":"f","Ru-98":"t","Ru-99":"t","Ru-100":"t","Ru-101":"t",
	"Ru-102":"t","Ru-103":"f","Ru-104":"t","Ru-105":"f","Ru-106":"f",
	"Rh-98":"f","Rh-99":"f","Rh-100":"f","Rh-101":"f","Rh-102":"f",
	"Rh-103":"t","Rh-104":"f","Rh-105":"f","Rh-106":"f","Rh-107":"f",
	"Rh-108":"f","Pd-99":"f","Pd-100":"f","Pd-101":"f","Pd-102":"t",
	"Pd-103":"f","Pd-104":"t","Pd-105":"t","Pd-106":"t","Pd-107":"l",
	"Pd-108":"t","Pd-109":"f","Pd-110":"t","Pd-111":"f","Pd-112":"f",
	"Ag-101":"f","Ag-102":"f","Ag-103":"f","Ag-104":"f","Ag-105":"f",
	"Ag-106":"f","Ag-107":"t","Ag-108":"f","Ag-109":"t","Ag-110":"f",
	"Ag-111":"f","Ag-112":"f","Ag-113":"f","Cd-102":"f","Cd-103":"f",
	"Cd-104":"f","Cd-105":"f","Cd-106":"t","Cd-107":"f","Cd-108":"t",
	"Cd-109":"f","Cd-110":"t","Cd-111":"t","Cd-112":"t","Cd-113":"t",
	"Cd-114":"t","Cd-115":"f","Cd-116":"t","Cd-117":"f","Cd-118":"f",
	"In-106":"f","In-107":"f","In-108":"f","In-109":"f","In-110":"f",
	"In-111":"f","In-112":"f","In-113":"t","In-114":"f","In-115":"t",
	"In-116":"f","In-117":"f","In-118":"f","In-119":"f","Sn-108":"f",
	"Sn-109":"f","Sn-110":"f","Sn-111":"f","Sn-112":"t","Sn-113":"f",
	"Sn-114":"t","Sn-115":"t","Sn-116":"t","Sn-117":"t","Sn-118":"t",
	"Sn-119":"t","Sn-120":"t","Sn-121":"f","Sn-122":"t","Sn-123":"f",
	"Sn-124":"t","Sn-125":"f","Sn-126":"f","Sn-127":"f","Sn-128":"f",
	"Sn-129":"f","Sn-130":"f","Sb-112":"f","Sb-113":"f","Sb-114":"f",
	"Sb-115":"f","Sb-116":"f","Sb-117":"f","Sb-118":"f","Sb-119":"f",
	"Sb-120":"f","Sb-121":"t","Sb-122":"f","Sb-123":"t","Sb-124":"f",
	"Sb-125":"f","Sb-126":"f","Sb-127":"f","Sb-128":"f","Sb-129":"f",
	"Sb-130":"f","Sb-131":"f","Sb-132":"f","Sb-133":"f","Te-114":"f",
	"Te-115":"f","Te-116":"f","Te-117":"f","Te-118":"f","Te-119":"f",
	"Te-120":"t","Te-121":"f","Te-122":"t","Te-123":"t","Te-124":"t",
	"Te-125":"t","Te-126":"t","Te-127":"f","Te-128":"t","Te-129":"f",
	"Te-130":"t","Te-131":"f","Te-132":"f","Te-133":"f","Te-134":"f",
	"I-117":"f","I-118":"f","I-119":"f","I-120":"f","I-121":"f",
	"I-122":"f","I-123":"f","I-124":"f","I-125":"f","I-126":"f",
	"I-127":"t","I-128":"f","I-129":"f","I-130":"f","I-131":"f",
	"I-132":"f","I-133":"f","I-134":"f","I-135":"f","Xe-118":"f",
	"Xe-119":"f","Xe-120":"f","Xe-121":"f","Xe-122":"f","Xe-123":"f",
	"Xe-124":"t","Xe-125":"f","Xe-126":"t","Xe-127":"f","Xe-128":"t",
	"Xe-129":"t","Xe-130":"t","Xe-131":"t","Xe-132":"t","Xe-133":"f",
	"Xe-134":"t","Xe-135":"f","Xe-136":"t","Xe-137":"f","Xe-138":"f",
	"Cs-123":"f","Cs-124":"f","Cs-125":"f","Cs-126":"f","Cs-127":"f",
	"Cs-128":"f","Cs-129":"f","Cs-130":"f","Cs-131":"f","Cs-132":"f",
	"Cs-133":"t","Cs-134":"f","Cs-135":"f","Cs-136":"f","Cs-137":"f",
	"Cs-138":"f","Cs-139":"f","Ba-124":"f","Ba-125":"f","Ba-126":"f",
	"Ba-127":"f","Ba-128":"f","Ba-129":"f","Ba-130":"t","Ba-131":"f",
	"Ba-132":"t","Ba-133":"f","Ba-134":"t","Ba-135":"t","Ba-136":"t",
	"Ba-137":"t","Ba-138":"t","Ba-139":"f","Ba-140":"f","Ba-141":"f",
	"Ba-142":"f","La-127":"f","La-128":"f","La-129":"f","La-130":"f",
	"La-131":"f","La-132":"f","La-133":"f","La-134":"f","La-135":"f",
	"La-136":"f","La-137":"f","La-138":"t","La-139":"t","La-140":"f",
	"La-141":"f","La-142":"f","La-143":"f","Ce-130":"f","Ce-131":"f",
	"Ce-132":"f","Ce-133":"f","Ce-134":"f","Ce-135":"f","Ce-136":"t",
	"Ce-137":"f","Ce-138":"t","Ce-139":"f","Ce-140":"t","Ce-141":"f",
	"Ce-142":"t","Ce-143":"f","Ce-144":"f","Ce-145":"f","Ce-146":"f",
	"Pr-133":"f","Pr-134":"f","Pr-135":"f","Pr-136":"f","Pr-137":"f",
	"Pr-138":"f","Pr-139":"f","Pr-140":"f","Pr-141":"t","Pr-142":"f",
	"Pr-143":"f","Pr-144":"f","Pr-145":"f","Pr-146":"f","Pr-147":"f",
	"Pr-148":"f","Pr-149":"f","Nd-134":"f","Nd-135":"f","Nd-136":"f",
	"Nd-137":"f","Nd-138":"f","Nd-139":"f","Nd-140":"f","Nd-141":"f",
	"Nd-142":"t","Nd-143":"t","Nd-144":"t","Nd-145":"t","Nd-146":"t",
	"Nd-147":"f","Nd-148":"t","Nd-149":"f","Nd-150":"t","Nd-151":"f",
	"Nd-152":"f","Pm-137":"f","Pm-138":"f","Pm-139":"f","Pm-140":"f",
	"Pm-141":"f","Pm-142":"f","Pm-143":"f","Pm-144":"f","Pm-145":"l",
	"Pm-146":"l","Pm-147":"l","Pm-148":"f","Pm-149":"f","Pm-150":"f",
	"Pm-151":"f","Pm-152":"f","Pm-153":"f","Pm-154":"f","Sm-140":"f",
	"Sm-141":"f","Sm-142":"f","Sm-143":"f","Sm-144":"t","Sm-145":"f",
	"Sm-146":"l","Sm-147":"t","Sm-148":"t","Sm-149":"t","Sm-150":"t",
	"Sm-151":"l","Sm-152":"t","Sm-153":"f","Sm-154":"t","Sm-155":"f",
	"Sm-156":"f","Sm-157":"f","Sm-158":"f","Eu-143":"f","Eu-144":"f",
	"Eu-145":"f","Eu-146":"f","Eu-147":"f","Eu-148":"f","Eu-149":"f",
	"Eu-150":"l","Eu-151":"t","Eu-152":"f","Eu-153":"t","Eu-154":"f",
	"Eu-155":"f","Eu-156":"f","Eu-157":"f","Eu-158":"f","Eu-159":"f",
	"Gd-144":"f","Gd-145":"f","Gd-146":"f","Gd-147":"f","Gd-148":"l",
	"Gd-149":"f","Gd-150":"l","Gd-151":"f","Gd-152":"t","Gd-153":"f",
	"Gd-154":"t","Gd-155":"t","Gd-156":"t","Gd-157":"t","Gd-158":"t",
	"Gd-159":"f","Gd-160":"t","Gd-161":"f","Gd-162":"f","Tb-147":"f",
	"Tb-148":"f","Tb-149":"f","Tb-150":"f","Tb-151":"f","Tb-152":"f",
	"Tb-153":"f","Tb-154":"f","Tb-155":"f","Tb-156":"f","Tb-157":"l",
	"Tb-158":"f","Tb-159":"t","Tb-160":"f","Tb-161":"f","Tb-162":"f",
	"Tb-163":"f","Tb-164":"f","Tb-165":"f","Dy-148":"f","Dy-149":"f",
	"Dy-150":"f","Dy-151":"f","Dy-152":"f","Dy-153":"f","Dy-154":"l",
	"Dy-155":"f","Dy-156":"t","Dy-157":"f","Dy-158":"t","Dy-159":"f",
	"Dy-160":"t","Dy-161":"t","Dy-162":"t","Dy-163":"t","Dy-164":"t",
	"Dy-165":"f","Dy-166":"f","Dy-167":"f","Dy-168":"f","Ho-153":"f",
	"Ho-154":"f","Ho-155":"f","Ho-156":"f","Ho-157":"f","Ho-158":"f",
	"Ho-159":"f","Ho-160":"f","Ho-161":"f","Ho-162":"f","Ho-163":"f",
	"Ho-164":"f","Ho-165":"t","Ho-166":"f","Ho-167":"f","Ho-168":"f",
	"Ho-169":"f","Er-154":"f","Er-155":"f","Er-156":"f","Er-157":"f",
	"Er-158":"f","Er-159":"f","Er-160":"f","Er-161":"f","Er-162":"t",
	"Er-163":"f","Er-164":"t","Er-165":"f","Er-166":"t","Er-167":"t",
	"Er-168":"t","Er-169":"f","Er-170":"t","Er-171":"f","Er-172":"f",
	"Er-173":"f","Er-174":"f","Er-175":"f","Tm-159":"f","Tm-160":"f",
	"Tm-161":"f","Tm-162":"f","Tm-163":"f","Tm-164":"f","Tm-165":"f",
	"Tm-166":"f","Tm-167":"f","Tm-168":"f","Tm-169":"t","Tm-170":"f",
	"Tm-171":"f","Tm-172":"f","Tm-173":"f","Tm-174":"f","Tm-175":"f",
	"Tm-176":"f","Yb-160":"f","Yb-161":"f","Yb-162":"f","Yb-163":"f",
	"Yb-164":"f","Yb-165":"f","Yb-166":"f","Yb-167":"f","Yb-168":"t",
	"Yb-169":"f","Yb-170":"t","Yb-171":"t","Yb-172":"t","Yb-173":"t",
	"Yb-174":"t","Yb-175":"f","Yb-176":"t","Yb-177":"f","Yb-178":"f",
	"Yb-179":"f","Yb-180":"f","Lu-165":"f","Lu-166":"f","Lu-167":"f",
	"Lu-168":"f","Lu-169":"f","Lu-170":"f","Lu-171":"f","Lu-172":"f",
	"Lu-173":"f","Lu-174":"f","Lu-175":"t","Lu-176":"t","Lu-177":"f",
	"Lu-178":"f","Lu-179":"f","Lu-180":"f","Lu-181":"f","Lu-182":"f",
	"Hf-166":"f","Hf-167":"f","Hf-168":"f","Hf-169":"f","Hf-170":"f",
	"Hf-171":"f","Hf-172":"f","Hf-173":"f","Hf-174":"t","Hf-175":"f",
	"Hf-176":"t","Hf-177":"t","Hf-178":"t","Hf-179":"t","Hf-180":"t",
	"Hf-181":"f","Hf-182":"f","Hf-183":"f","Hf-184":"f","Hf-185":"f",
	"Ta-169":"f","Ta-170":"f","Ta-171":"f","Ta-172":"f","Ta-173":"f",
	"Ta-174":"f","Ta-175":"f","Ta-176":"f","Ta-177":"f","Ta-178":"f",
	"Ta-179":"f","Ta-180":"t","Ta-181":"t","Ta-182":"f","Ta-183":"f",
	"Ta-184":"f","Ta-185":"f","Ta-186":"f","W-172":"f","W-173":"f",
	"W-174":"f","W-175":"f","W-176":"f","W-177":"f","W-178":"f",
	"W-179":"f","W-180":"t","W-181":"f","W-182":"t","W-183":"t",
	"W-184":"t","W-185":"f","W-186":"t","W-187":"f","W-188":"f",
	"W-189":"f","W-190":"f","Re-175":"f","Re-176":"f","Re-177":"f",
	"Re-178":"f","Re-179":"f","Re-180":"f","Re-181":"f","Re-182":"f",
	"Re-183":"f","Re-184":"f","Re-185":"t","Re-186":"l","Re-187":"t",
	"Re-188":"f","Re-189":"f","Re-190":"f","Re-191":"f","Os-179":"f",
	"Os-180":"f","Os-181":"f","Os-182":"f","Os-183":"f","Os-184":"t",
	"Os-185":"f","Os-186":"t","Os-187":"t","Os-188":"t","Os-189":"t",
	"Os-190":"t","Os-191":"f","Os-192":"t","Os-193":"f","Os-194":"f",
	"Os-195":"f","Os-196":"f","Ir-181":"f","Ir-182":"f","Ir-183":"f",
	"Ir-184":"f","Ir-185":"f","Ir-186":"f","Ir-187":"f","Ir-188":"f",
	"Ir-189":"f","Ir-190":"f","Ir-191":"t","Ir-192":"f","Ir-193":"t",
	"Ir-194":"f","Ir-195":"f","Ir-196":"f","Ir-197":"f","Pt-184":"f",
	"Pt-185":"f","Pt-186":"f","Pt-187":"f","Pt-188":"f","Pt-189":"f",
	"Pt-190":"t","Pt-191":"f","Pt-192":"t","Pt-193":"f","Pt-194":"t",
	"Pt-195":"t","Pt-196":"t","Pt-197":"f","Pt-198":"t","Pt-199":"f",
	"Pt-200":"f","Pt-201":"f","Pt-202":"f","Au-185":"f","Au-186":"f",
	"Au-187":"f","Au-188":"f","Au-189":"f","Au-190":"f","Au-191":"f",
	"Au-192":"f","Au-193":"f","Au-194":"f","Au-195":"f","Au-196":"f",
	"Au-197":"t","Au-198":"f","Au-199":"f","Au-200":"f","Au-201":"f",
	"Au-202":"f","Au-203":"f","Hg-189":"f","Hg-190":"f","Hg-191":"f",
	"Hg-192":"f","Hg-193":"f","Hg-194":"f","Hg-195":"f","Hg-196":"t",
	"Hg-197":"f","Hg-198":"t","Hg-199":"t","Hg-200":"t","Hg-201":"t",
	"Hg-202":"t","Hg-203":"f","Hg-204":"t","Hg-205":"f","Hg-206":"f",
	"Hg-207":"f","Hg-208":"f","Tl-192":"f","Tl-193":"f","Tl-194":"f",
	"Tl-195":"f","Tl-196":"f","Tl-197":"f","Tl-198":"f","Tl-199":"f",
	"Tl-200":"f","Tl-201":"f","Tl-202":"f","Tl-203":"t","Tl-204":"l",
	"Tl-205":"t","Tl-206":"f","Tl-207":"f","Tl-208":"f","Tl-209":"f",
	"Tl-210":"f","Pb-193":"f","Pb-194":"f","Pb-195":"f","Pb-196":"f",
	"Pb-197":"f","Pb-198":"f","Pb-199":"f","Pb-200":"f","Pb-201":"f",
	"Pb-202":"f","Pb-203":"f","Pb-204":"t","Pb-205":"l","Pb-208":"t",
	"Pb-209":"f","Pb-210":"f","Pb-211":"f","Bi-202":"f","Bi-203":"f",
	"Bi-204":"f","Bi-205":"f","Bi-206":"f","Bi-207":"f","Bi-208":"f",
	"Bi-209":"t","Bi-210":"f","Al-*26":"f","Kr-*85":"f","Cd-*15":"f"}        
	
	isotope="".join(isotope.split())
	for i in range(len(isotopes.keys())):
                #print "|"+isotope+stable_list[i][0]+"|"
		if isotope == isotopes.keys()[i]:
			return isotopes[isotope]
