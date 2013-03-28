'''
This script can be easily adapted to make many nugridse plots, which can
later be turned into a movie.
'''

# ~~~~~~~~~~~~~~~~ Start Declarations ~~~~~~~~~~~~~~~~ #

import os
import numpy
import pylab
import nugridse as ns
import matplotlib.pyplot as plt
from nvis import *

# ~~~~~~~~~~~~~~~~ End Declarations ~~~~~~~~~~~~~~~~ #




# ~~~~~~~~~~~~~~~~ Start global variables (edit to suit your needs) ~~~~~~~~~~~~~~~~ #

labelxtmp = 'Mass Coordinate' # x-Axis label
labelytmp = 'X' # y-Axis label
showtmp='False' 
atrixtmp = 'mass' 		# x-Axis attribute
limitstmp=[0,20,-9,0] 	# Plot window size. Format = [x1, x2, y1, y2]
logy_tf = True 		# True/False of whether the y-axis is log of the value
show_tf = False 	# whether to show each plot as it is generated (True slows the process down)
filedir= '/nfs/rpod3/fherwig/SEE/data/set1.1/ppd_wind/M20.0.Z1.0e-02.standard/H5_out' # The directory in which the H5_out files are

# Isotopes to be plotted (y-Axis attribute)
iso1='H-1'
iso2='He-4'
iso3='C-12'
iso4='O-16'
iso5='Ne-20
iso6='Si-28'
iso7='Fe-56'
iso8='Sr-86'
iso9='Ba-136'
iso10='Pb-206'

# ~~~~~~~~~~~~~~~~ End global variables ~~~~~~~~~~~~~~~~ #




# ~~~~~~~~~~~~~~~~ Start running code ~~~~~~~~~~~~~~~~ #

plt.ioff() 			# turns off the live plot so that a new plot doesn't open each time the code loops	
m20=ns.se(filedir) 	# loads the H5_out files using nugridse.py

i = 0 			# The array number of the cycle you want to start with. i.e. If the cycle numbers go up in 20, Cycle 20 is '0'
while i <= 4: 	# The array number of the cycle you want to end with. i.e. If the cycle numbers go up in 20, Cycle 100 is '4'
	
	
	
	
	# ~~~~~~~~~~~~~~~~ Start plot specific variables ~~~~~~~~~~~~~~~~ #
	cycletmp = m20.se.cycles[i].lstrip('0') 					# The current plots cycle number
	timetmp = '{:.2e}'.format(float(m20.se.ages[i])) 			# The current plots time stamp for YEARS
	#timetmp = '{:.2e}'.format(float(m20.se.ages[i])/31536000) 	# The current plot time stamp for SECONDS (uncomment one of the two)
	titletmp = 'Cycle = ' + cycletmp + ' Time= ' + str(timetmp) + ' years' # Plot title
	fnametmp = m20.se.cycles[i] 

	# ~~~~~~~~~~~~~~~~ End plot specific variables ~~~~~~~~~~~~~~~~ #
	
	
	
	
	m20.plot(atrix=atrixtmp, atriy=iso1, fname=fnametmp, legend=iso1,labelx=labelxtmp, labely=labelytmp ,shape='-', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)
	m20.plot(atrix=atrixtmp, atriy=iso2, fname=fnametmp, legend=iso2,labelx=labelxtmp, labely=labelytmp ,shape='--', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)
	m20.plot(atrix=atrixtmp, atriy=iso3, fname=fnametmp, legend=iso3,labelx=labelxtmp, labely=labelytmp ,shape='-.', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)
	m20.plot(atrix=atrixtmp, atriy=iso4, fname=fnametmp, legend=iso4,labelx=labelxtmp, labely=labelytmp ,shape=':', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)
	m20.plot(atrix=atrixtmp, atriy=iso5, fname=fnametmp, legend=iso5,labelx=labelxtmp, labely=labelytmp ,shape='o', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)
	m20.plot(atrix=atrixtmp, atriy=iso6, fname=fnametmp, legend=iso6,labelx=labelxtmp, labely=labelytmp ,shape='v', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)
	m20.plot(atrix=atrixtmp, atriy=iso7, fname=fnametmp, legend=iso7,labelx=labelxtmp, labely=labelytmp ,shape='^', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)	
	m20.plot(atrix=atrixtmp, atriy=iso8, fname=fnametmp, legend=iso8,labelx=labelxtmp, labely=labelytmp ,shape='<', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)
	m20.plot(atrix=atrixtmp, atriy=iso9, fname=fnametmp, legend=iso9,labelx=labelxtmp, labely=labelytmp ,shape='>', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)
	m20.plot(atrix=atrixtmp, atriy=iso10, fname=fnametmp, legend=iso10,labelx=labelxtmp, labely=labelytmp ,shape='.', title=titletmp, logy=logy_tf, show=show_tf, limits=limitstmp)

	os.chdir('/rpod2/home/nsbruce/M20') 		# The directory you wanted the files to be stored in
	plt.savefig('Plot.' + i +'.png', dpi=200) 	# the file names
	plt.clf() 			# clears the current plot so the code can loop and make it again
	os.chdir(filedir) 	# goes back to the directory with the files in it
	i=i+5 				# the sparsity, or how many plots to skip. If the cycle numbers go up in 20, '+5' plots every 100 cycles

# ~~~~~~~~~~~~~~~~ End running code ~~~~~~~~~~~~~~~~ #