########################################
# File to check for corrupt .h5 files
# using h5ls command (requirement)

# Adopt to path you want to look in
fol2scan = '/astro/reto/mnt/NuGrid/data/set1'

# How many cycles are saved in each file?
# or at least should be saved...
no_cycles = 1000.   # 1000 for Set1

# Then start npython
# run h5check_robot.py
# look at output file h5check_report.txt

# reto, 2012 
########################################

import os
import subprocess as sp

# message:
print 'Depending on the number of subfolders and files, this might take some time!'

# list all files in directory and add to allfiles list
allfiles = list()
for path, subdirs, files in os.walk(fol2scan):
    for name in files:
         allfiles.append(os.path.join(path, name))

# now go through and sort out h5 files
h5files = list()
for i in range(len(allfiles)):
   if allfiles[i][len(allfiles[i])-3:len(allfiles[i])] == '.h5':
      h5files.append(allfiles[i])

# check all the .h5 files and see if error
problemfilesh5ls = list()
problemfiles_nocyc = list()
for h5file in h5files:
   startcycle = float(h5file.split('.')[len(h5file.split('.'))-3])
   if startcycle != 0.:   # to avoid the 000000 problem we had before
      p = sp.Popen(['h5ls',h5file],stdout=sp.PIPE)
      out,err = p.communicate()
      if p.returncode != 0:
         problemfilesh5ls.append(h5file)
      else:   # check if every cycle is saved in output file
         cyclelist = list()
         out = out.split('\n')
         for i in range(len(out)):
            if out[i][0:5] == 'cycle':
               cycletmp = float(out[i].split()[0][6:len(out[i].split()[0])])
               cyclelist.append(cycletmp)
         tmp_nocycles = False
         tmp_misscyc  = False
         if len(cyclelist) != no_cycles:
            print len(cyclelist)
            tmp_nocycles = True
         for i in range(len(cyclelist)):
            if cyclelist[i] != startcycle + i:
               tmp_misscyc = startcycle + i

         if tmp_nocycles and tmp_misscyc == False:
            problemfiles_nocyc.append(h5file + ': Total number of cycles: ' + str(int(len(cyclelist))))
         elif tmp_nocycles and tmp_misscyc != False:
            problemfiles_nocyc.append(h5file + ': Total number of cycles: ' + str(int(len(cyclelist))) + '\tCycle ' + str(int(tmp_misscyc)) + ' is missing.')
         elif tmp_nocycles == False and tmp_misscyc != False:
            problemfiles_nocyc.append(h5file + ': Total number of cycles missing / shifted -> starting at cycle ' + str(int(tmp_misscyc)))
      

# write output file w/ problems, or return that everything is fine
print '\n\n\n'

if len(problemfilesh5ls) == 0 and len(problemfiles_nocyc) == 0:
   print 'All files are fine'
else:
   f_out = open('h5check_report.txt','w')
   f_out.writelines('Report\n')
   f_out.writelines('======\n\n')
   f_out.writelines('Using h5ls, problems in the following .h5 files were detected:\n')
   f_out.writelines('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
   for line in problemfilesh5ls:
      print 'Problems in file ' + line
      f_out.writelines(line + '\n')
   f_out.writelines('\n\n')
   f_out.writelines('The following problems with cycles in non-corrupt h5 files were detected:\n')
   f_out.writelines('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
   for line in problemfiles_nocyc:
      f_out.writelines(line + '\n')
   f_out.close()
