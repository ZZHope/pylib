########################################
# File to check for corrupt .h5 files
# using h5ls command (requirement)

# Adopt to path you want to look in
fol2scan = '/Users/reto/Desktop/h5check'

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
problemfiles = list()
for h5file in h5files:
   out = sp.call(['h5stat',h5file])
   if out != 0:
      problemfiles.append(h5file)

# write output file w/ problems, or return that everything is fine
print '\n\n\n'

if len(problemfiles) == 0:
   print 'All files are fine'
else:
   f_out = open('h5check_report.txt','w')
   f_out.writelines('Using h5ls, problems in the following .h5 files were detected:\n')
   for line in problemfiles:
      print 'Problems detected in file ' + line
      f_out.writelines(line + '\n')
   f_out.close()
