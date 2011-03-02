#  This script assumes that you are using python2.5, have all the essential libraries installed
#	and have the pyinstaller package.
#	The file organization should be as follows:
#	->pyinstaller (dir) ->->
#	->db.py
#	->this script
#	To run this script for the first time, run the following script after uncommenting the the following block of code
#	After you have run this once, you can recomment out the block and run it from there (if need be)

#echo downloading from svn
#svn co http://svn.pyinstaller.org/trunk pyinstaller
#echo installing
#python pyinstaller\source\linux\Make.py
#cd pyinstaller\source\linux\ 
#make
#cd ..
#cd ..
#cd ..

#	Do Not Comment Past Here

echo configuring

c:\python25\python.exe pyinstaller\Configure.py

echo building spec file.............

c:\python25\python.exe pyinstaller\Makespec.py -F h5Plot.py

echo compiling to binary.............

c:\python25\python.exe pyinstaller\Build.py h5Plot.spec

echo done
