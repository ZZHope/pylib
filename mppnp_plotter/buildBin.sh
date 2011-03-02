#!/bin/bash

#	To run this script for the first time, run the following script after uncommenting the the following block of code
#	After you have run this once, you can recomment out the block and run it from there (if need be)

#echo downloading from svn
#svn co http://svn.pyinstaller.org/trunk pyinstaller
#echo installing
#python pyinstaller/source/linux/Make.py
#cd pyinstaller/source/linux/ 
#make
#cd ..
#cd ..
#cd ..

#	Do Not Comment Past Here

echo configuring

python -O  pyinstaller/Configure.py

echo Starting compile

python -O pyinstaller/Makespec.py --onefile h5Plot.py 
python -O pyinstaller/Build.py h5Plot.spec 

echo all good
