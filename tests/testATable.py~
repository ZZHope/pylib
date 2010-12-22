from ascii_table import *
from asserts import *

def main():
	print '\nAscii Table test\n'
	a=AsciiTable('c12pg.dat','ascii')
	initTest(a)
	dataTest(a)
	'''
	keys=a.dcols.keys()
	for i in range(len(keys)):
		print keys[i]
		#print a.dcols[keys[i]]
	
	a.plot('upper','T9',plotType='ATable', shape='-')
	a.plot('upper','ado/CA88',plotType='ATable', shape='-')
	a.plot('CA88','ado/fit',plotType='ATable', shape='-')
	'''
	a=AsciiTable('c12.dat','ascii',sep=';')
	initTest(a)
	dataTest(a)
	'''
	keys=a.dcols.keys()
	for i in range(len(keys)):
		print keys[i]
		#print a.dcols[keys[i]]
	
	a.plot('upper','T9',plotType='ATable', shape='-')
	a.plot('upper','ado/CA88',plotType='ATable', shape='-')
	a.plot('CA88','ado/fit',plotType='ATable', shape='-')
	
	'''

def initTest(a):
	test='Init'
	assertEquals(len(a.hattrs),3,'Header length expected',test)
	
	assertEquals(a.hattrs[0],'1 12  6  1  1  1  0  0  0  1 13  7  0','Header expected',test)
	
	assertEquals(a.hattrs[1],'55   1.943','Header expected',test)
	
	assertEquals(a.hattrs[2],'c12pg','Header expected',test)
	
	assertEquals(len(a.dcols), 9,'Column length expected',test)
	
	assertIn(a.dcols, 'T9','Column expected',test)
	
	assertIn(a.dcols, 'ado/fit','Column expected',test)
	
def dataTest(a):
	test='Data'
	assertEquals(len(a.data['T9']),42 ,'Expected T9 Length', test)
	
	assertEquals(a.data['T9'][0],float(.006),'Expected T9 Data Value', test)
	
	assertEquals(a.data['T9'][41],float(1.250),'Expected T9 Data Value', test)
	
	assertEquals(len(a.data['fitted']),42 ,'Expected fitted Length', test)
	
	assertEquals(a.data['fitted'][0],float(1.21E-24),'Expected fitted Data Value', test)
	
	assertEquals(a.data['fitted'][41],float(1.42E+03),'Expected fitted Data Value', test)
	
	assertEquals(len(a.data['ado/fit']),42 ,'Expected ado/fit Length', test)
	
	assertEquals(a.data['ado/fit'][0],float(1.01),'Expected ado/fit Data Value', test)
	
	assertEquals(a.data['ado/fit'][41],float(1.01),'Expected ado/fit Data Value', test)
if __name__ == "__main__":
   	main()

