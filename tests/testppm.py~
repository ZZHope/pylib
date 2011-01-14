from ppm import *
from asserts import *
def main():
	print '\nPPM test\n'
	stddir='./YProfile'
	
	
	p=yprofile(stddir)
	headers= p.hattrs
	cols= p.dcols
	cyc=p.cattrs
	
	
	line=p.files[len(p.files)-1]
	keys=headers.keys()
	keys.sort()
	
	initTests(p)
	
	dataTests(p)
	'''
	#Data Test
	for i in range(len(cols)):
		print cols[i]
		j=p.get(cols[i],line,numType='file',resolution='a')
		print len(j)
	
	for i in range(len(cyc)):
		print cyc[i]
		j=p.get(cyc[i],line,numType='file',resolution='a')
		print len(j)
		
	for i in range(len(keys)):
		print keys[i]
		j=p.get(keys[i])
		print j
	
	#Num Type Tests
	print '\nInvalid NDump tests\n'
	print 'Negative NDUmp'
	j=p.get('j',-1) # should get warning
	print 'Large NDUmp'
	j=p.get('j',9000) # should get warning
	
	print '\nInvalid File tests\n'
	print 'File DNE'
	j=p.get('j','hello bobaaa',numType='file') # should get warning and get None returned
	print 'File is a int'
	j=p.get('j',45,numType='file') # should get warning and get None returned
	
	print '\nInvalid Time tests\n'
	print 'Negative Time'
	j=p.get('j',-10,numType='t')  # should get warning 	
	print 'Large Time'
	j=p.get('j',90000,numType='t')#	
	print 'Time is a string'
	j=p.get('j','abs',numType='t')
	'''
	'''
	#Plotting Tests
	p.plot('j','Y',) #Test with an Attri X with 2 resolutions and Attri Y with 2 Resolutions 
	p.plot('A','FVconv',) #Test with an Attri X with 2 resolutions and Attri Y with 1 Resolutions of length 1536
	p.plot('Y','RHOconv',)#Test with an Attri X with 2 resolutions and Attri Y with 1 Resolutions of length 768
		      #Feedback: Not equal length warning
	p.plot('UY H+He','Rho1',)#Test with an Attri X and Attri Y with 1 resolutions of length 768
	p.plot('FVconv','UY H+He',)#Test with an Attri X and Attri Y with 1 resolutions of length 1536
	p.plot('Ek H+He','Rho1',)#Test with an Attri X with 1 Resolutions of length 768 and Attri Y with 1 Resolutions of length 1536
				 ##Should produce feedback with no plot occuring
	p.plot('Rho1','Rho H+He',)#Test with an Attri X with 1 Resolutions of length 1536 and Attri Y with 1 Resolutions of length 768
				  #Should produce feedback with no plot occuring 
	p.plot('Ndump','EkXZHHeMax',)#Test with Two cycle Attributes
	
	p.plot('Time','courmx',)#Test with Two Top Attributes
	
	p.plot('Time','Ncycle') #Ncycle test for when it is = ****
	
		'''	

def initTests(p):
	test='Init'
	assertEquals(p.sldir,'./YProfile','Standard Directory expected',test)
	
	assertEquals(len (p.files),11,'Files length expected',test)

	assertEquals(p.files[5],'YProfile-01-0005.bobaaa','File name expected',test)
	
	assertEquals(len(p.hattrs),16,'Expected length of Header attributes',test)
	
	assertIn(p.hattrs,'Stellar Conv. Luminosity','Header Attributes',test)
	
	assertIn(p.hattrs,'At base of the convection zone p','Header Attributes',test)
	
	assertIn(p.hattrs,'Gravity turns off between radii High','Header Attributes',test)
	
	assertEquals(len (p.dcols),33,'Expected column attribute name list length',test)
	
	assertEquals(p.dcols[0],'j','Expected data attribute',test)
	
	assertEquals(p.dcols[20],'HUy','Expected data attribute',test)
	
	assertEquals(p.dcols[32],'Ek H+He','Expected data attribute',test)
	
	assertIn(p.dcols,'j','Data Columns',test)
	
	assertIn(p.dcols,'A','Data Columns',test)
	
	assertIn(p.dcols,'RHO1conv','Data Columns',test)
	
	assertIn(p.dcols,'HUy','Data Columns',test)
	
	
	assertEquals(len(p.cattrs),31,'Expected cycle attribute name list length',test)	
	
	assertEquals(p.cattrs[0],'Ndump','Expected cycle attribute',test)
	
	assertEquals(p.cattrs[27],'EkXZHHeMax','Expected cycle attribute',test)
	
	assertEquals(p.cattrs[30],'courmx','Expected cycle attribute',test)
	
	assertIn(p.cattrs,'Ndump','Cycle Attributes',test)
	
	assertIn(p.cattrs,'EkXZmax','Cycle Attributes',test)
	
	assertIn(p.cattrs,'HUyMax','Cycle Attributes',test)
	
	assertIn(p.cattrs,'EkXZHHeMax','Cycle Attributes',test)
	
	assertIn(p.cattrs,'courmx','Cycle Attributes',test)

def dataTests(p):
	test='Data'
	
	assertEquals(len(p.get('j', p.files[len(p.files)-1],numType='file', resolution ='a')),2,'Expected j Data length',test )
	
	assertEquals(len(p.get('j', p.files[len(p.files)-1],numType='file', resolution ='h')),1536,'Expected j Data length',test )
	
	assertEquals(len(p.get('j', p.files[len(p.files)-1],numType='file', resolution ='l')),768,'Expected j Data length',test )
	
	assertEquals(p.get('j', p.files[len(p.files)-1],numType='file', resolution ='h')[0],1536.0,'Expected j Data Value',test )
	
	assertEquals(float(p.get('j', p.files[len(p.files)-1],numType='file', resolution ='l')[len(p.get('j', p.files[len(p.files)-1],numType='file', resolution ='l'))-1]),1,'Expected j Data Value',test )
	
	assertEquals(len(p.get('Rho', p.files[len(p.files)-1],numType='file', resolution ='a')),1536,'Expected Rho Data length',test )
	
	assertEquals(len(p.get('Rho', p.files[len(p.files)-1],numType='file', resolution ='h')),1536,'Expected Rho Data length',test )
	
	assertEquals(len(p.get('Rho', p.files[len(p.files)-1],numType='file', resolution ='l')),1536,'Expected Rho Data length',test )
	
	assertEquals(p.get('Rho', p.files[len(p.files)-1],numType='file', resolution ='h')[0],float(9.17554E-05),'Expected Rho Data Value',test )
	
	assertEquals(float(p.get('Rho', p.files[len(p.files)-1],numType='file', resolution ='l')[len(p.get('Rho', p.files[len(p.files)-1],numType='file', resolution ='l'))-1]),float(3.45495E+01),'Expected Rho Data Value',test )

	assertEquals(len(p.get('Ek H+He', p.files[len(p.files)-1],numType='file', resolution ='a')),768,'Expected Ek H+He Data length',test )
	
	assertEquals(len(p.get('Ek H+He', p.files[len(p.files)-1],numType='file', resolution ='h')),768,'Expected Ek H+He Data length',test )
	
	assertEquals(len(p.get('Ek H+He', p.files[len(p.files)-1],numType='file', resolution ='l')),768,'Expected Ek H+He Data length',test )
	
	assertEquals(p.get('Ek H+He', p.files[len(p.files)-1],numType='file', resolution ='h')[0],0,'Expected Ek H+He Data Value',test )
	
	assertEquals(float(p.get('Ek H+He', p.files[len(p.files)-1],numType='file', resolution ='l')[len(p.get('Ek H+He', p.files[len(p.files)-1],numType='file', resolution ='l'))-1]),float(0.00000E+00),'Expected Ek H+He Data Value',test )
	
	assertEquals(len(p.get('Conv Ht', p.files[len(p.files)-1],numType='file', resolution ='l')),11,'Expected Conv Ht Data length',test )
	
	assertEquals(p.get('Conv Ht', p.files[len(p.files)-1],numType='file', resolution ='h')[0],float(3.02147E+01),'Expected Conv Ht Data Value',test )
	
	assertEquals(float(p.get('Conv Ht', p.files[len(p.files)-1],numType='file', resolution ='l')[len(p.get('Conv Ht', p.files[len(p.files)-1],numType='file', resolution ='l'))-1]),float(3.01692E+01),'Expected Conv Ht Data Value',test )

	assertEquals(len(p.get('Ncycle', p.files[len(p.files)-1],numType='file', resolution ='l')),11,'Expected Ncycle Data length',test )
	
	assertEquals(float(p.get('Ncycle', p.files[len(p.files)-1],numType='file', resolution ='h')[0]),0.0,'Expected Ncycle Data Value',test )
	
	assertEquals(p.get('Ncycle', p.files[len(p.files)-1],numType='file', resolution ='l')[len(p.get('Ncycle', p.files[len(p.files)-1],numType='file', resolution ='l'))-1],27320,'Expected Ncycle Data Value',test )
	
	assertEquals(p.get('gridY'),'1536', 'Expected GridY Data Value',test )
	
	assertEquals(p.get('At base of the convection zone R'),'9.50000E+00', 'Expected at base of the convection zone R Data Value',test )
	
	assertEquals(p.get('Stellar Conv. Luminosity'),'1.61400E-02 x 10^43 ergs,', 'Expected Stellar Conv. Luminosity Data Value',test )
if __name__ == "__main__":
   	main()

