from mppnp import *
from asserts import *
import os
	
def main():
	print '\nMPPNP test\n'
	try:
		import h5T
	except ImportError:
		print 'Error: Mppnp requires the module h5T'
		return None
	try:
		import h5py
	except ImportError:
		print 'Error: Mppnp requires the module h5py'
		return None
	try:
		import ascii_table
	except ImportError:
		print 'Error: Mppnp requires the module ascii_table'
		return None
	m= se('./mppnp')
	
	initTests(m)
	return None

def initTests(p):
	test='Init'
	
	assertEquals(p.sedir,'./mppnp',lineno()+': Standard Directory expected',test)
	
	assertEquals(p.pattern,'.h5',lineno()+': Pattern expected',test)
	
	assertEquals(len(p.sefiles),1,lineno()+': Number of files expected',test)
	
	assertEquals(len(p.se.cycles),50,lineno()+': Expected cycle length',test)
	
	assertEquals(len(p.se.ages),50,lineno()+': Expected age length',test)
	
	assertEquals(p.se.cycles[25], '0000000520',lineno()+': Expected cycle Value',test)
	
	assertEquals(p.se.ages[25],222057298101769.66,lineno()+': Expected age Value',test)
	
	assertEquals(len(p.se.hattrs),16,lineno()+': Expected hattrs length',test)
	
	assertEquals(p.se.hattrs[8],'overini',lineno()+': Expected hattrs Value',test)
	
	assertEquals(len(p.se.cattrs),4,lineno()+': Expected cycle length',test)
	
	assertEquals(p.se.cattrs[2],'deltat',lineno()+': Expected cycle Value',test)
	
	assertEquals(len(p.se.Tables),3,lineno()+': Expected Table length',test)
	
	assertEquals(len(p.se.Z),397,lineno()+': Expected cycle length',test)
	
	assertEquals(len(p.se.A),397,lineno()+': Expected cycle length',test)
	
	assertEquals(len(p.se.isomeric_states),397,lineno()+': Expected cycle length',test)
	
	assertEquals(len(p.se.isotopes),397,lineno()+': Expected cycle length',test)
	
	assertEquals(len(p.se.Tables[0]),len(p.se.Z),lineno()+': Expected Table-Z length',test)
	
	assertEquals(len(p.se.Tables[1]),len(p.se.A),lineno()+': Expected Table-A length',test)
	
	assertEquals(len(p.se.Tables[2]),len(p.se.isomeric_states),lineno()+': Expected Table-Isomer length',test)
	
	assertEquals(p.se.Z[100],30,lineno()+': Expected cycle length',test)
	
	assertEquals(p.se.A[100],64,lineno()+': Expected cycle length',test)
	
	assertEquals(p.se.isomeric_states[0],1,lineno()+': Expected cycle length',test)
	
	assertEquals(p.se.isotopes[100],'Zn-64',lineno()+': Expected cycle length',test)
	
	assertTrue(os.path.exists('./mppnp/h5Preproc.txt'),lineno()+': Preprocessor Was not Written',test)
	
def dataTests(p):
	test='Data'
	assertEquals(p.get(0,'age'),30065.513185346295,lineno()+': Expected age Value',test)
	
	assertEquals(p.get(0,'deltat'),5108.110334306285,lineno()+': Expected deltat Value',test)
	
	assertEquals(p.get(0,'shellnb'),425,lineno()+': Expected shellnb Value',test)
	
	assertEquals(p.get(0,'total_mass')[0],1.5,lineno()+': Expected total_mass Value',test)
	
	assertEquals(p.get(1000,'age'),74866845076710336.0,lineno()+': Expected age Value',test)
	
	assertEquals(p.get(1000,'deltat'),8409582784235.7109,lineno()+': Expected deltat Value',test)
	
	assertEquals(p.get(1000,'shellnb'),790,lineno()+': Expected shellnb Value',test)
	
	assertEquals(p.get(1000,'total_mass'),1.5,lineno()+': Expected total_mass Value',test)


if __name__ == "__main__":
   	main()

