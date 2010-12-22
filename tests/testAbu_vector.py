from ppn import *
from asserts import *
def main():
	print '\nabu_vector test\n'
	p=abu_vector('./ppn/')
	
	initTests(p)
	
	dataTests(p)
	
	return None
	
def initTests(p):
	test='Init'
	assertEquals(p.sldir,'./ppn/',lineno()+': Standard Directory expected',test)
	
	assertEquals(len(p.cattrs),8,lineno()+': Expected cattrs length',test)
	
	assertEquals(len(p.dcols),6,lineno()+': Expected data column length',test)
	
	assertEquals(len(p.files),2,lineno()+': Expected file list length',test)
	
	assertEquals(p.cattrs[0],'mod',lineno()+': Expected cattrs value',test)
	
	assertEquals(p.cattrs[-1],'densa',lineno()+': Expected cattrs value',test)
	
	assertEquals(p.dcols[0],'NUM',lineno()+': Expected Data column value',test)
	
	assertEquals(p.dcols[-1],'ISOTP',lineno()+': Expected Data column value',test)
	
def dataTests(p):
	test='Data'
	assertEquals(len(p.get('Z',0)),len(p.get('A',0)),lineno()+': Expected length Z-A',test)
	
	assertEquals(len(p.get('ABUNDNACE_MF',0)),1099,lineno()+': Expected data length',test)
	
	assertEquals(p.get('ABUNDNACE_MF',0)[1],7.28729000e-01,lineno()+': Expected abundence Value',test)
		
	assertEquals(p.get('ABUNDNACE_MF',0)[-1],0,lineno()+': Expected abundence Value',test)
	
	assertEquals(p.get('Z',0)[1],1,lineno()+': Expected data Value',test)
		
	assertEquals(p.get('Z',0)[-1],73,lineno()+': Expected data Value',test)
	
if __name__ == "__main__":
   	main()

