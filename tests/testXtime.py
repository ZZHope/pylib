from ppn import *
from asserts import *
def main():
	print '\nxtime test\n'
	p=xtime('./ppn/')
	
	initTests(p)
	
	dataTests(p)
	
	return None
	
def initTests(p):
	test='Init'
	assertEquals(p.sldir,'./ppn/',lineno()+': Standard Directory expected',test)
	
	assertEquals(len(p.data),39,lineno()+': Expected data length',test)
	
	assertEquals(len(p.cols),1103,lineno()+': Expected column length',test)
	
	assertEquals(p.cols[0],'t_y',lineno()+': Expected column value',test)
	
	assertEquals(p.cols[-1],'TAg80',lineno()+': Expected column value',test)
	
	assertEquals(len(p.col_num),1099,lineno()+': Expected column data length',test)

def dataTests(p):
	test='Data'
	assertEquals(len(p.get('H   2')),38,lineno()+': Expected data length',test)
		
	assertEquals(p.get('H   2')[0],1.41376e-05,lineno()+': Expected data Value',test)
		
	assertEquals(p.get('H   2')[-1],2.6127899999999998e-24,lineno()+': Expected data Value',test)
	
	assertEquals(len(p.get('BE  7')),38,lineno()+': Expected data length',test)
		
	assertEquals(p.get('BE  7')[0],0,lineno()+': Expected data Value',test)
		
	assertEquals(p.get('BE  7')[-1],1.5883899999999999e-20,lineno()+': Expected data Value',test)
if __name__ == "__main__":
   	main()

