from asserts import *
from mesa import *

def main():
	print '\nStar log test\n'
	
	p=star_log('LOGS') 
	
	initTest(p)
	
def initTest(p):
	test='init'
	assertEquals(len(p.header_attr),7,'Expected header_attr length',test)
	
	assertIn(p.header_attr,'initial_mass','Header', test)
	
	assertEquals(p.header_attr['initial_mass'],float(2.0000000000000000E+000),'Expected Header Value',test)
	
	assertIn(p.header_attr,'burn_min2','Header', test)
	
	assertEquals(p.header_attr['burn_min2'],float(1.0000000000000000E+003),'Expected Header Value',test)
	
	assertIn(p.header_attr,'he4_boundary_limit','Header', test)
	
	assertEquals(p.header_attr['he4_boundary_limit'],float(1.0000000000000000E-004),'Expected Header Value',test)
	
	assertEquals(len(p.cols),62,'Expected p.cols length',test)
	
	assertEquals(p.cols['model_number'],1,'Expected p.cols value',test)
	
	assertEquals(p.cols['star_age'],2,'Expected p.cols value',test)
	
	assertEquals(p.cols['h1_boundary_lgT'],25,'Expected p.cols value',test)
	
	assertEquals(p.cols['num_retries'],61,'Expected p.cols value',test)
	
	assertEquals(p.cols['num_backups'],62,'Expected p.cols value',test)
	
	assertEquals(len(p.data[0]),62,'Expected data length',test)
	
	assertEquals(p.data[0][0],114,'Expected Data Value',test)
	
	assertEquals(p.data[0][49],float(1.7123643547999989E-003),'Expected Data Value',test)
	
	assertEquals(p.data[99886][0],float(100000),'Expected Data Value',test)
	
	assertEquals(p.data[99886][59],float(5.0283367238775538E-003),'Expected Data Value',test)
if __name__ == "__main__":
   	main()

