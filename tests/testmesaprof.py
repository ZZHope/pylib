from asserts import *
from mesa2 import *

def main():
	print '\nprofile test\n'
	p=profile('LOGS',2,'log_num') 
	p.profiles_index()
	
	initTest(p)
	profiles_indexTest(p)
	dataTest(p)
	
	
def dataTest(p):
	test='Data'
	
	assertEquals(len(p.get('zone')),1866,'Expected zone Data length',test)
	
	assertEquals(p.get('zone')[0],1,'Expected model Data',test)
	
	assertEquals(len(p.get('logT')),1866,'Expected logT Data length',test)
	
	assertEquals(p.get('logT')[0],float(3.5040622934916339E+000),'Expected model Data',test)
	
	assertEquals(len(p.get('binding_energy')),1866,'Expected binding_energy Data length',test)
	
	assertEquals(p.get('binding_energy')[0],float(-2.2253068093770015E+012),'Expected model Data',test)
	
	assertEquals(len(p.get('cs_at_cell_bdy')),1866,'Expected  cs_at_cell_bdy Data length',test)
	
	assertEquals(p.get('cs_at_cell_bdy')[0],float(5.4369565672890819E+005),'Expected  cs_at_cell_bdy Data',test)
	
	assertEquals(p.get('cs_at_cell_bdy')[len(p.get('cs_at_cell_bdy'))-1],float(2.6145909546287167E+008),'Expected  cs_at_cell_bdy Data',test)
def profiles_indexTest(p):
	test='profiles_index'
	
	assertEquals(len(p.model),100,'Expected model length',test)
	
	assertEquals(p.model[0],2196,'Expected model Data',test)
	
	assertEquals(p.model[1],2701,'Expected model Data',test)
	
	assertEquals(p.model[98],99000,'Expected model Data',test)
	
	assertEquals(p.model[99],100000,'Expected model Data',test)
	
	assertEquals(len(p.log_ind),100,'Expected log_ind length',test)
	
	assertIn(p.log_ind, 64000, 'log_ind',test)
	
	assertEquals(p.log_ind[64000],98,'Expected log_ind Data',test)
	
	assertIn(p.log_ind, 94000, 'cols',test)
	
	assertEquals(p.log_ind[94000],1,'Expected log_ind Data',test)
	
def initTest(p):
	test= 'Prof init'
	
	assertEquals(len(p.cols),100,'Expected Column length',test)
	
	assertIn(p.cols, 'cs_at_cell_bdy', 'cols',test)
	
	assertEquals(p.cols['cs_at_cell_bdy'],100,'Expected Column Data',test)
	
	assertIn(p.cols, 'other', 'cols',test)
	
	assertEquals(p.cols['other'],99,'Expected Column Data',test)
	
	assertIn(p.cols, 'zone', 'cols',test)
	
	assertEquals(p.cols['zone'],1,'Expected Column Data',test)
	
	assertIn(p.cols, 'logT', 'cols',test)
	
	assertEquals(p.cols['logT'],2,'Expected Column Data',test)
	
	assertIn(p.cols, 'log_cdc', 'cols',test)
	
	assertEquals(p.cols['log_cdc'],43,'Expected Column Data',test)
	
	assertEquals(len(p.header_attr),45,'Expected Header length',test)
	
	assertIn(p.header_attr, 'model_number', 'Header',test)
	
	assertEquals(p.header_attr['model_number'],88000,'Expected Header Data',test)
	
	assertIn(p.header_attr, 'num_zones', 'Header',test)
	
	assertEquals(p.header_attr['num_zones'],1866,'Expected Header Data',test)
	
	assertIn(p.header_attr, 'center_eta', 'Header',test)
	
	assertEquals(p.header_attr['center_eta'],float(3.9258529326313237E+001),'Expected Header Data',test)
	
	assertIn(p.header_attr, 'burn_min1', 'Header',test)
	
	assertEquals(p.header_attr['burn_min1'],float(5.0000000000000000E+001),'Expected Header Data',test)
	
	assertIn(p.header_attr, 'burn_min2', 'Header',test)
	
	assertEquals(p.header_attr['burn_min2'],float(1.0000000000000000E+003),'Expected Header Data',test)
	
	assertEquals(len(p.data[0]),100,'Expected Data length',test)
	
	assertEquals(len(p.data[1000]),len(p.cols),'Expected Data length',test)
	
	assertEquals(p.data[0][0],1,'Expected Data Value',test)
	
	assertEquals(p.data[0][50],float(8.7444047193144983E-002),'Expected Data Value',test)
	
	assertEquals(p.data[0][99],float(5.4369565672890819E+005),'Expected Data Value',test)
	
	assertEquals(p.data[1865][0],float(1866),'Expected Data Value',test)
	
	assertEquals(p.data[1865][99],float(2.6145909546287167E+008),'Expected Data Value',test)
	
	
if __name__ == "__main__":
   	main()
