'''
asserts.py

a collection of assert statments and testing tools for testing a program

'''
import inspect


def assertEquals(actual,expected, msg='expected', test=''):
	
	if actual!=expected:
		print 'Assert true failed: '+msg+ ' '+str(expected)
		print 'Actual: ' +str(actual)
	else:
		if test=='':
			print 'Test passed'
		else:
			print test+ ' Test passed'
def assertIn(actual,expected, msg='', test=''):
	
	if expected not in actual:
		if msg=='':
			print 'Assert failed: '+str(expected)+ ' DNE'
		else:
			print 'Assert failed: '+str(expected)+ ' DNE in '+str(msg)
	else:
		if test=='':
			print 'Test passed'
		else:
			print test+ ' Test passed'

def assertTrue(boo,msg='', test=''):
	if boo:
		if test=='':
			print 'Test passed'
		else:
			print test+ ' Test passed'

	else:
		if msg=='':
			print 'Assert failed: False'
		else:
			print 'Assert failed: '+str(msg)
def lineno():
	return str(inspect.currentframe().f_back.f_lineno)
	
