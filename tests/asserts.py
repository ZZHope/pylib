'''
asserts.py

a collection of assert statmentd for testing a program

'''

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
