import numpy as np
from matplotlib import pyplot as plt

def check(sawarg):
	epstart = 0.01
	epend = 0.24
	enstart = 0.26
	enend = 0.49

	if ( sawarg < epstart ):
		value = (sawarg/epstart) * 0.5
	if ( sawarg >= epstart) and (sawarg <= epend ):
		value = 0.5
	if ( sawarg > epend) and (sawarg < enstart ):
		value = 0.5 - (sawarg-epend)/(enstart-epend)
	if ( sawarg >= enstart) and (sawarg <= enend ):
		value = -0.5
	if ( sawarg > enend) and (sawarg <= 0.5 ):
		value = -0.5 * (0.5-sawarg)/(0.5-enend)
	if ( sawarg > 0.5) and (sawarg < (1.0-enend) ):
		value = -0.5 * (sawarg-0.5)/(0.5-enend)
	if ( sawarg >= (1.0-enend)) and (sawarg <= (1.0-enstart) ):
		value = -0.5
	if ( sawarg > (1.0-enstart)) and (sawarg < (1.0-epend) ):
		value = 0.5 + (sawarg+epend-1.0)/(enstart-epend)
	if ( sawarg >= (1.0-epend)) and (sawarg <= (1.0-epstart) ):
		value = 0.5
	if ( sawarg > (1.0-epstart) ):
		value = 0.5 * (1.0-sawarg)/epstart

	return value

x = np.arange(0,1.00005,0.001)
y = np.array([check(i) for i in x])
