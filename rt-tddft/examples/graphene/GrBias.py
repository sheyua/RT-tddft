import sys; sys.path.append('../..')
from tools.tdpost import *

# point tdpost to the location of dump files
a = Rho('dump')									# ./ by default, ~ is okay, not ending with / is okay
a.load()									# load the charge dump files
a.mapcur()									# it takes one parameter half (True by default) to specify whether to hide the mirror image
#b = Vks('dump')
#b.load()									# load the potential dump files
