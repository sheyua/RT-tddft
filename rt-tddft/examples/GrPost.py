import sys; sys.path.append('..')
from tools.tdpost import *

# point tdpost to the location of dump files
a = Rho('GrBias/dump')
a.load()
a.comp_cur()
b = Rho('Gr/dump')
b.load()
# ground state [unbiased] charge density
b.plot_rho()
plt.title('Unbiased Ground State Rho [Half]', fontsize=25)
b.plot_rho(0,False)
plt.title('Unbiased Ground State Rho [Full]', fontsize=25)
# ground state charge density bias
#b.plot_rho()
#plt.title('Unbiased Ground State Charge Density', fontsize=25)
# first step current
a.plot_cur()
plt.title('First Timestep Current [Half]', fontsize=25)
a.plot_cur(1,False)
plt.title('First Timestep Current [Full]', fontsize=25)
# mid current between 2 layers
a.slice_hcur()
plt.title('Mid Cur Btw 2 Layers', fontsize=25)
