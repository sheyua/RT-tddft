import sys; sys.path.append('..')
from tools.tdpost import *
'''
Step1: point tdpost to the location of dump files
'''
a = Rho('GrBias/dump')
a.load()
b = Rho('Gr/dump')
b.load()
'''
Step2: ground state [unbiased] charge density
'''
b.plot_rho()
plt.title('Unbiased Ground State Rho', fontsize=25)
#b.plot_rho(0,False)
#plt.title('Unbiased Ground State Rho [Full]', fontsize=25)
'''
Step3: Vbias will introduce an initial charge imbalance
'''
gs_rho_bias = np.sum(a.hrho[:,0,:] - b.hrho[:,0,:], axis=0) # sum over spin
plt.figure()
plt.plot(np.arange(0,a.num_hmids,1)*a.dc,gs_rho_bias*100/a.dc, linewidth=2)
plt.xlim((0,a.c/2))
plt.xticks(fontsize=20)
plt.xlabel('z Position [A]', fontsize=25)
plt.yticks(fontsize=20)
plt.ylabel('Charge Density [0.01/A]', fontsize=25)
plt.title('Rho[Vbias=1eV] - Rho[Vbias=0eV]', fontsize=25)
plt.tight_layout()
'''
Step4: map out the charge flow through real-time propagation
'''
# mid current between 2 layers
a.slice_hcur()
plt.title('Mid Cur Btw 2 Layers', fontsize=25)
# map the real-time current
a.map_hcur()
plt.title('Cur to the R-Elec', fontsize=25)
