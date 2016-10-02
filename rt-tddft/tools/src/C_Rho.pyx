cimport numpy as cnp
import numpy as np
import os.path
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

cdef class Rho:
  # member variables
  cdef C_Rho *thisptr
  cdef public int nspin
  cdef public double dt
  cdef public double c
  cdef public double dc
  cdef public int init_step
  cdef public int num_step
  cdef public int num_mids
  cdef public int num_hmids
  cdef public cnp.ndarray vbias
  cdef public cnp.ndarray rho
  cdef public cnp.ndarray cur
  cdef public cnp.ndarray hrho
  cdef public cnp.ndarray hcur
  # constructor
  def __cinit__(self, str dump_dir):
    dump_dir = os.path.expanduser(dump_dir)
    dump_dir = os.path.abspath(dump_dir)
    cpp_string = dump_dir.encode('utf-8')
    self.thisptr = new C_Rho(cpp_string)
    # parse member variables
    self.nspin = self.thisptr.nspin
    self.dt = self.thisptr.dt
    self.c = self.thisptr.c
    self.init_step = self.thisptr.init_step
    self.num_step = self.thisptr.num_step
    self.num_mids = self.thisptr.num_mids
    self.num_hmids = np.ceil(self.num_mids/2).astype(int)
    self.dc = self.c / self.num_mids
  # destructor
  def __dealloc__(self):
    del self.thisptr
  # member methods
  def load(self):
    self.vbias = np.empty(shape=(self.num_step+1,), dtype=np.double, order='C')
    self.rho = np.empty(shape=(self.nspin, self.num_step+1,self.num_mids), dtype=np.double, order='C')
    # make cnp copies to the same memory location as python objects 
    cdef cnp.ndarray[double, ndim=1, mode='c'] vbias = self.vbias
    cdef cnp.ndarray[double, ndim=3, mode='c'] rho = self.rho
    self.thisptr.load(&vbias[0], &rho[0,0,0])
    self.hrho = 0.5*(self.rho[:,:,:self.num_hmids] + np.fliplr(self.rho[:,:,self.num_hmids:]))
  def comp_cur(self):
    if self.num_step == 0:
      print('There is no time propagation')
      return
    if self.rho == None:
      self.load();
    # compute current
    self.cur = np.empty(shape=(self.nspin, self.num_step,self.num_mids), dtype=np.double, order='C')
    cdef cnp.ndarray[double, ndim=3, mode='c'] cur = self.cur
    self.thisptr.comp_cur(&cur[0,0,0])
    self.hcur = 0.5*(self.cur - self.cur[:,:,::-1])[:,:,:self.num_hmids]

  def plot_rho(self, istep=0, half=True):
    # check istep
    if istep < self.init_step:
      print(istep, 'is less than the initial step', self.init_step)
      return
    elif istep > self.init_step+self.num_step:
      print(istep, 'is greater than the last step', self.init_step+self.num_step)
      return
    # make a plot
    plt.ion()
    plt.figure()
    if half:
      xcoor = np.arange(0.5,self.num_hmids,1)*self.dc
    else:
      xcoor = np.arange(0.5,self.num_mids,1)*self.dc
    for sdx in range(self.nspin):
      # for spin sdx
      if half:
        plt.plot(xcoor, self.hrho[sdx,istep,:]/self.dc, linewidth=2, label='spin '+str(sdx+1))
      else:
        plt.plot(xcoor, self.rho[sdx,istep,:]/self.dc, linewidth=2, label='spin '+str(sdx+1))
    if half:
      plt.xlim((0,self.c/2))
    else:
      plt.xlim((0,self.c))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('z Position [A]', fontsize=25)
    plt.ylabel('Charge Density [1/A]', fontsize=25)
    plt.legend(loc='best', fontsize=20)
    if istep==0:
      plt.title('Ground State Rho', fontsize=25)
    else:
      plt.title("Rho at %.0f atto" % (istep*self.dt), fontsize=25)
    plt.tight_layout()

  def plot_cur(self, istep=1, half=True):
    # check istep
    if istep < self.init_step+1:
      print(istep, 'is less than the first step', self.init_step+1)
      return
    elif istep > self.init_step+self.num_step:
      print(istep, 'is greater than the last step', self.init_step+self.num_step)
      return
    # check cur exist
    if self.cur is None:
      self.comp_cur()
    if self.cur is None:
      return
    # make a plot
    plt.ion()
    plt.figure()
    if half:
      xcoor = np.arange(0.5,self.num_hmids,1)*self.dc
    else:
      xcoor = np.arange(0.5,self.num_mids,1)*self.dc
    for sdx in range(self.nspin):
      # for spin sdx
      if half:
        plt.plot(xcoor, self.hcur[sdx,istep,:]*1e+6, linewidth=2, label='spin '+str(sdx+1))
      else:
        plt.plot(xcoor, self.cur[sdx,istep,:]*1e6, linewidth=2, label='spin '+str(sdx+1))
    if half:
      plt.xlim((0,self.c/2))
    else:
      plt.xlim((0,self.c))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('z Position [A]', fontsize=25)
    plt.ylabel('I [uAmp]', fontsize=25)
    plt.legend(loc='best', fontsize=20)
    plt.title("Current at %.0f atto" % (istep*self.dt), fontsize=25)
    plt.tight_layout()

  def slice_hcur(self, pos=0.5, with_vbias=True):
    # check pos
    if pos < 0 or pos >= 1:
      print('must slice at a fraction 0 <= pos < 1')
      return
    # make a plot
    pos = np.ceil(pos*self.num_hmids).astype(int)
    plt.ion()
    plt.figure()
    xcoor = np.arange(self.init_step+1,self.num_step+self.init_step+1,1)*self.dt
    for sdx in range(self.nspin):
      plt.plot(xcoor, self.hcur[sdx,:,pos]*1e+6, linewidth=2, label='spin '+str(sdx+1))
    plt.legend(loc='best', fontsize=20)
    plt.xticks(fontsize=20)
    plt.xlabel('Time [atto]', fontsize=25)
    plt.yticks(fontsize=20)
    plt.ylabel('I [uAmp]', fontsize=25)
    if with_vbias:
      ax2 = plt.gca().twinx()
      ax2.plot(xcoor, self.vbias[1:]*13.6057, linewidth=2, color='r')
    if with_vbias:
      ax2.set_ylabel('Vleft-Vright [eV]', fontsize=25, color='r')
      for t2 in ax2.get_yticklabels():
        t2.set_color('r')
        t2.set_fontsize(20)
    plt.title("Current at %.0f %%" % (pos*100/self.num_hmids), fontsize=25)
    plt.tight_layout()

#  def map_cur(self, half=True):
#    plt.ion()
#    plt.figure()
#    if half:
#      tdim = self.totHcur.shape[0]
#      zdim = self.totHcur.shape[1]
#    else:
#      tdim = self.cur.shape[1]
#      zdim = self.cur.shape[2]
#    # x-y meshgrid
#    dx = self.c/self.num_mids
#    xcoor = slice(0.5*dx, (zdim+0.5)*dx, dx)
#    dy = self.dt
#    ycoor = slice(dy, (tdim+1)*dy, dy)
#    U = np.mgrid[ycoor, xcoor]
#    # map totHcur
#    if half:
#      plt.pcolormesh(U[1], U[0], self.totHcur)
#    # map cur
#    else:
#      plt.pcolormesh(U[1], U[0], np.sum(self.cur, axis=0))
#    # x axis
#    plt.xticks(fontsize=20)
#    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#    plt.xlim( (0, zdim*dx) )
#    plt.xlabel('z Position [A]', fontsize=25)
#    # y axis
#    plt.yticks(fontsize=20)
#    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
#    plt.ylim( (self.init_step*dy, (tdim+self.init_step)*dy) )
#    plt.ylabel('Time [atto]', fontsize=25)
#    # color bar
#    cb = plt.colorbar()
#    cb.ax.tick_params(labelsize=20)
#    cb.ax.xaxis.set_label_position('top')
#    cb.ax.set_xlabel('I [uA]', fontsize=25)
#    # tighten layout
#    plt.tight_layout()
#    plt.plot(xcoor, rho/dx)
