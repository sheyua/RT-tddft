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
  cdef public int init_step
  cdef public int num_step
  cdef public int num_mids
  cdef public cnp.ndarray vbias
  cdef public cnp.ndarray rho
  cdef public cnp.ndarray cur
  cdef public cnp.ndarray totHcur
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
    # compute current
    self.cur = np.empty(shape=(self.nspin, self.num_step,self.num_mids), dtype=np.double, order='C')
    cdef cnp.ndarray[double, ndim=3, mode='c'] cur = self.cur
    self.thisptr.comp_cur(&cur[0,0,0])
    # sum over spin and do centeral estimation
    hcur = np.sum(self.cur, axis=0)
    hcur = 0.5*(hcur - hcur[:,::-1])[:,:np.ceil(self.num_mids/2).astype(int)]
    self.totHcur = 0.5*(hcur[:-1,:] + hcur[1:,:])
  def mapcur(self, half=True):
    plt.ion()
    plt.figure()
    if half:
      tdim = self.totHcur.shape[0]
      zdim = self.totHcur.shape[1]
    else:
      tdim = self.cur.shape[1]
      zdim = self.cur.shape[2]
    # x-y meshgrid
    dx = self.c/self.num_mids
    xcoor = slice(0.5*dx, (zdim+0.5)*dx, dx)
    dy = self.dt
    ycoor = slice(dy, (tdim+1)*dy, dy)
    U = np.mgrid[ycoor, xcoor]
    # map totHcur
    if half:
      plt.pcolormesh(U[1], U[0], self.totHcur)
    # map cur
    else:
      plt.pcolormesh(U[1], U[0], np.sum(self.cur, axis=0))
    # x axis
    plt.xticks(fontsize=20)
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.xlim( (0, zdim*dx) )
    plt.xlabel('z Position [A]', fontsize=25)
    # y axis
    plt.yticks(fontsize=20)
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    plt.ylim( (self.init_step*dy, (tdim+self.init_step)*dy) )
    plt.ylabel('Time [atto]', fontsize=25)
    # color bar
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=20)
    cb.ax.xaxis.set_label_position('top')
    cb.ax.set_xlabel('I [uA]', fontsize=25)
    # tighten layout
    plt.tight_layout()
  def plot_rho(self, istep=0, half=True):
    if istep < self.init_step:
      print(istep, 'is less than the initial step', self.init_step)
      return
