cimport numpy as cnp
import numpy as np
import os.path
from matplotlib import pyplot as plt

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
  def plot(self):
    plt.ion();
    plt.figure();
