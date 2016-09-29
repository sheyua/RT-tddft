import os.path

cdef class Rho:
  # member variables
  cdef C_Rho *thisptr
  cdef public int nspin
  cdef public double dt
  cdef public double c
  cdef public int num_step
  cdef public int init_step
  # constructor
  def __cinit__(self, str dump_dir):
    dump_dir = os.path.expanduser(dump_dir)
    dump_dir = os.path.abspath(dump_dir)
    cpp_string = dump_dir.encode('utf-8')
    self.thisptr = new C_Rho(cpp_string)
  # destructor
  def __dealloc__(self):
    del self.thisptr
  # member methods
  def load(self):
    self.thisptr.load()
    self.nspin = self.thisptr.nspin
    self.dt = self.thisptr.dt
    self.c = self.thisptr.c
    self.num_step = self.thisptr.num_step
    self.init_step = self.thisptr.init_step
  def plot(self):
    plt.ion();
    plt.figure();
