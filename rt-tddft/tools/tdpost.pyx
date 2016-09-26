import os.path

cdef class Rho:
  # member variables
  cdef C_Rho *thisptr
  # constructor
  def __cinit__(self, str dump_dir):
    dump_dir = os.path.expanduser(dump_dir)
    dump_dir = os.path.abspath(dump_dir)
    cpp_string = dump_dir.encode('utf-8')
    self.thisptr = new C_Rho(cpp_string)
  # destructor
  def __dealloc__(self):
    del self.thisptr

cdef class Vks:
  # member variables
  cdef C_Vks *thisptr
  # constructor
  def __cinit__(self, str dump_dir):
    dump_dir = os.path.expanduser(dump_dir)
    dump_dir = os.path.abspath(dump_dir)
    cpp_string = dump_dir.encode('utf-8')
    self.thisptr = new C_Vks(cpp_string)
  # destructor
  def __dealloc__(self):
    del self.thisptr
