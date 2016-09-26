from libcpp.string cimport string

cdef extern from "C_Rho.h" namespace "TD_Post":
  cdef cppclass C_Rho:
    # member variables
    int num_proc
    string dump_dir
    # constructor
    C_Rho(string) except+

cdef extern from "C_Vks.h" namespace "TD_Post":
  cdef cppclass C_Vks:
    # member variables
    int num_proc
    string dump_dir
    # constructor
    C_Vks(string) except+
