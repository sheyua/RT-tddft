from libcpp.string cimport string

cdef extern from "C_Rho.h" namespace "TD_Post":
  cdef cppclass C_Rho:
    # member variables
    string dump_dir
    int nspin
    double dt
    double c
    int num_step
    int init_step
    # constructor
    C_Rho(string) except+
    # methods
    void load()

cdef extern from "C_Vks.h" namespace "TD_Post":
  cdef cppclass C_Vks:
    # member variables
    string dump_dir
    int nspin
    double dt
    double c
    int num_step
    int init_step
    # constructor
    C_Vks(string) except+
    # methods
    void load()
