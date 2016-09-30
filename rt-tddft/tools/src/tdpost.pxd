from libcpp.string cimport string

cdef extern from "C_Rho.h" namespace "TD_Post":
  cdef cppclass C_Rho:
    # member variables
    string dump_dir
    int nspin
    double dt
    double c
    int init_step
    int num_step
    int num_mids
    # constructor
    C_Rho(string) except+
    # methods
    void load(double*, double*)

#cdef extern from "C_Vks.h" namespace "TD_Post":
#  cdef cppclass C_Vks:
#    # member variables
#    string dump_dir
#    int nspin
#    double dt
#    double c
#    int init_step
#    int num_step
#    # constructor
#    C_Vks(string) except+
#    # methods
#    void load()
