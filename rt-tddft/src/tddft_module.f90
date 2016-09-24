
!---
MODULE tddft_module
  !---
  ! This module contains the variables used for RT-tddft calculations
  USE kinds,    ONLY : dp
  implicit none
  save
  
  ! input parameters
  character(80) :: job                          ! 'transport'
  character(80) :: solver                       ! 'itsolver' or 'cgsolver'
  character(80) :: method                       ! 'CN', 'CN2', or 'CN-mid'
  real(dp)      :: conv_thr                     ! convergence threshold for the solver
  integer       :: max_iter                     ! maximum iteration of the solver
  logical       :: dump                         ! wether dump rho and vks or not
  character(80) :: dump_dir                     ! directory for dumping rho and vks
  real(dp)      :: dt                           ! time step length delta_t
  integer       :: num_step                     ! number of timesteps
  integer       :: init_step                    ! initial istep may be used for restart
  logical       :: e_mirror                     ! external bias use mirror image
  real(dp)      :: e_pstart                     ! bias voltage hits maximum
  real(dp)      :: e_pend                       ! bias voltage starts to drop
  real(dp)      :: e_nstart                     ! bias voltage drops to minimum
  real(dp)      :: e_nend                       ! bias voltage starts to increase again
  real(dp)      :: e_volt                       ! height of the bias voltage
  real(dp)      :: e_decay                      ! how fast the bias voltage decays

  ! shared global parameters
  integer, parameter :: iuntdwfc = 51           ! to save RT-tddft intermediate wfcs
  integer, parameter :: iuntdrho = 32           ! to save RT-tddft intermediate rhos
  integer, parameter :: iuntdvks = 33           ! to save RT-tddft intermediate vkss
  integer :: nwordtdwfc 
  integer, allocatable :: nbnd_occ(:)           ! occupied bands for each k-point
  integer :: max_nbnd_occ                       ! max number of occupied bands
  ! shared allocatable
  real(dp), allocatable :: rpos(:,:)            ! position operator in real space
  real(dp), allocatable :: num_elec(:)          ! total num_elecs for each spin 
  real(dp), allocatable :: dipole(:,:)          ! dipole moment for each spin
  complex(dp), allocatable :: tddft_psi(:,:)    ! time-propagated wvfcts

END MODULE tddft_module
!---
