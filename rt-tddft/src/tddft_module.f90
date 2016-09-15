
!---
MODULE tddft_module
  !---
  ! This module contains the variables used for RT-tddft calculations
  USE kinds, ONLY : DP
  
  IMPLICIT NONE
  SAVE
  
  ! ... Input parameters
  character(80) :: job             ! 'transport'
  real(dp) :: conv_threshold       ! convergence threshold for the linear solver
  real(dp) :: dt                   ! time step length delta_t
  integer  :: num_step             ! number of timesteps for real-time tddft
  integer  :: init_step            ! initial istep may be used for restart
  logical  :: e_mirror             ! external bias use mirror image
  real(dp) :: e_pstart             ! position where the bias voltage hits maximum
  real(dp) :: e_pend               ! position where the bias voltage starts to drop
  real(dp) :: e_nstart             ! position where the bias voltage drops to minimum
  real(dp) :: e_nend               ! position where the bias voltage starts to increase again
  real(dp) :: e_volt               ! height of the bias voltage
  real(dp) :: e_decay              ! how fast the bias voltage decays

  ! ... Shared global parameters
  complex(dp), parameter :: i_complex = (0.0_dp,1.0_dp)
  real(dp), allocatable :: r_pos(:,:)     ! position operator in real space
  real(dp), allocatable :: r_pos_s(:,:)   ! position operator in real space (smooth grid)
  integer, parameter :: iuntdwfc = 51     ! to save TDDFT intermediate wfcs
  integer :: nwordtdwfc 
  integer, allocatable :: nbnd_occ(:)     ! occupied bands for each k-point
  integer :: nbnd_occ_max                 ! max number of occupied bands
  real(dp) :: alpha_pv                    ! shift of conduction levels

  integer :: tddft_exit_code = 0

END MODULE tddft_module
!---
