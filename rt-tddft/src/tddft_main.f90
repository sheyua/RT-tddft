
!---
PROGRAM tddft_main
  !---
  ! This is the main driver of the real time TDDFT propagation.
  USE mp_global,       ONLY : mp_startup, nproc_pool_file, mp_global_end
  USE environment,     ONLY : environment_start, environment_end
  USE check_stop,      ONLY : check_stop_init
  USE control_flags,   ONLY : io_level, use_para_diag
  USE wvfct,           ONLY : nbnd
  USE tddft_module,    ONLY : job
  implicit none
  character(len=9)   :: code = 'RT-tddft'
  logical, external  :: check_para_diag

  ! initialize the environment
#ifdef __PARA
  call mp_startup(start_images=.true.)
#else
  call mp_startup(start_images=.false.)
#endif
  call environment_start(code)

  call tddft_read_input()
  call check_stop_init()

  ! read PW calculations setup
  io_level = 1
  call read_file

#ifdef __PARA
  use_para_diag = check_para_diag(nbnd)
#else
  use_para_diag = .false.
#endif

  ! read PW ground state wavefunctions
  call tddft_openfile()
  call support_check()

  ! summarize the calculation
  call tddft_welcome()
  call tddft_init()

  ! do time propagation
  select case(trim(job))
    case ('transport')
      call tddft_propagate() 
    case default
      call errore('tddft_main', 'RT-tddft cannot recognize this job type in input', 1)
  end select

  ! print timings and stop the code
  call tddft_closefile()
  call tddft_print_clocks()

  ! stop the program
#ifdef __PARA
  call mp_global_end()
#endif
  call environment_end(code)
  STOP

CONTAINS

  !---
  SUBROUTINE support_check()
    !---
    ! RT-tddft does not support the following features
    USE mp_bands,        ONLY : nbgrp
    USE ldaU,            ONLY : lda_plus_u
    USE uspp,            ONLY : okvan
    USE paw_variables,   ONLY : okpaw 
    USE noncollin_module,ONLY : noncolin
    USE ktetra,          ONLY : ltetra
    USE klist,           ONLY : two_fermi_energies
    USE control_flags,   ONLY : gamma_only
    USE mp_global,       ONLY : nproc_pool_file
    USE mp_pools,        ONLY : nproc_pool
    implicit none

    ! restrict RT-tddft to certain parallelization types
    if (nbgrp > 1) then
      call errore('tddft_main', 'RT-tddft does not support band-parallelization!', 1)
    endif

    ! restrict RT-tddft to certain xc functional types
    if (lda_plus_u) then
      call errore('tddft_main', 'RT-tddft does not support LDA plus U!',1)
    endif

    ! restrict RT-tddft to certain pseudo-potential types
    if (okvan) then
      call errore('tddft_main', 'RT-tddft does not support Vanderbilt-type potentials!',1)
    endif
    if (okpaw) then
      call errore('tddft_main', 'RT-tddft does not support PAW potentials!',1)
    endif

    ! restrict RT-tddft to certain collinear types
    if (noncolin) then
      call errore('tddft_main', 'RT-tddft does not support non-collinear spin polarization!',1)
    endif

    ! restrict RT-tddft to certain occupation methods
    if (ltetra) then
      call errore('tddft_main', 'RT-tddft does not support the tetrahedron method!',1)
    endif

    ! restrict RT-tddft to certain electronic distributions
    if (two_fermi_energies) then
      call errore('tddft_main', 'RT-tddft does not support two Fermi energies!',1)
    endif

    ! restrict RT-tddft to manually set gamma point
    if (gamma_only) then
      call errore('tddft_main', 'RT-tddft does not support gamma only calculation! &
          Please manually specify the Gamma point!',1)
    endif

    ! restrict RT-tddft to use the same number of proc as pwscf
    if ( nproc_pool_file /= nproc_pool ) then
      call errore('tddft_main', 'RT-tddft does not support different nproc from PWSCF!', 1)
    endif

  END SUBROUTINE support_check
  !---

END PROGRAM tddft_main
