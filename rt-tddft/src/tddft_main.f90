
!---
PROGRAM tddft_main
  !---
  ! This is the main driver of the real time TDDFT propagation.
  USE mp_global,       ONLY : mp_startup, nproc_pool_file
  USE environment,     ONLY : environment_start
  USE mp_bands,        ONLY : nbgrp
  USE check_stop,      ONLY : check_stop_init
  USE control_flags,   ONLY : io_level, use_para_diag
  USE uspp,            ONLY : okvan
  USE paw_variables,   ONLY : okpaw 
  USE noncollin_module,ONLY : noncolin
  USE wvfct,           ONLY : nbnd
!  USE tddft_module,    ONLY : job, tddft_exit_code
!  USE mp_images,       ONLY : nimage, my_image_id
!  USE mp_pools,        ONLY : nproc_pool
!  USE lsda_mod,        ONLY : nspin
!  USE fft_base,        ONLY : dffts
!  USE iotk_module  
!  USE xml_io_base
  USE io_global,       ONLY : stdout
  implicit none
  character (len=9)   :: code = 'QE'
  logical, external  :: check_para_diag

  ! initialize the environment
#ifdef __PARA
  call mp_startup(start_images=.true.)
#else
  call mp_startup(start_images=.false.)
#endif
  call environment_start(code)

  ! restrict RT-tddft to certain parallelization types
#ifndef __BANDS
  if (nbgrp > 1) then
    call errore('tddft_main', 'RT-tddft does not support band-parallelization ', 1)
  endif
#endif

  call start_clock('RT-tddft')
  call tddft_read_input()
  call check_stop_init()

  ! read PW calculations setup
  io_level = 1
  call read_file
  
  ! restrict RT-tddft to certain pseudo-potential types
  if (okvan) then
    call errore('tddft_main', 'RT-tddft does not support Vanderbilt-type potentials ',1)
  endif
  if (okpaw) then
    call errore('tddft_main', 'RT-tddft does not support PAW potentials ',1)
  endif
  ! restrict RT-tddft to certain collinear types
  if (noncolin) then
    call errore('tddft_main', 'RT-tddft does not support non-collinear spin polarization',1)
  endif


#ifdef __PARA
  use_para_diag = check_para_diag(nbnd)
#else
  use_para_diag = .false.
#endif

  ! read PW ground state wavefunctions
  call tddft_openfile()
!
!  if (gamma_only) then
!    call errore ('tdddft_main', 'Cannot run TDFFT with gamma_only == .true. ', 1)
!  endif
!if ((twfcollect .eqv. .false.)  .and. (nproc_pool_file /= nproc_pool)) &
!    call errore('tddft_main', 'Different number of CPU/pool. Set wf_collect=.true. in SCF', 1)
!#ifdef __BANDS
!  if (nbgrp > 1 .and. (twfcollect .eqv. .false.)) &
!    call errore('tddft_main', 'Cannot use band-parallelization without wf_collect in SCF', 1)
!#endif
!  if (noncolin) call errore('tdddft_main', 'non-collinear not supported yet', 1)
!
!  call tddft_allocate()
!  call tddft_setup()
!  call tddft_summarize()
!
!#ifdef __BANDS
!  call init_parallel_over_band(inter_bgrp_comm, nbnd)
!#endif
!
!  ! calculation
!  select case (trim(job))
!  case ('transport')
!     call molecule_optical_absorption
!  case default
!     call errore('tddft_main', 'RT-tddft cannot recognize this job type in input', 1)
!  end select
!  
!  ! print timings and stop the code
!  call tddft_closefil
   call print_clock_tddft()
!  call stop_run(tddft_exit_code)
!  call do_stop(tddft_exit_code)
!  
!  STOP
  
END PROGRAM tddft_main

