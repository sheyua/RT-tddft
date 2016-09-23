
!---
SUBROUTINE tddft_read_input()
  !---
  ! Read from stdin for the RT-tddft calculation parameters
  USE tddft_module
  USE io_files,         ONLY : prefix, tmp_dir  
  USE io_global,        ONLY : ionode, stdin, stdout
  USE constants,        ONLY : AU_SEC
  USE mp_images,        ONLY : my_image_id
  implicit none
  integer :: ios
  namelist /input_tddft/ job, prefix, tmp_dir, &
    solver, method, conv_threshold, max_iter, &
    dt, num_step, init_step, &
    e_mirror, e_pstart, e_pend, e_nstart, e_nend, e_volt, e_decay
  
  ! define input defult values
  job               = 'transport'
  prefix            = 'pwscf'
  tmp_dir           = './scratch/'

  solver            = 'itsolver'
  method            = 'CN-mid'
  conv_threshold    = 1.0d-16
  max_iter          = 200 

  dt                = 0.1d0
  num_step          = 1000
  init_step         = 1

  e_mirror          = .true.
  e_pstart          = 0.0d0
  e_pend            = 0.25d0
  e_nstart          = 0.25d0
  e_nend            = 0.5d0
  e_volt            = 0.0d0
  e_decay           = 0.0d0

  if ( ionode .and. my_image_id == 0) then
    call input_from_file()
    ! read input    
    read( stdin, input_tddft, err = 100, iostat = ios )
100 call errore('tddft_readin', 'reading input_tddft namelist', abs(ios))
  
    ! convert to atomic units
    dt = dt * 1.d-18 / ( 2.d0*AU_SEC)   ! au_sec using hartree
  endif

#ifdef __PARA
  ! broadcast input variables  
  call tddft_broadcast_input()
#endif

END SUBROUTINE tddft_read_input
!---

#ifdef __PARA
!---
SUBROUTINE tddft_broadcast_input()
  !--- 
  ! Broadcast tddft parameters to all processors 
  USE mp_world,     ONLY : world_comm
  USE mp,           ONLY : mp_bcast
  USE io_files,     ONLY : prefix, tmp_dir
  USE tddft_module
  implicit none
  integer, parameter :: root = 0    

  call mp_bcast(job, root, world_comm)
  call mp_bcast(prefix, root, world_comm)
  call mp_bcast(tmp_dir, root, world_comm)

  call mp_bcast(solver, root, world_comm)
  call mp_bcast(method, root, world_comm)
  call mp_bcast(conv_threshold, root, world_comm)
  call mp_bcast(max_iter, root, world_comm)

  call mp_bcast(dt, root, world_comm)
  call mp_bcast(num_step, root, world_comm)
  call mp_bcast(init_step, root, world_comm)

  call mp_bcast(e_mirror, root, world_comm)
  call mp_bcast(e_pstart, root, world_comm)
  call mp_bcast(e_pend, root, world_comm)
  call mp_bcast(e_nstart, root, world_comm)
  call mp_bcast(e_nend, root, world_comm)
  call mp_bcast(e_volt, root, world_comm)
  call mp_bcast(e_decay, root, world_comm)

END SUBROUTINE tddft_broadcast_input
!---
#endif

!---
SUBROUTINE tddft_openfile()
  ! --- 
  ! Open wavefunctions and k-points
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : npol
  USE buffers,          ONLY : open_buffer
  USE io_files,         ONLY : iunwfc, nwordwfc, iunigk, seqopn
  USE tddft_module,     ONLY : iuntdwfc, nwordtdwfc
  USE control_flags,    ONLY : io_level    
  implicit none 
  logical :: exst

  !
  ! ... nwordwfc is the record length (IN REAL WORDS)
  ! ... for the direct-access file containing wavefunctions
  ! ... io_level > 0 : open a file; io_level <= 0 : open a buffer
  !
  nwordwfc = nbnd*npwx*npol
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst )

  ! tddft_psi used for time propagation
  nwordtdwfc = nbnd*npwx*npol
  CALL open_buffer( iuntdwfc, 'tdwfc', nwordtdwfc, io_level, exst )

  ! iunigk contains the number of PW and the indices igk
  CALL seqopn( iunigk, 'igk', 'UNFORMATTED', exst )

END SUBROUTINE tddft_openfile
!---

!---
SUBROUTINE tddft_welcome()
  !---
  ! Print a short welcome summary of the calculation
  USE tddft_module
  USE io_global,    ONLY : stdout
  implicit none
  
  write(stdout,*)
  write(stdout,'(5X,''***************************'')')
  write(stdout,'(5X,''*** Welcome to RT-tddft ***'')')
  write(stdout,'(5X,''***************************'')')
  write(stdout,'(5X,''Calculation type      : '',A12)') trim(job)
  write(stdout,'(5X,''Solver                : '',A12)') trim(solver)
  write(stdout,'(5X,''Propagation method    : '',A12)') trim(method)
  write(stdout,'(5X,''Initial time step     : '',I12)') init_step
  write(stdout,'(5X,''Number or steps       : '',I12)') num_step
  write(stdout,'(5X,''Time step             : '',E12.4,'' Rydberg Atomic Time Unit'')') dt
  write(stdout,*)

  call flush_unit( stdout )

END SUBROUTINE tddft_welcome
!---

!---
SUBROUTINE tddft_closefile()
  !---
  ! Close wavefunction files
  !
  USE tddft_module,     ONLY : iuntdwfc
  USE io_files,         ONLY : iunwfc
  USE buffers,          ONLY : close_buffer
  implicit none

  call close_buffer( iunwfc, 'keep' )
  call close_buffer( iuntdwfc, 'keep' )

END SUBROUTINE tddft_closefile
!---

!---
SUBROUTINE tddft_print_clocks()
  !---
  ! Print clock information for tddft routines
  USE io_global,    ONLY : stdout
  USE tddft_module, ONLY : solver
  implicit none

  ! initialization
  write(stdout,*)
  write(stdout,'(5X,''***************************'')')
  write(stdout,'(5X,''*** Clock time RT-tddft ***'')')
  write(stdout,'(5X,''***************************'')')
  write(stdout,'(5X,''Initialization'')')
  call print_clock ('tddft_init')
  write(stdout,*)

  ! update
  write(stdout,'(5X,''Update rho and H'')')
  call print_clock ('tddft_update')
  write(stdout,*)
  
  ! time propagation
  write(stdout,'(5X,''Time propagation'')')
  call print_clock ('tddft_propagate')
  write(stdout,*)

  ! solver
  select case(trim(solver))
    case('cgsolver')
      write(stdout,'(5X,''CG square solver time'')')
      call print_clock ('cgsolver')
      write(stdout,*)
    case('itsolver')
      write(stdout,'(5X,''Iterative solver time'')')
      call print_clock ('itsolver')
      write(stdout,*)
    case default
  end select

  ! h_psi in iterations
  write(stdout,'(5X,''H|psi> in linear solver'')')
  call print_clock ('h_psi')
  write(stdout,*)

  ! compute num_elec and dipole
  write(stdout,'(5X,''Compute num_elec and dipole'')')
  call print_clock ('tddft_compute')
  write(stdout,*)

  call flush_unit( stdout )

END SUBROUTINE tddft_print_clocks
!---
