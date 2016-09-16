
!---
SUBROUTINE tddft_read_input()
  !---
  ! Read from stdin for the RT-tddft calculation parameters
  USE tddft_module
  USE io_files,         ONLY : prefix, tmp_dir  
  USE io_global,        ONLY : ionode, stdin, stdout
  USE constants,        ONLY : au_sec
  USE mp_images,        ONLY : my_image_id
  implicit none
  integer :: ios
  namelist /input_tddft/ job, prefix, tmp_dir, &
                        conv_threshold, dt, num_step, init_step, &
                        e_mirror, e_pstart, e_pend, e_nstart, e_nend, e_volt, e_decay
  
  ! define input defult values
  job          = 'transport'
  prefix       = 'pwscf'
  tmp_dir      = './scratch/'    
  conv_threshold = 1.0d-12                 ! convergence threshold    
  dt           = 2.d0                      ! time step (default: 2 attosecond)
  num_step     = 1000                      ! total time steps
  init_step    = 1                         ! initial step number
  e_mirror     = .true.
  e_pstart     = 0.0d0
  e_pend       = 0.25d0
  e_nstart     = 0.25d0
  e_nend       = 0.5d0
  e_volt       = 0.0d0
  e_decay      = 0.0d0

  if ( ionode .and. my_image_id == 0) then
    call input_from_file()
    ! read input    
    read( stdin, input_tddft, err = 100, iostat = ios )
100 call errore('tddft_readin', 'reading input_tddft namelist', abs(ios))
  
    ! convert to atomic units
    dt = dt * 1.d-18 / (2.d0*au_sec)           ! change from femto-second to a.u.
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
  call mp_bcast(conv_threshold, root, world_comm)
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
  USE io_files,         ONLY : iunwfc, nwordwfc, &
                               iunigk, seqopn
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
  write(stdout,'(5X,''Initial time step     : '',I12)') init_step
  write(stdout,'(5X,''Number or steps       : '',I12)') num_step
  write(stdout,'(5X,''Time step             : '',F12.4,'' Atomic Time Unit'')') dt
  write(stdout,*)

  call flush_unit( stdout )

END SUBROUTINE tddft_welcome
!---
  
!---
!SUBROUTINE tddft_init()
!END SUBROUTINE
!---




!-----------------------------------------------------------------------
SUBROUTINE tddft_closefil
  !-----------------------------------------------------------------------
  !
  ! ... Close files opened by TDDFT
  !
  USE tddft_module,     ONLY : iuntdwfc
  USE ldaU,             ONLY : lda_plus_U  
  USE io_files,         ONLY : prefix, iunhub, iunwfc
  USE buffers,          ONLY : close_buffer

  call close_buffer( iunwfc, 'keep' )
  call close_buffer( iuntdwfc, 'keep' )
  if ( lda_plus_u ) call close_buffer ( iunhub, status = 'keep' )

END SUBROUTINE tddft_closefil



!---
SUBROUTINE tddft_print_clocks()
  !---
  ! Print clock information for tddft routines
  USE io_global,  ONLY : stdout
  implicit none

  ! Initialization
  write(stdout,'(5X,''Initialization'')')
  call print_clock ('tddft_init')
  write(stdout,*)

!  !
!  write(stdout,*) '    Linear response'
!  call print_clock ('greenf')
!  call print_clock ('cgsolve')
!  call print_clock ('ch_psi')
!  call print_clock ('h_psi')
!  call print_clock ('s_psi')
!  write(stdout,*) '    Real time evolution'
!  call print_clock ('updateH')
!  call print_clock ('dipole')
!  write(stdout,*)
!  write(stdout,*) '    General routines'
!  call print_clock ('calbec')
!  call print_clock ('fft')
!  call print_clock ('ffts')
!  call print_clock ('fftw')
!  call print_clock ('cinterpolate')
!  call print_clock ('davcio')
!  call print_clock ('write_rec')
!  write(stdout,*)
!
!#ifdef __PARA
!  write(stdout,*) '    Parallel routines'
!  call print_clock ('reduce')  
!  call print_clock( 'fft_scatter' )
!  call print_clock( 'ALLTOALL' )
!  write(stdout,*)
!#endif
  call print_clock ('RT-tddft') 

END SUBROUTINE tddft_print_clocks
!---


!-----------------------------------------------------------------------
SUBROUTINE tddft_memory_report
  !-----------------------------------------------------------------------
  !
  ! ... Print estimated memory usage
  !
  USE io_global,                 ONLY : stdout
  USE noncollin_module,          ONLY : npol
  USE uspp,                      ONLY : okvan, nkb
  USE fft_base,                  ONLY : dffts
  USE pwcom
  IMPLICIT NONE
  integer, parameter :: Mb=1024*1024, complex_size=16, real_size=8

  ! the conversions to double prevent integer overflow in very large run
  write(stdout,'(5x,"Largest allocated arrays",5x,"est. size (Mb)",5x,"dimensions")')

  write(stdout,'(8x,"KS wavefunctions at k     ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     complex_size*nbnd*npol*DBLE(npwx)/Mb, npwx*npol,nbnd

  write(stdout,'(8x,"First-order wavefunctions ",f10.2," Mb",5x,"(",i8,",",i5,",",i3")")') &
     complex_size*nbnd*npol*DBLE(npwx)*10/Mb, npwx*npol,nbnd,10

  write(stdout,'(8x,"Charge/spin density       ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     real_size*dble(dffts%nnr)*nspin/Mb, dffts%nnr, nspin
  
  write(stdout,'(8x,"NL pseudopotentials       ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
     complex_size*nkb*DBLE(npwx)/Mb, npwx, nkb
  write(stdout,*)

END SUBROUTINE tddft_memory_report


