
!---
SUBROUTINE tddft_update(istep)
  !---
  ! Update the Hamiltonian
  USE kinds,                ONLY : dp
  USE klist,                ONLY : nks
  USE wavefunctions_module, ONLY : evc 
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE buffers,              ONLY : get_buffer
  USE tddft_module,         ONLY : e_mirror, e_pstart, e_pend, &
                                   e_nstart, e_nend, e_volt, e_decay
  USE extfield,             ONLY : emirror, epstart, epend, &
                                   enstart, enend, evolt
  USE scf,                  ONLY : vltot, rho, rho_core, rhog_core, v, vrs, kedtau
  USE fft_base,             ONLY : dfftp
  USE lsda_mod,             ONLY : nspin
  USE gvecs,                ONLY : doublegrid
  implicit none
  integer, intent(in) :: istep
  real(dp) :: dummy1, dummy2, dummy3, dummy4, dummy5, dummy6
  call start_clock('tddft_update')
  
  ! update dummy4 density frist
  ! for initialization step, read in evc for only one gamma point
  if( istep == 0 .and. nks == 1 ) then
    call get_buffer(evc, nwordwfc, iunwfc, 1)
  endif
  call sum_band()

  ! update external bias potential
  emirror   = e_mirror
  epstart   = e_pstart
  epend     = e_pend
  enstart   = e_nstart
  enend     = e_nend
  if ( istep == 0 ) then
    evolt = e_volt
    call add_efield(vltot, dummy1, rho%of_r, .true.)
  elseif ( abs(e_decay) < 1d-10 .or. istep <= 1.d0/e_decay ) then
    evolt = -e_volt * e_decay
    call add_efield(vltot, dummy1, rho%of_r, .true.)
    evolt = e_volt * (1.0d0 - istep*e_decay)
  else
    emirror = .false.
  endif
  
  ! calculate HXC-potential
  call v_of_rho( rho, rho_core, rhog_core, dummy2, dummy3, dummy4, dummy5, dummy1, dummy6, v )
  ! calculate total local potential (external + scf)
  call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    

  call stop_clock('tddft_update')
END SUBROUTINE tddft_update
!---
