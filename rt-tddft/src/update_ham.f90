!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE update_hamiltonian(istep)
  !-----------------------------------------------------------------------
  !
  ! ... Update the hamiltonian
  !
  USE kinds,         ONLY : dp
  USE ldaU,          ONLY : lda_plus_U
  USE scf,           ONLY : rho, rho_core, rhog_core, vltot, v, kedtau, vrs
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE io_global,     ONLY : stdout
  USE lsda_mod,      ONLY : nspin
  USE uspp,          ONLY : okvan, nkb
  USE dfunct,        ONLY : newd
  USE tddft_module,  ONLY : e_mirror, e_pstart, e_pend, e_nstart, e_nend, &
                            e_volt, e_decay
  USE becmod,        ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE wvfct,         ONLY : nbnd
  USE extfield,      ONLY : tefield, emirror, epstart, epend, enstart, enend, evolt
  implicit none
  integer, intent(in) :: istep
  real(dp) :: charge, ehart, etxc, vtxc, eth, etotefield

  call start_clock('updateH')
  
  ! calculate total charge density
  call deallocate_bec_type(becp)
  call sum_band()
  call allocate_bec_type(nkb, nbnd, becp)

  if (lda_plus_U) then
    call new_ns
  end if
    
  etotefield = 0.d0 ! etotefield is INOUT type
  epstart = e_pstart
  epend   = e_pend
  enstart = e_nstart
  enend   = e_nend
  ! manually add the external bias potential into vltot
  if ( istep == -1 ) then
    tefield = .true.
    emirror = e_mirror
    evolt = e_volt
    call add_efield(vltot, etotefield, rho%of_r, .true.)
  elseif ( istep <= 1.0d0/e_decay ) then
    tefield = .true.
    emirror = e_mirror
    evolt = -e_volt * e_decay
    call add_efield(vltot, etotefield, rho%of_r, .true.)
    evolt = e_volt * (1.0d0 - istep*e_decay)
  else
    tefield = .false.
    emirror = .false.
  endif
  write(stdout, *) 'GREP_TAG_after_add_efi', istep, evolt
  
  ! calculate HXC-potential
  call v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, etotefield, charge, v )

  ! calculate total local potential (external + scf)
  call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    
  write(stdout, *) 'GREP_TAG_after_set_vrs', istep, vrs(8019,1)
  
  ! calculate new D_nm matrix for ultrasoft pseudopotential
  if (okvan) then
  endif
  ! disable ultrasoft potentials
    
  call stop_clock('updateH')
    
END SUBROUTINE update_hamiltonian
