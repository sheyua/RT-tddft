!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine molecule_optical_absorption
  !----------------------------------------------------------------------
  !  ... Compute optical absorption spectrum by real-time TDDFT 
  !----------------------------------------------------------------------
  USE kinds,                       ONLY : dp
  USE io_global,                   ONLY : stdout, ionode
  USE io_files,                    ONLY : nwordwfc, iunwfc, iunigk
  USE ions_base,                   ONLY : nat, ntyp => nsp, ityp
  USE cell_base,                   ONLY : at, bg, omega, tpiba, tpiba2, alat
  USE wavefunctions_module,        ONLY : evc
  USE klist,                       ONLY : nks, nkstot, wk, xk, nelec, ngk
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k, ecutwfc
  USE lsda_mod,                    ONLY : current_spin, lsda, isk, nspin
  USE becmod,                      ONLY : becp
  USE mp_pools,                    ONLY : my_pool_id, inter_pool_comm, intra_pool_comm
  USE mp_global,                   ONLY : me_pool
  USE mp,                          ONLY : mp_sum, mp_barrier
  USE gvect,                       ONLY : ngm, g
  USE gvecs,                       ONLY : nls
  USE fft_base,                    ONLY : dfftp
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE fixed_occ,                   ONLY : tfixed_occ 
  USE uspp,                        ONLY : nkb, vkb, deeq
  USE ldaU,                        ONLY : lda_plus_U
  USE uspp_param,                  ONLY : nh
  USE scf,                         ONLY : rho, rho_core, rhog_core, vltot, v, vrs
  USE control_flags,               ONLY : tqr
  USE tddft_module

  IMPLICIT NONE

  !-- tddft variables ----------------------------------------------------
  complex(dp), allocatable :: b(:,:)
  complex(dp), allocatable :: tddft_hpsi(:,:), tddft_spsi(:,:)

  integer :: istep, lter, flag_global
  integer :: ik, is, ibnd
  complex(dp) :: ee                     ! i*dt/2
  real(dp) :: anorm, volRatio
  integer, external :: find_free_unit
  external tddft_ch_psi_all

  volRatio = omega / real(dfftp%nr1 * dfftp%nr2 * dfftp%nr3, dp) 

  ! TODO: gk_sort

  ! allocate memory
  call allocate_optical()

  ee = i_complex * dt / 2.d0  ! i*dt/2: do not change
  
  evc = cmplx(0.d0,0.d0)
  call tddft_cgsolver_initialize(npwx, nbnd_occ_max)
 
  ! print the legend
  if (ionode) call print_legend
  
  ! set wfc and ham to the most recent time
  if (nks > 1) rewind (iunigk)
  do ik = 1, nks
     current_k = ik
     current_spin = isk(ik)
     
     ! initialize at k-point k 
     call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
     g2kin = g2kin * tpiba2
     call init_us_2(npw, igk, xk(1,ik), vkb)
     
     ! read wfcs from file and compute becp
     evc = (0.d0, 0.d0)
     call get_buffer (evc, nwordwfc, iunwfc, ik)
  end do
  call update_hamiltonian(-1)
 
  call flush_unit(stdout)


  ! enter the main TDDFT loop 
  do istep = 1, num_step
     
    ! calculate dipole moment along x, y, and z direction
    call molecule_compute_dipole( charge, dipole )


    ! loop over k-points     
    if (nks > 1) rewind (iunigk)

    do ik = 1, nks
      current_k = ik
      current_spin = isk(ik)
        
      ! initialize at k-point k 
      call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      g2kin = g2kin * tpiba2
      call init_us_2(npw, igk, xk(1,ik), vkb)
        
      ! read wfcs from file and compute becp
      evc = (0.d0, 0.d0)
      call get_buffer (evc, nwordwfc, iunwfc, ik)
      if (istep > 1) then
        call get_buffer (tddft_psi, nwordtdwfc, iuntdwfc, ik)
      endif

      ! guess the wavefunction at the next timestep
      tddft_psi(:,:,1) = (0.d0,0.d0)
      do ibnd = 1, nbnd_occ(ik)
        tddft_psi(:,ibnd,1) = 2.d0*evc(:,ibnd) - tddft_psi(:,ibnd,2)
      enddo

      ! calculate H |psi_current>, S |psi_current>
      call h_psi(npwx, npw, nbnd_occ(ik), evc, tddft_hpsi)
      call s_psi(npwx, npw, nbnd_occ(ik), evc, tddft_spsi)
        
      ! calculate (S - H*dt*i/2) |\psi_current>
      b = (0.d0, 0.d0)
      b(1:npw, 1:nbnd_occ(ik)) = tddft_spsi(1:npw,1:nbnd_occ(ik)) - ee * tddft_hpsi(1:npw,1:nbnd_occ(ik))
        
      ! solve A * x = b
      call tddft_cgsolver(tddft_ch_psi_all, b, tddft_psi(:,:,1), npwx, npw, &
                          conv_threshold, ik, lter, flag_global, anorm, &
                          nbnd_occ(ik), ee)
        
      ! update the wavefunctions
      tddft_psi(:,:,2) = tddft_psi(:,:,1)
      evc(:,1:nbnd_occ(ik)) = tddft_psi(:,1:nbnd_occ(ik),1)

      ! save wavefunctions to disk
      ! for multiple k-points, sum_band will use the evc in file iunwfc to compute
      ! charge
      CALL save_buffer (evc, nwordwfc, iunwfc, ik)
      call save_buffer (tddft_psi, nwordtdwfc, iuntdwfc, ik)
        
    enddo ! ik

#ifdef __PARA
    ! reduce over k-points
    call mp_sum(charge, inter_pool_comm)
    call mp_sum(dipole, inter_pool_comm)
#endif

    ! print observables
    if (ionode) then
      do is = 1, nspin
        write(stdout,'(''CHARGE '',I1,1X,I6,3E16.6)') is, istep, charge(is)
        write(stdout,'(''DIP    '',I1,1X,I6,3E16.6)') is, istep, dipole(:,is)
      enddo
    endif
     
    ! update the hamiltonian (recompute charge and potential)
    call update_hamiltonian(istep)
     
    call flush_unit(stdout)
     
  enddo      ! end of TDDFT loop

  ! finish  
  call tddft_cgsolver_finalize()
  call deallocate_optical()
  
    
CONTAINS

  !====================================================================
  ! Print the legend key
  !====================================================================    
  SUBROUTINE print_legend
    write(stdout,'(5X,''Output quantities:'')')
    write(stdout,'(5X,''  CHARGE spin  istep  charge'')')
    write(stdout,'(5X,''  DIP    spin  istep  dipole(1:3)'')')
    write(stdout,*)
    call flush_unit(stdout)
  END SUBROUTINE print_legend

  
  !====================================================================
  ! Initialize and allocate memory
  !====================================================================    
  SUBROUTINE allocate_optical()
    IMPLICIT NONE
    
    allocate (tddft_hpsi(npwx,nbnd_occ_max))
    allocate (tddft_spsi(npwx,nbnd_occ_max))
    allocate (b(npwx,nbnd_occ_max))
    tddft_hpsi = (0.d0,0.d0)
    tddft_spsi = (0.d0,0.d0)
    b = (0.d0,0.d0)
    
  END SUBROUTINE allocate_optical
  
  
  !====================================================================
  ! Deallocate memory
  !====================================================================    
  SUBROUTINE deallocate_optical()
    IMPLICIT NONE

    deallocate (tddft_psi, tddft_hpsi, tddft_spsi, b)
    deallocate (charge, dipole)
    deallocate (r_pos, r_pos_s)
  END SUBROUTINE deallocate_optical
   
   
  !========================================================================
  ! output rho%of_r() and v%of_r()
  !========================================================================
!  SUBROUTINE report_rho_and_v(filename)
!    integer :: index0, idx, idy, idz, id
!    integer :: fzcur, fzfield
!    character(len=1024) :: filename
!    real(dp), allocatable :: charge(:), dipole(:,:), quadrupole(:,:,:), tmprho(:), tmpv(:)
!    allocate (charge(nspin), dipole(3,nspin), quadrupole(3,3,nspin), tmprho(nspin), tmpv(nspin))
!    if (istep == 1) then
!      index0 = 0
!#ifdef __PARA
!      write(filename, '(A8,I0)') "rho_",  me_pool
!        do idz = 1, me_pool
!          index0 = index0 + dfftp%npp(idz)
!        enddo
!#else
!      write(filename, '(A8)') "rho"
!#endif
!      fzcur = 32
!      open(unit=fzcur,file=trim(filename))
!      write(fzcur, *) "#Procid=", me_pool, nspin, dfftp%nnr, dfftp%nr1x, dfftp%nr2x
!    endif ! open charge density writer
!
!    if (istep == 1) then
!#ifdef __PARA
!      write(filename, '(A8,I0)') "v_",  me_pool
!#else
!      write(filename, '(A8)') "v"
!#endif
!      fzfield = 33
!      open(unit=fzfield,file=trim(filename))
!      write(fzfield, *) "#Procid=", me_pool, nspin, dfftp%nnr, dfftp%nr1x, dfftp%nr2x
!    endif ! open zKS potential writer
!
!    write(fzcur,*) "TimeStep: ", istep
!    write(fzfield,*) "TimeStep: ", istep
!    idz = 0
!    id = 0
!    do
!      if (id >= dfftp%nnr ) exit
!
!      tmprho = 0.d0
!      tmpv = 0.d0
!      do idy = 1, dfftp%nr2x
!      do idx = 1, dfftp%nr1x
!        do is = 1, nspin
!          tmprho(is) = tmprho(is) + rho%of_r( id+1, is )
!          tmpv(is) = tmpv(is) + v%of_r( id+1, is )
!        enddo
!        id = id + 1
!      enddo
!      enddo
!
!      do is = 1, nspin
!        write(fzcur,*) index0+idz, is, tmprho(is)*volRatio
!        write(fzfield,*) index0+idz, is, tmpv(is)/(dfftp%nr2x * dfftp%nr1x)
!      enddo
!      idz = idz + 1
!    enddo
!  END SUBROUTINE report_rho_and_v

END SUBROUTINE molecule_optical_absorption
