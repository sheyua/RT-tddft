
!---
SUBROUTINE tddft_init()
  !---
  ! Setup RT-tddft calculation
  !
  USE tddft_module
  USE kinds,            ONLY : dp
  USE wvfct,            ONLY : btype, nbndx
  USE klist,            ONLY : nkstot, lgauss, ngauss, degauss, nks, wk, xk
  USE io_global,        ONLY : stdout
  USE constants,        ONLY : pi
  USE pwcom,            ONLY : ef
  USE wvfct,            ONLY : nbnd, et, wg, ecutwfc, npw, npwx, igk, g2kin
  USE io_files,         ONLY : iunigk
  USE gvect,            ONLY : ngm, g
  USE cell_base,        ONLY : tpiba2
  USE fft_base,         ONLY : dfftp
  USE lsda_mod,         ONLY : nspin
  implicit none
  integer :: ik, ibnd
  real(dp) :: small, fac, xmax, emax
  call start_clock ('tddft_init')

  ! allocate memory needed by sum_band
  allocate(btype(nbndx,nkstot))
  btype = 1

  ! compute the numbers of occupied bands for each k points
  allocate(nbnd_occ(nks))
  if (lgauss) then
    write(stdout,*)
    write(stdout,'(5X,''Retrieve electronic occupation informations:'')')
    write(stdout,'(5X,''smearing method ngauss='',I4,2X,''degauss='',F8.4,'' Ry'')') ngauss, degauss
    write(stdout,*)

    ! discard conduction bands such that w0gauss(x,n) < small
    !small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
    !small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
    !small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
    small = 6.3491173359333d-8

    ! appropriate limit for different broadening
    if (ngauss == -99) then
      fac = 1.d0 / sqrt(small)
      xmax = 2.d0 * log(0.5d0*(fac + sqrt(fac*fac-4.d0)))
    else
      xmax = sqrt(-log(sqrt(pi)*small))
    endif
    emax = ef + xmax * degauss

    ! fill out number of bands occupied
    do ik = 1, nks
      do ibnd = 1, nbnd
         if (et(ibnd,ik) < emax) nbnd_occ(ik) = ibnd
      enddo
    enddo
  else
    do ik = 1, nks
      if ( ABS(wk(ik)) > 0.0d0 ) then
        do ibnd = 1, nbnd
          if ( ABS(wg(ibnd,ik)/wk(ik)) > 0.0d0 ) nbnd_occ(ik) = ibnd
        enddo
      endif
    enddo
  endif

  ! compute the maximum number of occupied bands
  nbnd_occ_max = 0
  do ik = 1, nks
    if (nbnd_occ(ik) > nbnd_occ_max) nbnd_occ_max = nbnd_occ(ik)
  enddo

  REWIND( iunigk )
  ! the following loop must NOT be called more than once in a run
  ! or else there will be problems with variable-cell calculations
  DO ik = 1, nks
     ! g2kin is used here as work space
     CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
     ! if there is only one k-point npw and igk stay in memory
     IF ( nks > 1 ) WRITE( iunigk ) igk
  END DO

  ! initialize allocatable variables
  allocate(tddft_psi (npwx,nbnd,2))
  allocate(charge(nspin), dipole(3,nspin))
  allocate(r_pos(3,dfftp%nnr), r_pos_s(3,dfftp%nnr))
  call allocate_memory()

  call stop_clock('tddft_init')

CONTAINS

  !---
  SUBROUTINE allocate_memory()
    implicit none
    tddft_psi = (0.d0,0.d0)
    charge = 0.d0
    dipole = 0.d0
    call molecule_setup_r
  END SUBROUTINE allocate_memory
  !---

END SUBROUTINE tddft_init
!---
