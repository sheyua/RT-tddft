
!---
SUBROUTINE tddft_init()
  !---
  ! Setup RT-tddft calculation
  !
  USE tddft_module
  implicit none
  call start_clock ('tddft_init')

  call gk_sort_local()

  call set_nbnd_occ()

  call allocate_sum_band()

  call set_tddft_allocatable()

  call stop_clock('tddft_init')
CONTAINS

  !---
  SUBROUTINE gk_sort_local()
    !---
    ! Local subroutine of gk_sort()
    USE io_files,   ONLY : iunigk
    USE klist,      ONLY : nks, xk
    USE gvect,      ONLY : ngm, g
    USE wvfct,      ONLY : ecutwfc, npw, igk, g2kin
    USE cell_base,  ONLY : tpiba2
    implicit none
    integer :: ik

    rewind( iunigk )
    ! the following loop must NOT be called more than once in a run
    ! or else there will be problems with variable-cell calculations
    do ik = 1, nks
       ! g2kin is used here as work space
       call gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
       ! if there is only one k-point npw and igk stay in memory
       if ( nks > 1 ) write( iunigk ) igk
    enddo

  END SUBROUTINE gk_sort_local
  !---

  !---
  SUBROUTINE set_nbnd_occ()
    !---
    ! Fill out the number of occupied bands for each k point
    USE kinds,      ONLY : dp
    USE klist,      ONLY : nks, lgauss, ngauss, degauss, wk
    USE io_global,  ONLY : stdout
    USE constants,  ONLY : pi
    USE pwcom,      ONLY : ef
    USE wvfct,      ONLY : nbnd, et, wg
    implicit none
    real(dp) :: small, fac, xmax, emax
    integer :: ik, ibnd

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
  
  END SUBROUTINE set_nbnd_occ
  !---

  !---
  SUBROUTINE allocate_sum_band()
    !---
    ! Allocate memory needed by sum_band
    USE wvfct,  ONLY : btype, nbndx
    USE klist,  ONLY : nkstot
    USE becmod, ONLY : allocate_bec_type, becp
    USE wvfct,  ONLY : nbnd
    USE uspp,   ONLY : nkb
    implicit none
    
    allocate(btype(nbndx,nkstot))
    btype = 1
    call allocate_bec_type(nkb, nbnd, becp)
  
  END SUBROUTINE allocate_sum_band
  !---


  !---
  SUBROUTINE set_tddft_allocatable()
    !---
    ! Allocate and initialize all array variables
    USE wvfct,      ONLY : nbnd, npwx
    USE lsda_mod,   ONLY : nspin
    USE fft_base,   ONLY : dfftp
    implicit none
    
    ! allocate variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( tddft_psi (npwx,nbnd,2))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( charge(nspin)         ) 
    allocate( dipole(3,nspin)       )
    allocate( r_pos(3,dfftp%nnr)    )
    allocate( r_pos_s(3,dfftp%nnr)  )
    
    ! initialize variables
    tddft_psi = cmplx(0.d0,0.d0)
    charge = 0.d0
    dipole = 0.d0
    r_pos = 0.d0
    r_pos_s = 0.d0
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call molecule_setup_r
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END SUBROUTINE set_tddft_allocatable
  !---

END SUBROUTINE tddft_init
!---
