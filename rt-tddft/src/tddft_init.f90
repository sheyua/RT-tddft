
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
  call set_rpos()
  call tddft_update(0, 0)
  call tddft_compute(0)

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
    max_nbnd_occ = 0
    do ik = 1, nks
      if (nbnd_occ(ik) > max_nbnd_occ) max_nbnd_occ = nbnd_occ(ik)
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
    allocate( tddft_psi(npwx,nbnd) )
    allocate( num_elec(nspin) ) 
    allocate( dipole(3,nspin) )
    allocate( rpos(3,dfftp%nnr) )
    
    ! initialize variables
    tddft_psi = (0.d0, 0.d0)
    num_elec = 0.d0
    dipole = 0.d0
    rpos = 0.d0

  END SUBROUTINE set_tddft_allocatable
  !---

  !---
  SUBROUTINE set_rpos()
    !---
    ! Setup the position operator in real space. The origin is set to center of ionic charges. (r is in units of alat)
    USE kinds,        ONLY : dp
    USE mp_global,    ONLY : me_pool
    USE fft_base,     ONLY : dfftp
    USE ions_base,    ONLY : nat, tau, ityp, zv
    USE cell_base,    ONLY : at, bg
    USE tddft_module, ONLY : rpos
    implicit none
    real(dp) :: zvtot, x0(3), rtmp(3), inv_nr(3)
    integer :: ia, i, j, k, index, index0, ir, ipol
  
    ! calculate the center of ionic charges
    zvtot = 0.d0
    x0 = 0.d0
    do ia = 1, nat
       zvtot = zvtot + zv(ityp(ia))
       x0(:) = x0(:) + tau(:,ia)*zv(ityp(ia))
    enddo
    x0 = x0 / zvtot
  
    ! grid density
    inv_nr(1) = 1.d0 / real(dfftp%nr1,dp)
    inv_nr(2) = 1.d0 / real(dfftp%nr2,dp)
    inv_nr(3) = 1.d0 / real(dfftp%nr3,dp)
  
    index0 = 0
#ifdef __MPI
    do i = 1, me_pool
      index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
    enddo
#endif
  
    ! loop over real space grid
    do ir = 1, dfftp%nnr
      index = index0 + ir - 1
      k     = index / (dfftp%nr1x*dfftp%nr2x)
      index = index - (dfftp%nr1x*dfftp%nr2x)*k
      j     = index / dfftp%nr1x
      index = index - dfftp%nr1x*j
      i     = index
  
      do ipol = 1, 3
        rtmp(ipol) = real(i,dp)*inv_nr(1)*at(ipol,1) + &
            real(j,dp)*inv_nr(2)*at(ipol,2) + &
            real(k,dp)*inv_nr(3)*at(ipol,3)
      enddo
  
      ! minimum image convenction
      rtmp = rtmp - x0
      call cryst_to_cart( 1, rtmp, bg, -1 )
      rtmp = rtmp - anint(rtmp)
      call cryst_to_cart( 1, rtmp, at, 1 )
      
      rpos(:,ir) = rtmp(:)
    enddo
  
  END SUBROUTINE set_rpos
  !---

END SUBROUTINE tddft_init
!---
