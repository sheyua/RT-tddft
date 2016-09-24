
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
#ifdef __PARA
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
