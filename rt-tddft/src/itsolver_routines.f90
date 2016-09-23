
!---
SUBROUTINE itsolver_initialize(max_nbnd_occ)
  !---
  ! Allocate memory for the solver
  USE wvfct,            ONLY : npwx
  USE itsolver_module
  implicit none
  integer, intent(in) :: max_nbnd_occ
 
  allocate( y(npwx, max_nbnd_occ), y_swp(npwx, max_nbnd_occ) )

END SUBROUTINE itsolver_initialize 
!---
 
!---
SUBROUTINE itsolver_finalize()
  !---
  ! Deallocate memory
  USE itsolver_module
  implicit none

  deallocate( y, y_swp )

END SUBROUTINE itsolver_finalize
!---

!---
SUBROUTINE itsolver (b, x, tol, nbnd, delta_t, max_iter)
  !---
  ! Solves [ I + idtH ]x = b
  USE kinds,            ONLY : dp
  USE wvfct,            ONLY : npwx, npw
  USE mp_pools,         ONLY : intra_pool_comm
  USE mp,               ONLY : mp_max
  USE itsolver_module
  implicit none
  integer, intent(in) :: max_iter
  integer, intent(in) :: nbnd
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: delta_t
  complex(dp), intent(in) :: b(npwx, nbnd)
  complex(dp), intent(inout) :: x(npwx, nbnd)
  call start_clock ('itsolver')

  e_coef = (0.d0, 1.d0) * delta_t

  ! initial step
  ! compute y(0)
  call h_psi(npwx, npw, nbnd, x, y_swp)
  y(1:npw, 1:nbnd) = b(1:npw, 1:nbnd) - &
    e_coef*y_swp(1:npw, 1:nbnd) - x(1:npw, 1:nbnd)
  ! compute x(1)
  x(1:npw, 1:nbnd) = x(1:npw, 1:nbnd) + y(1:npw, 1:nbnd)

  ! start the iterative loop
  do iter = 1, max_iter
    ! compute y(iter)
    call h_psi(npwx, npw, nbnd, y, y_swp)
    y(1:npw, 1:nbnd) = -e_coef*y_swp(1:npw, 1:nbnd)

    ! compute x(iter+1)
    x(1:npw, 1:nbnd) = x(1:npw, 1:nbnd) + y(1:npw, 1:nbnd)

    ! compute residual(iter)
    residual = maxval(abs(y(1:npw, 1:nbnd)))
#ifdef __PARA
    mp_max(residual, intra_pool_comm)
#endif

    ! exit_flag
    if(residual < tol) then
      exit
    endif
  enddo

  if(residual >= tol) then
    call errore('itsolver', 'RT-tddft::itsolver cannot achieve convergence',1)
  endif

  call stop_clock ('itsolver')
  RETURN
END SUBROUTINE itsolver
!---
