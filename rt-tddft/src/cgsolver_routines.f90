! This CG square routine is created by Xiaofeng Qian @ MIT
! 

!---
SUBROUTINE cgsolver_initialize()
  !---
  ! Allocate memory for the solver
  USE wvfct,            ONLY : npwx
  USE cgsolver_module
  implicit none
 
  allocate (r(npwx), rt(npwx), Ax(npwx), u(npwx), p(npwx), q(npwx), &
            qh(npwx), uh(npwx), vh(npwx) )

END SUBROUTINE cgsolver_initialize 
!---
 
!---
SUBROUTINE cgsolver_finalize()
  !---
  ! Deallocate memory
  USE cgsolver_module
  implicit none

  deallocate( r, rt, Ax, u, p, q, qh, uh, vh )

END SUBROUTINE cgsolver_finalize
!---


!---
SUBROUTINE cgsolver (A, b, x, tol, nbnd, dt, max_iter)
  !---
  ! ... Conjugate-Gradient Square method for solving:   A * x = b
  ! ... where: A*x is evaluated by subroutine 'A', and 'A' is implicit
  ! ... general square-matrix.
  USE kinds,            ONLY : dp
  USE wvfct,            ONLY : npwx, npw
  USE mp_pools,         ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum
  USE cgsolver_module
  implicit none
  integer, intent(in) :: max_iter           ! the number of bands
  integer, intent(in) :: nbnd               ! the number of bands
  real(dp), intent(in) :: tol               ! threshold for convergence
  real(dp), intent(in) :: dt                ! i*dt/2
  complex(dp), intent(in) :: b(npwx,nbnd)   ! input: the known term
  complex(dp), intent(out) :: x(npwx,nbnd)  ! the solution of the linear system
  external A                                ! the subroutine computing A*x
  !----------------------------------------------------------------------
  complex(dp), external :: zdotc
  real(dp), external    :: ddot
  integer :: imin, stag, i, ibnd
  call start_clock ('cgsolver')

  ! initialize module variables
  tolb = 0.d0
  n2b  = 0.d0
  relres = 0.d0
  normr= 0.d0
  rho  = 0.d0
  rho1 = 0.d0
  r    = (0.d0, 0.d0)
  Ax   = (0.d0, 0.d0)
  
  flag           = 1
  imin           = 0
  
  !----------------------------------------------------------------------
  ! loop over bands
  !----------------------------------------------------------------------
  do ibnd = 1, nbnd

     n2b = dble(zdotc(npw, b(1,ibnd), 1, b(1,ibnd), 1))
#ifdef __MPI
     call mp_sum(n2b, intra_pool_comm)
#endif
     n2b = dsqrt(n2b)
     tolb = tol * n2b

     call A(x(1,ibnd), Ax, dt, 1)
     
     r = b(:,ibnd) - Ax
     normr = dble(zdotc( npw, r, 1, r, 1))
#ifdef __MPI
     call mp_sum(normr, intra_pool_comm)
#endif  
     normr = dsqrt(normr)
     
     if (normr < tolb) then
        flag = 0
        cycle
     endif
     
     rt = r
     normrmin = normr
     stag = 0
     rho  = cmplx(1.d0, 0.d0)
     
     ! CG iteration
     do i = 1, max_iter
        rho1 = rho
        rho = zdotc(npw, rt, 1, r, 1)
#ifdef __MPI
        call mp_sum(rho, intra_pool_comm)
#endif
        
        if (rho == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag  = 2
           return
        endif

        if (i == 1) then
           u = r
           p = u
        else
           beta = rho / rho1
           if ( beta == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
              flag  = 3
              return
           endif
           u = r + beta * q
           p = u + beta * (q + beta * p)
        endif
        
        call A(p, vh, dt, 1)
        
        rtvh = zdotc(npw, rt, 1, vh, 1)
#ifdef __MPI
        call mp_sum(rtvh, intra_pool_comm)
#endif
        
        if (rtvh == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag = 4
           return
        else
           alpha = rho / rtvh
        endif
        
        if (alpha == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag = 5
           stag = 1
           return
        endif
        
        q  = u - alpha * vh
        
        uh = u + q
        
        x(:,ibnd)  = x(:,ibnd) + alpha * uh
        
        call A(x(1,ibnd), Ax, dt, 1)
        
        normr = sum(conjg(b(:,ibnd) - Ax) * (b(:,ibnd) - Ax))
#ifdef __MPI
        call mp_sum(normr, intra_pool_comm)
#endif
        
        if (normr <= tolb) then
           flag = 0
           exit
        endif
        
        if (stag == 1) then
           flag  = 5
           return
        endif
        
        if (normr < normrmin)  then
           normrmin = normr
           imin = i
        endif
        
        call A(uh, qh, dt, 1)
        
        r  = r - alpha * qh
        
     enddo ! i

  end do ! ibnd
  !----------------------------------------------------------------------
  ! end of the loop over bands
  !----------------------------------------------------------------------
  
  if (flag > 0) then
     call errore('cgsolver', 'RT-tddft::cgsolver cannot achieve convergence', flag)
     stop
  end if

  call stop_clock ('cgsolver')
  return
END SUBROUTINE cgsolver
