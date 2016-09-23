
!---
MODULE cgsolver_module
  !---
  ! Variables used in the linear cg solver
  USE kinds,    ONLY : dp
  implicit none
  save

  integer     :: flag
  real(dp)    :: tolb, normr, relres, n2b, normrmin
  complex(dp) :: rho, rho1, alpha, beta, rtvh
  complex(dp), allocatable :: r(:), Ax(:), rt(:), vh(:), u(:), &
                              uh(:), q(:), qh(:), p(:)

END MODULE cgsolver_module
!---
