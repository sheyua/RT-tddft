
!---
MODULE itsolver_module
  !---
  ! Variables used in the linear iterative solver
  USE kinds,    ONLY : dp
  implicit none
  save

  complex(dp), allocatable :: y(:,:), y_swp(:,:)
  complex(dp) :: e_coef
  integer :: iter
  real(dp) :: residual

END MODULE itsolver_module
!---
