
!---
SUBROUTINE propgator_Euler(x, prop_x, delta_t, dim_m)
  !---
  ! Propagate array (or vector x) to be [I - i*Hdt]x 
  USE kinds,        ONLY : dp
  USE wvfct,        ONLY : npwx, npw
  implicit none
  integer, intent(in) :: dim_m   ! product will be performed on npw * dim_m
  real(dp), intent(in) :: delta_t
  complex(dp), intent(in) :: x(npwx, dim_m)
  complex(dp), intent(inout) :: prop_x(npwx, dim_m)
  complex(dp) :: e_coef

  ! compute the product of the hamiltonian H and the vector x
  e_coef = (0.d0, 1.d0) * delta_t
  
  ! call h_psi
  prop_x = (0.d0, 0.d0)
  call h_psi (npwx, npw, dim_m, x, prop_x)

  ! compute prop_x
  prop_x(1:npw, 1:dim_m) = x(1:npw, 1:dim_m) - e_coef * prop_x(1:npw, 1:dim_m)

  RETURN
END SUBROUTINE propgator_Euler
!---
