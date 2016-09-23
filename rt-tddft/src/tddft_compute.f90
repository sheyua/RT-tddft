
!---
SUBROUTINE tddft_compute(istep)
  USE kinds,        ONLY : dp
  USE fft_base,     ONLY : dfftp
  USE lsda_mod,     ONLY : nspin
  USE scf,          ONLY : rho
  USE cell_base,    ONLY : omega, alat
  USE mp_global,    ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE io_global,    ONLY : stdout
  USE tddft_module, ONLY : rpos, num_elec, dipole
  implicit none
  integer, intent(in) :: istep
  integer :: is, ipol
  real(dp) :: volRatio
  call start_clock('tddft_compute')
  
  ! grid volume ratio
  volRatio = omega / real(dfftp%nr1 * dfftp%nr2 * dfftp%nr3, dp)
  
  ! compute charge and dipole
  do is = 1, nspin
    num_elec(is) = sum(rho%of_r(:,is))
    do ipol = 1, 3
      dipole(ipol, is) = sum(rpos(ipol,:)*rho%of_r(:,is))
    enddo
  enddo

  num_elec = num_elec * volRatio
  dipole = dipole * alat * volRatio

#ifdef __PARA
  call mp_sum(num_elec, intra_pool_comm)
  call mp_sum(dipole, intra_pool_comm)
#endif

    ! print observables
    write(stdout, *) 'CHARGE', 1, istep, num_elec(1)
    write(stdout,'(''DIP    '',I1,1X,I6,3E16.6)') 1, istep, dipole(:,1)

  call stop_clock('tddft_compute')
END SUBROUTINE tddft_compute
!---
