
!---
SUBROUTINE tddft_compute(istep)
  !---
  ! Compute charge density, dipole moment, dump charge cdf and Kohn-Sham
  ! potential
  USE kinds,        ONLY : dp
  USE fft_base,     ONLY : dfftp
  USE lsda_mod,     ONLY : nspin
  USE scf,          ONLY : rho
  USE cell_base,    ONLY : omega, alat
  USE mp_global,    ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE io_global,    ONLY : stdout
  USE extfield,     ONLY : evolt
  USE tddft_module, ONLY : rpos, num_elec, dipole, dt, dump, init_step
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

#ifdef __MPI
  call mp_sum(num_elec, intra_pool_comm)
  call mp_sum(dipole, intra_pool_comm)
#endif

  ! print legend
  if(istep == 0) then
    write(stdout,'(5X,A4,1X,A6,1X,A16,1X,A16,1X,A16)') 'Spin', 'istep', 'Vbias', 'Charge', 'dipole' 
    ! print observables
    do is = 1, nspin
      write(stdout,'(5X,I4,1X,I6,1X,ES16.8,1X,ES16.8,1X,3ES16.8)') is, init_step-1, evolt, num_elec(is), dipole(:,is)
    enddo
  else
    ! print observables
    do is = 1, nspin
      write(stdout,'(5X,I4,1X,I6,1X,ES16.8,1X,ES16.8,1X,3ES16.8)') is, istep, evolt, num_elec(is), dipole(:,is)
    enddo
  endif
  call flush_unit( stdout )

  ! dump rho and vks
  if(dump) then
    call dump_rho_vks()
  endif

  call stop_clock('tddft_compute')
CONTAINS

  !---
  SUBROUTINE dump_rho_vks()
    !---
    ! dump rho vks
    USE mp_global,      ONLY : me_pool
    USE tddft_module,   ONLY : iuntdrho, iuntdvks
    USE constants,      ONLY : AU_SEC, BOHR_RADIUS_ANGS
    USE scf,            ONLY : v
    USE cell_base,      ONLY : at
    implicit none
    integer :: index0, idx, idy, idz, id_xy
    real(dp) :: tmprho(2), tmpvks(2)

    index0 = 0
#ifdef __MPI
    do idz = 1, me_pool
      index0 = index0 + dfftp%npp(idz)
    enddo
#endif

    ! print legend
    if(istep == 0) then
      write(iuntdrho, '(''#ProcID= '',I6,1X,I1,1X,ES16.8,1X,I6,1X,I6,1X,ES16.8)') me_pool, nspin, dt*(2.d18*AU_SEC), &
        index0, dfftp%nr3, sqrt(sum(at(3,:)*at(3,:)))*alat*BOHR_RADIUS_ANGS
      write(iuntdvks, '(''#ProcID= '',I6,1X,I1,1X,ES16.8,1X,I6,1X,I6,1X,ES16.8)') me_pool, nspin, dt*(2.d18*AU_SEC), &
        index0, dfftp%nr3, sqrt(sum(at(3,:)*at(3,:)))*alat*BOHR_RADIUS_ANGS
      write(iuntdrho,'(''TimeStep: '',I8,1X,ES16.8)') init_step-1, evolt
      write(iuntdvks,'(''TimeStep: '',I8,1X,ES16.8)') init_step-1, evolt
    else
      write(iuntdrho,'(''TimeStep: '',I8,1X,ES16.8)') istep, evolt
      write(iuntdvks,'(''TimeStep: '',I8,1X,ES16.8)') istep, evolt
    endif

    ! compute rho and vks
    idz = 0
    id_xy = 0
    do 
      if(id_xy >= dfftp%nnr) exit

      tmprho = 0.d0
      tmpvks = 0.d0
      do idy = 1, dfftp%nr2x
      do idx = 1, dfftp%nr1x
        id_xy = id_xy + 1
        do is = 1, nspin
          tmprho(is) = tmprho(is) + rho%of_r(id_xy, is)
          tmpvks(is) = tmpvks(is) + v%of_r(id_xy, is)
        enddo
      enddo
      enddo

      do is = 1, nspin
        write(iuntdrho,'(I1,1X,ES16.8,1X,I0)') is, tmprho(is)*volRatio, index0+idz
        write(iuntdvks,'(I1,1X,ES16.8,1X,I0)') is, tmpvks(is)/(dfftp%nr2*dfftp%nr1), index0+idz
      enddo
      idz = idz + 1

    enddo

  END SUBROUTINE dump_rho_vks
  !---

END SUBROUTINE tddft_compute
!---
