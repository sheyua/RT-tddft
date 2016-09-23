
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
  USE extfield,     ONLY : evolt
  USE tddft_module, ONLY : rpos, num_elec, dipole, dt
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

  ! print legend
  if(istep == 0) then
    write(stdout,'(5X,A4,1X,A6,1X,A16,1X,A16,1X,A16)') 'Spin', 'istep', 'Vbias', 'Charge', 'dipole' 
  endif

  ! print observables
  do is = 1, nspin
    write(stdout,'(5X,I4,1X,I6,1X,ES16.8,1X,ES16.8,1X,3ES16.8)') is, istep, evolt, num_elec(is), dipole(:,is)
  enddo

  call flush_unit( stdout )

  call stop_clock('tddft_compute')
END SUBROUTINE tddft_compute
!---


!    !========================================================================
!    ! output rho%of_r() and v%of_r()
!    !========================================================================
!    if (istep == 1) then
!      index0 = 0
!#ifdef __PARA
!      write(filename, '(A8,I0)') "rho_",  me_pool
!        do idz = 1, me_pool
!          index0 = index0 + dfftp%npp(idz)
!        enddo
!#else
!      write(filename, '(A8)') "rho"
!#endif
!      fzcur = 32
!      open(unit=fzcur,file=trim(filename))
!      write(fzcur, *) "#Procid=", me_pool, nspin, dfftp%nnr, dfftp%nr1x, dfftp%nr2x
!    endif ! open charge density writer
!
!    if (istep == 1) then
!#ifdef __PARA
!      write(filename, '(A8,I0)') "v_",  me_pool
!#else
!      write(filename, '(A8)') "v_"
!#endif
!      fzfield = 33
!      open(unit=fzfield,file=trim(filename))
!      write(fzfield, *) "#Procid=", me_pool, nspin, dfftp%nnr, dfftp%nr1x, dfftp%nr2x
!    endif ! open zKS potential writer
!
!    write(fzcur,*) "TimeStep: ", istep
!    write(fzfield,*) "TimeStep: ", istep
!    idz = 0
!    id = 0
!    do
!      if (id >= dfftp%nnr ) exit
!
!      tmprho = 0.d0
!      tmpv = 0.d0
!      do idy = 1, dfftp%nr2x
!      do idx = 1, dfftp%nr1x
!        do is = 1, nspin
!          tmprho(is) = tmprho(is) + rho%of_r( id+1, is )
!          tmpv(is) = tmpv(is) + v%of_r( id+1, is )
!        enddo
!        id = id + 1
!      enddo
!      enddo
!
!      do is = 1, nspin
!        write(fzcur,*) index0+idz, is, tmprho(is)*volRatio
!        write(fzfield,*) index0+idz, is, tmpv(is)/(dfftp%nr2x * dfftp%nr1x)
!      enddo
!      idz = idz + 1
!    enddo
!    !========================================================================
!    ! output rho%of_r() and v%of_r()
!    !========================================================================
