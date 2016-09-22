
!---
SUBROUTINE tddft_propagate()
  !---
  ! Real-time propagation of the electronic states
  USE kinds,                ONLY : dp
  USE wvfct,                ONLY : npwx, nbnd, npw
  USE io_files,             ONLY : iunigk
  USE klist,                ONLY : nks
  USE wavefunctions_module, ONLY : evc
  USE tddft_module
  implicit none
  complex(dp) :: dt_cmplx
  integer :: istep, ik
  complex(dp), allocatable :: b(:,:)
  external propgator_Euler
  call start_clock ('tddft_propagate')

  ! setup parameters
  call cgsolver_initialize()
  allocate(b(npwx, max_nbnd_occ))

  ! main loop
  do istep = init_step, (init_step+num_step-1)

!    ! CN itegrator
!    if(nks > 1) rewind(iunigk)
!    do ik = 1, nks
! 
!      ! init potential inorder to call h_psi
!      call init_k()
!
!      ! prepare b
!      ! b = [I - iH(t)dt/hbar]*psi(t-dt)
!      call propgator_Euler(evc, b, 0.5d0*dt, nbnd_occ(ik))
!
!      ! solve A * x = b
!      tddft_psi(1:npw,1:nbnd_occ(ik)) = 2.d0*evc(1:npw,1:nbnd_occ(ik)) - &
!        tddft_psi(1:npw,1:nbnd_occ(ik))
!      call cgsolver(propgator_Euler, b, tddft_psi, conv_threshold, nbnd_occ(ik), &
!        -0.5d0*dt, max_iter)
!
!      !call iterative_solver()
!      call save_k()
!    enddo
!    call tddft_update(istep, 1)
!    call tddft_compute(istep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! second half integrator
    if(nks > 1) rewind(iunigk)
    do ik = 1, nks
 
      ! init potential inorder to call h_psi
      call init_k()

      ! prepare b
      ! b = [I - iH(t)0.25dt/hbar]*psi(t)
      call propgator_Euler(evc, b, 0.25d0*dt, nbnd_occ(ik))
      !b = (0.d0, 0.d0)
      !call h_psi(npwx, npw, nbnd_occ(ik), evc, b)
      !b(1:npw, 1:nbnd_occ(ik)) = evc(1:npw, 1:nbnd_occ(ik)) - &
      !  (0.d0, 1.d0)*dt*0.25d0*b(1:npw, 1:nbnd_occ(ik))

      ! compute psi(t+0.5dt)
      tddft_psi(1:npw,1:nbnd_occ(ik)) = 2.d0*evc(1:npw,1:nbnd_occ(ik)) - &
        tddft_psi(1:npw,1:nbnd_occ(ik))
      call cgsolver(propgator_Euler, b, tddft_psi, conv_threshold, nbnd_occ(ik), &
        -0.25d0*dt, max_iter)

      call save_k()
    enddo

    call tddft_update(istep, 2)

    ! second half integrator
    if(nks > 1) rewind(iunigk)
    do ik = 1, nks
 
      ! init potential inorder to call h_psi
      call init_k()

      ! prepare b
      ! b = [I - iH(t+0.5dt)dt/hbar]*psi(t)
      call propgator_Euler(tddft_psi, b, 0.5d0*dt, nbnd_occ(ik))
      !b = (0.d0, 0.d0)
      !call h_psi(npwx, npw, nbnd_occ(ik), tddft_psi, b)
      !b(1:npw, 1:nbnd_occ(ik)) = tddft_psi(1:npw, 1:nbnd_occ(ik)) - &
      !  (0.d0, 1.d0)*dt*0.5d0*b(1:npw, 1:nbnd_occ(ik))

      ! compute psi(t+dt)
      tddft_psi(1:npw,1:nbnd_occ(ik)) = 2.d0*evc(1:npw,1:nbnd_occ(ik)) - &
        tddft_psi(1:npw,1:nbnd_occ(ik))
      call cgsolver(propgator_Euler, b, tddft_psi, conv_threshold, nbnd_occ(ik), &
        -0.5d0*dt, max_iter)

      call save_k()
    enddo

    call tddft_update(istep, 3)
    call tddft_compute(istep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  enddo

  ! deallocate  parameters
  call cgsolver_finalize()
  deallocate(b)

  call stop_clock('tddft_propagate')
CONTAINS

  !---
  SUBROUTINE init_k()
    !---
    ! Initialize at k-point k
    USE klist,                  ONLY : xk
    USE gvect,                  ONLY : ngm, g
    USE wvfct,                  ONLY : ecutwfc, npw, igk, g2kin
    USE cell_base,              ONLY : tpiba2
    USE uspp,                   ONLY : vkb
    USE buffers,                ONLY : get_buffer
    USE wavefunctions_module,   ONLY : evc
    USE io_files,               ONLY : nwordwfc, iunwfc
    USE tddft_module,           ONLY : tddft_psi, nwordtdwfc, iuntdwfc
    implicit none

    ! initialize potential
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin = g2kin * tpiba2
    call init_us_2(npw, igk, xk(1,ik), vkb)

    ! initialize wavefunctions
    call get_buffer(evc, nwordwfc, iunwfc, ik)
    if (istep == 1) then
      call get_buffer(tddft_psi, nwordwfc, iunwfc, ik)
    elseif (istep > 1) then
      call get_buffer(tddft_psi, nwordtdwfc, iuntdwfc, ik)
    endif

  END SUBROUTINE init_k
  !---

  !---
  SUBROUTINE save_k()
    !---
    ! Save evc and tddft_psi into the right buffer
    USE buffers,                ONLY : save_buffer
    USE io_files,               ONLY : nwordwfc, iunwfc
    USE wavefunctions_module,   ONLY : evc
    USE klist,                  ONLY : nks
    implicit none

    ! evc(t+1)
    call save_buffer(tddft_psi, nwordwfc, iunwfc, ik)
    ! tddft_psi(t+1) = evc(t)
    call save_buffer(evc, nwordtdwfc, iuntdwfc, ik)
    
    ! needed for one gamma point sum_band
    if(nks == 1) then
        evc(:,:) = tddft_psi(:,:)
    endif

  END SUBROUTINE save_k
  !---

END SUBROUTINE tddft_propagate
!---



!---
!  SUBROUTINE iterative_solver()
!    !---
!    ! Solve psi(t+dt) from H(t) and b
!    ! [I + iH(t)dt/hbar]*psi(t+dt) = b
!    USE wvfct,                  ONLY : npwx, npw
!    USE wavefunctions_module,   ONLY : evc
!    USE io_global,              ONLY : stdout
!    implicit none
!    integer :: num_iter_ik
!
!    ! prepare b
!    ! b = [I - iH(t)dt/hbar]*psi(t-dt)
!    call h_psi(npwx, npw, nbnd_occ(ik), tddft_psi, tddft_hpsi)
!    b(1:npw, 1:nbnd_occ(ik)) = tddft_psi(1:npw, 1:nbnd_occ(ik)) - &
!        dt_cmplx * tddft_hpsi(1:npw, 1:nbnd_occ(ik))
!
!    ! guess wavefunction psi(t+dt)
!    call h_psi(npwx, npw, nbnd_occ(ik), evc, tddft_hpsi)
!    tddft_psi(1:npw, 1:nbnd_occ(ik)) = b(1:npw, 1:nbnd_occ(ik)) - &
!        dt_cmplx * tddft_hpsi(1:npw, 1:nbnd_occ(ik))
!
!    ! compute increment y
!    y(1:npw, 1:nbnd_occ(ik)) = tddft_psi(1:npw, 1:nbnd_occ(ik)) - &
!        evc(1:npw, 1:nbnd_occ(ik))
!
!    ! start iterations
!    num_iter_ik = 0
!    do while(maxval(abs(y(1:npw, 1:nbnd_occ(ik)))) > conv_threshold)
!      ! compute diff(psi)
!      call h_psi(npwx, npw, nbnd_occ(ik), y, tddft_hpsi)
!      y(1:npw, 1:nbnd_occ(ik)) = -dt_cmplx * tddft_hpsi(1:npw, 1:nbnd_occ(ik))
!      
!      ! update tddft_psi
!      tddft_psi(1:npw, 1:nbnd_occ(ik)) = tddft_psi(1:npw, 1:nbnd_occ(ik)) + &
!        y(1:npw, 1:nbnd_occ(ik))
!      num_iter_ik = num_iter_ik + 1
!      
!      ! falls out of maximum iterations
!      if(num_iter_ik > max_iter) then
!        ! debug
!        write(stdout, *) 'max|diff(psi)|:', maxval( abs(y(1:npw,1:nbnd_occ(ik))) )
!        call errore('tddft_propagate::iterative_solver', 'Iterative solver &
!        does not converge after maximum number of steps',1)
!      endif
!    enddo
!    if(num_iter_ik > num_iter) num_iter = num_iter_ik
!
!  END SUBROUTINE iterative_solver
!  !---
