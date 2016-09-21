
!---
SUBROUTINE tddft_propagate()
  !---
  ! Real-time propagation of the electronic states
  USE kinds,        ONLY : dp
  USE wvfct,        ONLY : npwx, nbnd
  USE io_files,     ONLY : iunigk
  USE klist,        ONLY : nks
  USE tddft_module
  implicit none
  complex(dp) :: e_coef
  integer :: istep, ik
  complex(dp), allocatable :: b(:,:)
  complex(dp), allocatable :: tddft_hpsi(:,:)
  complex(dp), allocatable :: y(:,:)
  call start_clock ('tddft_propagate')

  ! setup parameters
  e_coef = (0.0_dp, 1.0_dp) * dt
  allocate(b(npwx, max_nbnd_occ))
  b = (0.d0, 0.d0)
  allocate(tddft_hpsi(npwx, max_nbnd_occ))
  tddft_hpsi = (0.d0, 0.d0)
  allocate(y(npwx, nbnd))
  y = (0.d0, 0.d0)

  ! main loop
  do istep = init_step, (init_step+num_step-1)
    num_iter = 0
    ! loop over k-points
    if(nks > 1) rewind(iunigk)
    do ik = 1, nks
      call init_k()
      call iterative_solver()
      call save_k()
    enddo

    ! update the hamiltonian (recompute charge and potential)
    call tddft_update(istep)
    call tddft_compute(istep)
  enddo

  ! deallocate  parameters
  deallocate(b, tddft_hpsi, y)

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
  SUBROUTINE iterative_solver()
    !---
    ! Solve psi(t+dt) from H(t) and b
    ! [I + iH(t)dt/hbar]*psi(t+dt) = b
    USE wvfct,                  ONLY : npwx, npw
    USE wavefunctions_module,   ONLY : evc
    USE io_global,              ONLY : stdout
    implicit none
    integer :: num_iter_ik

    ! prepare b
    ! b = [I - iH(t)dt/hbar]*psi(t-dt)
    call h_psi(npwx, npw, nbnd_occ(ik), tddft_psi, tddft_hpsi)
    b(1:npw, 1:nbnd_occ(ik)) = tddft_psi(1:npw, 1:nbnd_occ(ik)) - &
        e_coef * tddft_hpsi(1:npw, 1:nbnd_occ(ik))

    ! guess wavefunction psi(t+dt)
    call h_psi(npwx, npw, nbnd_occ(ik), evc, tddft_hpsi)
    tddft_psi(1:npw, 1:nbnd_occ(ik)) = b(1:npw, 1:nbnd_occ(ik)) - &
        e_coef * tddft_hpsi(1:npw, 1:nbnd_occ(ik))

    ! compute increment y
    y(1:npw, 1:nbnd_occ(ik)) = tddft_psi(1:npw, 1:nbnd_occ(ik)) - &
        evc(1:npw, 1:nbnd_occ(ik))

    ! start iterations
    num_iter_ik = 0
    do while(maxval(abs(y(1:npw, 1:nbnd_occ(ik)))) > conv_threshold)
      ! compute diff(psi)
      call h_psi(npwx, npw, nbnd_occ(ik), y, tddft_hpsi)
      y(1:npw, 1:nbnd_occ(ik)) = -e_coef * tddft_hpsi(1:npw, 1:nbnd_occ(ik))
      
      ! update tddft_psi
      tddft_psi(1:npw, 1:nbnd_occ(ik)) = tddft_psi(1:npw, 1:nbnd_occ(ik)) + &
        y(1:npw, 1:nbnd_occ(ik))
      num_iter_ik = num_iter_ik + 1
      
      ! falls out of maximum iterations
      if(num_iter_ik > max_iter) then
        ! debug
        write(stdout, *) 'max|diff(psi)|:', maxval( abs(y(1:npw,1:nbnd_occ(ik))) )
        call errore('tddft_propagate::iterative_solver', 'Iterative solver &
        does not converge after maximum number of steps',1)
      endif
    enddo
    if(num_iter_ik > num_iter) num_iter = num_iter_ik

  END SUBROUTINE iterative_solver
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
