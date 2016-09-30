
!---
SUBROUTINE tddft_propagate()
  !---
  ! Real-time propagation of the electronic states
  USE kinds,                ONLY : dp
  USE wvfct,                ONLY : npwx
  USE tddft_module
  implicit none
  complex(dp) :: dt_cmplx
  integer :: istep, ik
  complex(dp), allocatable :: b(:,:)
  external propgator_Euler
  call start_clock ('tddft_propagate')

  ! setup solver and parameters
  select case(trim(solver))
    case('cgsolver')
      call cgsolver_initialize()
    case('itsolver')
      call itsolver_initialize(max_nbnd_occ)
    case default
      call errore('tddft_propagate', 'RT-tddft::tddft_progapate cannot recognize this solver in input', 1)
  end select
  allocate(b(npwx, max_nbnd_occ))

  ! main loop
  do istep = (init_step+1), (init_step+num_step)

    select case(trim(method))
      case('CN')
        call method_CN()
      case('CN-mid')
        call method_CN_mid()
      case('CN2')
        call method_CN2()
      case default
        call errore('tddft_propagate', 'RT-tddft::tddft_progapate cannot recognize this propagation method in input', 1)
    end select

  enddo

  ! finalize solver and parameters
  select case(trim(solver))
    case('cgsolver')
      call cgsolver_finalize()
    case('itsolver')
      call itsolver_finalize()
    case default
      call errore('tddft_propagate', 'RT-tddft::tddft_progapate cannot recognize this solver in input', 1)
  end select
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

  !---
  SUBROUTINE method_CN()
    !---
    ! CN method, o(dt^2) locally, unconditionally stable
    USE klist,                ONLY : nks
    USE io_files,             ONLY : iunigk
    USE wavefunctions_module, ONLY : evc
    USE wvfct,                ONLY : npw
    implicit none

    if(nks > 1) rewind(iunigk)
    do ik = 1, nks
 
      ! init potential inorder to call h_psi
      call init_k()

      ! prepare b
      ! b = [I - iH(t)0.5dt/hbar]*psi(t)
      call propgator_Euler(evc, b, 0.5d0*dt, nbnd_occ(ik))

      ! solve A * x = b
      tddft_psi(1:npw,1:nbnd_occ(ik)) = 2.d0*evc(1:npw,1:nbnd_occ(ik)) - &
        tddft_psi(1:npw,1:nbnd_occ(ik))
      select case(trim(solver))
        case('cgsolver')
          call cgsolver(propgator_Euler, b, tddft_psi, conv_thr, nbnd_occ(ik), &
            -0.5d0*dt, max_iter)
        case('itsolver')
          call itsolver(b, tddft_psi, conv_thr, nbnd_occ(ik), &
            0.5d0*dt, max_iter)
        case default
          call errore('tddft_propagate', 'RT-tddft::tddft_progapate cannot recognize this solver in input', 1)
      end select

      call save_k()
    enddo
    call tddft_update(istep, 1)
    call tddft_compute(istep)

  END SUBROUTINE method_CN
  !---

  !---
  SUBROUTINE method_CN_mid()
    !---
    ! mid-point CN method, o(dt^3) locally, unconditionally stable
    USE klist,                ONLY : nks
    USE io_files,             ONLY : iunigk
    USE wavefunctions_module, ONLY : evc
    USE wvfct,                ONLY : npw
    implicit none

    ! first half integrator
    if(nks > 1) rewind(iunigk)
    do ik = 1, nks
 
      ! init potential inorder to call h_psi
      call init_k()

      ! prepare b
      ! b = [I - iH(t)0.25dt/hbar]*psi(t)
      call propgator_Euler(evc, b, 0.25d0*dt, nbnd_occ(ik))

      ! guess psi(t+0.5dt)
      tddft_psi(1:npw,1:nbnd_occ(ik)) = 2.d0*evc(1:npw,1:nbnd_occ(ik)) - &
        tddft_psi(1:npw,1:nbnd_occ(ik))
      select case(trim(solver))
        case('cgsolver')
          call cgsolver(propgator_Euler, b, tddft_psi, conv_thr, nbnd_occ(ik), &
            -0.25d0*dt, max_iter)
        case('itsolver')
          call itsolver(b, tddft_psi, conv_thr, nbnd_occ(ik), &
            0.25d0*dt, max_iter)
        case default
          call errore('tddft_propagate', 'RT-tddft::tddft_progapate cannot recognize this solver in input', 1)
      end select

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

      ! compute psi(t+dt)
      tddft_psi(1:npw,1:nbnd_occ(ik)) = 2.d0*evc(1:npw,1:nbnd_occ(ik)) - &
        tddft_psi(1:npw,1:nbnd_occ(ik))
      select case(trim(solver))
        case('cgsolver')
          call cgsolver(propgator_Euler, b, tddft_psi, conv_thr, nbnd_occ(ik), &
            -0.5d0*dt, max_iter)
        case('itsolver')
          call itsolver(b, tddft_psi, conv_thr, nbnd_occ(ik), &
            0.5d0*dt, max_iter)
        case default
          call errore('tddft_propagate', 'RT-tddft::tddft_progapate cannot recognize this solver in input', 1)
      end select

      call save_k()
    enddo
    call tddft_update(istep, 3)
    call tddft_compute(istep)

  END SUBROUTINE method_CN_mid
  !---

  !---
  SUBROUTINE method_CN2()
    !---
    ! second order CN method, o(dt^3) locally
    USE klist,                ONLY : nks
    USE io_files,             ONLY : iunigk
    USE wavefunctions_module, ONLY : evc
    USE wvfct,                ONLY : npw
    implicit none

    if(nks > 1) rewind(iunigk)
    do ik = 1, nks
 
      ! init potential inorder to call h_psi
      call init_k()

      ! prepare b
      ! b = [I - iH(t)dt/hbar]*psi(t-dt)
      call propgator_Euler(tddft_psi, b, dt, nbnd_occ(ik))

      ! solve A * x = b
      tddft_psi(1:npw,1:nbnd_occ(ik)) = 2.d0*evc(1:npw,1:nbnd_occ(ik)) - &
        tddft_psi(1:npw,1:nbnd_occ(ik))
      select case(trim(solver))
        case('cgsolver')
          call cgsolver(propgator_Euler, b, tddft_psi, conv_thr, nbnd_occ(ik), &
            -dt, max_iter)
        case('itsolver')
          call itsolver(b, tddft_psi, conv_thr, nbnd_occ(ik), &
            dt, max_iter)
        case default
          call errore('tddft_propagate', 'RT-tddft::tddft_progapate cannot recognize this solver in input', 1)
      end select

      call save_k()
    enddo
    call tddft_update(istep, 1)
    call tddft_compute(istep)

  END SUBROUTINE method_CN2
  !---

END SUBROUTINE tddft_propagate
!---
