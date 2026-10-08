!Copyright (c) 2026, Anjia Ying and Lin Fu.
!Additional source for the IOIM-S inlet in a derivative of Xcompact3d
!(https://github.com/xcompact3d/Incompact3d).
!SPDX-License-Identifier: BSD 3-Clause

module spline_interp_mod_fast_3d
  implicit none
  private
  public :: spline_solver_type
  public :: init_spline_solver_fast
  public :: interp_mid_natural_spline_fast_3d

  type :: spline_solver_type
     integer :: width_half = 0
     integer :: skpt = 0
     integer :: width_all = 0
     real(8), allocatable :: bfac(:)   ! factorised main diagonal
     real(8), allocatable :: gamma(:)  ! forward-elimination coefficient
     real(8), allocatable :: alpha(:)  ! interpolation weight
     real(8), allocatable :: beta(:)   ! interpolation weight
  end type spline_solver_type

contains

  !===========================================================
  ! Precompute factors for a fixed window and skip.
  !===========================================================
  subroutine init_spline_solver_fast(width_all, skpt, solvern)
    implicit none
    integer, intent(in) :: width_all, skpt
    type(spline_solver_type), intent(out) :: solvern

    integer :: width_half, i, j
    real(8) :: xloc

    width_half = width_all/2

    solvern%width_half = width_half
    solvern%skpt       = skpt
    solvern%width_all  = width_all

    allocate(solvern%bfac(width_all))
    allocate(solvern%gamma(width_all))
    allocate(solvern%alpha(skpt))
    allocate(solvern%beta(skpt))

    solvern%bfac  = 0.0d0
    solvern%gamma = 0.0d0

    if (width_all > 2) then
       solvern%bfac(2) = 4.0d0
       do i = 3, width_all-1
          solvern%gamma(i) = 1.0d0 / solvern%bfac(i-1)
          solvern%bfac(i)  = 4.0d0 - solvern%gamma(i)
       end do
    end if

    do j = 1, skpt
       xloc = dble(j-1) / dble(skpt)
       solvern%alpha(j) = 1.0d0 - xloc
       solvern%beta(j)  = xloc
    end do

  end subroutine init_spline_solver_fast


  !===========================================================
  ! Batched natural cubic spline over the logical time window.
  !
  ! u_sq(nyn,nzn,3,width_all): physical storage
  ! sqN(a)=b: logical sample a is stored in physical slot b
  ! y_new: samples interpolated on [width_half, width_half+1)
  !===========================================================
  subroutine interp_mid_natural_spline_fast_3d(u_sq, nyn, nzn, sqN, solvern, y_new)
    implicit none
    integer, intent(in) :: nyn, nzn
    type(spline_solver_type), intent(in) :: solvern
    real(8), intent(in)  :: u_sq(nyn, nzn, 3, solvern%width_all)
    integer, intent(in)  :: sqN(solvern%width_all)
    real(8), intent(out) :: y_new(nyn, nzn, 3, solvern%skpt)

    integer :: width_all, width_half, skpt
    integer :: iy, iz, ic, i, j
    real(8) :: rhs(solvern%width_all), M(solvern%width_all)
    real(8) :: tmp
    real(8) :: uL, uR

    width_all  = solvern%width_all
    width_half = solvern%width_half
    skpt       = solvern%skpt

    !---------------------------------------------------------
    ! Check that sqN maps each logical index onto a physical slot.
    ! logical index a -> physical slot b
    !---------------------------------------------------------
    do i = 1, width_all
       if (sqN(i) < 1 .or. sqN(i) > width_all) then
          error stop "interp_mid_natural_spline_fast_3d: sqN out of range"
       end if
    end do

    !---------------------------------------------------------
    ! Spline each (iy,iz,ic) in logical order 1..width_all.
    !---------------------------------------------------------
    do iz = 1, nzn
       do iy = 1, nyn
          do ic = 1, 3

             M   = 0.0d0
             rhs = 0.0d0

             if (width_all > 2) then
                do i = 2, width_all-1
                   rhs(i) = 6.0d0 * ( u_sq(iy,iz,ic,sqN(i+1)) &
                                - 2.0d0*u_sq(iy,iz,ic,sqN(i))   &
                                +        u_sq(iy,iz,ic,sqN(i-1)) )
                end do

                do i = 3, width_all-1
                   rhs(i) = rhs(i) - solvern%gamma(i) * rhs(i-1)
                end do

                M(width_all-1) = rhs(width_all-1) / solvern%bfac(width_all-1)
                do i = width_all-2, 2, -1
                   M(i) = (rhs(i) - M(i+1)) / solvern%bfac(i)
                end do
             end if

             uL = u_sq(iy,iz,ic,sqN(width_half))
             uR = u_sq(iy,iz,ic,sqN(width_half+1))

             do j = 1, skpt
                tmp = ((solvern%alpha(j)**3 - solvern%alpha(j)) * M(width_half)   &
                     +  (solvern%beta(j)**3  - solvern%beta(j))  * M(width_half+1)) / 6.0d0

                y_new(iy,iz,ic,j) = solvern%alpha(j) * uL &
                                   + solvern%beta(j)  * uR &
                                   + tmp
             end do

          end do
       end do
    end do

  end subroutine interp_mid_natural_spline_fast_3d

end module spline_interp_mod_fast_3d