!Copyright (c) 2026, Anjia Ying and Lin Fu.
!Additional source for the IOIM-S inlet in a derivative of Xcompact3d
!(https://github.com/xcompact3d/Incompact3d).
!SPDX-License-Identifier: BSD 3-Clause

module ioim_fftw_plans
  use, intrinsic :: iso_c_binding
  implicit none

  logical :: ioim_plans_initialized = .false.

  type(c_ptr) :: plan_z_fwd1 = c_null_ptr
  type(c_ptr) :: plan_z_fwd2 = c_null_ptr
  type(c_ptr) :: plan_z_fwd3 = c_null_ptr

  type(c_ptr) :: plan_z_bwd1 = c_null_ptr
  type(c_ptr) :: plan_z_bwd2 = c_null_ptr
  type(c_ptr) :: plan_z_bwd3 = c_null_ptr

  type(c_ptr) :: plan_t_fwd  = c_null_ptr

contains

  subroutine ioim_create_fftw_plans(nzloc, ntloc, &
                                    input_Z_realm1, input_Z_realm2, input_Z_realm3, &
                                    output_Z_cplxm1, output_Z_cplxm2, output_Z_cplxm3, &
                                    input_X_cplx, output_X_cplx)
    include 'fftw3.f03'

    integer, intent(in) :: nzloc, ntloc
    real(8), intent(inout) :: input_Z_realm1(nzloc), input_Z_realm2(nzloc), input_Z_realm3(nzloc)
    complex(8), intent(inout) :: output_Z_cplxm1(nzloc/2+1), output_Z_cplxm2(nzloc/2+1), output_Z_cplxm3(nzloc/2+1)
    complex(8), intent(inout) :: input_X_cplx(ntloc), output_X_cplx(ntloc)

    if (.not. ioim_plans_initialized) then
      plan_z_fwd1 = fftw_plan_dft_r2c_1d(nzloc, input_Z_realm1, output_Z_cplxm1, FFTW_MEASURE)
      plan_z_fwd2 = fftw_plan_dft_r2c_1d(nzloc, input_Z_realm2, output_Z_cplxm2, FFTW_MEASURE)
      plan_z_fwd3 = fftw_plan_dft_r2c_1d(nzloc, input_Z_realm3, output_Z_cplxm3, FFTW_MEASURE)

      plan_z_bwd1 = fftw_plan_dft_c2r_1d(nzloc, output_Z_cplxm1, input_Z_realm1, FFTW_MEASURE)
      plan_z_bwd2 = fftw_plan_dft_c2r_1d(nzloc, output_Z_cplxm2, input_Z_realm2, FFTW_MEASURE)
      plan_z_bwd3 = fftw_plan_dft_c2r_1d(nzloc, output_Z_cplxm3, input_Z_realm3, FFTW_MEASURE)

      plan_t_fwd  = fftw_plan_dft_1d(ntloc, input_X_cplx, output_X_cplx, FFTW_FORWARD, FFTW_MEASURE)

      ioim_plans_initialized = .true.
    endif
  end subroutine ioim_create_fftw_plans


  subroutine ioim_destroy_fftw_plans()
    include 'fftw3.f03'

    if (ioim_plans_initialized) then
      if (c_associated(plan_z_fwd1)) call fftw_destroy_plan(plan_z_fwd1)
      if (c_associated(plan_z_fwd2)) call fftw_destroy_plan(plan_z_fwd2)
      if (c_associated(plan_z_fwd3)) call fftw_destroy_plan(plan_z_fwd3)

      if (c_associated(plan_z_bwd1)) call fftw_destroy_plan(plan_z_bwd1)
      if (c_associated(plan_z_bwd2)) call fftw_destroy_plan(plan_z_bwd2)
      if (c_associated(plan_z_bwd3)) call fftw_destroy_plan(plan_z_bwd3)

      if (c_associated(plan_t_fwd))  call fftw_destroy_plan(plan_t_fwd)

      call fftw_cleanup()
      ioim_plans_initialized = .false.
    endif
  end subroutine ioim_destroy_fftw_plans

end module ioim_fftw_plans