!Copyright (c) 2026, Anjia Ying and Lin Fu.
!Additional source for the IOIM-S inlet in a derivative of Xcompact3d
!(https://github.com/xcompact3d/Incompact3d).
!SPDX-License-Identifier: BSD 3-Clause

module ioim_checkpoint_binary
  implicit none
contains

  subroutine ioim_checkpoint_save_binary(filename, icyc10, UVW_TYZ_lag, syy, cyt, UP2, &
                                         N_I,N_Gamma, ny, nz, nv, width_interp)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: icyc10, N_I,N_Gamma, ny, nz, nv, width_interp
    complex(8), intent(in) :: UVW_TYZ_lag(N_Gamma,ny,nz/2,nv)
    complex(8), intent(in) :: syy(N_I,nz/2,nv)
    complex(8), intent(in) :: cyt(N_I,ny,nz/2,nv)
    real(8),    intent(in) :: UP2(ny,nz,nv,width_interp)

    integer :: iu
    integer :: version

    version = 1

    open(newunit=iu, file=filename, form='unformatted', access='stream', &
         status='replace', action='write')

    write(iu) version
    write(iu) icyc10
    write(iu) N_I,N_Gamma, ny, nz, nv, width_interp

    write(iu) UVW_TYZ_lag
    write(iu) syy
    write(iu) cyt
    write(iu) UP2

    close(iu)
  end subroutine ioim_checkpoint_save_binary


  subroutine ioim_checkpoint_load_binary(filename, file_exists, file_match, icyc10_rec, &
                                         UVW_TYZ_lag, syy, cyt, UP2, &
                                         N_I,N_Gamma, ny, nz, nv, width_interp)
    implicit none
    character(len=*), intent(in) :: filename
    logical, intent(out) :: file_exists, file_match
    integer, intent(out) :: icyc10_rec
    integer, intent(in) :: N_I,N_Gamma, ny, nz, nv, width_interp
    complex(8), intent(out) :: UVW_TYZ_lag(N_Gamma,ny,nz/2,nv)
    complex(8), intent(out) :: syy(N_I,nz/2,nv)
    complex(8), intent(out) :: cyt(N_I,ny,nz/2,nv)
    real(8),    intent(out) :: UP2(ny,nz,nv,width_interp)
    
    real(8), allocatable, dimension(:,:,:,:) :: UP2_read

    integer :: iu
    integer :: version
    integer :: N_I_f, N_Gamma_f, ny_f, nz_f, nv_f, width_interp_f
    integer :: iiw,iiwf

    inquire(file=filename, exist=file_exists)

    file_match = .false.
    icyc10_rec = 0

    if (.not. file_exists) return

    open(newunit=iu, file=filename, form='unformatted', access='stream', &
         status='old', action='read')

    read(iu) version
    read(iu) icyc10_rec
    read(iu) N_I_f,N_Gamma_f, ny_f, nz_f, nv_f, width_interp_f
    
    if (N_I_f /= N_I .or. N_Gamma_f/=N_Gamma .or. ny_f /= ny .or. nz_f /= nz .or. &
        nv_f /= nv) then
      close(iu)
      return
    endif
        
    if ( (width_interp_f .ne. width_interp) ) then
        write(*,*) '>>> *** Width of the sliding cubic-spline window not match! I would periodically assign its value. *** <<<'
    endif

    allocate(UP2_read( ny , nz , nv , width_interp_f ))

    read(iu) UVW_TYZ_lag
    read(iu) syy
    read(iu) cyt
    read(iu) UP2_read
    
    ! Fill UP2 element by element from the checkpoint buffer.
    do iiw = 1,width_interp
        
        iiwf = mod((iiw-1),width_interp_f)+1
        
        UP2(:,:,:,iiw) = UP2_read(:,:,:,iiwf)
        
    enddo
    
    deallocate(UP2_read)
    
    close(iu)
    file_match = .true.

  end subroutine ioim_checkpoint_load_binary

end module ioim_checkpoint_binary