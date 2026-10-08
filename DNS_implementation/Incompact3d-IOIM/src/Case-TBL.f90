!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause
!
!Modifications for the IOIM-S, LWS-S, and LWS-R turbulent boundary-layer inlet.
!Copyright (c) 2026, Anjia Ying and Lin Fu.
!Distributed under the same BSD 3-Clause License. The original copyright
!notice above is retained.

module tbl

  use decomp_2d
  use decomp_2d_constants
  use decomp_2d_mpi
  use variables
  use param
  use spline_interp_mod_fast_3d
  use ioim_fftw_plans
  use ioim_checkpoint_binary

  implicit none

  integer :: fs
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_tbl, boundary_conditions_tbl, postprocess_tbl, visu_tbl, visu_tbl_init, &
            momentum_forcing_sponge

contains

  subroutine init_tbl (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d_io
    use param , only : zptwofive
    use var, only : Umx_XY, Vmx_XY, xp, yp, zp
    use MPI


    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    integer :: k,j,i,ierror,ii,is,it,code
    integer :: ig,jg,kg,l,ierr
    integer, dimension (:), allocatable :: seed
    real(mytype) :: randu,randv,randw
    real(mytype) :: xloc,yloc,zloc,delta_i,delta_x,utau_i,nu
    real(mytype) :: ys,ufy,vfy,arg_sin,arg_cos,dup_dx,dvfy_dy,dvp_dy
    real(mytype) :: up,vp,wp
    real(mytype), dimension(4) :: yscl,uamp,vamp,freq,uconv,zlam,phase

    if (iscalar==1) then

       phi1(:,:,:,:) = zptwofive !change as much as you want
          if ((nclyS1==2).and.(xstart(2)==1)) then
             !! Generate a hot patch on bottom boundary
             phi1(:,1,:,:) = one
          endif
          if ((nclySn==2).and.(xend(2)==ny)) THEN
             phi1(:,xsize(2),:,:) = zptwofive
          endif

    endif
    ux1=zero;uy1=zero;uz1=zero

    if (iin.eq.4) then
       ! Cold start for recycling inflow. Mean field is the MATLAB profile.
       ! Fluctuations follow STREAmS add_synthetic_perturbations: four outer
       ! Fourier modes, with w recovered from continuity, plus a weak random
       ! seed inside the local boundary layer. Mode 1 in STREAmS is prepared
       ! but not applied; the same is done here.
       nu = one/re
       delta_i = yp(ny)
       do j = 2, ny
          if (Umx_XY(1,j).ge.0.99_mytype) then
             delta_i = yp(j)
             exit
          endif
       enddo
       utau_i = sqrt(max(nu*Umx_XY(1,2)/max(yp(2),1.0e-12_mytype),1.0e-12_mytype))

       yscl  = (/0.30_mytype, 0.45_mytype, 0.60_mytype, 0.50_mytype/)*delta_i
       uamp  = (/0.30_mytype, 0.20_mytype, 0.08_mytype, 0.04_mytype/)
       vamp  = (/-0.06_mytype, -0.05_mytype, -0.04_mytype, -0.03_mytype/)
       freq  = (/1.20_mytype, 0.60_mytype, 0.40_mytype, 0.20_mytype/)/delta_i
       uconv = (/0.90_mytype, 0.90_mytype, 0.90_mytype, 0.90_mytype/)
       zlam  = (/0.333_mytype, 0.25_mytype, 0.20_mytype, 0.166_mytype/)*zlz
       phase = zero
       if (nrank.eq.0) then
          call random_seed(size=ii)
          allocate(seed(ii))
          do i = 1, ii
             seed(i) = 63946*(i-1)+17
          enddo
          call random_seed(put=seed)
          deallocate(seed)
          do l = 1, 4
             call random_number(phase(l))
          enddo
          phase = two*pi*phase
       endif
       call MPI_BCAST(phase, 4, real_type, 0, MPI_COMM_WORLD, ierr)

       do k=1,xsize(3)
          kg = k+xstart(3)-1
          zloc = zp(kg)
          do j=1,xsize(2)
             jg = j+xstart(2)-1
             yloc = yp(jg)
             do i=1,xsize(1)
                ig = i+xstart(1)-1
                xloc = xp(ig)
                ux1(i,j,k) = Umx_XY(ig,jg)
                uy1(i,j,k) = Vmx_XY(ig,jg)
                uz1(i,j,k) = zero

                delta_x = delta_i
                do ii = 2, ny
                   if (Umx_XY(ig,ii).ge.0.99_mytype) then
                      delta_x = yp(ii)
                      exit
                   endif
                enddo

                up = zero
                vp = zero
                wp = zero
                do l = 1, 4
                   arg_sin = freq(l)*xloc/uconv(l)
                   arg_cos = two*pi*zloc/zlam(l)+phase(l)
                   ys = yloc/yscl(l)
                   ufy = ys*exp(-ys)
                   vfy = ys*ys*exp(-(ys*ys))
                   up = up+uamp(l)*ufy*sin(arg_sin)*cos(arg_cos)
                   vp = vp+vamp(l)*vfy*sin(arg_sin)*cos(arg_cos)
                   dup_dx = uamp(l)*ufy*freq(l)/uconv(l)*cos(arg_sin)
                   dvfy_dy = two*ys/yscl(l)*exp(-(ys*ys))*(one-ys*ys)
                   dvp_dy = vamp(l)*dvfy_dy*sin(arg_sin)
                   wp = wp-(dup_dx+dvp_dy)*sin(arg_cos)*zlam(l)/(two*pi)
                enddo
                if ((yloc.lt.delta_x).and.(init_noise.ne.zero)) then
                   call random_number(randu)
                   call random_number(randv)
                   call random_number(randw)
                   up = up+0.03_mytype*(randu-half)*(yloc/delta_i)
                   vp = vp+0.03_mytype*(randv-half)*(yloc/delta_i)
                   wp = wp+0.03_mytype*(randw-half)*(yloc/delta_i)
                endif
                ux1(i,j,k) = ux1(i,j,k)+up
                uy1(i,j,k) = uy1(i,j,k)+vp
                uz1(i,j,k) = uz1(i,j,k)+wp
             enddo
          enddo
       enddo

       if (xstart(2).eq.1) then
          ux1(:,1,:) = zero
          uy1(:,1,:) = zero
          uz1(:,1,:) = zero
       endif

       if (nrank.eq.0) write(*,*) '# Start with predicted mean + synthetic-mode fluctuation as initial field, iin=4, delta_in=', delta_i, ' u_tau=', utau_i
    else
       ! Original TBL cold start: Blasius inlet duplicated over the domain.
       call blasius()
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
                uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
                uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
             enddo
          enddo
       enddo
    endif

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_tbl
  !********************************************************************
  subroutine boundary_conditions_tbl (ux,uy,uz,phi)

    use navier, only : tbl_flrt
    use param , only : zero, zptwofive, irestart
    USE var, only: yp
    ! USE variables, only : bxx1,bxy1,bxz1
    USE var, only: ta1,tb1,tc1,ta2,tb2,tc2,ta3,tb3,tc3
    USE var, only: td1,te1,tf1,td2,te2,tf2,td3,te3,tf3
    USE var, only: tg1,tg2,tg3
    USE var, only: Um_inlet,Vm_inlet,Um_recyc,Vm_recyc
    USE var, only: length_xR
    USE var, only: eps_syy
    USE var, only: indU_p1,indU_p2
    USE var, only: UVW_record_TYZ
    USE var, only: UVW_record_TYZ_lagAll,UVW_record_TYZ_lagEnd
    USE var, only: syy_cpu,cyt_cpu,syy1_cpu,cyt1_cpu
    USE var, only: tyt_cpu
    USE var, only: wcalcR_cpu,up_cpu
    USE var, only: wcalc_i,wcalc_o
    USE var, only: up_cpuZ
    USE var, only: wcalc_phy,up_phy
    USE var, only: u_sq,u_interp
    USE var, only: sq_ft
    USE var, only: nv
    USE var, only: dt
    USE var, only: dp,cdp
    USE var, only: niu,threshold_bt,U_infty,u_tau_R,u_tau_I,l_tau_R,l_tau_I,Re_tau_R,Re_tau_I
    USE var, only: gama,gama8,gamaU
    USE var, only: jref_I,jref_R,jref
    USE var, only: ind_Inner1,ind_Inner2,ind_Outer1,ind_Outer2
    USE var, only: wt_Inner1, wt_Inner2, wt_Outer1, wt_Outer2, W_yita, W_o, W_i
    USE var, only: icyc10_rec
    USE var, only: bdamp,indrec
    USE var, only: UmT_recyc_LWS , VmT_recyc_LWS
    
    USE var, only: i_inlet,xp,m_proc,n_mrank
    
    USE var, only: ioim_local_init_done,solvern
    
    USE var, only: comm_inlet, color_io, key_io, rank_inlet
    
    USE param, only: re
    
    use variables , only : method_rr,delta_in,x_recyc
    use variables , only : N_I,n_sp,width_interp,N_Upsilon,N_Gamma,i_zshift
    
    use MPI
    
    use, intrinsic :: iso_c_binding

    implicit none
    
    include 'fftw3.f03'

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxr,uyr,uzr
    real(mytype),dimension(zsize(2)) :: u_m1,v_m1,w_m1
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: u_shift,v_shift,w_shift
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    
    real(mytype),dimension(ysize(3)) :: Ut_I_inner_Z,Vt_I_inner_Z,Wt_I_inner_Z
    real(mytype),dimension(ysize(3)) :: Ut_I_outer_Z,Vt_I_outer_Z,Wt_I_outer_Z
    real(mytype),dimension(ysize(2),ysize(3)) :: Ut_I_YZ,Vt_I_YZ,Wt_I_YZ
    
    real(mytype),dimension(ysize(3)) :: Un_I_Z,Vn_I_Z,Wn_I_Z
    real(mytype),dimension(ysize(2),ysize(3)) :: Un_I_YZ,Vn_I_YZ,Wn_I_YZ
    
    real(mytype) :: Ut_I_inner,Vt_I_inner,Wt_I_inner
    real(mytype) :: Ut_I_outer,Vt_I_outer,Wt_I_outer
    real(mytype),dimension(ysize(2)) :: Ut_I,Vt_I,Wt_I
    
    real(mytype),dimension(ysize(2)) :: Um_recyc_LWS,Vm_recyc_LWS
    real(mytype),dimension(ysize(2),ysize(3)) :: Um_inlet_LWS,Vm_inlet_LWS
    real(mytype),dimension(ysize(3)) :: Um_inner_Z_LWS,Vm_inner_Z_LWS,Um_outer_Z_LWS,Vm_outer_Z_LWS
    real(mytype),dimension(ysize(3)) :: Um_outer_Z_LWS_end
    real(mytype) :: T_LWS
    
    real(mytype) :: x, y, z, um, wt_1, wt_2, wt_3, wt_4, flux_scaling
    real(mytype) :: yita
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,yI_inner,yI_outer,yR_inner,yR_outer

    integer :: i, j,jj, k, is,inda,indb,indc,indd,n_shift,k_shift
    integer :: code,ierror
    
    real(mytype) :: alp, bet, delta_I_099,delta_R_099
    
    integer :: icyc1,icyc00,icyc01,icyc10,icyc11,icyc1_now,icyc10_now
    integer :: iitr,iitp
    integer :: iit,it1,it2
    integer :: iiy,iiz,iiv,iz1
    integer :: izg,iiz_max,iiy_max
    integer :: iermpi
    integer :: ii_p,infoa
    integer :: iunit,m,m1,m2
    integer :: n_ioim_cal,t_ioim_cal
    real(mytype) :: wt1,wt2
    character(len=100) :: filename
    
    integer:: i_inlet_recv,ystart3_recv,yend3_recv
    
    real(mytype) :: alpha_t
    
    real(mytype) :: weight_1,weight_2
    
    !!! FFT
    real(8), dimension(zsize(3)) :: input_Z_realm1,input_Z_realm2,input_Z_realm3
    complex(cdp), dimension(zsize(3)/2+1) :: output_Z_cplxm1,output_Z_cplxm2,output_Z_cplxm3
    
    complex(cdp), dimension(N_I) :: input_X_cplx
    complex(cdp), dimension(N_I) :: output_X_cplx
    real(mytype) :: scalez
    
    complex(cdp) :: zimg,tmp_ioim
    
    integer :: imm,iuu
    integer :: i_rank
    complex(cdp), dimension(N_I) :: expz
    
    
    complex(cdp), allocatable, dimension(:,:,:) :: syy,syy_half
    complex(cdp), allocatable, dimension(:,:,:,:) :: cyt, UVW_TYZ_lag,cyt_half, UVW_TYZ_lag_half
    real(dp), allocatable, dimension(:,:,:,:) :: UP2
    
    real(kind=dp), allocatable :: buff_r(:)
    complex(kind=cdp), allocatable :: buff_c(:)
    
    integer :: nnz
    integer :: statusx,ierr
    integer :: ios,nlines,ireadC
    logical :: file_exists,file_match
    character(len=1000) :: testlengthx
    
    integer :: nyn, nzn, ind_store
    
    integer,dimension(N_I) :: sq_start
    integer,dimension(N_I) :: sq_end
    
    integer,dimension(width_interp) :: sqN
    
    real(dp)    :: deni_loc(ysize(2),nv), deno_loc(ysize(2),nv)
    real(dp)    :: deni_glb(ysize(2),nv), deno_glb(ysize(2),nv)
    real(dp)    :: gamma_io(ysize(2),nv)
    complex(cdp) :: num_loc(ysize(2),nv), num_glb(ysize(2),nv)
    real(dp)    :: eps_gamma
    
    nyn=ysize(2)
    nzn=ysize(3)

    ! Recycling inflow is the only path that uses the IOIM pool.
    if (iin.eq.4) then

    if (.not. ioim_local_init_done) then
        call init_spline_solver_fast(width_interp, n_sp, solvern)

        call ioim_create_fftw_plans(zsize(3), N_I, &
             input_Z_realm1, input_Z_realm2, input_Z_realm3, &
             output_Z_cplxm1, output_Z_cplxm2, output_Z_cplxm3, &
             input_X_cplx, output_X_cplx)

        ioim_local_init_done = .true.
    endif
    
    file_exists = .false.
    file_match  = .false.
    
    zimg = cmplx(0.0_cdp, 1.0_cdp, kind=cdp)
    
    !!! Recyling-rescaling methods for inflow turbulence
    if (iin.eq.4) then
       
       ta1=zero
       tb1=zero
       tc1=zero
       ta2=zero
       tb2=zero
       tc2=zero
       ta3=zero
       tb3=zero
       tc3=zero
       
       td2=zero
       te2=zero
       tf2=zero
       td3=zero
       te3=zero
       tf3=zero
       
       n_shift = zsize(3)/2
       
       uxr(1,:,:) = ux(indrec,:,:)
       uyr(1,:,:) = uy(indrec,:,:)
       uzr(1,:,:) = uz(indrec,:,:)

       call transpose_x_to_y(uxr,ta2)
       call transpose_x_to_y(uyr,tb2)
       call transpose_x_to_y(uzr,tc2)
       
       ! Recycling-plane extraction. IOIM-S also stores the spanwise spectrum.
           
           call transpose_y_to_z(ta2,ta3)
           call transpose_y_to_z(tb2,tb3)
           call transpose_y_to_z(tc2,tc3)
           
           if (method_rr .eq. 1) then
               
               do imm = 1,zsize(2)

                  input_Z_realm1 = ta3(1,imm,1:zsize(3))
                  input_Z_realm2 = tb3(1,imm,1:zsize(3))
                  input_Z_realm3 = tc3(1,imm,1:zsize(3))

                  call fftw_execute_dft_r2c(plan_z_fwd1, input_Z_realm1, output_Z_cplxm1)
                  call fftw_execute_dft_r2c(plan_z_fwd2, input_Z_realm2, output_Z_cplxm2)
                  call fftw_execute_dft_r2c(plan_z_fwd3, input_Z_realm3, output_Z_cplxm3)

                  ta3(2,imm,1:zsize(3):2) = real(output_Z_cplxm1(2:zsize(3)/2+1))
                  tb3(2,imm,1:zsize(3):2) = real(output_Z_cplxm2(2:zsize(3)/2+1))
                  tc3(2,imm,1:zsize(3):2) = real(output_Z_cplxm3(2:zsize(3)/2+1))

                  ta3(3,imm,1:zsize(3):2) = aimag(output_Z_cplxm1(2:zsize(3)/2+1))
                  tb3(3,imm,1:zsize(3):2) = aimag(output_Z_cplxm2(2:zsize(3)/2+1))
                  tc3(3,imm,1:zsize(3):2) = aimag(output_Z_cplxm3(2:zsize(3)/2+1))

               enddo
           
               !-------------------------------------------------------
           endif
       
           ! Obtain the mean velocity profile
           if(zstart(1)==1) then
              do j=1,zsize(2)
                 u_m1(j) = sum(ta3(1,j,:))/real(zsize(3))
                 v_m1(j) = sum(tb3(1,j,:))/real(zsize(3))
                 w_m1(j) = sum(tc3(1,j,:))/real(zsize(3))
              enddo
              do k=1,zsize(3)
                 td3(1,:,k) = u_m1
                 te3(1,:,k) = v_m1
                 tf3(1,:,k) = w_m1
              enddo
           endif
           call transpose_z_to_y(td3,td2)
           call transpose_z_to_y(te3,te2)
           call transpose_z_to_y(tf3,tf2)
       
           call transpose_z_to_y(ta3,ta2)
           call transpose_z_to_y(tb3,tb2)
           call transpose_z_to_y(tc3,tc2)
       
          ta2(1,:,:) = ta2(1,:,:) - td2(1,:,:)
          tb2(1,:,:) = tb2(1,:,:) - te2(1,:,:)
          tc2(1,:,:) = tc2(1,:,:) - tf2(1,:,:)
          
          Um_recyc_LWS = td2(1,:,1)
          Vm_recyc_LWS = te2(1,:,1)

       ! method_rr = 1: IOIM-S; 2: LWS-S; 3: LWS-R
       if (method_rr .eq. 1) then
       ! IOIM-S: time-separated pool, footprint removal and inner-outer blending.
       
       if (itime .eq. 1) then
           
           i_inlet = 0
           if (ystart(1) .eq. 1) then
               
                i_inlet = 1
        
                allocate(UVW_record_TYZ( N_I , ysize(2) , ysize(3) , nv))
                UVW_record_TYZ=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(UVW_record_TYZ_lagAll( N_Gamma , ysize(2) , ysize(3) , nv))
                UVW_record_TYZ_lagAll=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(UVW_record_TYZ_lagEnd( ysize(2) , ysize(3) , nv))
                UVW_record_TYZ_lagEnd=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                eps_syy = 1e-12
        
                !!!
                allocate(syy_cpu( N_I , ysize(3) , nv ))
                syy_cpu=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(syy1_cpu( N_I , nv ))
                syy1_cpu=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(cyt_cpu( N_I , ny , ysize(3) , nv ))
                cyt_cpu=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(cyt1_cpu( N_I , ny , nv ))
                cyt1_cpu=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(tyt_cpu( N_I , nv ))
                tyt_cpu=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                !!!
                allocate(wcalc_i( N_I , ny , ysize(3) , nv))
                wcalc_i=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(wcalc_o( N_I , ny , ysize(3) , nv))
                wcalc_o=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(wcalcR_cpu( N_I , ysize(3) , nv))
                wcalcR_cpu=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                allocate(up_cpu( N_I , ny , ysize(3) , nv ))
                up_cpu=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                !!!
                allocate(up_cpuZ( ny , ysize(3) , nv ))
                up_cpuZ=cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
        
                !!!
                allocate(wcalc_phy( ny , ysize(3) , nv ))
                wcalc_phy=zero
        
                allocate(up_phy( ny , ysize(3) , nv ))
                up_phy=zero
        
                allocate(u_sq( ysize(2),ysize(3),nv,width_interp ))
                u_sq=zero
        
                allocate(u_interp( ysize(2),ysize(3),nv,n_sp ))
                u_interp=zero
            
                allocate(sq_ft(N_I))
                sq_ft = zero
                do i = 1, N_I
                    sq_ft(i) = 2.0*pi*real(i-1)/real(N_I)
                end do
        
           endif
           
           call MPI_Barrier(MPI_COMM_WORLD, ierr)
           
           if (nrank .ne. 0) then
                
                call MPI_Send(i_inlet, 1,  MPI_INTEGER, 0, &
                & 101, MPI_COMM_WORLD, iermpi)
                call MPI_Send(ystart(3), 1,  MPI_INTEGER, 0, &
                & 102, MPI_COMM_WORLD, iermpi)
                call MPI_Send(yend(3), 1,  MPI_INTEGER, 0, &
                & 103, MPI_COMM_WORLD, iermpi)
                
           endif
           
           if (nrank .eq. 0) then
               
                n_mrank = 0
               
                do i = 0,nproc-1
                    
                    if (i .ne. 0) then
                        
                        call MPI_Recv(i_inlet_recv, 1,  MPI_INTEGER, &
                        & i, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                        call MPI_Recv(ystart3_recv, 1,  MPI_INTEGER, &
                        & i, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                        call MPI_Recv(yend3_recv, 1,  MPI_INTEGER, &
                        & i, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                        
                    else
                        i_inlet_recv = i_inlet
                        ystart3_recv = ystart(3)
                        yend3_recv   = yend(3)
                    endif
                    
                    if (i_inlet_recv .eq. 1) then
                        n_mrank           = n_mrank + 1
                        m_proc(n_mrank,1) = i
                        m_proc(n_mrank,2) = ystart3_recv
                        m_proc(n_mrank,3) = yend3_recv
                    endif
                    
                enddo
               
           endif
           
           if (i_inlet .eq. 1) then
                color_io = 1
                key_io   = nrank
           else
                color_io = MPI_UNDEFINED
                key_io   = nrank
           endif

           call MPI_COMM_SPLIT(MPI_COMM_WORLD, color_io, key_io, comm_inlet, ierr)
           
           if (ystart(1) .eq. 1) then
               niu = 1.0/re
               threshold_bt = 0.99   ! threshlod to determine the boundary thickness
               U_infty = 1.0
           
               u_tau_R = sqrt(niu*Um_recyc(2)/yp(2))
               u_tau_I = sqrt(niu*Um_inlet(2)/yp(2))
           
               l_tau_R = niu/u_tau_R
               l_tau_I = niu/u_tau_I
           
               do j=1,ysize(2)
                   if(Um_inlet(j).ge.threshold_bt*U_infty) then
                      jref_I = j
                      exit 
                   endif
               enddo
               do j=1,ysize(2)
                   if(Um_recyc(j).ge.threshold_bt*U_infty) then
                      jref_R = j
                      exit 
                   endif
               enddo
           
               Re_tau_R = yp(jref_R)/l_tau_R
               Re_tau_I = yp(jref_I)/l_tau_I
               
               delta_I_099 = yp(jref_I)
               delta_R_099 = yp(jref_R)
           
               do j=1,ysize(2)
                   if(yp(j)/l_tau_R.ge.Re_tau_I) then
                      jref = j
                      exit 
                   endif
               enddo
               
               if (nrank.eq.0) write(*,*) 'Largest length of skipt time instances for IOIM sampling:',N_Upsilon
               
               alp = 4.0
               bet = 0.2
               
               bdamp = 0.5*(1.-tanh(4.*(yp/delta_I_099-2.)))
               
               !!! Rescaling
               gama8 = delta_R_099/delta_I_099
               gama  = l_tau_R/l_tau_I
               gamaU = u_tau_R/u_tau_I
               
               call set_lund_maps(gama, gama8, delta_I_099, alp, bet)

               length_xR = xp(indU_p2)-xp(indU_p1)
               
           endif
           
           
           
       endif
       
       if (itime .eq. 1) then
           
           if (nrank .eq. 0) then
               if (i_zshift .ne.0) then
                   write(*,*) '********************** Spanwise shift called! **********************'
               else
                   write(*,*) '********************** Spanwise shift NOT called! **********************'
               endif
           endif
           
            if (nrank .eq. 0) then

                allocate(UVW_TYZ_lag_half( N_Gamma , ny , nz/2 , nv ))
                allocate(syy_half( N_I , nz/2 , nv ))
                allocate(cyt_half( N_I , ny , nz/2 , nv ))
                
                allocate(UVW_TYZ_lag( N_Gamma , ny , nz , nv ))
                allocate(syy( N_I , nz , nv ))
                allocate(cyt( N_I , ny , nz , nv ))
                allocate(UP2( ny , nz , nv , width_interp ))

                UVW_TYZ_lag_half = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                syy_half     = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                cyt_half     = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                
                UVW_TYZ_lag = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                syy     = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                cyt     = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                UP2     = zero

                file_exists = .false.
                file_match = .false.
                icyc10_rec = 0
                ! A cold start must not inherit an inlet pool from a previous run.
                if (irestart .ne. 0) then
                   call ioim_checkpoint_load_binary("data/inlet/ioim_checkpoint.bin", &
                        file_exists, file_match, icyc10_rec, &
                        UVW_TYZ_lag_half, syy_half, cyt_half, UP2, &
                        N_I,N_Gamma, ny, nz, nv, width_interp)
                   if (file_exists .and. file_match) then
                      UVW_TYZ_lag(:,:,1:nz:2,:) = UVW_TYZ_lag_half
                      syy(:,1:nz:2,:) = syy_half
                      cyt(:,:,1:nz:2,:) = cyt_half
                   else
                      icyc10_rec = 0
                   endif
                else
                   write(*,*) 'Cold start: ioim_checkpoint.bin is not read.'
                endif
                
            endif
            
            call MPI_Bcast(file_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, iermpi)
            call MPI_Bcast(file_match, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, iermpi)
            call MPI_Bcast(icyc10_rec, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, iermpi)
              
            if (nrank .eq. 0) then
                
                !write(*,*) nrank,file_exists , file_match

                if (file_exists .and. file_match) then

                    do i = 1,n_mrank

                        i_rank       = m_proc(i,1)
                        ystart3_recv = m_proc(i,2)
                        yend3_recv   = m_proc(i,3)
                        nnz = yend3_recv-ystart3_recv+1

                        if (i_rank.ne.0) then
                            call MPI_Send(UVW_TYZ_lag(:,:,ystart3_recv:yend3_recv,:), (N_Gamma * ny * nnz * nv), &
                                    MPI_DOUBLE_COMPLEX, i_rank, 101, MPI_COMM_WORLD, iermpi)
                            call MPI_Send(syy(:,ystart3_recv:yend3_recv,:), (N_I * nnz * nv), &
                                    MPI_DOUBLE_COMPLEX, i_rank, 103, MPI_COMM_WORLD, iermpi)
                            call MPI_Send(cyt(:,:,ystart3_recv:yend3_recv,:), (N_I * ny * nnz * nv), &
                                    MPI_DOUBLE_COMPLEX, i_rank, 104, MPI_COMM_WORLD, iermpi)
                            call MPI_Send(UP2(:,ystart3_recv:yend3_recv,:,:), (ny * nnz * nv * width_interp), &
                                    MPI_DOUBLE_PRECISION, i_rank, 105, MPI_COMM_WORLD, iermpi)
                        else
                            UVW_record_TYZ_lagAll = UVW_TYZ_lag(:,:,ystart(3):yend(3),:)
                            syy_cpu       = syy(:,ystart(3):yend(3),:)
                            cyt_cpu       = cyt(:,:,ystart(3):yend(3),:)
                            u_sq           = UP2(:,ystart(3):yend(3),:,:)
                        endif
                    enddo

                    write(*,'(A,I0,A)') '>>>>>>>>>> ***** Binary checkpoint exists with ICYC10_REC = ', &
                                            icyc10_rec, '. I would read it. ***** <<<<<<<<<<'

                else

                    print *, '>>>>>>>>>> ***** Binary checkpoint does not exist or size mismatches. I would initialize it. ***** <<<<<<<<<<'
                    call execute_command_line('mkdir -p data/inlet')

                endif

                deallocate(UVW_TYZ_lag, syy, cyt, UVW_TYZ_lag_half, syy_half, cyt_half, UP2)

            endif
              
            if (file_exists .and. file_match) then
                
                if ( (nrank .ne. 0) .and. (i_inlet .eq. 1) ) then
                    
                    call MPI_Recv(UVW_record_TYZ_lagAll, (N_Gamma * ny * ysize(3) * nv) ,  &
                            & MPI_DOUBLE_COMPLEX, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                            
                    call MPI_Recv(syy_cpu, (N_I * ysize(3) * nv) ,  &
                            & MPI_DOUBLE_COMPLEX, 0, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                            
                    call MPI_Recv(cyt_cpu, (N_I * ny * ysize(3) * nv) ,  &
                            & MPI_DOUBLE_COMPLEX, 0, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                            
                    call MPI_Recv(u_sq, (ny * ysize(3) * nv * width_interp) ,  &
                            & MPI_DOUBLE_PRECISION, 0, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                    
                endif
                
            endif
          
       endif
       
       if (ystart(1) .eq. 1) then
           
           icyc1_now  = itime
           icyc10_now = icyc1_now/n_sp
           
           ! --------------------------------------------------------------------------------------------------------------
           icyc1  = icyc1_now + icyc10_rec*n_sp
           icyc10 = icyc1/n_sp
           icyc11 = icyc10 - (N_I - 1)
           icyc00 = mod(icyc1,n_sp)
           
           ! ------------------------------------------------------------------------------------
           if( icyc00.eq.0 ) then
              
              iitr    = mod((icyc10-1), N_Gamma )+1
              UVW_record_TYZ_lagAll(iitr,:,:,1) = ta2(2,:,:) + zimg*ta2(3,:,:)
              UVW_record_TYZ_lagAll(iitr,:,:,2) = tb2(2,:,:) + zimg*tb2(3,:,:)
              UVW_record_TYZ_lagAll(iitr,:,:,3) = tc2(2,:,:) + zimg*tc2(3,:,:)
              
              if(icyc10.ge.N_I) then
                  
                  ind_store = mod((icyc10-1),N_Gamma)+1
                  
                  do i = 1, N_I
                      sq_start(i) = mod(ind_store - N_I + i - 1 + N_Gamma, N_Gamma) + 1
                  end do

                  ! Oldest N_I samples in the circular pool, used to build the IOIM block.
                  do i = 1, N_I
                      sq_end(i) = mod(ind_store + i - 1, N_Gamma) + 1
                  end do
                  
                  UVW_record_TYZ = UVW_record_TYZ_lagAll(sq_end,:,:,:)
                  
                  if (mod(ystart(3),2).eq.1) then
                      iz1 = 1
                  else
                      iz1 = 2
                  endif
                  
                  do iiv = 1, nv
                      do iiz = iz1, ysize(3), 2
                          do iiy = 1, jref
                              input_X_cplx  = UVW_record_TYZ(:,iiy,iiz,iiv)
                              call fftw_execute_dft(plan_t_fwd, input_X_cplx, output_X_cplx)
                              wcalc_o(:,iiy,iiz,iiv) = output_X_cplx
                          end do
                      end do
                  end do
                  
                  wcalcR_cpu = wcalc_o(:,jref,:,:)
                  
                  if (icyc11 .ge. N_Upsilon) then
                    alpha_t = 1.0/real(N_Upsilon)
                  else
                    alpha_t = 1.0/real(icyc11)
                  endif
                  
                  do iiz = iz1,ysize(3),2
                    
                     do iiy = 1,jref
                        cyt1_cpu(:,iiy,:) = wcalc_o(:,iiy,iiz,:) * conjg( wcalcR_cpu(:,iiz,:) )
                     enddo
                    
                     syy1_cpu = wcalcR_cpu(:,iiz,:) * conjg( wcalcR_cpu(:,iiz,:) )
                                
                     cyt_cpu(:,:,iiz,:) = (1.0-alpha_t)*cyt_cpu(:,:,iiz,:) + alpha_t*cyt1_cpu
                     syy_cpu(:,iiz,:)   = (1.0-alpha_t)*syy_cpu(:,iiz,:)   + alpha_t*syy1_cpu
                  
                  enddo
                  
              endif
                  
              ! Stage 2/3: recycling pool is full, so the oldest N_I block can be time-separated.
              ! Stage 2 (icyc10 < N_Gamma+N_Upsilon): provisional inner signal is the
              ! undecomposed recycling fluctuation. Stage 3: remove the IOIM footprint.
              if( icyc10 .ge. N_Gamma ) then
              
                if( icyc10 .ge. N_Gamma + N_Upsilon ) then
                  do iiz = iz1, ysize(3), 2
                      do iiy = 1, jref
                          do iiv = 1, nv
                              do iit = 1, N_I
                                  tmp_ioim = (cyt_cpu(iit,iiy,iiz,iiv) / conjg(syy_cpu(iit,iiz,iiv))+eps_syy) * wcalc_o(iit,jref,iiz,iiv)
                                  wcalc_i(iit,iiy,iiz,iiv) = wcalc_o(iit,iiy,iiz,iiv) - tmp_ioim
                              end do
                          end do
                      end do
                  end do
                else
                  wcalc_i = wcalc_o
                endif
                
                do j=1,ysize(2)
         
                   wcalc_i(:,j,:,:) = (wt_Inner1(j)*wcalc_i(:,ind_Inner1(j),:,:) + wt_Inner2(j)*wcalc_i(:,ind_Inner2(j),:,:))  ! N_I,ny,ysize(3),3
                   wcalc_o(:,j,:,:) = (wt_Outer1(j)*wcalc_o(:,ind_Outer1(j),:,:) + wt_Outer2(j)*wcalc_o(:,ind_Outer2(j),:,:))
                   
                enddo
                
                eps_gamma = 1.0d-20
                num_loc  = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                deni_loc = 0.0
                deno_loc = 0.0
                
                do iiv = 1, nv
                    do iiz = 1, ysize(3)
                        do j = 1, ysize(2)
                            do iit = 1, N_I
                                num_loc(j,iiv)  = num_loc(j,iiv)  + conjg(wcalc_i(iit,j,iiz,iiv)) * wcalc_o(iit,j,iiz,iiv)
                                deni_loc(j,iiv) = deni_loc(j,iiv) + real(conjg(wcalc_i(iit,j,iiz,iiv)) * wcalc_i(iit,j,iiz,iiv), kind=mytype)
                                deno_loc(j,iiv) = deno_loc(j,iiv) + real(conjg(wcalc_o(iit,j,iiz,iiv)) * wcalc_o(iit,j,iiz,iiv), kind=mytype)
                            end do
                        end do
                    end do
                end do
                
                if (i_inlet .eq. 1) then
                    call MPI_ALLREDUCE(num_loc,  num_glb,  ysize(2)*nv, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_inlet, ierr)
                    call MPI_ALLREDUCE(deni_loc, deni_glb, ysize(2)*nv, MPI_DOUBLE_PRECISION, MPI_SUM, comm_inlet, ierr)
                    call MPI_ALLREDUCE(deno_loc, deno_glb, ysize(2)*nv, MPI_DOUBLE_PRECISION, MPI_SUM, comm_inlet, ierr)
                endif
                
                do iiv = 1, nv
                    do j = 1, ysize(2)
                        
                        gamma_io(j,iiv) = abs(num_glb(j,iiv)) / sqrt(deni_glb(j,iiv) * deno_glb(j,iiv) + eps_gamma)
                        
                        W_i(j,iiv) = sqrt(1 - (1 - (gamma_io(j,iiv)**2.0))*(W_o(j)**2.0) ) - gamma_io(j,iiv)*W_o(j)
                        
                        wcalc_i(:,j,:,iiv) = bdamp(j)*( W_i(j,iiv)*wcalc_i(:,j,:,iiv)/gamaU + W_o(j)*wcalc_o(:,j,:,iiv) )
                        
                    end do
                end do
                
                iitp = (N_I+1)/2

                expz(:) = exp(zimg * sq_ft(:) * real(iitp-1)) / real(N_I)
                
                do iiz = iz1, ysize(3), 2
                    do iiy = 1, ysize(2)
                        do iiv = 1, nv
                            
                            tmp_ioim = (0.0d0, 0.0d0)
                            
                            do iit = 1, N_I
                                tmp_ioim = tmp_ioim + wcalc_i(iit,iiy,iiz,iiv) * expz(iit)
                            end do
                            
                            up_cpuZ(iiy,iiz,iiv) = tmp_ioim
                            
                        end do
                    end do
                end do
                
                td2(1,:,:) =  real(up_cpuZ(:,:,1))
                td2(2,:,:) = aimag(up_cpuZ(:,:,1))
                      
                te2(1,:,:) =  real(up_cpuZ(:,:,2))
                te2(2,:,:) = aimag(up_cpuZ(:,:,2))
                      
                tf2(1,:,:) =  real(up_cpuZ(:,:,3))
                tf2(2,:,:) = aimag(up_cpuZ(:,:,3))
                      
                call transpose_y_to_z(td2,td3)
                call transpose_y_to_z(te2,te3)
                call transpose_y_to_z(tf2,tf3)
                
                scalez = 1.0_mytype / real(zsize(3),mytype)
                
                iiy_max = min(zsize(2), jref - zstart(2) + 1)

                do iiy = 1, iiy_max

                    output_Z_cplxm1(1) = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                    output_Z_cplxm2(1) = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)
                    output_Z_cplxm3(1) = cmplx(0.0_cdp, 0.0_cdp, kind=cdp)

                    do k = 1, zsize(3)/2
                        output_Z_cplxm1(k+1) = cmplx(td3(1,iiy,2*k-1), td3(2,iiy,2*k-1), kind=cdp)
                        output_Z_cplxm2(k+1) = cmplx(te3(1,iiy,2*k-1), te3(2,iiy,2*k-1), kind=cdp)
                        output_Z_cplxm3(k+1) = cmplx(tf3(1,iiy,2*k-1), tf3(2,iiy,2*k-1), kind=cdp)
                    end do

                    call fftw_execute_dft_c2r(plan_z_bwd1, output_Z_cplxm1, input_Z_realm1)
                    call fftw_execute_dft_c2r(plan_z_bwd2, output_Z_cplxm2, input_Z_realm2)
                    call fftw_execute_dft_c2r(plan_z_bwd3, output_Z_cplxm3, input_Z_realm3)

                    td3(1,iiy,:) = input_Z_realm1 * scalez
                    te3(1,iiy,:) = input_Z_realm2 * scalez
                    tf3(1,iiy,:) = input_Z_realm3 * scalez

                end do

                call transpose_z_to_y(td3,td2)
                call transpose_z_to_y(te3,te2)
                call transpose_z_to_y(tf3,tf2)
                          
                up_phy(1:jref,:,1) = td2(1,1:jref,:)
                up_phy(1:jref,:,2) = te2(1,1:jref,:)
                up_phy(1:jref,:,3) = tf2(1,1:jref,:)
                
                !-------------------------------------------------------
                
                ind_store = mod((icyc10-1),width_interp)+1
                
                do i = 1, width_interp
                   sqN(i) = mod(ind_store+i-1, width_interp) + 1
                end do
                
                u_sq(:,:,:,ind_store) = up_phy
                
                if (icyc10 .ge. N_Gamma + width_interp - 1) then
                    call interp_mid_natural_spline_fast_3d(u_sq, nyn, nzn, sqN, solvern, u_interp)
                endif
                
              endif
           
           endif
           
           if (icyc1_now.eq.1 .and. file_exists .and. file_match) then
               
               ind_store = mod((icyc10-1),width_interp)+1
               
               do i = 1, width_interp
                  sqN(i) = mod(ind_store+i-1, width_interp) + 1
               end do
               
               call interp_mid_natural_spline_fast_3d(u_sq, nyn, nzn, sqN, solvern, u_interp)
               
           endif
               
           ! UVW_TYZ, syy_cpu, cyt_cpu
           if( (mod(itime,icheckpoint).eq.0) ) then
                          
            if (nrank .eq. 0) then
               write(*,'(A,I0,A)') &
               & 'Output the checkpoint for inlet turbulence! Round:', icyc10/N_I,'.'
                            
               allocate(UVW_TYZ_lag( N_Gamma , ny , nz , nv ))
               allocate(syy( N_I , nz , nv ))
               allocate(cyt( N_I , ny , nz , nv ))
               allocate(UP2( ny , nz , nv , width_interp ))
                              
               UVW_TYZ_lag   = cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
               syy       = cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
               cyt       = cmplx(0.0_cdp,  0.0_cdp, kind=cdp)
               UP2       = zero
            endif
                          
            if (nrank .ne. 0) then
                
                call MPI_Send(UVW_record_TYZ_lagAll, (N_Gamma * ny * ysize(3) * nv) ,  &
                        & MPI_DOUBLE_COMPLEX, 0, 101, MPI_COMM_WORLD, iermpi)
                            
                call MPI_Send(syy_cpu, (N_I * ysize(3) * nv) ,  &
                        & MPI_DOUBLE_COMPLEX, 0, 103, MPI_COMM_WORLD, iermpi)
                            
                call MPI_Send(cyt_cpu, (N_I * ny * ysize(3) * nv) ,  &
                        & MPI_DOUBLE_COMPLEX, 0, 104, MPI_COMM_WORLD, iermpi)
                            
                call MPI_Send(u_sq, (ny * ysize(3) * nv * width_interp) ,  &
                        & MPI_DOUBLE_PRECISION, 0, 105, MPI_COMM_WORLD, iermpi)
                            
            else
                UVW_TYZ_lag(:,:,ystart(3):yend(3),:) = UVW_record_TYZ_lagAll
                syy(:,ystart(3):yend(3),:)       = syy_cpu
                cyt(:,:,ystart(3):yend(3),:)     = cyt_cpu
                UP2(:,ystart(3):yend(3),:,:)     = u_sq
            endif
                          
            if (nrank .eq. 0) then
                              
                do i = 1,n_mrank
                                
                    i_rank       = m_proc(i,1)
                    ystart3_recv = m_proc(i,2)
                    yend3_recv   = m_proc(i,3)
                    nnz = yend3_recv-ystart3_recv+1
                                
                    if (i_rank .ne. 0) then
                        allocate(buff_c(N_Gamma * ny * nnz * nv))
                        call MPI_Recv(buff_c, (N_Gamma * ny * nnz * nv) ,  &
                                & MPI_DOUBLE_COMPLEX, i_rank, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                        UVW_TYZ_lag(:,:,ystart3_recv:yend3_recv,:) = reshape(buff_c, [N_Gamma , ny , nnz , nv])
                        deallocate(buff_c)
                                
                        allocate(buff_c(N_I * nnz * nv))
                        call MPI_Recv(buff_c, (N_I * nnz * nv) ,  &
                                & MPI_DOUBLE_COMPLEX, i_rank, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                        syy(:,ystart3_recv:yend3_recv,:) = reshape(buff_c, [N_I , nnz , nv])
                        deallocate(buff_c)
                                
                        allocate(buff_c(N_I * ny * nnz * nv))
                        call MPI_Recv(buff_c, (N_I * ny * nnz * nv) ,  &
                                & MPI_DOUBLE_COMPLEX, i_rank, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                        cyt(:,:,ystart3_recv:yend3_recv,:) = reshape(buff_c, [N_I , ny , nnz , nv])
                        deallocate(buff_c)
                                
                        allocate(buff_r(ny * nnz * nv * width_interp))
                        call MPI_Recv(buff_r, (ny * nnz * nv * width_interp) ,  &
                                & MPI_DOUBLE_PRECISION, i_rank, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
                        UP2(:,ystart3_recv:yend3_recv,:,:) = reshape(buff_r, [ny , nnz , nv , width_interp])
                        deallocate(buff_r)
                    endif
                                
                enddo
                          
                call execute_command_line('mkdir -p data/inlet')

                call ioim_checkpoint_save_binary("data/inlet/ioim_checkpoint.bin", &
                     icyc10, UVW_TYZ_lag(:,:,1:nz:2,:), syy(:,1:nz:2,:), cyt(:,:,1:nz:2,:), UP2, &
                     N_I,N_Gamma, ny, nz, nv, width_interp)

                write(*,*) '~~~~~~~~~~~~~~~~~~~ Finish writing binary checkpoint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

                deallocate(UVW_TYZ_lag, UP2, syy, cyt)
                              
            endif
                          
           endif
           ! ------------------------------------------------------------------------------------
           
           ! ........................................... IOIM ....................................................
           ! Stage 1: icyc10 < N_Gamma, instantaneous LWS.
           ! Stage 2: pool full, provisional inner signal, spline once the inlet window is full.
           ! Stage 3: icyc10 >= N_Gamma+N_Upsilon, full IOIM decomposition.
           t_ioim_cal = N_Gamma + N_Upsilon
           
           if ( (mod(itime,ilist).eq.0 .or. itime.eq.1) .and. (nrank .eq. 0) ) then
              if (icyc10 .lt. N_Gamma) then
                 write(*,'(A,I0,A,I0,A)') &
                 & '>>>>>>>>>> ***** Inflow stage 1/3 (LWS): ', icyc10,'/',N_Gamma,'. ***** <<<<<<<<<<'
              else if (icyc10 .lt. t_ioim_cal) then
                 write(*,'(A,I0,A,I0,A)') &
                 & '>>>>>>>>>> ***** Inflow stage 2/3 (time-separated, undecomposed inner): ', icyc10,'/',t_ioim_cal,'. ***** <<<<<<<<<<'
              else
                 write(*,'(A,I0,A)') &
                 & '>>>>>>>>>> ***** Inflow stage 3/3 (IOIM): sample ', icyc10,'. ***** <<<<<<<<<<'
              endif
           endif
           
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
           if ( icyc10 .ge. N_Gamma + width_interp - 1 ) then
               
               icyc01 = mod(icyc1,n_sp)+1
               
               Un_I_YZ = u_interp(:,:,1,icyc01)
               Vn_I_YZ = u_interp(:,:,2,icyc01)
               Wn_I_YZ = u_interp(:,:,3,icyc01)
               
               if(i_zshift .ne. 0) then
               
                   td2(1,:,:) = Un_I_YZ
                   te2(1,:,:) = Vn_I_YZ
                   tf2(1,:,:) = Wn_I_YZ
               
                   call half_span_shift(td2, te2, tf2, ta3, tb3, tc3, u_shift, v_shift, w_shift, ta2, tb2, tc2, n_shift)
               
               else
                   
                   ta2(1,:,:) = Un_I_YZ
                   tb2(1,:,:) = Vn_I_YZ
                   tc2(1,:,:) = Wn_I_YZ
               
               endif
               
               do k = 1,ysize(3)
                  ta2(1,:,k) = ta2(1,:,k) + Um_inlet
                  tb2(1,:,k) = tb2(1,:,k) + Vm_inlet
               enddo
               
           else
               
               ! ........................................... LWS ....................................................
               call blend_lund_fluctuations(ta2(1,:,:), tb2(1,:,:), tc2(1,:,:), Ut_I_YZ, Vt_I_YZ, Wt_I_YZ)

               if (i_zshift .ne. 0) then
                  td2(1,:,:) = Ut_I_YZ
                  te2(1,:,:) = Vt_I_YZ
                  tf2(1,:,:) = Wt_I_YZ
                  call half_span_shift(td2, te2, tf2, ta3, tb3, tc3, u_shift, v_shift, w_shift, ta2, tb2, tc2, n_shift)
               else
                  ta2(1,:,:) = Ut_I_YZ
                  tb2(1,:,:) = Vt_I_YZ
                  tc2(1,:,:) = Wt_I_YZ
               endif
               
               do k = 1,ysize(3)
                  ta2(1,:,k) = ta2(1,:,k) + Um_inlet
                  tb2(1,:,k) = tb2(1,:,k) + Vm_inlet
               enddo
               
           endif
           
       endif
       
       !-------------------------------------------------------------------------------------------------------------
       
       elseif (method_rr .eq. 2) then
       ! LWS-S: Lund rescaling of the instantaneous recycling plane.
       if (itime .eq. 1) then
           
          niu = 1.0/re
          threshold_bt = 0.99   ! threshlod to determine the boundary thickness
          U_infty = 1.0
           
          u_tau_R = sqrt(niu*Um_recyc(2)/yp(2))
          u_tau_I = sqrt(niu*Um_inlet(2)/yp(2))
           
          l_tau_R = niu/u_tau_R
          l_tau_I = niu/u_tau_I
           
          do j=1,ysize(2)
              if(Um_inlet(j).ge.threshold_bt*U_infty) then
              jref_I = j
              exit 
              endif
          enddo
          do j=1,ysize(2)
              if(Um_recyc(j).ge.threshold_bt*U_infty) then
              jref_R = j
              exit 
              endif
          enddo
           
          Re_tau_R = yp(jref_R)/l_tau_R
          Re_tau_I = yp(jref_I)/l_tau_I
               
          delta_I_099 = yp(jref_I)
          delta_R_099 = yp(jref_R)
           
          do j=1,ysize(2)
              if(yp(j)/l_tau_R.ge.Re_tau_I) then
              jref = j
              exit 
              endif
          enddo
           
          alp = 4.0
          bet = 0.2
               
          bdamp = 0.5*(1.-tanh(4.*(yp/delta_I_099-2.)))
               
          !!! Rescaling
          gama8 = delta_R_099/delta_I_099
          gama  = l_tau_R/l_tau_I
          gamaU = u_tau_R/u_tau_I
               
          if (nrank .eq. 0) then
              write(*,'(A,F8.3,A,F8.3,A)') '>>>***** Inlet of LWS called! Outer scaling factor: ',gama8,', Inner scaling factor: ',gama,'. ***** <<<'
          endif
               
          call set_lund_maps(gama, gama8, delta_I_099, alp, bet)

       endif
       
       if (ystart(1) .eq. 1) then

            ! ........................................... LWS ....................................................
            call blend_lund_fluctuations(ta2(1,:,:), tb2(1,:,:), tc2(1,:,:), Ut_I_YZ, Vt_I_YZ, Wt_I_YZ)
            
            !ta2(1,:,:) = Ut_I_YZ
            !tb2(1,:,:) = Vt_I_YZ
            !tc2(1,:,:) = Wt_I_YZ
            
            if(i_zshift .ne. 0) then
               
                td2(1,:,:) = Ut_I_YZ
                te2(1,:,:) = Vt_I_YZ
                tf2(1,:,:) = Wt_I_YZ
               
                call half_span_shift(td2, te2, tf2, ta3, tb3, tc3, u_shift, v_shift, w_shift, ta2, tb2, tc2, n_shift)
               
            else
                   
                ta2(1,:,:) = Ut_I_YZ
                tb2(1,:,:) = Vt_I_YZ
                tc2(1,:,:) = Wt_I_YZ
               
            endif
            
            do k = 1,ysize(3)
               ta2(1,:,k) = ta2(1,:,k) + Um_inlet
               tb2(1,:,k) = tb2(1,:,k) + Vm_inlet
            enddo
               
       endif
       
       ! ---------------------------------------------------------------------------------------------------
       
       elseif (method_rr .eq. 3) then
       ! LWS-R: Lund rescaling with a running-mean recycling profile.
           
       
       !write(*,*) 'ITIME!!!!!!!!!!!!!',itime,icheckpoint
       
       if (itime .eq. 1) then
              
              inquire(file="data/inlet/V_LWS_recyc.txt", exist=file_exists)
              
              nlines = 0
              filename = "data/inlet/V_LWS_recyc.txt"
              open(newunit=iunit, file=filename, status='old', action='read', iostat=ios)
              do
                  read(iunit, '(A)', iostat=ios) testlengthx
                  if (ios /= 0) exit
                  nlines = nlines + 1
              enddo
              close(iunit)
                
              if (nlines .eq. ny) then
                  file_match = .true.
              else
                  file_match = .false.
       endif
              
              icyc10_rec  = 0
          
              if ((file_exists .eqv. .true.) .and. (file_match .eqv. .true.)) then
                                    
                 filename = "data/inlet/U_LWS_recyc.txt"
                 open(newunit=iunit, file=filename, status='old', action='read', iostat=ios)
                 read(iunit, '(I10)', iostat=ios) icyc10_rec
                 do j = 1,ny
                    read(iunit, '(1(E30.15))') UmT_recyc_LWS(j)
                 enddo
                 close(iunit)
                    
                 filename = "data/inlet/V_LWS_recyc.txt"
                 open(newunit=iunit, file=filename, status='old', action='read', iostat=ios)
                 do j = 1,ny
                    read(iunit, '(1(E30.15))') VmT_recyc_LWS(j)
                 enddo
                 close(iunit)
                 
              else
                 
                 if (nrank .eq. 0) then
                     print *, '>>>>>>>>>> ***** Historical mean data does not exist or match the size. I would initialize it. ***** <<<<<<<<<<'
                     call execute_command_line('mkdir -p data/inlet')
                 endif
                
                 icyc10_rec  = 0
                
              endif
          
       endif
       
        niu = 1.0/re
        threshold_bt = 0.99   ! threshlod to determine the boundary thickness
        U_infty = 1.0
               
        delta_I_099 = delta_in
        
        !!!
           
        alp = 4.0
        bet = 0.2
               
        bdamp = 0.5*(1.-tanh(4.*(yp/delta_I_099-2.)))
       
        icyc10 = itime+icyc10_rec
        
        T_LWS = min(10.0*delta_I_099/U_infty , real(icyc10)*dt)
               
        UmT_recyc_LWS = (Um_recyc_LWS*dt + UmT_recyc_LWS*(T_LWS-dt))/T_LWS
        VmT_recyc_LWS = (Vm_recyc_LWS*dt + VmT_recyc_LWS*(T_LWS-dt))/T_LWS
           
        do j=1,ysize(2)
            if(UmT_recyc_LWS(j).ge.threshold_bt*U_infty) then
            jref_R = j
            exit 
            endif
        enddo
        
        wt1 = (UmT_recyc_LWS(jref_R) - threshold_bt*U_infty)/(UmT_recyc_LWS(jref_R) - UmT_recyc_LWS(jref_R-1))
        wt2 = 1.0-wt1
        delta_R_099 = wt1*yp(jref_R-1) + wt2*yp(jref_R)
        
        if( (mod(itime,icheckpoint).eq.0) .and. (nrank.eq.0) ) then
            
           filename = "data/inlet/U_LWS_recyc.txt"
           open(newunit=iunit, file=filename, status='replace', action='write')
           write(iunit, '(I10)') icyc10
           do j = 1,ny
              write(iunit, '(1(E30.15))') UmT_recyc_LWS(j)
           enddo
           close(iunit)
           
           filename = "data/inlet/V_LWS_recyc.txt"
           open(newunit=iunit, file=filename, status='replace', action='write')
           do j = 1,ny
              write(iunit, '(1(E30.15))') VmT_recyc_LWS(j)
           enddo
           close(iunit)
           
        endif
               
        !!! Rescaling
        gama8 = delta_R_099/delta_I_099
        gama  = gama8**0.125
        gamaU = 1.0/gama
               
        if ( (mod(itime,ilist).eq.0 .or. itime.eq.1) .and. (nrank .eq. 0) ) then
            write(*,'(A,F8.3,A,F8.3,A,F8.3,A)') &
            & '>>>>>>>>>> ***** Original LWS. The dynamic parameters: gama8: ', gama8,', T_LWS: ', &
            & T_LWS,'/',10.0*delta_I_099/U_infty,'. ***** <<<<<<<<<<'
        endif
               
        call set_lund_maps(gama, gama8, delta_I_099, alp, bet)

       if (ystart(1) .eq. 1) then

            ! ........................................... LWS ....................................................
            
            call blend_lund_fluctuations(ta2(1,:,:), tb2(1,:,:), tc2(1,:,:), Ut_I_YZ, Vt_I_YZ, Wt_I_YZ)
            do j=1,ysize(2)
                
                Um_inner_Z_LWS = (wt_Inner1(j)*UmT_recyc_LWS(ind_Inner1(j)) + wt_Inner2(j)*UmT_recyc_LWS(ind_Inner2(j)))
                Um_inner_Z_LWS = Um_inner_Z_LWS/gamaU
                
                Um_outer_Z_LWS = (wt_Outer1(j)*UmT_recyc_LWS(ind_Outer1(j)) + wt_Outer2(j)*UmT_recyc_LWS(ind_Outer2(j)))
                Um_outer_Z_LWS = U_infty + bdamp(j)*(Um_outer_Z_LWS - U_infty)/gamaU
                
                Vm_inner_Z_LWS = (wt_Inner1(j)*VmT_recyc_LWS(ind_Inner1(j)) + wt_Inner2(j)*VmT_recyc_LWS(ind_Inner2(j)))
                Vm_outer_Z_LWS = (wt_Outer1(j)*VmT_recyc_LWS(ind_Outer1(j)) + wt_Outer2(j)*VmT_recyc_LWS(ind_Outer2(j)))
                
                Um_inlet_LWS(j,:) = (1.0-W_yita(j))*Um_inner_Z_LWS + W_yita(j)*Um_outer_Z_LWS
                Vm_inlet_LWS(j,:) = (1.0-W_yita(j))*Vm_inner_Z_LWS + W_yita(j)*Vm_outer_Z_LWS
             
            enddo
               
            !ta2(1,:,:) = Ut_I_YZ + Um_inlet_LWS
            !tb2(1,:,:) = Vt_I_YZ + Vm_inlet_LWS
            !tc2(1,:,:) = Wt_I_YZ
            
            if(i_zshift .ne. 0) then
               
                td2(1,:,:) = Ut_I_YZ + Um_inlet_LWS
                te2(1,:,:) = Vt_I_YZ + Vm_inlet_LWS
                tf2(1,:,:) = Wt_I_YZ
               
                call half_span_shift(td2, te2, tf2, ta3, tb3, tc3, u_shift, v_shift, w_shift, ta2, tb2, tc2, n_shift)
               
            else
                   
                ta2(1,:,:) = Ut_I_YZ + Um_inlet_LWS
                tb2(1,:,:) = Vt_I_YZ + Vm_inlet_LWS
                tc2(1,:,:) = Wt_I_YZ
               
            endif
           
       endif
       
       ! ---------------------------------------------------------------------------------------------------
       
       endif
       ! ***********************************************************************************************************************
       
       call transpose_y_to_x(ta2,uxr)
       call transpose_y_to_x(tb2,uyr)
       call transpose_y_to_x(tc2,uzr)
       
    endif
    endif
    !-------------------------------------------------------
    
    ! Inlet. iin=4 uses the recycling plane; every other value keeps the original Blasius inlet.
    if (iin.eq.4) then
    
        do k=1,xsize(3)
          do j=1,xsize(2)
            bxx1(j,k)=uxr(1,j,k)
            bxy1(j,k)=uyr(1,j,k)
            bxz1(j,k)=uzr(1,j,k)
          enddo
        enddo
    
    else

        call blasius()
        if (iscalar==1) then
           do k=1,xsize(3)
              do j=1,xsize(2)
                 phi(1,:,:,:)=zptwofive
                 if ((xstart(2)==1)) then
                    phi(:,1,:,:) = one
                 endif
                 if ((xend(2)==ny)) THEN
                    phi(:,xsize(2),:,:) = zptwofive
                 endif
              enddo
           enddo
        endif

    endif
    

    !OUTFLOW based on a 1D convection equation

    udx=one/dx
    udy=one/dy
    udz=one/dz
    uddx=half/dx
    uddy=half/dy
    uddz=half/dz

    do k=1,xsize(3)
       do j=1,xsize(2)

          cx=ux(nx,j,k)*gdt(itr)*udx

          if (cx<zero) cx=zero
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
          if (iscalar==1) phi(nx,:,:,:) =  phi(nx,:,:,:) - cx*(phi(nx,:,:,:)-phi(nx-1,:,:,:))
          enddo
    enddo

    !! Bottom Boundary
    if (ncly1 == 2) then
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i, k) = zero
          byy1(i, k) = zero
          byz1(i, k) = zero
        enddo
      enddo
      
      elseif (ncly1 == 1) then
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byx1(i, k) = ux(i, 2, k)
             byy1(i, k) = zero
             byz1(i, k) = uz(i, 2, k)
          enddo
       enddo
       
    endif
    !! Top Boundary
    if (nclyn == 2) then
    
    !write(*,*) 'Top boundary',xsize(1),xsize(2),xsize(3),nrank
    
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = uy(i, xsize(2) - 1, k)
             byzn(i, k) = uz(i, xsize(2) - 1, k)
          enddo
       enddo
       
      
      elseif (nclyn == 1) then
       do k = 1, xsize(3)
          do i = 1, xsize(1)
             byxn(i, k) = ux(i, xsize(2) - 1, k)
             byyn(i, k) = zero
             byzn(i, k) = uz(i, xsize(2) - 1, k)
          enddo
       enddo
       
    endif
    
    !!! Side Boundaries
    !! Side Boundary - 1
    if (nclz1 == 2) then
      do k = 1, xsize(2)
        do i = 1, xsize(1)
          bzx1(i, k) = zero
          bzy1(i, k) = zero
          bzz1(i, k) = zero
        enddo
      enddo
      
      elseif (nclz1 == 1) then
       do k = 1, xsize(2)
          do i = 1, xsize(1)
             bzx1(i, k) = ux(i, k, 2)
             bzy1(i, k) = uy(i, k, 2)
             bzz1(i, k) = zero
          enddo
       enddo
       
    endif
    !! SideBoundary - 2
    if (nclzn == 2) then
       do k = 1, xsize(2)
          do i = 1, xsize(1)
             bzxn(i, k) = ux(i, k, xsize(3) - 1)
             bzyn(i, k) = uy(i, k, xsize(3) - 1)
             bzzn(i, k) = uz(i, k, xsize(3) - 1)
          enddo
       enddo
       
      
      elseif (nclzn == 1) then
       do k = 1, xsize(2)
          do i = 1, xsize(1)
             bzxn(i, k) = ux(i, k, xsize(3) - 1)
             bzyn(i, k) = uy(i, k, xsize(3) - 1)
             bzzn(i, k) = zero
          enddo
       enddo
       
    endif
    !!!!!!!!!!!!!!!!!!

    !SCALAR   
    if (itimescheme/=7) then
    if (iscalar/=0) then
          if ((nclyS1==2).and.(xstart(2)==1)) then
             !! Generate a hot patch on bottom boundary
             phi(1,1,:,:) = one
          endif
          if ((nclySn==2).and.(xend(2)==ny)) THEN
             phi(1,xsize(2),:,:) = phi(1,xsize(2)-1,:,:)
          endif
    endif
    endif

    !update of the flow rate (what is coming in the domain is getting out)
    call tbl_flrt(ux,uy,uz)

    return
  end subroutine boundary_conditions_tbl
  
  !********************************************************************
  !********************************************************************

  subroutine inflow (ux,uy,uz,phi)
  !
  !*******************************************************************************
  
    USE param
    USE variables
    USE MPI
    USE var, only: ux_inflow, uy_inflow, uz_inflow
    !USE var, only: ux_rescale, uy_rescale, uz_rescale
  
    implicit none
  
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    real(mytype) :: um
    integer :: i,j,k,itime_input
    
    um=zpfive*(u1+u2)
    do k=1,xsize(3)
    do j=1,xsize(2)
      bxx1(j,k)=um
      bxy1(j,k)=zero
      bxz1(j,k)=zero
    enddo
    enddo

    ! if(nrank.eq.1) write(*,*) '2222222222222222222222222'
    
    if (iin.eq.1.or.iin.eq.2) then
    
      call random_number(bxo)
      call random_number(byo)
      call random_number(bzo)

      do k=1,xsize(3)
      do j=1,xsize(2)
        bxx1(j,k)=bxx1(j,k)+(two*bxo(j,k)-one)*inflow_noise*um
        bxy1(j,k)=bxy1(j,k)+(two*byo(j,k)-one)*inflow_noise*um
        bxz1(j,k)=bxz1(j,k)+(two*bzo(j,k)-one)*inflow_noise*um
        if (iscalar.eq.1) then
          phi(1,j,k,:)=one
        endif
      enddo
      enddo
    else if (iin.eq.3) then
      ! Reading from files (when precursor simulations exist)
      itime_input=mod(itime,ntimesteps)
      if (itime_input==0) itime_input=ntimesteps
      if (mod(itime,ilist)==0.and.nrank==0) print *,'Reading inflow from a file, time step: ', itime_input
      do k=1,xsize(3)
      do j=1,xsize(2)
        ! Case 1: Inflow is turbulence added to mean flow profile
        !bxx1(j,k)=bxx1(j,k)+ux_inflow(itime_input,j,k)
        !bxy1(j,k)=bxy1(j,k)+uy_inflow(itime_input,j,k)
        !bxz1(j,k)=bxz1(j,k)+uz_inflow(itime_input,j,k)
        ! Case 2: Inflow is full velocity field
        bxx1(j,k)=ux_inflow(itime_input,j,k)
        bxy1(j,k)=uy_inflow(itime_input,j,k)
        bxz1(j,k)=uz_inflow(itime_input,j,k)
      enddo
      enddo
    !else if (iin.eq.4) then
    !  ! Reading from files (when precursor simulations exist)
    !  itime_input=mod(itime,ntimesteps)
    !  if (itime_input==0) itime_input=ntimesteps
    !  if (mod(itime,ilist)==0.and.nrank==0) print *,'Reading inflow from a file, time step: ', itime_input
    !  do k=1,xsize(3)
    !  do j=1,xsize(2)
    !    bxx1(j,k)=ux_rescale(itime_input,j,k)
    !    bxy1(j,k)=uy_rescale(itime_input,j,k)
    !    bxz1(j,k)=uz_rescale(itime_input,j,k)
    !  enddo
    !  enddo
      
    endif

    return
  end subroutine inflow 
  
  !*******************************************************************************
  !
  subroutine outflow (ux,uy,uz,phi)
  !
  !*******************************************************************************

    USE param
    USE variables
    USE MPI

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609._mytype
    uxmin=1609._mytype
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
        if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
      enddo
    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE,uxmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uxmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    cx=zpfive*(uxmax+uxmin)*gdt(itr)*udx
    do k=1,xsize(3)
      do j=1,xsize(2)
        bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
        bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
        bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
        if (iscalar.eq.1) then
          phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
        endif
      enddo
    enddo

    return
  end subroutine outflow 

  !********************************************************************
  !********************************************************************
  subroutine blasius()

    use decomp_2d_io
    use MPI
    use param, only : zero, zptwo, zpeight, one, nine

    implicit none

    real(mytype) :: eta_bl, f_bl, g_bl, x_bl,h_bl
    real(mytype) :: delta_int, delta_eta, eps_eta


    real(mytype) :: x, y, z
    integer :: i, j, k, is

    do k=1,xsize(3)
       do j=1,xsize(2)
          if (istret.eq.0) y=real(j+xstart(2)-1-1,mytype)*dy
          if (istret.ne.0) y=yp(j+xstart(2)-1)

          eta_bl=y*real(4.91,mytype)/nine

          !OLD POLYNOMIAL FITTING

          delta_eta=zero
          eps_eta=zero
          delta_int=zptwo

          if (eta_bl>=(real(7.5,mytype)/nine)) then
             delta_eta=eta_bl-real(7.5,mytype)/nine
             eta_bl=real(7.5,mytype)/nine
             eps_eta=real(0.00015,mytype)
          end if

          f_bl=1678.64209592595000_mytype*eta_bl**14-11089.69250174290_mytype*eta_bl**13 &
               +31996.435014067000_mytype*eta_bl**12-52671.52497797990_mytype*eta_bl**11 &
               +54176.169116766700_mytype*eta_bl**10-35842.82047060970_mytype*eta_bl**9  &
               +15201.308887124000_mytype*eta_bl**8 -4080.171379356480_mytype*eta_bl**7  &
               +702.12963452810300_mytype*eta_bl**6 -56.20639258053180_mytype*eta_bl**5  &
               -17.018112827391400_mytype*eta_bl**4 +0.819582894357566_mytype*eta_bl**3  &
               -0.0601348202321954_mytype*eta_bl**2 +2.989739912704050_mytype*eta_bl**1

          f_bl=f_bl+(1-exp(-delta_eta/delta_int))*eps_eta



          if (eta_bl >= (7.15_mytype/nine)) then
             delta_int=zpeight
             delta_eta=eta_bl-7.15_mytype/nine
             eta_bl   =       7.15_mytype/nine
             eps_eta  =     0.0005_mytype
          end if

          g_bl=4924.052847797540_mytype*eta_bl**14-34686.2970972733000_mytype*eta_bl**13 &
               +108130.253843618_mytype*eta_bl**12-195823.099139525000_mytype*eta_bl**11 &
               +227305.908339065_mytype*eta_bl**10-176106.001047617000_mytype*eta_bl**9  &
               +92234.5885895112_mytype*eta_bl**8 -32700.3687158807000_mytype*eta_bl**7  &
               +7923.51008739107_mytype*eta_bl**6 -1331.09245288739000_mytype*eta_bl**5  &
               +130.109496961069_mytype*eta_bl**4 -7.64507811014497000_mytype*eta_bl**3  &
               +6.94303207046209_mytype*eta_bl**2 -0.00209716712558639_mytype*eta_bl**1 ! &

          g_bl=g_bl+(1-exp(-delta_eta/delta_int))*eps_eta



          x_bl=one/(4.91_mytype**2*xnu)

          bxx1(j,k)=f_bl/1.0002014996204402_mytype/1.0000000359138641_mytype !To assure 1.0 in infinity
          bxy1(j,k)=g_bl*sqrt(xnu/x_bl)/1.000546554_mytype
          bxz1(j,k)=zero
       enddo
    enddo

    !STORE VALUE F_BL_INF G_BL_INF (ONLY ONE MORE TIME)------------------

    y=yly
    eta_bl=y*4.91_mytype/nine  !The 9 is due to interpolation

    delta_eta=zero
    eps_eta=zero
    delta_int=zptwo

    if (eta_bl>=(7.5_mytype/nine)) then
       delta_eta=eta_bl-7.5_mytype/nine
       eta_bl   =       7.5_mytype/nine
       eps_eta  =   0.00015_mytype
    end if

    !To assure 1.0 in infinity
    f_bl_inf=1678.6420959259500_mytype*eta_bl**14-11089.69250174290_mytype*eta_bl**13 &
            +31996.435014067000_mytype*eta_bl**12-52671.52497797990_mytype*eta_bl**11 &
            +54176.169116766700_mytype*eta_bl**10-35842.82047060970_mytype*eta_bl**9  &
            +15201.308887124000_mytype*eta_bl**8 -4080.171379356480_mytype*eta_bl**7  &
            +702.12963452810300_mytype*eta_bl**6 -56.20639258053180_mytype*eta_bl**5  &
            -17.018112827391400_mytype*eta_bl**4 +0.819582894357566_mytype*eta_bl**3  &
            -0.0601348202321954_mytype*eta_bl**2 +2.989739912704050_mytype*eta_bl**1


    f_bl_inf=f_bl_inf+(1-exp(-delta_eta/delta_int))*eps_eta
    f_bl_inf=f_bl_inf/1.0002014996204402_mytype/1.0000000359138641_mytype !To assure 1.0 in infinity

#ifdef DEBG
    if (nrank == 0) write(*,*)'f_bl_inf ', f_bl_inf
#endif

    if (eta_bl>= (7.15_mytype/nine)) then
       delta_int=zpeight
       delta_eta=eta_bl-7.15_mytype/nine
       eta_bl   =       7.15_mytype/nine
       eps_eta  =     0.0005_mytype
    end if

    g_bl_inf=4924.05284779754_mytype*eta_bl**14-34686.2970972733000_mytype*eta_bl**13 &
            +108130.253843618_mytype*eta_bl**12-195823.099139525000_mytype*eta_bl**11 &
            +227305.908339065_mytype*eta_bl**10-176106.001047617000_mytype*eta_bl**9  &
            +92234.5885895112_mytype*eta_bl**8 -32700.3687158807000_mytype*eta_bl**7  &
            +7923.51008739107_mytype*eta_bl**6 -1331.09245288739000_mytype*eta_bl**5  &
            +130.109496961069_mytype*eta_bl**4 -7.64507811014497000_mytype*eta_bl**3  &
            +6.94303207046209_mytype*eta_bl**2 -0.00209716712558639_mytype*eta_bl**1


    g_bl_inf=g_bl_inf+(1-exp(-delta_eta/delta_int))*eps_eta
    g_bl_inf=g_bl_inf/1.000546554_mytype
#ifdef DEBG
    if (nrank == 0) write(*,*)'g_bl_inf ', g_bl_inf
#endif

    return
  end subroutine blasius

  !############################################################################
  subroutine postprocess_tbl(ux1,uy1,uz1,pp3,ep1)

    USE MPI
    USE decomp_2d_io
    USE var, only : umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, only : nzmsize,npress
    USE variables
    
    !USE param, only: ivisu, ioutput, itime
    USE param
    USE tools
    
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3
    character(len=30) :: filename
    integer :: out_cntr,inum
    
    if (iin.ne.3) then
      ! Record outflow 
      if (ioutflow.eq.1) then
        out_cntr=mod(itime,ntimesteps)
        if (out_cntr==0) out_cntr=ntimesteps
        call append_outflow(ux1,uy1,uz1,out_cntr)
        if (mod(itime,ntimesteps)==0) then
          call write_outflow(itime/ntimesteps)
        endif
      endif
    endif
    
    ! For statistics !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    out_cntr=mod(itime,ntimesteps)
    if (out_cntr==0) out_cntr=ntimesteps
    call append_outflow_XY(ux1,uy1,uz1,pp3,out_cntr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine postprocess_tbl

  subroutine visu_tbl_init (visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)

    visu_initialised = .true.

  end subroutine visu_tbl_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_tbl
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs TBL-specific visualization
  !!
  !############################################################################
  subroutine visu_tbl(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use var, only: Um_recXY, Vm_recXY, Wm_recXY, UU_recXY, UV_recXY, UW_recXY
    use var, only: VV_recXY, VW_recXY, WW_recXY, Pm_recXY, PN_recXY
    
    use var, only: Um_XY, Vm_XY, Wm_XY, UU_XY, UV_XY, UW_XY
    use var, only: VV_XY, VW_XY, WW_XY, Pm_XY, PN_XY
    use var, only: UmA_XY, VmA_XY, WmA_XY, UUA_XY, UVA_XY, UWA_XY
    use var, only: VVA_XY, VWA_XY, WWA_XY, PmA_XY, PNA_XY
    USE var, only: dp
    USE var, only: num_now,num_rec
    
    use var, only : pp1
    use var, only : pp2, ppi2, dip2, ph2
    use var, only : ppi3, dip3, ph3
    use var, only : npress
    use tools, only : rescale_pressure
    
    USE param, only: ioutput
    
    use visu, only : write_field
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    
    real(kind=dp), allocatable :: buff_r(:)
    integer,dimension(3) :: zstart_recv,zend_recv,zsize_recv
    integer :: i,j
    integer :: iermpi
    character(len=30) :: filename
    
    integer, intent(in) :: num
    
    integer :: ios_rec, pos_rec
    logical :: exists_rec
    character(len=256) :: line_rec, key_rec, value_rec
    
    !integer :: dim11,dim12,dim13,dim21,dim22,dim23

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    
    !dim11 = size(ta3,1)
    !dim12 = size(ta3,2)
    !dim13 = size(ta3,3)
    !dim21 = size(td2,1)
    !dim22 = size(td2,2)
    !dim23 = size(td2,3)
    !write(*,*) 'Rank-',nrank,' ta3',dim11,dim12,dim13,' td2',dim21,dim22,dim23
    
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
    
    
    !if (num.eq.1) then
    !   inquire(file='restart.info', exist=exists_rec)
    !   if (.not. exists_rec) then
    !       num_rec = 0
    !       if (nrank.eq.0) print *, 'The restart.info not found, let num_rec = ', num_rec
    !   else
    !       open(unit=10, file='restart.info', status='old', action='read', iostat=ios_rec)

    !       if (ios_rec == 0) then
    !           do
    !               read(10, '(A)', iostat=ios_rec) line_rec
    !               if (ios_rec /= 0) exit

    !               line_rec = adjustl(line_rec)
    !               if (len_trim(line_rec) == 0) cycle
    !              if (line_rec(1:1) == '!') cycle

    !              pos_rec = index(line_rec, '=')
    !              if (pos_rec <= 0) cycle

     !              key_rec   = adjustl(trim(line_rec(:pos_rec-1)))
     !              value_rec = adjustl(trim(line_rec(pos_rec+1:)))

    !               if (trim(key_rec) == 'num') then
    !                   read(value_rec, *, iostat=ios_rec) num_rec
    !                   if (ios_rec /= 0) num_rec = 0
    !                  exit
    !               end if
    !           end do

    !           close(10)
    !       else
    !           num_rec = 0
    !       end if
    !   end if
    !endif
    
    num_now = itime/ioutput + num_rec
    if (nrank.eq.0) print *, 'The recorded number of output data found. num_rec = ', num_now

    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)
    call write_field(di1, ".", "vort",num_now, flush=.true.) ! Reusing temporary array, force flush
    
    
    !!! Um_recXY, Vm_recXY, Wm_recXY
    ta1 = Um_recXY
    tb1 = Vm_recXY
    tc1 = Wm_recXY
    
    call transpose_x_to_y(ta1,ta2)
    call transpose_x_to_y(tb1,tb2)
    call transpose_x_to_y(tc1,tc2)
    
    call transpose_y_to_z(ta2,ta3)
    call transpose_y_to_z(tb2,tb3)
    call transpose_y_to_z(tc2,tc3)
    
    !do i=1,zsize(1)
    !    do j=1,zsize(2)
    !        Um_XY(i,j) = sum(ta3(i,j,:))/real(zsize(3))
    !        Vm_XY(i,j) = sum(tb3(i,j,:))/real(zsize(3))
    !        Wm_XY(i,j) = sum(tc3(i,j,:))/real(zsize(3))
    !    enddo
    !enddo
    
    Um_XY = sum(ta3,dim = 3)/real(zsize(3))/real(ioutput)
    Vm_XY = sum(tb3,dim = 3)/real(zsize(3))/real(ioutput)
    Wm_XY = sum(tc3,dim = 3)/real(zsize(3))/real(ioutput)
    
    !!! UU_recXY, UV_recXY, UW_recXY
    ta1 = UU_recXY
    tb1 = UV_recXY
    tc1 = UW_recXY
    
    call transpose_x_to_y(ta1,ta2)
    call transpose_x_to_y(tb1,tb2)
    call transpose_x_to_y(tc1,tc2)
    
    call transpose_y_to_z(ta2,ta3)
    call transpose_y_to_z(tb2,tb3)
    call transpose_y_to_z(tc2,tc3)
    
    !do i=1,zsize(1)
    !    do j=1,zsize(2)
    !        UU_XY(i,j) = sum(ta3(i,j,:))/real(zsize(3))
    !        UV_XY(i,j) = sum(tb3(i,j,:))/real(zsize(3))
    !        UW_XY(i,j) = sum(tc3(i,j,:))/real(zsize(3))
    !    enddo
    !enddo
    
    UU_XY = sum(ta3,dim = 3)/real(zsize(3))/real(ioutput)
    UV_XY = sum(tb3,dim = 3)/real(zsize(3))/real(ioutput)
    UW_XY = sum(tc3,dim = 3)/real(zsize(3))/real(ioutput)
    
    !!! VV_recXY, VW_recXY, WW_recXY
    ta1 = VV_recXY
    tb1 = VW_recXY
    tc1 = WW_recXY
    
    call transpose_x_to_y(ta1,ta2)
    call transpose_x_to_y(tb1,tb2)
    call transpose_x_to_y(tc1,tc2)
    
    call transpose_y_to_z(ta2,ta3)
    call transpose_y_to_z(tb2,tb3)
    call transpose_y_to_z(tc2,tc3)
    
    !do i=1,zsize(1)
    !    do j=1,zsize(2)
    !        VV_XY(i,j) = sum(ta3(i,j,:))/real(zsize(3))
    !        VW_XY(i,j) = sum(tb3(i,j,:))/real(zsize(3))
    !        WW_XY(i,j) = sum(tc3(i,j,:))/real(zsize(3))
    !    enddo
    !enddo
    
    VV_XY = sum(ta3,dim = 3)/real(zsize(3))/real(ioutput)
    VW_XY = sum(tb3,dim = 3)/real(zsize(3))/real(ioutput)
    WW_XY = sum(tc3,dim = 3)/real(zsize(3))/real(ioutput)
    
    !!! Pm_recXY, PN_recXY
    ta1 = Pm_recXY
    tb1 = PN_recXY
    
    call transpose_x_to_y(ta1,ta2)
    call transpose_x_to_y(tb1,tb2)
    
    call transpose_y_to_z(ta2,ta3)
    call transpose_y_to_z(tb2,tb3)
    
    !do i=1,zsize(1)
    !    do j=1,zsize(2)
    !        Pm_XY(i,j) = sum(ta3(i,j,:))/real(zsize(3))
    !        PN_XY(i,j) = sum(tb3(i,j,:))/real(zsize(3))
    !    enddo
    !enddo
    
    Pm_XY = sum(ta3,dim = 3)/real(zsize(3))/real(ioutput)
    PN_XY = sum(tb3,dim = 3)/real(zsize(3))/real(ioutput)
    
    if (nrank .ne. 0) then
        
        call MPI_Send(zstart, 3 , MPI_INTEGER, 0, 112, MPI_COMM_WORLD, iermpi)
        call MPI_Send(zend, 3 ,  MPI_INTEGER, 0, 113, MPI_COMM_WORLD, iermpi)
        
        call MPI_Send(Um_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 101, MPI_COMM_WORLD, iermpi)
        call MPI_Send(Vm_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 102, MPI_COMM_WORLD, iermpi)
        call MPI_Send(Wm_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 103, MPI_COMM_WORLD, iermpi)
        call MPI_Send(UU_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 104, MPI_COMM_WORLD, iermpi)
        call MPI_Send(UV_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 105, MPI_COMM_WORLD, iermpi)
        call MPI_Send(UW_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 106, MPI_COMM_WORLD, iermpi)
        call MPI_Send(VV_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 107, MPI_COMM_WORLD, iermpi)
        call MPI_Send(VW_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 108, MPI_COMM_WORLD, iermpi)
        call MPI_Send(WW_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 109, MPI_COMM_WORLD, iermpi)
        call MPI_Send(Pm_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 110, MPI_COMM_WORLD, iermpi)
        call MPI_Send(PN_XY, (zsize(1) * zsize(2)) ,  &
             & MPI_DOUBLE_PRECISION, 0, 111, MPI_COMM_WORLD, iermpi)
    else
        UmA_XY(zstart(1):zend(1),zstart(2):zend(2)) = Um_XY
        VmA_XY(zstart(1):zend(1),zstart(2):zend(2)) = Vm_XY
        WmA_XY(zstart(1):zend(1),zstart(2):zend(2)) = Wm_XY
        UUA_XY(zstart(1):zend(1),zstart(2):zend(2)) = UU_XY
        UVA_XY(zstart(1):zend(1),zstart(2):zend(2)) = UV_XY
        UWA_XY(zstart(1):zend(1),zstart(2):zend(2)) = UW_XY
        VVA_XY(zstart(1):zend(1),zstart(2):zend(2)) = VV_XY
        VWA_XY(zstart(1):zend(1),zstart(2):zend(2)) = VW_XY
        WWA_XY(zstart(1):zend(1),zstart(2):zend(2)) = WW_XY
        PmA_XY(zstart(1):zend(1),zstart(2):zend(2)) = Pm_XY
        PNA_XY(zstart(1):zend(1),zstart(2):zend(2)) = PN_XY
    endif
    
    if (nrank.eq.0) then
        
        do i = 1,nproc-1
            
            call MPI_Recv(zstart_recv, 3 , MPI_INTEGER, i, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            call MPI_Recv(zend_recv, 3 ,  MPI_INTEGER, i, 113, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            
            zsize_recv = zend_recv - zstart_recv + 1
            
            allocate(buff_r(zsize_recv(1) * zsize_recv(2)))
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            UmA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            VmA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            WmA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            UUA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 105, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            UVA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 106, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            UWA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 107, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            VVA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 108, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            VWA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 109, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            WWA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            PmA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            call MPI_Recv(buff_r, (zsize_recv(1) * zsize_recv(2)) ,  &
                 & MPI_DOUBLE_PRECISION, i, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE, iermpi)
            PNA_XY(zstart_recv(1):zend_recv(1),zstart_recv(2):zend_recv(2)) = &
            & reshape(buff_r,[zsize_recv(1) , zsize_recv(2)])
            
            deallocate(buff_r)
            
        enddo
    
        write(filename, '("data/Um-",I0,".bin")') num_now
        open(unit=201, file=filename, access='stream', form='unformatted', action='write', status='replace')
  
        write(filename, '("data/Vm-",I0,".bin")') num_now
        open(unit=202, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/Wm-",I0,".bin")') num_now
        open(unit=203, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/UU-",I0,".bin")') num_now
        open(unit=204, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/UV-",I0,".bin")') num_now
        open(unit=205, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/UW-",I0,".bin")') num_now
        open(unit=206, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/VV-",I0,".bin")') num_now
        open(unit=207, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/VW-",I0,".bin")') num_now
        open(unit=208, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/WW-",I0,".bin")') num_now
        open(unit=209, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/Pm-",I0,".bin")') num_now
        open(unit=210, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(filename, '("data/PN-",I0,".bin")') num_now
        open(unit=211, file=filename, access='stream', form='unformatted', action='write', status='replace')
    
        write(201) UmA_XY
        write(202) VmA_XY
        write(203) WmA_XY
        write(204) UUA_XY
        write(205) UVA_XY
        write(206) UWA_XY
        write(207) VVA_XY
        write(208) VWA_XY
        write(209) WWA_XY
        write(210) PmA_XY
        write(211) PNA_XY
    
        close(201)
        close(202)
        close(203)
        close(204)
        close(205)
        close(206)
        close(207)
        close(208)
        close(209)
        close(210)
        close(211)
    
    endif
    
    Um_recXY = zero
    Vm_recXY = zero
    Wm_recXY = zero
    UU_recXY = zero
    UV_recXY = zero
    UW_recXY = zero
    VV_recXY = zero
    VW_recXY = zero
    WW_recXY = zero
    Pm_recXY = zero
    PN_recXY = zero

  end subroutine visu_tbl
  
  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: momentum_forcing_sponge
  !!      AUTHOR: Anjia Ying
  !!
  !############################################################################
  subroutine momentum_forcing_sponge(dux1, duy1, duz1, ux1, uy1, uz1)
  
    USE var, only : ta1,tb1,tc1
    USE var, only : ta2,tb2,tc2
    USE var, only : ta3,tb3,tc3
    USE var, only : td3,te3,tf3
    USE var, only: yp,xp,xlx
    
    USE var, only: Umx_XY,Vmx_XY,ind_xs,sigx_s
    
    use variables, only : x_sponge,y_sponge, in_xsponge, in_ysponge

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
    
    real(mytype) :: As_s, Ns_s, Lx_s, xs_s,dx_s, Ly_s, ys_s,dy_s
    integer :: i,j,k,pd_xfirst
    
    As_s = 4._mytype
    Ns_s = 3._mytype
    
    Lx_s = xlx
    xs_s = x_sponge
    
    Ly_s = yly
    ys_s = y_sponge
    
    if (itime .eq. 1) then
       pd_xfirst = 0
       if (x_sponge .ge. 0) then
          do i = 1,xsize(1)
           if(xp(xstart(1)+i-1).ge.xs_s) then
            
               dx_s = xp(xstart(1)+i-1) - xs_s
               sigx_s(i) = As_s*(Ns_s+1.)*(Ns_s+2.)*(dx_s**Ns_s*(Lx_s-xp(xstart(1)+i-1)) / (Lx_s-xs_s)**(Ns_s+2.))
            
               if(pd_xfirst.eq.0) ind_xs = i
               pd_xfirst = 1
            
           endif
          enddo
       endif
    endif
    
    do i = 1,xsize(3)
       ta1(ind_xs:xsize(1),:,i) = ux1(ind_xs:xsize(1),:,i) - Umx_XY(ind_xs:xsize(1),xstart(2):xend(2))
       tb1(ind_xs:xsize(1),:,i) = uy1(ind_xs:xsize(1),:,i) - Vmx_XY(ind_xs:xsize(1),xstart(2):xend(2))
    enddo
    
    !if (nrank.eq.0) then
    !do i = 1,xsize(1)
    !    write(*,*) xp(i),sigx_s(i),ind_xs
    !enddo
    !endif
    
    if (x_sponge .ge. 0) then
       do i = ind_xs,xsize(1)
          dux1(i,:,:,1) = dux1(i,:,:,1) - in_xsponge*sigx_s(i) * ta1(i,:,:)
          duy1(i,:,:,1) = duy1(i,:,:,1) - in_xsponge*sigx_s(i) * tb1(i,:,:)
          duz1(i,:,:,1) = duz1(i,:,:,1) - in_xsponge*sigx_s(i) * uz1(i,:,:)
       enddo
    endif

  end subroutine momentum_forcing_sponge


  subroutine set_lund_maps(gama_in, gama8_in, delta_i, alp, bet)
    ! Inner/outer interpolation maps and the Lund blending weight.
    ! The search is monotone, so each station starts from the previous index.
    use var, only: yp, yly, ysize
    use var, only: ind_Inner1, ind_Inner2, ind_Outer1, ind_Outer2
    use var, only: wt_Inner1, wt_Inner2, wt_Outer1, wt_Outer2, W_yita, W_o
    real(mytype), intent(in) :: gama_in, gama8_in, delta_i, alp, bet
    integer :: j, jj
    real(mytype) :: yR_inner, yR_outer, yita

    do j = 1, ysize(2)
       if (j > 1) then
          ind_Inner2(j) = ind_Inner2(j-1)
          ind_Outer2(j) = ind_Outer2(j-1)
       else
          ind_Inner2(j) = 1
          ind_Outer2(j) = 1
       endif
       yR_inner = yp(j)*gama_in
       yR_outer = yp(j)*gama8_in

       if (yR_inner <= yly) then
          do jj = ind_Inner2(j), ysize(2)
             if (yp(jj) >= yR_inner) then
                ind_Inner2(j) = jj
                exit
             endif
          enddo
          if (ind_Inner2(j) /= 1) then
             ind_Inner1(j) = ind_Inner2(j)-1
             wt_Inner2(j) = (yR_inner-yp(ind_Inner1(j)))/(yp(ind_Inner2(j))-yp(ind_Inner1(j)))
             wt_Inner1(j) = 1.0-wt_Inner2(j)
          else
             ind_Inner1(j) = 1
             wt_Inner1(j) = 0.5
             wt_Inner2(j) = 0.5
          endif
       else
          ind_Inner1(j) = ysize(2)
          ind_Inner2(j) = ysize(2)
          wt_Inner1(j) = 0.5
          wt_Inner2(j) = 0.5
       endif

       if (yR_outer <= yly) then
          do jj = ind_Outer2(j), ysize(2)
             if (yp(jj) >= yR_outer) then
                ind_Outer2(j) = jj
                exit
             endif
          enddo
          if (ind_Outer2(j) /= 1) then
             ind_Outer1(j) = ind_Outer2(j)-1
             wt_Outer2(j) = (yR_outer-yp(ind_Outer1(j)))/(yp(ind_Outer2(j))-yp(ind_Outer1(j)))
             wt_Outer1(j) = 1.0-wt_Outer2(j)
          else
             ind_Outer1(j) = 1
             wt_Outer1(j) = 0.5
             wt_Outer2(j) = 0.5
          endif
       else
          ind_Outer1(j) = ysize(2)
          ind_Outer2(j) = ysize(2)
          wt_Outer1(j) = 0.5
          wt_Outer2(j) = 0.5
       endif

       yita = yp(j)/delta_i
       W_yita(j) = 0.5 * ( 1.0 + tanh( ( alp*(yita-bet) ) / ( (1.0-2.0*bet)*yita+bet ) ) / tanh(alp) )
       W_o(j) = max(1d-12, min((1.0_mytype-1d-12), W_yita(j)))
    enddo
  end subroutine set_lund_maps

  subroutine blend_lund_fluctuations(u_src, v_src, w_src, u_out, v_out, w_out)
    ! Lund inner/outer blend of a recycling-plane fluctuation field.
    use var, only: ysize, bdamp, gamaU, W_yita
    use var, only: ind_Inner1, ind_Inner2, ind_Outer1, ind_Outer2
    use var, only: wt_Inner1, wt_Inner2, wt_Outer1, wt_Outer2
    real(mytype), intent(in) :: u_src(:,:), v_src(:,:), w_src(:,:)
    real(mytype), intent(out) :: u_out(:,:), v_out(:,:), w_out(:,:)
    integer :: j
    real(mytype), dimension(size(u_src,2)) :: u_in, u_out_sc, v_in, v_out_sc, w_in, w_out_sc

    do j = 1, ysize(2)
       u_in = wt_Inner1(j)*u_src(ind_Inner1(j),:) + wt_Inner2(j)*u_src(ind_Inner2(j),:)
       u_out_sc = wt_Outer1(j)*u_src(ind_Outer1(j),:) + wt_Outer2(j)*u_src(ind_Outer2(j),:)
       v_in = wt_Inner1(j)*v_src(ind_Inner1(j),:) + wt_Inner2(j)*v_src(ind_Inner2(j),:)
       v_out_sc = wt_Outer1(j)*v_src(ind_Outer1(j),:) + wt_Outer2(j)*v_src(ind_Outer2(j),:)
       w_in = wt_Inner1(j)*w_src(ind_Inner1(j),:) + wt_Inner2(j)*w_src(ind_Inner2(j),:)
       w_out_sc = wt_Outer1(j)*w_src(ind_Outer1(j),:) + wt_Outer2(j)*w_src(ind_Outer2(j),:)
       u_out(j,:) = bdamp(j)*((1.0-W_yita(j))*u_in + W_yita(j)*u_out_sc)/gamaU
       v_out(j,:) = bdamp(j)*((1.0-W_yita(j))*v_in + W_yita(j)*v_out_sc)/gamaU
       w_out(j,:) = bdamp(j)*((1.0-W_yita(j))*w_in + W_yita(j)*w_out_sc)/gamaU
    enddo
  end subroutine blend_lund_fluctuations

  subroutine half_span_shift(uyz, vyz, wyz, ua, va, wa, us, vs, ws, uo, vo, wo, n_shift)
    ! Shift a y-pencil plane by half the span and return it in uo/vo/wo.
    real(mytype), intent(in) :: uyz(:,:,:), vyz(:,:,:), wyz(:,:,:)
    real(mytype), intent(inout) :: ua(:,:,:), va(:,:,:), wa(:,:,:)
    real(mytype), intent(inout) :: us(:,:,:), vs(:,:,:), ws(:,:,:)
    real(mytype), intent(inout) :: uo(:,:,:), vo(:,:,:), wo(:,:,:)
    integer, intent(in) :: n_shift
    integer :: k, k_shift

    call transpose_y_to_z(uyz, ua)
    call transpose_y_to_z(vyz, va)
    call transpose_y_to_z(wyz, wa)
    do k = 1, zsize(3)
       k_shift = mod(k+n_shift-1, zsize(3)) + 1
       us(1,:,k) = ua(1,:,k_shift)
       vs(1,:,k) = va(1,:,k_shift)
       ws(1,:,k) = wa(1,:,k_shift)
    enddo
    call transpose_z_to_y(us, uo)
    call transpose_z_to_y(vs, vo)
    call transpose_z_to_y(ws, wo)
  end subroutine half_span_shift

end module tbl
