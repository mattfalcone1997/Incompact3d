!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module stats
  use param, only : mytype
  use decomp_2d, only : DECOMP_INFO
#ifdef HAVE_FFTW
  include "fftw3.f"
#endif
  implicit none

  character(len=*), parameter :: io_statistics = "statistics-io", &
       stat_dir = "statistics"

  integer :: stats_time
  real(mytype), dimension(:,:,:), allocatable ::  uvwp_mean
  real(mytype), dimension(:,:,:), allocatable ::  uu_mean
  real(mytype), dimension(:,:,:), allocatable ::  pu_mean
  real(mytype), dimension(:,:,:), allocatable ::  uuu_mean
  real(mytype), dimension(:,:,:), allocatable ::  pdudx_mean
  real(mytype), dimension(:,:,:), allocatable ::  dudx_mean
  real(mytype), dimension(:,:,:), allocatable ::  dudxdudx_mean

  real(mytype), dimension(:,:,:), allocatable ::  pdvdy_q_mean

  real(mytype), dimension(:,:,:), allocatable :: uv_quadrant_mean
  real(mytype), dimension(:,:,:), allocatable :: uuuu_mean

  real(mytype), dimension(:,:,:), allocatable :: spectra_x_mean
  real(mytype), dimension(:,:,:,:), allocatable :: spectra_z_mean

  type(DECOMP_INFO) :: uvwp_info, uu_info, pu_info
  type(DECOMP_INFO) :: uuu_info, pdudx_info, dudx_info
  type(DECOMP_INFO) :: dudxdudx_info, pdvdy_q_info
  type(DECOMP_INFO) :: uv_quadrant_info, uuuu_info
  real(mytype), dimension(:,:,:), allocatable :: lambda2mean, lambda22mean

  
  real(mytype), dimension(:,:), allocatable :: coefy
  real(mytype) :: coefy_bc1(3), coefy_bcn(3)
  real(mytype) :: coefx_bc1(3), coefx_bcn(3)
  real(mytype) :: coefz_bc1(3), coefz_bcn(3)

  real(mytype) :: coefx_b, coefx_f
  real(mytype) :: coefz_b, coefz_f
  real(mytype), dimension(:), allocatable :: h_quads
  integer :: spectra_level, nspectra, spectra_nlocs, spectra_size
  real(mytype), dimension(:), allocatable :: spectra_xlocs
  integer, allocatable, dimension(:) :: spectra_x_indices
  type(DECOMP_INFO) :: spectra_x_info, spectra_z_info
  type(DECOMP_INFO) :: dstat_plane


  private
  public overall_statistic, h_quads, spectra_level

contains

  subroutine init_statistic_adios2

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    use param, only : istatbudget, istatpstrain, istatlambda2, istatquadrant
    use param, only : istatflatness
    use var, only : numscalar
    
    implicit none

    integer :: ierror
    character(len=30) :: varname
    integer :: is
    logical, save :: initialised = .false.

    if (.not. initialised) then
       call decomp_2d_init_io(io_statistics)
    
       call decomp_2d_register_variable(io_statistics, "uvwp_mean", 3, 0, 0, mytype,opt_decomp=uvwp_info)
       call decomp_2d_register_variable(io_statistics, "uu_mean", 3, 0, 0, mytype,opt_decomp=uu_info)
              
       if (istatbudget) then
        call decomp_2d_register_variable(io_statistics, "pu_mean", 3, 0, 0, mytype,opt_decomp=pu_info)
       call decomp_2d_register_variable(io_statistics, "uuu_mean", 3, 0, 0, mytype,opt_decomp=uuu_info)
       call decomp_2d_register_variable(io_statistics, "pdudx_mean", 3, 0, 0, mytype,opt_decomp=pdudx_info)
       call decomp_2d_register_variable(io_statistics, "dudx_mean", 3, 0, 0, mytype,opt_decomp=pdudx_info)
       call decomp_2d_register_variable(io_statistics, "dudxdudx_mean", 3, 0, 0, mytype,opt_decomp=dudxdudx_info)

      endif

       if (istatpstrain) then
        call decomp_2d_register_variable(io_statistics, "pdvdy_q_mean", 3, 0, 0, mytype, opt_decomp=pdvdy_q_info)
       endif

       if (istatlambda2) then
        call decomp_2d_register_variable(io_statistics, "lambda2mean", 3, 1, 3, mytype, opt_nplanes=1)
        call decomp_2d_register_variable(io_statistics, "lambda22mean", 3, 1, 3, mytype, opt_nplanes=1)
       endif

       if (istatflatness) then
        call decomp_2d_register_variable(io_statistics, "uuuu_mean", 3, 0, 0, mytype, opt_decomp=uuuu_info)
       endif
       if (istatquadrant) then
        call decomp_2d_register_variable(io_statistics, "uv_quadrant_mean", 3, 0, 0, mytype,opt_decomp=uv_quadrant_info)
       endif

       do is=1, numscalar
          write(varname,"('phi',I2.2)") is
          call decomp_2d_register_variable(io_statistics, varname, 3, 1, 3, mytype, opt_nplanes=1)
          write(varname,"('phiphi',I2.2)") is
          call decomp_2d_register_variable(io_statistics, varname, 3, 1, 3, mytype, opt_nplanes=1)
       enddo

       initialised = .true.
    endif
       
  end subroutine init_statistic_adios2
  
  !
  ! Initialize to zero all statistics
  !
  subroutine init_statistic

    use param, only : zero, iscalar, istatbudget, istatpstrain
    use param, only : istatquadrant, istatlambda2, nquads, istatspectra
    use param, only : istatflatness
    use variables, only : nx, ny, nz
    use decomp_2d, only : zsize, decomp_info_init
    use MPI

    implicit none
    integer :: code

    allocate(uvwp_mean(zsize(1),zsize(2),4))
    allocate(uu_mean(zsize(1),zsize(2),6))

    uvwp_mean = zero
    uu_mean = zero

    if (istatbudget) then
      allocate(pu_mean(zsize(1),zsize(2),3))
      allocate(uuu_mean(zsize(1),zsize(2),10))
      allocate(pdudx_mean(zsize(1),zsize(2),3))
      allocate(dudx_mean(zsize(1),zsize(2),9))
      allocate(dudxdudx_mean(zsize(1),zsize(2),12))

      pu_mean = zero
      uuu_mean = zero
      pdudx_mean = zero
      dudx_mean = zero
      dudxdudx_mean = zero
      
    endif

    if (istatflatness) then
      allocate(uuuu_mean(zsize(1),zsize(2),3))
      uuuu_mean = zero
    endif

    if (istatpstrain) then
      if (.not.istatbudget) then
        write(*,*) "If p-strain conditional averaging used, istatbudget must be T"
        call MPI_Abort(MPI_COMM_WORLD,1,code)
      endif

      allocate(pdvdy_q_mean(zsize(1),zsize(2),4))

      pdvdy_q_mean = zero
    endif

    if (istatquadrant) then
      allocate(uv_quadrant_mean(zsize(1),zsize(2),4*nquads))

      uv_quadrant_mean = zero
    endif
    if (istatlambda2) then
#ifndef HAVE_LAPACK 
      write(*,*) "Cannot calaculate mean lambda2 without LAPACK"
      call MPI_Abort(MPI_COMM_WORLD,1,code)
#endif  

      allocate(lambda2mean(zsize(1),zsize(2),1))
      allocate(lambda22mean(zsize(1),zsize(2),1))

      lambda2mean = zero
      lambda22mean = zero
    endif

    if (istatspectra) then
#ifndef HAVE_FFTW
      write(*,*) "Cannot use spectra unless built with FFTw"
      call MPI_Abort(MPI_COMM_WORLD,1,code)
#endif      
      call init_spectra
    endif
  
    call decomp_info_init(nx, ny, 4, uvwp_info)
    call decomp_info_init(nx, ny, 6, uu_info)
    call decomp_info_init(nx, ny, 3, pu_info)
    call decomp_info_init(nx, ny, 10, uuu_info)
    call decomp_info_init(nx, ny, 3, pdudx_info)
    call decomp_info_init(nx, ny, 9, dudx_info)
    call decomp_info_init(nx, ny, 12, dudxdudx_info)

    call decomp_info_init(nx, ny, 4*nquads, uv_quadrant_info)

    call decomp_info_init(nx, ny, 4, pdvdy_q_info)
    call decomp_info_init(nx, ny, 3, uuuu_info)

    call init_statistic_adios2

    call grad_init

  end subroutine init_statistic

  subroutine init_spectra
    use param
    use decomp_2d
    use dbg_schemes, only : abs_prec
    real(mytype), dimension(:,:,:), allocatable :: spectra_x_in, spectra_z_in
    real(mytype), dimension(:,:,:), allocatable :: spectra_x_out, spectra_z_out
    real(mytype) :: x
#ifdef HAVE_FFTW
    integer :: plane_type = FFTW_MEASURE

    if (spectra_level == 1) then
      nspectra = 4
    else if (spectra_level == 2) then
      nspectra = 5
    else
      write(*,*) "Invalid spectra level"
    call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif

    if (spectra_nlocs > 0) then
      allocate(spectra_x_indices(spectra_nlocs))
      do j = 1, spectra_nlocs
        do i = 1, nx
          x = real(i-1,kind=mytype)*dx
          if (abs_prec(spectra_xlocs(i)-x) <= half*dx) then
            spectra_x_indices(j) = i
            exit
          endif
        enddo
      enddo
    else
      spectra_nlocs = 1
    endif

    allocate(spectra_x_mean(nx/2, zsize(2),nspectra))
    allocate(spectra_z_mean(zsize(2),nz/2,spectra_nlocs,nspectra))

    allocate(spectra_x_out(nx/2,xsize(2),xsize(3)))
    allocate(spectra_z_out(zsize(1),zsize(2),nz/2))

    allocate(spectra_x_in(nx,xsize(2),xsize(3)))
    allocate(spectra_z_in(zsize(1),zsize(2),nz))

#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft_r2c(plan_x, 1, xsize(1), &
         xsize(2)*xsize(3), spectra_x_in, xsize(1), 1, &
         xsize(1), spectra_x_out, spectra_x_info%xsz(1), 1, spectra_x_info%xsz(1), &
         plan_type)
#else
    call sfftw_plan_many_dft_r2c(plan_x, 1, xsize(1), &
         xsize(2)*xsize(3), spectra_x_in, xsize(1), 1, &
         xsize(1), spectra_x_out, spectra_x_info%xsz(1), 1, spectra_x_info%xsz(1), &
         plan_type)
#endif

#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft_r2c(plan1, 1, zsize(3), &
         zsize(1)*zsize(2), spectra_z_in, zsize(3), &
         zsize(1)*zsize(2), 1, spectra_z_out, nz/2, &
         zsize(1)*zsize(2), 1, plan_type)
#else
    call sfftw_plan_many_dft_r2c(plan1, 1, zsize(3), &
         zsize(1)*zsize(2), spectra_z_in, zsize(3), &
         zsize(1)*zsize(2), 1, spectra_z_out, nz/2, &
         zsize(1)*zsize(2), 1, plan_type)
#endif
   
    deallocate(spectra_x_in,spectra_x_out)
    deallocate(spectra_z_in,spectra_z_out)
#endif    
  end subroutine
  !
  ! Read all statistics if possible
  !
  subroutine restart_statistic

    use param, only : initstat, irestart, ifirst, zero
    use variables, only : nstat
    use var, only : tmean

    implicit none

    ! No reading for statistics when nstat > 1 or no restart
    call init_statistic()

    if (nstat.gt.1 .or. irestart.eq.0) then
       initstat = ifirst
       return
    else
       call init_statistic_adios2
    endif


    ! Temporary array
    tmean = zero

    ! Read all statistics
    call read_or_write_all_stats(.true.)

  end subroutine restart_statistic

  function gen_statname(stat) result(newname)

    implicit none
    
    character(len=*), intent(in) :: stat
    character(len=30) :: newname
    
#ifndef ADIOS2
    write(newname, "(A,'.dat',I7.7)") stat, stats_time
#else
    write(newname, *) stat
#endif
    
  end function gen_statname

  !
  ! Statistics: perform all IO
  !
  subroutine read_or_write_all_stats(flag_read)

    use param, only : iscalar, itime, istatbudget, istatpstrain
    use param, only : initstat2, istatlambda2, istatquadrant, istatflatness
    use variables, only : numscalar
    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_write_mode, decomp_2d_read_mode, &
         decomp_2d_open_io, decomp_2d_close_io, decomp_2d_start_io, decomp_2d_end_io

    implicit none

    ! Argument
    logical, intent(in) :: flag_read

    ! Local variables
    integer :: is, it
    character(len=30) :: filename
    integer :: io_mode
    integer :: ierror
    

    ! File ID to read or write
    if (flag_read) then
        it = itime - 1
    else
        it = itime
     endif
     stats_time = it

    if (nrank==0) then
      print *,'==========================================================='
      if (flag_read) then
        print *,'Reading stat file', stats_time
      else
        print *,'Writing stat file', stats_time
      endif
    endif

    if (flag_read) then
       io_mode = decomp_2d_read_mode
    else
       io_mode = decomp_2d_write_mode
    endif
#ifdef ADIOS2
    call decomp_2d_open_io(io_statistics, stat_dir, io_mode)
    call decomp_2d_start_io(io_statistics, stat_dir)
#endif
    
    call read_or_write_one_stat(flag_read, gen_statname("uvwp_mean"), uvwp_mean,uvwp_info)
    call read_or_write_one_stat(flag_read, gen_statname("uu_mean"), uu_mean,uu_info)
   
    if (istatbudget) then
      call read_or_write_one_stat(flag_read, gen_statname("uuu_mean"), uuu_mean, uuu_info)
      call read_or_write_one_stat(flag_read, gen_statname("pdudx_mean"), pdudx_mean, pdudx_info)
      call read_or_write_one_stat(flag_read, gen_statname("dudx_mean"), dudx_mean,dudx_info)
      call read_or_write_one_stat(flag_read, gen_statname("pu_mean"), pu_mean,pu_info)
      call read_or_write_one_stat(flag_read, gen_statname("dudxdudx_mean"), dudxdudx_mean,dudxdudx_info)

    endif

    if (istatflatness) then
      call read_or_write_one_stat(flag_read, gen_statname("uuuu_mean"), uuuu_mean, uuuu_info)
    endif

    if (istatpstrain .and. itime > initstat2) then
      call read_or_write_one_stat(flag_read, gen_statname("pdvdy_q_mean"), pdvdy_q_mean, pdvdy_q_info)
    endif

    if (istatquadrant .and. itime>initstat2) then
      call read_or_write_one_stat(flag_read, gen_statname("uv_quadrant_mean"), uv_quadrant_mean, uv_quadrant_info)
    endif
    if (istatlambda2) then
      call read_or_write_one_stat(flag_read, gen_statname("lambda2mean"), lambda2mean,dstat_plane)
      call read_or_write_one_stat(flag_read, gen_statname("lambda22mean"), lambda22mean,dstat_plane)
    endif
#ifdef ADIOS2
    call decomp_2d_end_io(io_statistics, stat_dir)
    call decomp_2d_close_io(io_statistics, stat_dir)
#endif
    
    if (nrank==0) then
      if (flag_read) then
        print *,'Read stat done!'
      else
        print *,'Write stat done!'
      endif
      print *,'==========================================================='
    endif

  end subroutine read_or_write_all_stats

  !
  ! Statistics: perform one IO
  !
  subroutine read_or_write_one_stat(flag_read, filename, array,opt_decomp)

    use decomp_2d, only : mytype, xstS, xenS, zsize
    use decomp_2d_io, only : decomp_2d_read_one, decomp_2d_write_one

    implicit none

    ! Arguments
    logical, intent(in) :: flag_read
    character(len=*), intent(in) :: filename
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: array
    type(DECOMP_INFO), intent(in) :: opt_decomp
    integer :: ierror

    if (flag_read) then
       ! There was a check for nvisu = 1 before
       call decomp_2d_read_one(3, array, stat_dir, filename, io_statistics, reduce_prec=.false.,opt_decomp=opt_decomp)
    else
       call decomp_2d_write_one(3, array, stat_dir, filename, 0, io_statistics, reduce_prec=.false.,opt_decomp=opt_decomp)
    endif

  end subroutine read_or_write_one_stat

  subroutine read_write_spectra_x(flag_read, filename)
    use decomp_2d
    use MPI
    logical, intent(in) :: flag_read
    character(len=*), intent(in) :: filename

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: color, key
    integer, dimension(2) :: dims, coords
    logical, dimension(2) :: periods 
    integer :: split_comm_y, split_comm_z, code
    integer :: arr_type, comm_rank, fh
    integer :: i, j, k


    sizes = [xsize(1)/2,ysize(2),nspectra]
    subsizes = [xsize(1)/2,xsize(2),nspectra]
    starts = [0, xstart(2)-1,0]

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, code)
    key = coords(1)
    color = coords(2)
 
    call MPI_Comm_split(DECOMP_2D_COMM_CART_X, color, key, split_comm_y,code)
    call MPI_Comm_split(DECOMP_2D_COMM_CART_X, key, color, split_comm_z,code)

    call MPI_Type_create_subarray(3,sizes,subsizes,&
                                   starts,MPI_ORDER_FORTRAN,&
                                   complex_type,arr_type,code)
    call MPI_Type_commit(arr_type,code)
    call MPI_Comm_rank(split_comm_z,comm_rank,code)

    if (flag_Read) then
      call MPI_File_open(split_comm_y,filename,MPI_MODE_RDONLY,&
                          MPI_INFO_NULL,fh, code)

      call MPI_File_read_all(fh,spectra_x_mean,1,arr_type,&
                              MPI_STATUS_IGNORE,code)                          
      call MPI_File_close(fh, code)
    else
      if (comm_rank == 0) then
        call MPI_File_open(split_comm_y,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,&
                          MPI_INFO_NULL,fh, code)
        
        call MPI_File_write_all(fh,spectra_x_mean,1,arr_type,&
                                MPI_STATUS_IGNORE,code)                          
        call MPI_File_close(fh, code)
  
      endif
    endif
    call MPI_Type_free(arr_type,code)

  end subroutine

  subroutine read_write_spectra_z(flag_read, filename)
    logical, intent(in) :: flag_read
    character(len=*), intent(in) :: filename
  end subroutine

  !
  ! Statistics : Intialize, update and perform IO
  !
  subroutine overall_statistic(ux1,uy1,uz1,phi1,pp3,ep1)

    use param
    use variables
    use decomp_2d
    use decomp_2d_io
    use tools, only : rescale_pressure
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use ibm_param, only : ubcx, ubcy, ubcz
    use var, only : ux2, ux3, uy2, uy3, uz2, uz3
    use var, only : ta2
    use var, only : nxmsize, nymsize, nzmsize
    use var, only : ppi3, dip3
    use var, only : pp2, ppi2, dip2
    use var, only : pp1, ta1, di1
    implicit none

    !! Inputs
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar),intent(in) :: phi1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

    !! Locals

    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: dudy, dvdy, dwdy
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: dudx, dvdx, dwdx

    integer :: is
    character(len=30) :: filename

    if (itime.lt.initstat) then
       return
    elseif (itime.eq.initstat) then
       call init_statistic()

       if (istatout < 1) istatout = icheckpoint
    elseif (itime.eq.ifirst) then
       if (itempaccel.eq.1) then
          call init_statistic()
       else
          call restart_statistic()
       endif
       if (istatout < 1) istatout = icheckpoint
    endif
    
    if (mod(itime,istatcalc) /= 0) return
    if (itempaccel == 1 .and. mod(itime,istatout)/=0) return 

    !! Mean pressure
    !WORK Z-PENCILS
    call interzpv(ppi3,pp3(:,:,:,1),dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(ppi3,pp2,ph3) !nxm nym nz
    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call interxpv(td1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    ! Convert to physical pressure
    call rescale_pressure(td1)

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)
    call transpose_x_to_y(td1,td2)

    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)
    call transpose_y_to_z(td2,td3)

    call update_average_scalar(uvwp_mean(:,:,4), td3, ep1)

    !! Mean velocity
    call update_average_vector(uvwp_mean(:,:,1),uvwp_mean(:,:,2), uvwp_mean(:,:,3), &
                               ux3, uy3, uz3, td3)

    !! Second-order velocity moments
    call update_variance_vector(uu_mean(:,:,1),uu_mean(:,:,2), uu_mean(:,:,3),&
                                uu_mean(:,:,4), uu_mean(:,:,5), uu_mean(:,:,6), &
                                ux3, uy3, uz3, td3)

    if (istatbudget) then
      call grad_x(ux1, ta1)
      call grad_x(uy1, tb1)
      call grad_x(uz1, tc1)
      
      !y-derivatives
      call grad_y(ux2,ta2)
      call grad_y(uy2,tb2)
      call grad_y(uz2,tc2)

      !!z-derivatives

      call grad_z(ux3,ta3)
      call grad_z(uy3,tb3)
      call grad_z(uz3,tc3)

      call transpose_y_to_z(ta2,dudy)
      call transpose_y_to_z(tb2,dvdy)
      call transpose_y_to_z(tc2,dwdy)

      call transpose_x_to_y(ta1,ta2)
      call transpose_x_to_y(tb1,tb2)
      call transpose_x_to_y(tc1,tc2)

      call transpose_y_to_z(ta2,dudx)
      call transpose_y_to_z(tb2,dvdx)
      call transpose_y_to_z(tc2,dwdx)

    call update_skewness_tensor(uuu_mean(:,:,1), uuu_mean(:,:,2), uuu_mean(:,:,3),&
                               uuu_mean(:,:,4),uuu_mean(:,:,5),uuu_mean(:,:,6),&
                               uuu_mean(:,:,7),uuu_mean(:,:,8),uuu_mean(:,:,9),&
                               uuu_mean(:,:,10), ux3, uy3, uz3, td3)  

    call update_pvelograd_vector(pdudx_mean(:,:,1), pdudx_mean(:,:,2), pdudx_mean(:,:,3),&
                                dudx, dvdy, tc3, td3)

    call update_pvelo_vector(pu_mean(:,:,1), pu_mean(:,:,2), pu_mean(:,:,3), ux3, uy3, uz3,td3)

    call update_velograd_tensor(dudx_mean, dudx,dudy,ta3,dvdx,dvdy,&
                                tb3,dwdx,dwdy,tc3,td3)

    call update_velograd2_tensor(dudxdudx_mean,&
                                 dudx,dvdx,dwdx,dudy,dvdy,dwdy,td3)
                      
    endif                             

    if (istatflatness) then
      call update_flatness(uuuu_mean,ux3, uy3, uz3,td3)
    endif
    if (istatpstrain .and. itime>initstat2) then
      call update_pstrain_cond_avg(pdvdy_q_mean(:,:,1),pdvdy_q_mean(:,:,2),&
                                  pdvdy_q_mean(:,:,3),pdvdy_q_mean(:,:,4),&
                                  td3,dvdy,uvwp_mean(:,:,4),dudx_mean(:,:,5),td3)
    endif

    if (istatquadrant .and. itime>initstat2) then
      call update_uv_quad_avg(uv_quadrant_mean, uu_mean(:,:,1), uu_mean(:,:,2),&
                              uvwp_mean(:,:,1), uvwp_mean(:,:,2),ux3,uy3, td3)
    endif

    if (istatspectra .and. itime>initstat2) then
      call update_spectra_avg(spectra_x_mean,spectra_z_mean,ux1,uy1,uz1,ux3,uy3,uz3)
    endif
    if (istatlambda2) then
      call update_lambda2_avg(lambda2mean,lambda22mean,&
                                dudx,dudy,ta3,dvdx,dvdy,&
                                tb3,dwdx,dwdy,tc3,td3)
    endif

    ! Write all statistics
    if (mod(itime,istatout)==0) then
       call read_or_write_all_stats(.false.)
    endif

  end subroutine overall_statistic

  !
  ! Basic function, can be applied to arrays
  !
  elemental real(mytype) function one_minus_ep1(var, ep1)

    use decomp_2d, only : mytype
    use param, only : iibm, one

    implicit none

    ! inputs
    real(mytype), intent(in) :: var, ep1

    if (iibm==2) then
      one_minus_ep1 = (one - ep1) * var
    else
      one_minus_ep1 = var
    endif

  end function one_minus_ep1

  !
  ! Update um, the average of ux
  !
  subroutine update_average_scalar(um, ux, ep,mask,istat2)

    use decomp_2d, only : mytype, xsize, zsize
    use param, only : itime, initstat, istatcalc, itempaccel, initstat2
    use var, only : di1, tmean

    implicit none

    ! inputs
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: um
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux
    logical, dimension(zsize(1),zsize(2),zsize(3)), optional, intent(in) :: mask
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ep
    logical, optional, intent(in) :: istat2

    real(mytype), dimension(zsize(1), zsize(2)) :: stat_z
    real(mytype) :: stat_inc
    integer :: i, j
    
    if (present(mask)) then
      stat_z = sum(ux,dim=3,mask=mask) / real(zsize(3),kind=mytype)

    else
      stat_z = sum(ux,dim=3) / real(zsize(3),kind=mytype)
    endif

    if (itempaccel==1) then
      um(:,:) = stat_z(:,:)
      return
    endif

    if (present(istat2)) then
      if (istat2) then
        stat_inc = real((itime-initstat2)/istatcalc+1, kind=mytype)
      else
        stat_inc = real((itime-initstat)/istatcalc+1, kind=mytype)
      endif
    else
      stat_inc = real((itime-initstat)/istatcalc+1, kind=mytype)
    endif

    do j = 1, zsize(2)
      do i = 1, zsize(1)
        um(i,j) = um(i,j) + (stat_z(i,j) - um(i,j))/ stat_inc
      enddo
    enddo

  end subroutine update_average_scalar

  ! subroutine average_z(array, stat_xz,mask)
  !   use decomp_2d
  !   use param, only : zero

  !   real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: array
  !   real(mytype), dimension(zsize(1),zsize(2),1) :: stat_xz
  !   logical, dimension(zsize(1),zsize(2),zsize(3)), optional, intent(in) :: mask
  !   real(mytype) :: avg_z
  !   integer :: i, j, k

    
  !   stat_xz = sum(array,dim=3,mask=mask) / real(zsize(3),kind=mytype)


  !   do i = 1, zsize(1)
  !     do j = 1, zsize(2)
  !       avg_z = sum(array)
  !       stat_xz(i,j,1) = stat_xz(i,j,1) + array(i,j,k)/zsize(3)
  !     enddo
  !   enddo


  !  end subroutine
  !
  ! Update (um, vm, wm), the average of (ux, uy, uz)
  !
  subroutine update_average_vector(um, vm, wm, ux, uy, uz, ep)

    use decomp_2d, only : mytype, xsize, zsize

    implicit none

    ! inputs
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: um, vm, wm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(um, ux, ep)
    call update_average_scalar(vm, uy, ep)
    call update_average_scalar(wm, uz, ep)

  end subroutine update_average_vector

  !
  ! Update (uum, vvm, wwm, uvm, uwm, vwm), the variance of the vector (ux, uy, uz)
  !
  subroutine update_variance_vector(uum, vvm, wwm, uvm, uwm, vwm, ux, uy, uz, ep)

    use decomp_2d, only : mytype, xsize, zsize

    implicit none

    ! inputs
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: uum, vvm, wwm, uvm, uwm, vwm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(uum, ux*ux, ep)
    call update_average_scalar(vvm, uy*uy, ep)
    call update_average_scalar(wwm, uz*uz, ep)
    call update_average_scalar(uvm, ux*uy, ep)
    call update_average_scalar(uwm, ux*uz, ep)
    call update_average_scalar(vwm, uy*uz, ep)

  end subroutine update_variance_vector

  subroutine update_skewness_tensor(uuum,uuvm,uuwm,uvvm,uvwm,&
                                    uwwm,vvvm,vvwm,vwwm,wwwm,&
                                    ux, uy, uz, ep)
  use decomp_2d, only : mytype, xsize, zsize
  real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: uuum,uuvm,uuwm,uvvm,uvwm,&
                                                                 uwwm,vvvm,vvwm,vwwm,wwwm
  real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

  call update_average_scalar(uuum, ux*ux*ux, ep)
  call update_average_scalar(uuvm, ux*ux*uy, ep)
  call update_average_scalar(uuwm, ux*ux*uz, ep)
  call update_average_scalar(uvvm, ux*uy*uy, ep)
  call update_average_scalar(uvwm, ux*uy*uz, ep)
  call update_average_scalar(uwwm, ux*uz*uz, ep)
  call update_average_scalar(vvvm, uy*uy*uy, ep)
  call update_average_scalar(vvwm, uy*uy*uz, ep)
  call update_average_scalar(vwwm, uy*uz*uz, ep)
  call update_average_scalar(wwwm, uz*uz*uz, ep)
                                    
  end subroutine update_skewness_tensor

  subroutine update_pvelograd_vector(pdudxm,pdvdym,pdwdzm,dudx,dvdy,dwdz,ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: pdudxm,pdvdym,pdwdzm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: dudx, dvdy, dwdz, ep

    call update_average_scalar(pdudxm, ep*dudx, ep)
    call update_average_scalar(pdvdym, ep*dvdy, ep)
    call update_average_scalar(pdwdzm, ep*dwdz, ep)
  end subroutine

  subroutine update_pvelo_vector(pum, pvm, pwm, ux, uy, uz, ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: pum, pvm, pwm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(pum, ep*ux, ep)
    call update_average_scalar(pvm, ep*uy, ep)
    call update_average_scalar(pwm, ep*uz, ep)

  end subroutine update_pvelo_vector

  subroutine update_velograd_tensor(dudxm, dudx, dudy, dudz, dvdx,&
                                    dvdy, dvdz, dwdx, dwdy, dwdz, ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2),9), intent(inout) :: dudxm


    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dudy, dudz, ep
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dvdx, dvdy, dvdz   
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dwdx, dwdy, dwdz

    call update_average_scalar(dudxm(:,:,1), dudx, ep)
    call update_average_scalar(dudxm(:,:,2), dudy, ep)
    call update_average_scalar(dudxm(:,:,3), dudz, ep)
    call update_average_scalar(dudxm(:,:,4), dvdx, ep)
    call update_average_scalar(dudxm(:,:,5), dvdy, ep)
    call update_average_scalar(dudxm(:,:,6), dvdz, ep)
    call update_average_scalar(dudxm(:,:,7), dudx, ep)
    call update_average_scalar(dudxm(:,:,8), dwdy, ep)
    call update_average_scalar(dudxm(:,:,9), dwdz, ep)


  end subroutine update_velograd_tensor

  subroutine update_velograd2_tensor(dudxdudxm,&
                                     dudx,dvdx,dwdx,dudy,dvdy,dwdy,ep)
    use decomp_2d, only : mytype, xsize, zsize

    real(mytype), dimension(zsize(1),zsize(2),12), intent(inout) :: dudxdudxm

    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dvdx, dwdx, ep
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudy, dvdy, dwdy

    call update_average_scalar(dudxdudxm(:,:,1), dudx*dudx, ep)
    call update_average_scalar(dudxdudxm(:,:,2), dudx*dvdx, ep)
    call update_average_scalar(dudxdudxm(:,:,3), dudx*dwdx, ep)
    call update_average_scalar(dudxdudxm(:,:,4), dvdx*dvdx, ep)
    call update_average_scalar(dudxdudxm(:,:,5), dvdx*dwdx, ep)
    call update_average_scalar(dudxdudxm(:,:,6), dwdx*dwdx, ep)

    call update_average_scalar(dudxdudxm(:,:,7), dudy*dudy, ep)
    call update_average_scalar(dudxdudxm(:,:,8), dudy*dvdy, ep)
    call update_average_scalar(dudxdudxm(:,:,9), dudy*dwdy, ep)
    call update_average_scalar(dudxdudxm(:,:,10), dvdy*dvdy, ep)
    call update_average_scalar(dudxdudxm(:,:,11), dvdy*dwdy, ep)
    call update_average_scalar(dudxdudxm(:,:,12), dwdy*dwdy, ep)

  end subroutine update_velograd2_tensor
  
  subroutine update_flatness(uuuum, ux, uy, uz, ep)
    use decomp_2d, only : mytype, xsize, zsize
    real(mytype), dimension(zsize(1),zsize(2),3), intent(inout) :: uuuum
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, uz, ep
  
    call update_average_scalar(uuuum(:,:,1), ux*ux*ux*ux, ep)
    call update_average_scalar(uuuum(:,:,2), uy*uy*uy*uy, ep)
    call update_average_scalar(uuuum(:,:,3), uz*uz*uz*uz, ep)

  end subroutine update_flatness

  subroutine update_uv_quad_avg(uvqm, uum, vvm, um, vm, ux, uy, ep)
    use dbg_schemes, only : sqrt_prec, abs_prec
    use decomp_2d, only : mytype, xsize, zsize
    use param, only: nquads, zero
    use MPI
    real(mytype), dimension(zsize(1),zsize(2),4*nquads), intent(inout) :: uvqm
    real(mytype), dimension(zsize(1),zsize(2)), intent(in) :: um, vm
    real(mytype), dimension(zsize(1),zsize(2)), intent(in) :: uum, vvm
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: ux, uy, ep

    logical, dimension(:,:,:,:,:) , allocatable :: mask
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: uv_fluct

    logical :: quad(4), thresh
    integer :: i, j, k, l, ierr
    real(mytype) :: u_rms, v_rms, u_fluct, v_fluct

    allocate(mask(zsize(1),zsize(2),zsize(3),nquads,4))
    uv_fluct = zero
      do k = 1, zsize(3)
        do j = 1, zsize(2)
          do i = 1, zsize(1)
            u_rms = sqrt_prec(uum(i,j) - um(i,j)*um(i,j))
            v_rms = sqrt_prec(vvm(i,j) - vm(i,j)*vm(i,j))

            u_fluct = ux(i,j,k) - um(i,j)
            v_fluct = uy(i,j,k) - vm(i,j)

            uv_fluct(i,j,k) = u_fluct*v_fluct

            quad(1) =  u_fluct .gt. zero .and. v_fluct .gt. zero
            quad(2) =  u_fluct .lt. zero .and. v_fluct .gt. zero
            quad(3) =  u_fluct .lt. zero .and. v_fluct .lt. zero
            quad(4) =  u_fluct .gt. zero .and. v_fluct .lt. zero

            do l = 1, nquads

              thresh = abs_prec(uv_fluct(i,j,k)) > h_quads(l)*u_rms*v_rms
              

              mask(i,j,k,l,1) =  quad(1) .and. thresh
              mask(i,j,k,l,2) =  quad(2) .and. thresh
              mask(i,j,k,l,3) =  quad(3) .and. thresh
              mask(i,j,k,l,4) =  quad(4) .and. thresh

          enddo
        enddo
      enddo
    enddo

    if (any(isnan(uv_fluct))) then
      write(*,*) "NaN detected in quadrant analysis"
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif 
    do j = 1, nquads
      do i = 1, 4
        call update_average_scalar(uvqm(:,:,(j-1)*4+i),uv_fluct,&
                                  ep,mask=mask(:,:,:,j,i),istat2=.true.)
      enddo
    enddo

    deallocate(mask)

  end subroutine update_uv_quad_avg
  subroutine update_pstrain_cond_avg(pdvdy_q1m,pdvdy_q2m,&
                                     pdvdy_q3m,pdvdy_q4m,&
                                     p,dvdy, pm, dvdym,ep)
    use decomp_2d, only : mytype, xsize, zsize
    use param, only : zero, zpfive
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: pdvdy_q1m, pdvdy_q2m
    real(mytype), dimension(zsize(1),zsize(2)), intent(inout) :: pdvdy_q3m, pdvdy_q4m
    real(mytype), dimension(zsize(1),zsize(2)), intent(in) :: pm, dvdym
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: p,dvdy, ep

    logical, dimension(zsize(1),zsize(2),zsize(3)) :: mask1, mask2, mask3, mask4
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: p_fluct, dvdy_fluct
    integer :: i,j,k

    do k = 1, zsize(3)
      do j = 1, zsize(2)
        do i = 1, zsize(1)
          p_fluct(i,j,k) = p(i,j,k) - pm(i,j)
          dvdy_fluct(i,j,k) = dvdy(i,j,k) - dvdym(i,j)
          mask1(i,j,k) = p_fluct(i,j,k) .gt. zero .and.  dvdy_fluct(i,j,k) .gt. zero
          mask2(i,j,k) = p_fluct(i,j,k) .lt. zero .and.  dvdy_fluct(i,j,k) .gt. zero
          mask3(i,j,k) = p_fluct(i,j,k) .lt. zero .and.  dvdy_fluct(i,j,k) .lt. zero
          mask4(i,j,k) = p_fluct(i,j,k) .gt. zero .and.  dvdy_fluct(i,j,k) .lt. zero
        enddo
      enddo
    enddo

    call update_average_scalar(pdvdy_q1m, p_fluct*dvdy_fluct, ep, mask=mask1,istat2=.true.)
    call update_average_scalar(pdvdy_q2m, p_fluct*dvdy_fluct, ep, mask=mask2,istat2=.true.)
    call update_average_scalar(pdvdy_q3m, p_fluct*dvdy_fluct, ep, mask=mask3,istat2=.true.)
    call update_average_scalar(pdvdy_q4m, p_fluct*dvdy_fluct, ep, mask=mask4,istat2=.true.)

  end subroutine update_pstrain_cond_avg

  subroutine update_lambda2_avg(lambda2m,lambda22m,&
                              dudx,dudy,dudz,dvdx,dvdy,dvdz,&
                              dwdx,dwdy,dwdz,ep)
    use decomp_2d, only : mytype, xsize, zsize
    use MPI
    use param

    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: lambda2m,lambda22m
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dudx, dudy, dudz
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dvdx, dvdy, dvdz
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: dwdx, dwdy, dwdz, ep

    real(mytype), dimension(3,3) :: S, Omega, S2O2, U_mat
    integer :: lwork, info, info_sum, nb
    real(mytype), dimension(3) :: lambda
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: lambda2
    real(mytype), dimension(:), allocatable :: work
    integer :: i,j,k,code

#ifdef HAVE_LAPACK
    interface
      integer function ilaenv(ISPEC,NAME,OPTS, N1, N2, N3, N4)
        integer :: ispec
        character(*) :: NAME, OPTS
        integer :: N1, N2, N3, N4
      end function ilaenv
    end interface

    nb = ilaenv(1,'DSYTRD','NU',3,3,-1,-1)
    lwork = (nb+2)*3
    allocate(work(lwork))

    info_sum = 0
    do k = 1, zsize(3)
      do j = 1, zsize(2)
        do i = 1, zsize(1)
          S(1,1) = dudx(i,j,k)
          S(2,2) = dvdy(i,j,k)
          S(3,3) = dwdz(i,j,k)

          S(1,2) = zpfive* ( dudy(i,j,k) + dvdx(i,j,k))
          S(2,1) = S(1,2)

          S(1,3) = zpfive* ( dudz(i,j,k) + dwdx(i,j,k))
          S(3,1) = S(1,3)

          S(3,2) = zpfive* ( dwdy(i,j,k) + dvdz(i,j,k))
          S(2,3) = S(3,2)

          Omega(1,1) = zero
          Omega(2,2) = zero
          Omega(3,3) = zero

          Omega(1,2) = zpfive* ( dudy(i,j,k) - dvdx(i,j,k))
          Omega(2,1) = -Omega(1,2)

          Omega(1,3) = zpfive* ( dudz(i,j,k) - dwdx(i,j,k))
          Omega(3,1) = -Omega(1,3)

          Omega(3,2) = zpfive* ( dwdy(i,j,k) - dvdz(i,j,k))
          Omega(2,3) = -Omega(3,2)

          S2O2 = matmul(S,S) + matmul(Omega,Omega)

          call dsyev('N','U',3,S2O2,U_mat,lambda,work,lwork,info)
          lambda2(i,j,k) = lambda(2)
          info_sum = info_sum + info
        enddo
      enddo
    enddo

    if (info_sum/=0) then
      write(*,*) "Problems calculating eigenvalules"
      call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif

    call update_average_scalar(lambda2m, lambda2, ep)
    call update_average_scalar(lambda22m, lambda2*lambda2, ep)

#endif    
  end subroutine    
  
  subroutine update_spectra_avg(spec_xm,spec_zm, ux1, uy1, uz1, ux3, uy3, uz3)
    use MPI
    use decomp_2d
    use variables, only : nx, nz

    real(mytype), dimension(nx/2,xsize(2),nspectra) :: spec_xm
    real(mytype), dimension(zsize(1),zsize(2),spectra_size) :: spec_zm

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)) :: ux3, uy3, uz3
    real(mytype), dimension(zsize(1),zsize(2),4) :: um
#ifdef HAVE_FFTW

    complex(mytype), dimension(nx/2,xsize(2),xsize(3)) :: u_spec_x, vspec_x, w_spec_x
    complex(mytype), dimension(zsize(1),zsize(2),nz/2) :: u_spec_z, vspec_z, w_spec_z

    complex(mytype), dimension(nx/2,xsize(2),nspectra) :: spect_x_avg_z, spect_x_avg_z_global
    complex(mytype), dimension(xsize(2),nz/2,spectra_nlocs,nspectra) :: spect_z_avg_z, spect_z_avg_z_global
    
    integer :: i, j, k, code, arr_type, comm_rank, fh
    integer :: color, key
    integer, dimension(2) :: dims, coords
    logical, dimension(2) :: periods 
    integer :: split_comm_y, split_comm_z
    real(mytype) :: stat_inc
    
    spect_x_avg_z = zero
    spect_z_avg_z = zero

    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, periods, coords, ierr)
    key = coords(1)
    color = coords(2)
 
    call MPI_Comm_split(DECOMP_2D_COMM_CART_X, color, key, split_comm_y,ierr)
    call MPI_Comm_split(DECOMP_2D_COMM_CART_X, key, color, split_comm_z,ierr)


      
    call dfftw_execute_r2c(plan_x,ux1,u_spec_x)
    call dfftw_execute_r2c(plan_x,uy1,v_spec_x)
    call dfftw_execute_r2c(plan_x,uz1,w_spec_x)

    do k = 1, xsize(3)
      do j =1, xsize(2)
        do i = 1, nx/2
          spect_x_avg_z(i,j,1) = spect_x_avg_z(i,j,1) &
                  + u_spec_x(i,j,k)*conjg(u_spec_x(i,j,k))
          spect_x_avg_z(i,j,2) = spect_x_avg_z(i,j,2) &
                  + v_spec_x(i,j,k)*conjg(v_spec_x(i,j,k))                    
          spect_x_avg_z(i,j,3) = spect_x_avg_z(i,j,3) &
                  + w_spec_x(i,j,k)*conjg(w_spec_x(i,j,k))                    
          spect_x_avg_z(i,j,4) = spect_x_avg_z(i,j,4) &
                  + u_spec_x(i,j,k)*conjg(v_spec_x(i,j,k))                    
        enddo
      enddo
    enddo      

    call dfftw_execute_r2c(plan_z,ux3,u_spec_z)
    call dfftw_execute_r2c(plan_z,uy3,v_spec_z)
    call dfftw_execute_r2c(plan_z,uz3,w_spec_z)

    spect_z_avg_z = zero
    if (spectra_nlocs == 0) then
        do k = 1, nz/2
          do j = 1, zsize(2)
            do i = 1, zsize(1)

            spect_z_avg_x(j,k,1,1) = spect_z_avg_z(j,k,1,1) &
                      +  u_spec_z(i,j,k)*conjg(u_spec_z(i,j,k))/xsize(1)
            spect_z_avg_x(j,k,1,2) = spect_z_avg_z(j,k,1,2) &
                      +  v_spec_z(i,j,k)*conjg(v_spec_z(i,j,k))/xsize(1)             
            spect_z_avg_x(j,k,1,3) = spect_z_avg_z(j,k,1,3) &
                      +  w_spec_z(i,j,k)*conjg(w_spec_z(i,j,k))/xsize(1)
            spect_z_avg_x(j,k,1,4) = spect_z_avg_z(j,k,1,4) &
                      +  u_spec_z(i,j,k)*conjg(v_spec_z(i,j,k))/xsize(1)
          enddo
        enddo
      enddo

      call MPI_Allreduce(spect_z_avg_x,spect_z_avg_x_global,size(spect_z_avg_x),&
                         complex_type,MPI_SUM,split_comm_y,code)    else

    endif

    if (spectra_level == 2) then
      write(*,*) "Not implemented yet"
      call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif
    spect_x_avg_z = spect_x_avg_z / real(zsize(3),kind=mytype)
    
    call MPI_Allreduce(spect_x_avg_z,spect_x_avg_z_global,size(spect_x_avg_z),&
                        complex_type,MPI_SUM,split_comm_y,code)

    stat_inc = real((itime-initstat2)/istatcalc+1, kind=mytype)
    spectra_x_mean = spectra_x_mean &
                    + (spect_x_avg_z_global - spectra_x_mean) &
                    / stat_inc


#endif    
  end subroutine
  subroutine grad_init
    use param, only : zpfive, one, three,onepfive, two
    use var, only : yp, ny, dx, dz

    real(mytype) :: ddy1, ddy2
    integer :: i

    allocate(coefy(2:ny-1,3))

    do i = 2, ny-1
      ddy1 = yp(i) - yp(i-1)
      ddy2 = yp(i+1) - yp(i)

      coefy(i,1) = -(ddy2)/(ddy1 * (ddy1 + ddy2))
      coefy(i,2) = (ddy2 - ddy1) / (ddy1 * ddy2)
      coefy(i,3) = ddy1 / (ddy2 * (ddy1 + ddy2))
    enddo

    coefx_b = - one / dx
    coefx_f = one /dx

    coefz_b = - one/dz
    coefz_f = one/dz

    ddy1 = yp(2) - yp(1)
    ddy2 = yp(3) - yp(2)

    coefy_bc1(1) = -( two*ddy1 + ddy2)/(ddy1*(ddy1 + ddy2))
    coefy_bc1(2) = (ddy1 + ddy2)/(ddy1*ddy2)
    coefy_bc1(3) = -ddy1/(ddy2*(ddy1 + ddy2))

    ddy1 = yp(ny-1) - yp(ny-2)
    ddy2 = yp(ny-2) - yp(ny-3)

    coefy_bcn(1) = ddy2/(ddy1*(ddy1 + ddy2))
    coefy_bcn(2) = -(ddy1 + ddy2)/(ddy1*ddy2)
    coefy_bcn(3) = (two*ddy2 + ddy1)/(ddy2*(ddy1 + ddy2))

    coefx_bc1(1) = - onepfive/ dx
    coefx_bc1(2) = two / dx
    coefx_bc1(3) = - zpfive / dx

    coefx_bcn(1) = zpfive / dx
    coefx_bcn(2) = -two / dx
    coefx_bcn(3) = onepfive / dx

    coefz_bc1(1) = - onepfive/ dz
    coefz_bc1(2) = two / dz
    coefz_bc1(3) = - zpfive / dz

    coefz_bcn(1) = zpfive / dz
    coefz_bcn(2) = -two / dz
    coefz_bcn(3) = onepfive / dz
  end subroutine

  subroutine grad_x(array, grad)
    use decomp_2d, only : xsize
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: array
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: grad

    integer :: i, j, k


    do k = 1, xsize(3)
      do j = 1, xsize(2)

        grad(1,j,k) = coefx_bc1(1)*array(1,j,k) + coefx_bc1(2)*array(2,j,k) &
                       + coefx_bc1(3)*array(3,j,k)

        do i = 2,xsize(1)-1
          grad(i,j,k) = coefx_b*array(i-1,j,k) + coefx_f*array(i+1,j,k)
        enddo

        grad(xsize(1),j,k) = coefx_bcn(1)*array(xsize(1)-2,j,k) + &
                      coefx_bcn(2)*array(xsize(1)-1,j,k) + &
                      coefx_bcn(3)*array(xsize(1),j,k)
      enddo
    enddo

  end subroutine grad_x

  subroutine grad_y(array, grad)
    use decomp_2d, only : ysize
    real(mytype), dimension(ysize(1),ysize(2),ysize(3)), intent(in) :: array
    real(mytype), dimension(ysize(1),ysize(2),ysize(3)), intent(out) :: grad
    integer :: i, j, k

    do k = 1, ysize(3)
      do j = 2, ysize(2)-1
        do i = 1,ysize(1)
          grad(i,j,k) = coefy(j,1)*array(i,j-1,k) &
                        + coefy(j,2)*array(i,j,k) &
                        + coefy(j,3)*array(i,j+1,k)
        enddo
      enddo

      do i = 1, ysize(1)
        grad(i,1,k) = coefy_bc1(1)*array(i,1,k) + coefy_bc1(2)*array(i,2,k) &
                        + coefy_bc1(3)*array(i,3,k)

        grad(i,ysize(2),k) = coefy_bcn(1)*array(i,ysize(2)-2,k) + &
                              coefy_bcn(2)*array(i,ysize(2)-1,k) + &
                              coefy_bcn(3)*array(i,ysize(2),k)
      enddo
    enddo
  end subroutine grad_y

  subroutine grad_z(array, grad)
    use decomp_2d, only : zsize
    
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(in) :: array
    real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(out) :: grad
    integer :: i, j, k

    do k = 2, zsize(3)-1
      do j = 1, zsize(2)
        do i = 1,zsize(1)
          grad(i,j,k) = coefz_b*array(i,j,k-1) + coefz_f*array(i,j,k+1)
        enddo
      enddo
    enddo

    do j = 1, zsize(2)
      do i = 1,zsize(1)
        grad(i,j,1) = coefz_bc1(1)*array(i,j,1) + coefz_bc1(2)*array(i,j,2) &
                        + coefz_bc1(3)*array(i,j,3)       

        grad(i,j,zsize(3)) = coefz_bcn(1)*array(i,j,zsize(3)-2) + &
                        coefz_bcn(2)*array(i,j,zsize(3)-1) + &
                        coefz_bcn(3)*array(i,j,zsize(3))
      enddo
    enddo
  end subroutine grad_z
endmodule stats

