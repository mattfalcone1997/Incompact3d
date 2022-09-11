!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module stats
  use param, only : mytype
  use decomp_2d, only : DECOMP_INFO
  implicit none

  character(len=*), parameter :: io_statistics = "statistics-io", &
       stat_dir = "statistics"

  integer :: stats_time

  real(mytype), dimension(:,:,:), allocatable ::  pmean
  real(mytype), dimension(:,:,:), allocatable ::  umean, uumean
  real(mytype), dimension(:,:,:), allocatable ::  vmean, vvmean
  real(mytype), dimension(:,:,:), allocatable ::  wmean, wwmean
  real(mytype), dimension(:,:,:), allocatable ::  uvmean, uwmean
  real(mytype), dimension(:,:,:), allocatable ::  vwmean
  real(mytype), dimension(:,:,:), allocatable ::  phimean, phiphimean
  type(DECOMP_INFO) :: dstat_plane

  private
  public overall_statistic

contains

  subroutine init_statistic_adios2

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io

    use var, only : numscalar
    
    implicit none

    integer :: ierror
    character(len=30) :: varname
    integer :: is
    logical, save :: initialised = .false.

    if (.not. initialised) then
       call decomp_2d_init_io(io_statistics)
    
       call decomp_2d_register_variable(io_statistics, "umean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "vmean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "wmean", 3, 1, 3, mytype, opt_nplanes=1)
       
       call decomp_2d_register_variable(io_statistics, "pmean", 3, 1, 3, mytype, opt_nplanes=1)
       
       call decomp_2d_register_variable(io_statistics, "uumean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "vvmean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "wwmean", 3, 1, 3, mytype, opt_nplanes=1)

       call decomp_2d_register_variable(io_statistics, "uvmean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "uwmean", 3, 1, 3, mytype, opt_nplanes=1)
       call decomp_2d_register_variable(io_statistics, "vwmean", 3, 1, 3, mytype, opt_nplanes=1)

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

    use param, only : zero, iscalar
    use variables, only : nx, ny, nz
    use decomp_2d, only : zsize, decomp_info_init
    

    implicit none

    allocate(umean(zsize(1),zsize(2),1))
    allocate(vmean(zsize(1),zsize(2),1))
    allocate(wmean(zsize(1),zsize(2),1))
    allocate(pmean(zsize(1),zsize(2),1))
    allocate(uumean(zsize(1),zsize(2),1))
    allocate(vvmean(zsize(1),zsize(2),1))
    allocate(wwmean(zsize(1),zsize(2),1))
    allocate(uvmean(zsize(1),zsize(2),1))
    allocate(uwmean(zsize(1),zsize(2),1))
    allocate(vwmean(zsize(1),zsize(2),1))
    allocate(phimean(zsize(1),zsize(2),1))
    allocate(phiphimean(zsize(1),zsize(2),1))

    
    pmean = zero
    umean = zero
    uumean = zero
    vmean = zero
    vvmean = zero
    wmean = zero
    wwmean = zero
    uvmean = zero
    uwmean = zero
    vwmean = zero

    call decomp_info_init(nx, ny, 1,dstat_plane)

    write(*,*) dstat_plane%xsz
    write(*,*) dstat_plane%ysz
    write(*,*) dstat_plane%zsz
    call init_statistic_adios2

  end subroutine init_statistic

  !
  ! Read all statistics if possible
  !
  subroutine restart_statistic

    use param, only : initstat, irestart, ifirst, zero
    use variables, only : nstat
    use var, only : tmean

    implicit none

    ! No reading for statistics when nstat > 1 or no restart
    if (nstat.gt.1 .or. irestart.eq.0) then
       call init_statistic()
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

    use param, only : iscalar, itime
    use variables, only : numscalar
    use decomp_2d, only : nrank
    use decomp_2d_io, only : decomp_2d_write_mode, decomp_2d_read_mode, &
         decomp_2d_open_io, decomp_2d_close_io, decomp_2d_start_io, decomp_2d_end_io
    use var, only : pmean
    use var, only : umean, uumean
    use var, only : vmean, vvmean
    use var, only : wmean, wwmean
    use var, only : uvmean, uwmean
    use var, only : vwmean
    use var, only : phimean, phiphimean

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
    
    call read_or_write_one_stat(flag_read, gen_statname("pmean"), pmean)
    call read_or_write_one_stat(flag_read, gen_statname("umean"), umean)
    call read_or_write_one_stat(flag_read, gen_statname("vmean"), vmean)
    call read_or_write_one_stat(flag_read, gen_statname("wmean"), wmean)

    call read_or_write_one_stat(flag_read, gen_statname("uumean"), uumean)
    call read_or_write_one_stat(flag_read, gen_statname("vvmean"), vvmean)
    call read_or_write_one_stat(flag_read, gen_statname("wwmean"), wwmean)

    call read_or_write_one_stat(flag_read, gen_statname("uvmean"), uvmean)
    call read_or_write_one_stat(flag_read, gen_statname("uwmean"), uwmean)
    call read_or_write_one_stat(flag_read, gen_statname("vwmean"), vwmean)

    if (iscalar==1) then
       do is=1, numscalar
          write(filename,"('phi',I2.2)") is
          call read_or_write_one_stat(flag_read, gen_statname(trim(filename)), phimean(:,:,:,is))
          write(filename,"('phiphi',I2.2)") is
          call read_or_write_one_stat(flag_read, gen_statname(trim(filename)), phiphimean(:,:,:,is))
       enddo
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
  subroutine read_or_write_one_stat(flag_read, filename, array)

    use decomp_2d, only : mytype, xstS, xenS, zsize
    use decomp_2d_io, only : decomp_2d_read_one, decomp_2d_write_one

    implicit none

    ! Arguments
    logical, intent(in) :: flag_read
    character(len=*), intent(in) :: filename
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: array
    integer :: ierror

    if (flag_read) then
       ! There was a check for nvisu = 1 before
       call decomp_2d_read_one(3, array, stat_dir, filename, io_statistics, reduce_prec=.false.,opt_decomp=dstat_plane)
    else
       call decomp_2d_write_one(3, array, stat_dir, filename, 0, io_statistics, reduce_prec=.false.,opt_decomp=dstat_plane)
    endif

  end subroutine read_or_write_one_stat

  !
  ! Statistics : Intialize, update and perform IO
  !
  subroutine overall_statistic(ux1,uy1,uz1,phi1,pp3,ep1)

    use param
    use variables
    use decomp_2d
    use decomp_2d_io
    use tools, only : rescale_pressure

    use var, only : nxmsize, nymsize, nzmsize
    use var, only : ppi3, dip3
    use var, only : pp2, ppi2, dip2
    use var, only : pp1, ta1, di1

    use var, only : tmean
    use var, only : pmean
    use var, only : umean, uumean
    use var, only : vmean, vvmean
    use var, only : wmean, wwmean
    use var, only : uvmean, uwmean
    use var, only : vwmean
    use var, only : phimean, phiphimean

    implicit none

    !! Inputs
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar),intent(in) :: phi1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

    !! Locals
    integer :: is
    character(len=30) :: filename

    if (itime.lt.initstat) then
       return
    elseif (itime.eq.initstat) then
       call init_statistic()
    elseif (itime.eq.ifirst) then
       call restart_statistic()
    endif

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
    call interxpv(ta1,pp1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    ! Convert to physical pressure
    call rescale_pressure(ta1)

    call update_average_scalar(pmean, ta1, ep1)

    !! Mean velocity
    call update_average_vector(umean, vmean, wmean, &
                               ux1, uy1, uz1, ep1)

    !! Second-order velocity moments
    call update_variance_vector(uumean, vvmean, wwmean, uvmean, uwmean, vwmean, &
                                ux1, uy1, uz1, ep1)

    ! Write all statistics
    if (mod(itime,icheckpoint)==0) then
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
  subroutine update_average_scalar(um, ux, ep)

    use decomp_2d, only : mytype, xsize, zsize
    use param, only : itime, initstat
    use var, only : di1, tmean

    implicit none

    ! inputs
    real(mytype), dimension(zsize(1),zsize(2),1), intent(inout) :: um
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, ep

    real(mytype), dimension(zsize(1), zsize(2), 1) :: stat_z

    call average_z(ux, stat_z)
    um = um + (stat_z - um) / real(itime-initstat+1, kind=mytype)

  end subroutine update_average_scalar

  subroutine average_z(array, stat_xz)
    use decomp_2d
    use param, only : zero

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: array
    real(mytype), dimension(zsize(1),zsize(2),1) :: stat_xz

    real(mytype), dimension(:,:,:), allocatable :: array_ypencil
    real(mytype), dimension(:,:,:), allocatable :: array_zpencil
    integer :: i, j, k
    allocate(array_ypencil(ysize(1),ysize(2),ysize(3)))
    allocate(array_zpencil(zsize(1),zsize(2),zsize(3)))

    call transpose_x_to_y(array, array_ypencil)
    call transpose_y_to_z(array_ypencil, array_zpencil)

    stat_xz(:,:,:) = zero
    do i = 1, zsize(1)
      do j = 1, zsize(2)
        do k = 1, zsize(3)
          stat_xz(i,j,1) = stat_xz(i,j,1) + array_zpencil(i,j,k)/zsize(3)
        enddo
      enddo
    enddo


   end subroutine
  !
  ! Update (um, vm, wm), the average of (ux, uy, uz)
  !
  subroutine update_average_vector(um, vm, wm, ux, uy, uz, ep)

    use decomp_2d, only : mytype, xsize, xstS, xenS

    implicit none

    ! inputs
    real(mytype), dimension(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)), intent(inout) :: um, vm, wm
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(um, ux, ep)
    call update_average_scalar(vm, uy, ep)
    call update_average_scalar(wm, uz, ep)

  end subroutine update_average_vector

  !
  ! Update (uum, vvm, wwm, uvm, uwm, vwm), the variance of the vector (ux, uy, uz)
  !
  subroutine update_variance_vector(uum, vvm, wwm, uvm, uwm, vwm, ux, uy, uz, ep)

    use decomp_2d, only : mytype, xsize, xstS, xenS

    implicit none

    ! inputs
    real(mytype), dimension(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)), intent(inout) :: uum, vvm, wwm, uvm, uwm, vwm
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, uy, uz, ep

    call update_average_scalar(uum, ux*ux, ep)
    call update_average_scalar(vvm, uy*uy, ep)
    call update_average_scalar(wwm, uz*uz, ep)
    call update_average_scalar(uvm, ux*uy, ep)
    call update_average_scalar(uwm, ux*uz, ep)
    call update_average_scalar(vwm, uy*uz, ep)

  end subroutine update_variance_vector

endmodule stats

